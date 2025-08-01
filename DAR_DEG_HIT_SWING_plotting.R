options(scipen=99999)
options(stringsAsFactors=F)
options(digits=7)

library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
library(patchwork)


outdir="DARDEG"

#chrOrder=paste0("chr", c(1:19, "X"))
chrOrder=paste0("chr", c(1:19))



###################
# get HITS for those regions

g.background = read_tsv("./data/genome_windows_30000.bed", col_names = c("chrom", "start", "end")) %>% 
  group_by(chrom) %>% 
  summarise(glen = plyr::round_any(max(end), 30000, floor))
g.background %>% count(chrom) %>% print(n=22)



makesym = function(lst){
  lst2=lst %>% rename(chrom1=chrom2, start1=start2, end1=end2, 
                      chrom2=chrom1, start2=start1, end2=end1)
  return(bind_rows(lst, lst2))
}

makesym.win = function(lst){
  lst2=lst %>% rename(window1=window2, 
                      window2=window1)
  return(bind_rows(lst, lst2))
}

# extract random shifts from file names
randids = sort(as.numeric((tibble(fl=list.files(sprintf("./results/%s/randomregions/", outdir), pattern = "*.tsv")) %>% separate(fl, into=c("src", "rand", "ext"),remove = F) )$rand))

# read counts after shifting
get_hit_counts = function(fn, shifted_window_mapping.sel){
  return(read_tsv(fn, show_col_types = FALSE) %>%
           #add ID so we remove symmetric connections within same regions (encoded by contact row number)
           mutate(contactID=row_number()) %>% 
           #get region information    
           makesym.win() %>%
           left_join(shifted_window_mapping.sel, by=join_by(window1==window)) %>% 
           left_join(shifted_window_mapping.sel, by=join_by(window2==window)) %>% 
           # counting
           #remove duplicated connections within same group
           distinct(grpstr.x, grpstr.y, contactID, .keep_all = T) %>% #
           group_by(grpstr.x) %>% 
           summarize(count=sum(!is.na(NPMI.y)))  )
}


# overlap contact counts with DAR/DEG regions
contact_count_pergroup = c()
for(i in 1:length(randids)){
  k=randids[[i]]
  print(paste(i, k))
  
  shifted_window_mapping.sel = read_tsv(sprintf("./results/%s/randomregions/dardegs_%d.tsv", outdir, k)) %>% 
    mutate(window = sprintf("%s:%s-%s", chrom, start.adj, end.adj)) %>% 
    select(window, grpstr, grp)
  
  sd_hit_counts = get_hit_counts(sprintf("./intermediate/contacts_in_regions_20240115/sd_dardegs_%s.txt", k), shifted_window_mapping.sel) %>% 
    rename(count.sd=count)
  hc_hit_counts = get_hit_counts(sprintf("./intermediate/contacts_in_regions_20240115/wt_dardegs_%s.txt", k), shifted_window_mapping.sel) %>% 
    rename(count.hc=count)
  
  combined = shifted_window_mapping.sel %>% 
    distinct(grpstr, grp) %>% 
    left_join(sd_hit_counts, by=join_by(grpstr==grpstr.x)) %>% 
    left_join(hc_hit_counts, by=join_by(grpstr==grpstr.x)) %>% 
    mutate(count.sd=replace_na(count.sd, 0)) %>% 
    mutate(count.hc=replace_na(count.hc, 0)) %>%
    ungroup() %>% 
    mutate(
      rand=k,
      runid=i
    )
  
  contact_count_pergroup = contact_count_pergroup %>% bind_rows(combined) 
}

#rewrite location strings
contact_count_pergroup = contact_count_pergroup %>% separate(grpstr, sep="[:-]", into=c("chrom", "start", "end"), remove = F)
contact_count_pergroup = contact_count_pergroup %>% mutate(start=as.numeric(start),
                                                           end=as.numeric(end))


contact_count_pergroup %>% write_tsv(sprintf("./results/%s/fullist.tsv.gz", outdir))



###################

MIN.COUNT=10
DEC.CUT =0.8


compar=function(a,b){ return(a/pmax(b,1))}


hit_entrichments = contact_count_pergroup %>% mutate(sd_hc_ratio=compar(count.sd, count.hc),
                                                     hc_sd_ratio=compar(count.hc, count.sd),
                                                     count_both = count.hc+count.sd)

# get decile values from counts in randomized sets
thrs = hit_entrichments %>%
  filter(!rand==0) %>% 
  group_by(chrom, grp, grpstr, start) %>%
  summarize(thrs.sd_hc = quantile(sd_hc_ratio, DEC.CUT),
            thrs.hc_sd = quantile(hc_sd_ratio, DEC.CUT),
            thrs.swingcounts = quantile(count_both, DEC.CUT))

thrs %>% summary()

# get counts from observed data - rand 0 for 0 bp shifting
hit_entrichments_obs = hit_entrichments %>% 
  filter(rand==0) %>% 
  left_join(thrs) %>% 
  mutate(sd_hit_region = sd_hc_ratio>thrs.sd_hc & count_both>=MIN.COUNT,
         hc_hit_region = hc_sd_ratio>thrs.sd_hc & count_both>=MIN.COUNT,
         swing_region = (count_both>=thrs.swingcounts)& (count_both>MIN.COUNT),
         swing_only_region = swing_region & !sd_hit_region &!hc_hit_region
  )

#hit_entrichments_obs = hit_entrichments_obs %>% separate(grpstr, sep = "[:-]", into = c("chrom","start","stop"), remove=F)
hit_entrichments_obs %>% write_tsv(sprintf("./results/%s/region_annotations_%s.tsv.gz", outdir, DEC.CUT))


hit_entrichments_obs %>% 
  filter(sd_hit_region) %>% select(chrom, start, end, sd_hc_ratio) %>% 
  write_tsv(sprintf("./results/%s/SD_hit_regions_p%s.bed", outdir, DEC.CUT))
hit_entrichments_obs %>% 
  filter(hc_hit_region) %>% select(chrom, start, end, hc_sd_ratio) %>% 
  write_tsv(sprintf("./results/%s/HC_hit_regions_p%s.bed", outdir, DEC.CUT))
hit_entrichments_obs %>% 
  filter(swing_only_region) %>% select(chrom, start, end, count_both) %>% 
  write_tsv(sprintf("./results/%s/SWING_regions_p%s.bed", outdir, DEC.CUT))



print(hit_entrichments_obs %>% count(swing_region))

print(hit_entrichments_obs %>% count(hc_hit_region, sd_hit_region, swing_only_region))
print(hit_entrichments_obs %>% count(hc_hit_region, sd_hit_region))
print(hit_entrichments_obs %>% filter(hc_hit_region==sd_hit_region, hc_hit_region==T))




###################
#
# Plot data
#


# get summed-up chromosome lengths for manhattan-like plot
g.background.pos = read_tsv("./data/genome_windows_30000.bed", col_names = c("chrom", "start", "end")) %>% group_by(chrom) %>% summarise(gstart = max(end)) %>% 
  filter(chrom %in% chrOrder) %>%
  arrange(factor(chrom, levels=chrOrder)) %>%
  mutate(gstart.sum = lag(cumsum(gstart), default=0))

#hit_entrichments_obs = read_tsv(sprintf("./results/%s/region_annotations_%s.tsv.gz", outdir, k))


pltdata=
  hit_entrichments_obs %>% 
  left_join(g.background.pos, by = join_by(chrom==chrom)) %>%
  mutate(relstart = start+gstart.sum,
         relstop = end+gstart.sum  )


plot_selected_swing = ggplot()+
  #geom_point(aes(x=relstart, y=count_both, color=paste(sd_hit_region, hc_hit_region, swing_only_region), size=ifelse(swing_region>0, 20,0)))+
  geom_point(data=pltdata %>% filter(swing_region==0), aes(x=relstart, y=count_both), color="#aaaaaa99", size=2, alpha=.7)+
  geom_point(data=pltdata %>% filter(sd_hit_region |hc_hit_region | swing_region), aes(x=relstart, y=count_both, color=paste(sd_hit_region, hc_hit_region, swing_only_region)), size=4, alpha=.7)+
  geom_vline(xintercept = g.background.pos$gstart.sum, color="grey", linetype=2) +
  #  geom_text(data=g.background, aes(label=c(1:19), x=gstart.sum+50000000, y=-15)) +
  #  ggtitle("Position and scores of selected regions")+
  theme_bw()+
  scale_color_manual(values=list(#"FALSE FALSE FALSE"="#aaaaaa99", 
                                 "FALSE TRUE FALSE"="#fbb040",
                                 "TRUE FALSE FALSE"="#2a9d8f",
                                 "FALSE FALSE TRUE"="#9e1f63"
  ),
  #labels=c("not selected","SWING HIT", "HC HIT region", "SD HIT region"))+
  labels=c("SWING HIT", "HC HIT region", "SD HIT region"))+
  #geom_text(data=g.background.pos, aes(label=paste0("Chr", c(1:19)), x=gstart.sum+50000000, y=-30), size=5)+
  geom_text(data=g.background.pos, aes(label=c(1:19), x=gstart.sum+50000000, y=-30), size=5)+
  scale_size_binned_area(
    limits = c(0, 100),
    breaks = seq(0,100, 10)
  )+
  geom_hline(yintercept = 0)+
  guides(size="none")+
  theme(legend.position = "right",
        text=element_text(size=20))+
  ylab("hc+sd count")+
  xlab("Genomic position")+
  labs(color=element_text("Category"))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()
        )


plot_selected_swing
ggsave(filename = sprintf("./results/%s/region_annotations_%s_sum_250201.pdf", outdir, DEC.CUT), height = 8, width=15)




pltdata2=
  hit_entrichments_obs %>% 
  mutate(hc_sd_ratio.s = ifelse(hc_sd_ratio>sd_hc_ratio, hc_sd_ratio, NA),
         sd_hc_ratio.s = ifelse(hc_sd_ratio<sd_hc_ratio, sd_hc_ratio, NA)) %>% 
  left_join(g.background.pos, by = join_by(chrom==chrom)) %>%
  mutate(relstart = start+gstart.sum,
         relstop = end+gstart.sum
  ) 


plot_selected_hits =   ggplot()+
  geom_point(data=pltdata2 %>% filter(!sd_hit_region, !hc_hit_region, !swing_only_region, !is.na(hc_sd_ratio.s)) , aes(x=relstart, y=hc_sd_ratio.s), color="#aaaaaa99", size=2, alpha=.7)+
  geom_point(data=pltdata2 %>% filter(!sd_hit_region, !hc_hit_region, !swing_only_region, !is.na(sd_hc_ratio.s)) , aes(x=relstart, y=-sd_hc_ratio.s), color="#aaaaaa99", size=2, alpha=.7)+
  geom_point(data=pltdata2 %>% filter(sd_hit_region+hc_hit_region+swing_only_region>0, !is.na(hc_sd_ratio.s)), aes(x=relstart, y=hc_sd_ratio.s, color=paste(sd_hit_region, hc_hit_region, swing_only_region)), size=4, alpha=.7)+
  geom_point(data=pltdata2 %>% filter(sd_hit_region+hc_hit_region+swing_only_region>0, !is.na(sd_hc_ratio.s)), aes(x=relstart, y=-sd_hc_ratio.s, color=paste(sd_hit_region, hc_hit_region, swing_only_region)), size=4, alpha=.7)+
  geom_vline(xintercept = g.background.pos$gstart.sum, color="grey", linetype=2) +
  theme_bw()+
  scale_color_manual(values=list(#"FALSE FALSE FALSE"="#aaaaaa99", 
    "FALSE TRUE FALSE"="#fbb040",
    "TRUE FALSE FALSE"="#2a9d8f",
    "FALSE FALSE TRUE"="#9e1f63"
  ),
  labels=c("SWING HIT", "HC HIT region", "SD HIT region"))+ #"not selected",
  ylim(c(-50,110))+ 
  geom_text(data=g.background.pos, aes(label=c(1:19), x=gstart.sum+50000000, y=-50))+
#  scale_size_binned_area(
#    limits = c(0, 110),
#    breaks = seq(0,110, 10)
#  )+
  guides(size="none")+
  theme(legend.position = "right",
        text=element_text(size=20))+
  ylab("HC/SD ratio")+
  xlab("Genomic position")+
  labs(color=element_text("Category"))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()
  )


plot_selected_hits
ggsave(filename = sprintf("./results/%s/region_annotations_%s_ratio_250201.pdf", outdir, DEC.CUT), height = 8, width=15)


