rm(list=ls())

library(tidyverse)
library(qtl2)
grav2 <- read_cross2("./data/lt_qtl_3-10-21.yaml")

# 1 
map <- insert_pseudomarkers(grav2$gmap, step=1)

# 2 
pr <- calc_genoprob(grav2, map, error_prob=0.002)

# 3

#kinship <- calc_kinship(pr)

# 4 - set up pseudo markers and kinship

grid <- calc_grid(grav2$gmap, step=1)

pr_grid <- probs_to_grid(pr, grid)
kinship_grid <- calc_kinship(pr_grid)

# 5 - run scan

out <- scan1(pr_grid, grav2$pheno, kinship = kinship_grid)



# 6 - look for peaks

peaks <- find_peaks(out, map, threshold=3, peakdrop=3, drop=2)
peaks2.9 <- find_peaks(out, map, threshold=2.9, peakdrop=3, drop=2)


peaks95 <- find_peaks(out, map, prob = 0.95, peakdrop = 2, threshold = 3)



# 7 p-value

operm2 <- scan1perm(pr, grav2$pheno, n_perm=100, kinship = kinship_grid)
sig <- summary(operm, alpha=c(0.2, 0.053))


# 8 Plot output

# par(mar=c(5.1, 4.1, 1.1, 1.1))
# ymx <- maxlod(out) # overall maximum LOD score
# plot(out, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
# plot(out, map, lodcolumn=2, col="violetred", add=TRUE)
# legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out), bg="gray90")

map_df <- data.frame(unlist(map)) %>% 
  add_rownames(var="LG.locus") %>% 
  separate(LG.locus, c("LG", "locus"), extra = "merge", fill ="left")

out_df <- as.data.frame(out) %>% rownames_to_column(var = "locus")

out_map <- left_join(out_df, map_df, by = "locus")


## reorder df to go from LG1 to LG42

out_map$LG <-fct_reorder(out_map$LG, as.numeric(out_map$LG))
colors <- rep(c("darkgray","steelblue"), 21)

### RDA cands

RDA_cands <- read.csv("./data/cands_1-8-21_(3-11edit).csv")

map_3.10.21 <- read.csv("./data/map_3-10-21.csv") #file made using rqtl_file_formater


CM_3.10.21 <- map_3.10.21 %>% filter(LocusID %in% RDA_cands$SNPID) %>% select(MappedChromosome ,LocusID, pred_cm)
RDA_cands <- left_join(RDA_cands, CM_3.10.21, by =c("SNPID" = "LocusID"))


RDA_cands$LG <-  RDA_cands$LakeTroutLG
RDA_cands$unlist.map. <- RDA_cands$CM

RDA_LIPID_LOCI <- map_df %>% filter(locus %in% RDA_cands$SNPID[1:13])
RDA_COND_LOCI <- map_df %>% filter(locus %in% RDA_cands$SNPID[14:16])




## length
# plot_len <- ggplot(out_map, aes(x=unlist.map., y = meanLen, color = LG))+
#   geom_line()+
#   geom_hline(yintercept = 3, color = "purple")+
#   geom_hline(yintercept = sig[4], color = "red", linetype = "dashed")+
#   facet_grid(~LG, scales = "free")+
#   scale_color_manual(values =colors)+
#   theme_classic()+
#   labs(x="cM", y = "LOD", title="Length")+
#   theme(panel.spacing = unit(0, "lines"),
#         legend.position = "none",
#         axis.text.x = element_text(size = 0))

## fat
plot_fat <- ggplot(out_map, aes(x=unlist.map., y = meanFat, color = LG))+
 #geom_vline(data = RDA_LIPID_LOCI, aes(xintercept = unlist.map.), size = 3, alpha = .25, color ="black")+
  geom_line(alpha = 0.8)+
  geom_hline(yintercept = 3, color = "purple")+
  geom_hline(yintercept = sig[2], color = "red", linetype = "dashed")+
  facet_grid(~factor(LG), scales = "free")+
  scale_color_manual(values =colors)+
  theme_classic()+
  labs(x="cM", y = "LOD", title="Lipid content")+
  theme(panel.spacing = unit(0, "lines"), 
        legend.position = "none",
        axis.text.x = element_text(size = 0))

## condition
plot_condition <- ggplot(out_map, aes(x=unlist.map., y = meanK, color = LG))+
  geom_line(alpha = 0.8)+
  geom_hline(yintercept = 3, color = "purple")+
  geom_hline(yintercept = sig[6], color = "red", linetype = "dashed")+
  facet_grid(~factor(LG), scales = "free")+
  scale_color_manual(values =colors)+
  theme_classic()+
  labs(x="cM", y = "LOD", title = "Condition")+
  theme(panel.spacing = unit(0, "lines"), 
        legend.position = "none",
        axis.text.x = element_text(size = 0))


plot_fat+ geom_point(data = RDA_LIPID_LOCI, aes(x = unlist.map., y = 3), size = 3)

 
 plot_condition + geom_point(data = RDA_COND_LOCI, aes(x = unlist.map., y = 3), size = 3)+
   
plot_fat<- plot_fat+ geom_vline(data = RDA_LIPID_LOCI, aes(xintercept = unlist.map.), size = 3, alpha = .25, color ="gray")+ geom_point(data = RDA_LIPID_LOCI, aes(x = unlist.map., y = 3),position = "jitter" )

plot_condition <- plot_condition+ geom_vline(data = RDA_COND_LOCI, aes(xintercept = unlist.map.),size =3, alpha = .25, color ="gray")+ geom_point(data = RDA_COND_LOCI, aes(x = unlist.map., y = 3),position = "jitter") 




ggpubr::ggarrange(plot_fat, plot_condition, ncol =1, nrow=2)
ggsave("lt_qtl_plot.png")



########### publication figures

chr <- out_map$LG

pos <- out_map$unlist.map


continuous_pos <- function(chr, pos){
  df <- data.frame(LG = chr, cm=pos)
  
  pos_summary <- df %>% group_by(LG) %>% summarise(max(cm)) %>% arrange(as.numeric(LG))
  
  cont_vals= NULL
  lg=NULL
  
  for(i in unique(df$LG)){
    
    if(i == 1){
      new_vals=df[df$LG ==i,"cm"]
    } else {
      vals = df[df$LG ==i,"cm"]
      new_vals = vals+sum(pos_summary[1:i-1,2])
    }
    
    
    cont_vals = c(cont_vals, new_vals)
    lg = c(lg, rep(i, length(new_vals)))
    
  }
  
  data.frame(lg, cont_vals)
  
}

continuous_map <- continuous_pos(out_map$LG, out_map$unlist.map)

out_map_cont <- bind_cols(out_map, continuous_map)
colorDF <- data.frame(LG=as.character(1:42), color = colors)
out_map_cont <- left_join(out_map_cont, colorDF, by = "LG")






### data saved for figures:
#  out_map_cont.RData, sig.RData,


out_map_cont <- read.csv("./data/out_map_3-10-21.csv")


## LG labels

labs <- out_map_cont %>% group_by(LG) %>% summarise(LG_lab= max(cont_vals)-((max(cont_vals)-min(cont_vals))/2)) %>% arrange(as.numeric(LG))

labs$LG_b <- c("1", "", "3", "","5", "","7", "","9", "","11", "","13", "","15", "","17", "","19", "","21", "","23", "","25", "","27", "","29", "","31", "","33", "","35", "","37", "","39", "","41", "")



colors <- rep(c("darkgray","steelblue"), 21)

map_df_cont <- continuous_pos(map_df$LG, map_df$unlist.map.)
map_df_cont <- bind_cols(map_df,map_df_cont)
map_df_cont <- left_join(map_df_cont, colorDF, by = "LG")


RDA_LIPID_LOCI <- map_df_cont %>% filter(locus %in% RDA_cands$SNPID[1:13])
RDA_COND_LOCI <- map_df_cont %>% filter(locus %in% RDA_cands$SNPID[14:16])
## ADD LOD
RDA_LIPID_LOCI <- left_join(RDA_LIPID_LOCI, RDA_cands[1:13,c("SNPID", "LOD")], by = c("locus"="SNPID"))
RDA_COND_LOCI <- left_join(RDA_COND_LOCI, RDA_cands[14:16,c("SNPID", "LOD")], by = c("locus"="SNPID"))

RDA_LIPID_LOCI_labs <- RDA_LIPID_LOCI %>% filter(locus %in% c("CLocus_42626", "CLocus_91657" ,"CLocus_23471" , "CLocus_47288" ))




plot_fat <- ggplot(out_map_cont %>% filter(lg ==38), aes(x=cont_vals, y = meanFat))+
  #geom_vline(data = RDA_LIPID_LOCI, aes(xintercept = unlist.map.), size = 3, alpha = .25, color ="black")+
  geom_line(alpha = 0.8, color = out_map_cont$color)+
  #scale_color_manual(values = c("darkgray","steelblue"))+
  geom_hline(yintercept = 3, color = "purple")+
  geom_hline(yintercept = sig[2], color = "red", linetype = "dashed")+
  #facet_grid(~factor(LG), scales = "free")+
  scale_x_continuous(breaks = labs$LG_lab, labels = labs$LG_b)+
  theme_classic()+
  labs(x="", y = "LOD", title="")+
  theme(panel.spacing = unit(0, "lines"), 
        legend.position = "none",
        axis.text.x = element_text(size = 6))



RDA_LIPID_LOCI_labs$cont_vals_adjusted =c((1577.85-(51.16-43.9)), (1620.27-(41.42-35)),  (1620.27-(41.42-36)), (553.68-(101-96.5)))
RDA_LIPID_LOCI_labs$LOD = c(3.76, 2.95, 2.95,3.94)
#plot_fat <- plot_fat+ggrepel::geom_label_repel(data = RDA_LIPID_LOCI_labs, aes(label = locus, x= cont_vals_adjusted, y = LOD),segment.size = 0.5, nudge_y = 1, nudge_x = 300, segment.color="black", min.segment.length = 1, direction="y", size = 2)

#plot_fat <- plot_fat+ggrepel::geom_label_repel(data = RDA_LIPID_LOCI, aes(label = locus, x= cont_vals, y = LOD),segment.size = 0.5, nudge_y = -.5, nudge_x = 200, segment.color="black", min.segment.length = .2, direction="y", size = 2)
plot_fat <- plot_fat+geom_point(data = RDA_LIPID_LOCI, aes(x= cont_vals, y = 2.9), size = .5, shape = 21, color = RDA_LIPID_LOCI$color)


### plot fat zoomed

ggplot(out_map_cont %>% filter(LG ==38), aes(x=unlist.map., y = meanFat))+
  #geom_vline(data = RDA_LIPID_LOCI, aes(xintercept = unlist.map.), size = 3, alpha = .25, color ="black")+
  geom_line(alpha = 0.8)+
  #scale_color_manual(values = c("darkgray","steelblue"))+
  geom_hline(yintercept = 3, color = "purple")+
  geom_hline(yintercept = sig[2], color = "red", linetype = "dashed")+
  #facet_grid(~factor(LG), scales = "free")+
  scale_x_continuous(breaks = labs$LG_lab, labels = labs$LG_b)+
  theme_classic()+
  labs(x="", y = "LOD", title="")+
  theme(panel.spacing = unit(0, "lines"), 
        legend.position = "none",
        axis.text.x = element_text(size = 6))



## condition
plot_condition <- ggplot(out_map_cont, aes(x=cont_vals, y = meanK, color = LG))+
  geom_line(alpha = 0.8, color = out_map_cont$color)+
  geom_hline(yintercept = 3, color = "purple")+
  geom_hline(yintercept = sig[6], color = "red", linetype = "dashed")+
  #facet_grid(~factor(LG), scales = "free")+
  scale_color_manual(values =colors)+
  scale_x_continuous(breaks = labs$LG_lab, labels = labs$LG_b)+
  theme_classic()+
  labs(x="Linkage Group", y = "LOD", title = "")+
  theme(panel.spacing = unit(0, "lines"), 
        legend.position = "none",
        axis.text.x = element_text(size = 6))

#plot_condition <- plot_condition+ggrepel::geom_label_repel(data = RDA_COND_LOCI, aes(label = locus, x= cont_vals, y = LOD),segment.size = 0.5, nudge_y = -.5, nudge_x = 200, segment.color="black", min.segment.length = .2, direction="y", size = 2)
plot_condition <- plot_condition+geom_point(data = RDA_COND_LOCI, aes(x= cont_vals, y = 2.5), size = .5, shape = 21, color = RDA_COND_LOCI$color)


ggpubr::ggarrange(plot_fat, plot_condition, ncol =1, nrow=2, labels = c("A", "B"))

save(out_map_cont, file = "qtl/out_map_cont.RData")

