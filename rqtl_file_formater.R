### build QTL dataset
library(tidyverse)
alns <- read.table("./data/align_dat_qtl.txt")


library(readxl)
linkMap <- read_excel("./data/Supplementary Material - Table S1 - Linkage Map Information.xlsx", 
                                                                      skip = 1)
head(alns)
head(linkMap)

## filter alingment to only LGs and good alignment
alns2 <- alns %>% filter(V2 %in% unique(grep("LG", alns$V2, value=T)), V4==60) %>% group_by(V1) %>% mutate(N=n()) %>% 
filter(N==1) %>% select(-N)



# join to existing map data
mappedAlns <- left_join(alns2, linkMap[,1:6], by = c("V1"="RAD_LOCUS_NAME"))
colnames(mappedAlns) <- c("LocusID", "MappedChromosome", "PosInBP", "qual", colnames(mappedAlns[5:9]))
mappedAlns$MappedChromosome <- as.numeric(gsub("LG","",mappedAlns$MappedChromosome))

wrong_chrom <- mappedAlns %>% filter(LocusID %in% grep("qtl", mappedAlns$LocusID, value=T), MappedChromosome != LG)

## remove all loci that were mapped to a different LG in the LT genome than in Seth's paper
mappedAlns <- mappedAlns %>% filter(!LocusID %in% wrong_chrom$LocusID)
   
   
# view relationahip between position in BP and centi-morgans
ggplot(mappedAlns, aes(x = PosInBP , y = FEMALE_POS, color = MappedChromosome ))+
  geom_point()+
  facet_wrap(~MappedChromosome)+
  labs(x = "Position (BP) draft lake trout genome", y = "genetic pos (cM)")+
  theme_bw()

## use a loess relationship to predict positions in cm of markers aligned to the same genome that were not included in original map.

predict_cm <- function(x){
  m_total <- NULL
  for(lg in unique(x$MappedChromosome)){
    d <- x[x$MappedChromosome == lg,]
    d_map <- d %>% drop_na()
    m <- predict(loess(loess(d_map$FEMALE_POS~d_map$PosInBP)), d$PosInBP)
    m_dat <- data.frame(m, lg)
    m_total <- bind_rows(m_total, m_dat)
    
    
  }
  m_total
}


pred_cm <- predict_cm(mappedAlns)

# join predicted cm to original dataset
mappedAlns$pred_cm <- pred_cm$m


# plot to compare how well predicted values (black) relate to known values (gray)
ggplot(mappedAlns, aes(x = PosInBP , y = pred_cm))+
  geom_point(data=mappedAlns, aes(x =PosInBP, y = FEMALE_POS), color = "gray", alpha=0.5)+
  #geom_point(, alpha = 0.05)+
  geom_smooth(method = "loess")+
  facet_wrap(~MappedChromosome)+
  labs(x = "Position (BP) draft lake trout genome", y = "genetic pos (cM)")+
  theme_bw()


# Drop loci that have no distance in cm - these are loci that are were not flanked by qtl loci and so we would have to extrapolate

mappedAlns <- mappedAlns[!is.na(mappedAlns$pred_cm), ]



## create QTL map file

map <- mappedAlns %>% filter(!LocusID%in% grep("qtlCLocus", mappedAlns$LocusID, value=T)) %>% 
  select(LocusID, MappedChromosome, pred_cm) %>% mutate(pred_cm = round(pred_cm, 2))



#write.csv(map, "./data/map.csv", row.names = F)


#--------------------------------------------------
### Create genotypes file
#--------------------------------------------------

dat95 <- read.table("./data/LT_90_parsed_genotype_calls.txt", header=F)

#dat95$V2 <- paste(dat95$V1, dat95$V2, sep = "_")
dat95_t <- t(dat95)
colnames(dat95_t) <- c("SampleID", dat95_t[3,2:ncol(dat95_t)])
dat95_t <- as.data.frame(dat95_t[8:nrow(dat95_t),])

# change column names to match locus names
colnames(dat95_t) <- c("SampleID", paste0("CLocus_", gsub("_.*","",colnames(dat95_t)[2:length(colnames(dat95_t))])))

dat95_t <- select(dat95_t, c("SampleID", map$LocusID))


#--------------------------------------------------
### Create phenotypes file
#--------------------------------------------------

# customize bio dat to get single mean values for wt, length, fat, and rankFat
bio_data <- read.csv("./data/master_BioData_ranks_poSAL.csv")

bio_data <- bio_data %>% filter(SampleID %in% dat95_t$SampleID) %>% 
  mutate(K = 100000*weight_g/length_mm^3) %>% group_by(SampleID) %>% 
  summarise(cross = first(cross), 
            family =paste(first(dam_col), first(sire_col), sep = "_"), 
            sire=first(sire_col),
            dam=first(dam_col), 
            meanRank = first(meanRank), 
            meanFat=first(meanFat), 
            meanWT = mean(weight_g, na.rm=T),
            meanLen=mean(length_mm, na.rm = T),
            meanK=mean(K, na.rm=T), 
            sex = first(sex[which(sex !="?")])) %>% arrange(cross, family) %>%
  filter(!is.na(sex)) %>% 
  select(SampleID, sex, meanFat, meanLen, meanK, family, cross)


## Remove all individuals with family sizes < 
bio_data %>% group_by(cross, family) %>% summarise(N =n()) %>% arrange(-N) 

small_fams <- bio_data %>% group_by(cross, family) %>% summarise(N =n()) %>% filter(N<6) 

bio_data <- bio_data %>% filter(!family%in% small_fams$family)

phenotypes <- bio_data %>% select(SampleID, meanFat, meanLen, meanK)
covars <- bio_data %>% select(SampleID, sex, cross, family)

write.csv(phenotypes, "./data/phenotypes.csv", row.names = F)
write.csv(covars, "./data/covariates.csv", row.names = F)


## arrange genotypes file based on phenotypes file

pheno_and_geno <- left_join(bio_data, dat95_t, by = "SampleID") %>% 
  select(-sex, -meanFat, -meanLen, -meanK, -family, -cross)

pheno_and_geno[pheno_and_geno==0] <- "A"
pheno_and_geno[pheno_and_geno==1] <- "H"
pheno_and_geno[pheno_and_geno==2] <- "B"
pheno_and_geno[pheno_and_geno==9] <- "-"

#write.csv(pheno_and_geno, "./data/genotypes.csv", row.names = F)

