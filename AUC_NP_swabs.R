library(tidyr) # spread()
library(MESS) # AUC

# METADATA
meta <- read.table("DataIN/20200902 metadati PRNT.txt", row.names = 1, sep = "\t", header = T, check.names=F, stringsAsFactors = F)
meta <- meta[,c(4,5,6,9)]
meta$`Date of start symptoms` <- as.Date(meta$`Date of start symptoms`, format="%d/%m/%Y")
meta$`Date of hospital admission` <- as.Date(meta$`Date of hospital admission`, format="%d/%m/%Y")


# CT DATASET
CT <- read.table("DataIN/CT - NP.txt", sep = "\t", header = T, check.names=F, stringsAsFactors = F)
CT$DateT <- as.Date(CT$DateT, format="%d/%m/%Y")


# ET (40 - CT) calculation
CT$E <- 40 - CT$E
CT$RdRp <- 40 - CT$RdRp
CT$N <- 40 - CT$N


# Calculation of the average CT values for missing dates
ID <- unique(CT$ID)
CTmeanE <- data.frame()
CTmeanR <- data.frame()
CTmeanN <- data.frame()

for(i in 1:length(ID)){
  m <- CT[CT$ID %in% ID[i],]
  m <- m[order(m$DateT), ]
  mE <- m[,c(1,2,5)]
  mR <- m[,c(1,3,5)]
  mN <- m[,c(1,3,5)]
  mE <- mE[complete.cases(mE),]
  mR <- mR[complete.cases(mR),]
  mN <- mN[complete.cases(mN),]
  
  # GENE E
  if(nrow(mE) > 1){
    for(j in 1:(nrow(mE)-1)){
      diff <- as.numeric(mE[j + 1, 3] - mE[j, 3])
      if(diff > 1) CTmeanE <- rbind(CTmeanE, cbind(ID = ID[i], E = mean(c(mE[j + 1, 2], mE[j, 2])), DateT = as.character(seq(from = (mE[j, 3] + 1), to = (mE[j + 1, 3] - 1), by = 1))))
    }
  }
  
  # GENE R
  if(nrow(mR) > 1){
    for(j in 1:(nrow(mR)-1)){
      diff <- as.numeric(mR[j + 1, 3] - mR[j, 3])
      if(diff > 1) CTmeanR <- rbind(CTmeanR, cbind(ID = ID[i], RdRp = mean(c(mR[j + 1, 2], mR[j, 2])), DateT = as.character(seq(from = (mR[j, 3] + 1), to = (mR[j + 1, 3] - 1), by = 1))))
    }
  }
  
  # GENE N
  if(nrow(mN) > 1){
    for(j in 1:(nrow(mN)-1)){
      diff <- as.numeric(mN[j + 1, 3] - mN[j, 3])
      if(diff > 1) CTmeanN <- rbind(CTmeanN, cbind(ID = ID[i], N = mean(c(mN[j + 1, 2], mN[j, 2])), DateT = as.character(seq(from = (mN[j, 3] + 1), to = (mN[j + 1, 3] - 1), by = 1))))
    }
  }
}

CTmean <- merge(CTmeanE, CTmeanR, by = c("ID", "DateT"), all = T)
CTmean <- merge(CTmean, CTmeanN, by = c("ID", "DateT"), all = T)
CTmean <- cbind(ID = CTmean$ID, CTmean[,c(3:5)], DateT = CTmean$DateT)
rm(mE, mR, mN, diff, i, ID, j, m, CTmeanE, CTmeanN, CTmeanR)


# Clearance calculation (for symptomatic patients from symptoms onset date, from asymptomatic ones from the hospital admission date)
CT3 <- rbind(CT, CTmean)
CT3 <- merge(CT3, meta, by.x = "ID", by.y = "row.names")
CT3 <- cbind(CT3, Clearance = CT3$DateT - CT3$`Date of start symptoms`)
CT3[CT3$Symptoms == 0, colnames(CT3) %in% "Clearance"] <- CT3[CT3$Symptoms == 0, colnames(CT3) %in% "DateT"] - CT3[CT3$Symptoms == 0, colnames(CT3) %in% "Date of hospital admission"]


# Loop for each gene
labGene <- c("E", "RdRp", "N")
for(i in 1:length(labGene)){
  
  # count the time points for each patient 
  CT4 <- CT[, c(1, 1+i, 5)]
  CT4 <- CT4[complete.cases(CT4),]
  id <- table(CT4$ID)
  
  # select patients with at least 2 time points
  id <- id[id > 1]
  CT5 <- CT3[CT3$ID %in% names(id), ]
  CT5 <- CT5[,c(1, 1+i, 10)]
  CT5 <- CT5[complete.cases(CT5),]
  CT5[,2] <- as.numeric(as.character(CT5[,2]))
  names(CT5)[2] <- "CT"
  
  # AUC calculation
  Auc <- data.frame()
  ID <- unique(CT5$ID)
  for(w in 1:length(ID)){
    m <- CT5[CT5$ID %in% ID[w],]
    Auc <- rbind(Auc, cbind(ID = ID[w], AUC = auc(m$Clearance, m$CT)))
  }
  assign(paste("Auc_", labGene[i], sep = ""), Auc)
  
}


# Merge AUC of RdRp, E and N genes
Auc <- merge(Auc_E, Auc_RdRp, by = "ID", all = T)
colnames(Auc) <- c("ID", "AUC of E Gene", "AUC of RdRp Gene")
Auc <- merge(Auc, Auc_N, by = "ID", all = T)
names(Auc)[4] <- "AUC of N Gene"
write.table(Auc, "DataIN/AUC values.txt", sep="\t", row.names = F, col.names = T,  quote = F)


rm(m, i, ID, CT3, CT4, CT5, CTmean, id, labGene)




