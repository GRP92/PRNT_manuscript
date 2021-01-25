library(ggplot2)
library(corrplot)
library(RColorBrewer) 
library("readxl") # read excel


# ================================================== CorHP FUNCTION ==================================================

CorHP <- function(m, lab3, limR, Dir){
  
  # remove col with n° of values < 5
  sm <- colSums(!is.na(m))
  m <- m[,sm > 4]
  
  # remove col with n° of zeros = n° rows
  sm <- colSums(m == 0, na.rm = T)
  ir <- which(sm == nrow(m))
  if(length(ir) > 0) m <- m[,-ir]
  
  # save data matrix
  write.table(cbind(ID = rownames(m), m), paste0(Dir, "/Correlation Matrix - ", lab3, ".txt"), sep="\t", row.names = F, col.names = T,  quote = F)
  
  # correlation and pv
  c <- cor(m, method = "spearman", use="pairwise")
  pv <- apply(m, 2, function(y){apply(m, 2, function(x){cor.test(x, y, method="spearman", na.action = "na.omit")$p.value})})
  # c <- rcorr(m, type = "spearman")[[1]]
  # pv <- rcorr(m, type = "spearman")[[3]]
  
  # filtering
  pv[is.na(pv)] <- 1
  c[which(pv > 0.05)] <- 0
  c[is.na(c)] <- 0
  c[abs(c) < limR] <- 0
  s <- round(colSums(c), 4)
  c <- c[s != 1, s != 1]
  write.table(cbind(ID = rownames(c), c), paste0(Dir, "/Correlation Matrix (R values) - ", lab3, ".txt"), sep="\t", row.names = F, col.names = T,  quote = F)
  
  # plot
  sz <- round(5 + (ncol(m)/5.5), 1)
  pdf(paste0(Dir, "/CorHP - ", lab3, ".pdf"), width = sz, height = sz)
  corrplot(c, cl.cex = 1+(sz/25), type="upper", order="original", tl.col = "black", col = rev(brewer.pal(n=6, name="RdYlBu")), title = "")
  dev.off()
}



# ================================================== RegLin FUNCTION ==================================================

RegLin <- function(plt, labA, labB, Sym, reg, colBool, colr, r){
  
  p <- ggplot(data=plt, aes(x=B, y=A)) + ylab(labA) + xlab(labB) + theme_linedraw(base_size = 19) 
  if(colBool) p <- p + geom_point(size = 5, alpha = 0.7, aes(colour = Sym)) + scale_color_manual(values=mycols2)
  if(!colBool) p <- p + geom_point(size = 5, alpha = 0.7, col = colr)
  if(reg) p <- p + geom_smooth(method="lm", se=T, level = 0.95, color = "black", alpha = 0.2)
  p <- p + annotate(geom = 'text', label = r, x = Inf, y = Inf, hjust = 1.3, vjust = 1.5, parse=TRUE, size = 6)
  p <- p + theme(axis.text.y = element_text(colour = "black", family = "sans", size = 16),
                 axis.text.x = element_text(colour = "black", family = "sans", size = 16),
                 axis.title.y = element_text(vjust= 4, family = "sans"),
                 legend.title = element_text(colour = "white"),
                 plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
                 panel.grid = element_line(color = "gray87")) 
  p <- p + scale_y_continuous(labels = function(x) format(x, scientific = FALSE))
  
  p
  
}




# ================================================== META MATRIX ==================================================

meta <- read.table("DataIN/20200902 metadati PRNT.txt", row.names = 1, sep = "\t", header = T, check.names=F, stringsAsFactors = F)
meta <- meta[,-c(1,4,5,12,14)]
meta$Sex <- gsub("M","0", meta$Sex)
meta$Sex <- gsub("F","1", meta$Sex)
meta <- as.matrix(meta)
class(meta) <- "numeric"
meta <- as.data.frame(meta)



# ================================================== PROTEOMICS MATRIX ==================================================

PRT <- read.table("DataIN/Olink_KD.csv", sep = ",", header = T, check.names=F, stringsAsFactors = F)
PRT$ID <- gsub(" ","", PRT$ID)
PRTcov <- PRT[grepl("CACTUS", PRT$ID2, fixed = TRUE),]
rownames(PRTcov) <- PRTcov$ID
PRTcov <- PRTcov[,-c(1:2)]
PRTcov <- as.matrix(PRTcov)
class(PRTcov) <- "numeric"
rm(PRT)


# ================================================== FACS MATRIX ==================================================

CD40LposST <- read_excel(paste0("DataIN/MATRIX/CD40L - CovPosST.xlsx"))
CD40LposUS <- read_excel(paste0("DataIN/MATRIX/CD40L - CovPosUS.xlsx"))
BposUS <- read_excel(paste0("DataIN/MATRIX/Bcells - CovPosUS.xlsx"))

CD40LposST <- t(CD40LposST)
colnames(CD40LposST) <- CD40LposST[1,]
CD40LposST <- CD40LposST[-1,]

CD40LposUS <- t(CD40LposUS)
colnames(CD40LposUS) <- CD40LposUS[1,]
CD40LposUS <- CD40LposUS[-1,]

BposUS <- t(BposUS)
colnames(BposUS) <- BposUS[1,]
BposUS <- BposUS[-1,]

FACS <- merge(CD40LposUS, CD40LposST, by = "row.names")
rownames(FACS) <- FACS[,1]
FACS <- FACS[,-1]

FACS <- merge(FACS, BposUS, by = "row.names")
rownames(FACS) <- FACS[,1]
FACS <- FACS[,-1]

rm(CD40LposUS)




# ================================================== CORRELATIONS HEATMAP ==================================================

m <- merge(meta, FACS, by = "row.names")
rownames(m) <- m$Row.names
m <- as.matrix(m[,-1])

m <- merge(m, PRTcov, by = "row.names")
rownames(m) <- m$Row.names
m <- as.matrix(m[,-1])
class(m) <- "numeric"
CorHP(m, lab3 = "META vs FACS vs PROTEOMICS", limR = 0, Dir = "Plots")




# ================================================== CORRELATIONS PLOTS ==================================================

ComparisonList <- list()

# FIGURE 1 ddPCR Correlations
mx1 <- read.table("DataIN/Dataset Time-Matched 1.txt", sep = "\t", header = T, row.names = 1, check.names=F, stringsAsFactors = F)
ComparisonList[[1]] <- as.data.frame(mx1[, c(1, 2)]) # Fig 1b
ComparisonList[[2]] <- as.data.frame(mx1[, c(1, 3)]) # Fig 1a
ComparisonList[[3]] <- as.data.frame(mx1[, c(4, 2)]) # Fig 1d
ComparisonList[[4]] <- as.data.frame(mx1[, c(1, 4)])
ComparisonList[[5]] <- as.data.frame(mx1[, c(4, 3)]) # Fig 1c


# Supplementary Fig. 1 IgG - Neutralization Correlation
mx1 <- read.table("DataIN/Dataset Time-Matched 2.txt", sep = "\t", header = T, row.names = 1, check.names=F, stringsAsFactors = F)
ComparisonList[[6]] <- as.data.frame(mx1[, c(2, 1)]) # Fig S1


# B-cells Correlations
mx1 <- merge(meta[, c(6, 7)], BposUS[, 4], by = "row.names")
names(mx1)[4] <- "SARS-CoV-2-specific B cells"
rownames(mx1) <- mx1[,1]
mx1 <- as.matrix(mx1[,-1])
class(mx1) <- "numeric"
ComparisonList[[7]] <- as.data.frame(mx1[, c(2, 3)])
ComparisonList[[8]] <- as.data.frame(mx1[, c(1, 3)])


# CD40L Correlations
mx1 <- merge(meta[, c(6, 7)], CD40LposST[, 1], by = "row.names")
names(mx1)[4] <- "SARS-CoV-2-specific CD4 T-cells"
rownames(mx1) <- mx1[,1]
mx1 <- as.matrix(mx1[,-1])
class(mx1) <- "numeric"
ComparisonList[[9]] <- as.data.frame(mx1[, c(2, 3)])
ComparisonList[[10]] <- as.data.frame(mx1[, c(1, 3)])


# FIGURE 4 SLAMF1 Correlations
mx1 <- merge(meta[, c(6, 7)], PRTcov[, colnames(PRTcov) %in% "SLAMF1"], by = "row.names")
names(mx1)[4] <- "SLAMF1"
rownames(mx1) <- mx1[,1]
mx1 <- as.matrix(mx1[,-1])
class(mx1) <- "numeric"
ComparisonList[[11]] <- as.data.frame(mx1[, c(2, 3)]) # Fig 4f
ComparisonList[[12]] <- as.data.frame(mx1[, c(1, 3)]) # Fig 4e


# SLAMF1 - Bcells Correlation
mx1 <- merge(BposUS[, 4], PRTcov[, colnames(PRTcov) %in% "SLAMF1"], by = "row.names")
names(mx1)[2] <- "SARS-CoV-2-specific B cells"
names(mx1)[3] <- "SLAMF1"
rownames(mx1) <- mx1[,1]
mx1 <- as.matrix(mx1[,-1])
class(mx1) <- "numeric"
ComparisonList[[13]] <- as.data.frame(mx1[, c(1, 2)])


# FIGURE 1 AUC - IGG Correlation
#mx1 <- read.table("DataIN/Dataset Time-Matched 3.txt", sep = "\t", header = T, row.names = 1, check.names=F, stringsAsFactors = F)
ComparisonList[[14]] <- as.data.frame(meta[, c(8, 7)]) # Fig 1e



for(j in 1:length(ComparisonList)){
  
  plt <- ComparisonList[[j]]
  labA = names(plt)[1]
  labB = names(plt)[2]
  colnames(plt) <- c("A", "B")
  plt <- plt[complete.cases(plt),]
  
  # Correlation
  cs <- cor(plt$A, plt$B, method = "spearman")
  pvs <- cor.test(plt$A, plt$B, method="spearman")$p.value
  
  if(pvs < 0.0001) {
    ann = paste0("atop('rho' == ", round(cs, 3), ", p < 0.0001)")
  } else {
    ann = paste0("atop('rho' == ", round(cs, 3), ", p == ", round(pvs, 5), ")")
  }
  
  mxy <- max(abs(plt$A))
  wtd <- round(1.8 - nchar(mxy)/7, 1)
  
  p <- ggplot(data=plt, aes(x=B, y=A)) + ylab(labA) + xlab(labB) + theme_linedraw(base_size = 19) + 
    geom_point(size = 4, alpha = 0.6, col = "black") + geom_smooth(method="lm", se=T, level = 0.95, color = "black", alpha = 0.2)
  
  if(j == 6 | j == 11 | j == 12) {
    p <- p + annotate(geom = 'text', label = ann, x = Inf, y = Inf, hjust = 3.5, vjust = 1.5, parse=TRUE, size = 6)
  } else {
    p <- p + annotate(geom = 'text', label = ann, x = Inf, y = Inf, hjust = 1.3, vjust = 1.5, parse=TRUE, size = 6)
  }
  
  p <- p + theme(axis.text.y = element_text(colour = "black", family = "sans", size = 14),
                 axis.text.x = element_text(colour = "black", family = "sans", size = 14),
                 axis.title.y = element_text(vjust= 2, family = "sans"),
                 legend.title = element_text(colour = "white"),
                 plot.margin = unit(c(0.5, 0.5, 0.5, wtd), "cm"),
                 panel.grid = element_line(color = "gray87")) 
  
  
  if(j == 1 | j == 2  | j == 4 | j == 14) p <- p + scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", format(x, scientific = T)))))})
  
  p
  ggsave(paste0("Plots/Linear Regression - ", gsub("/","|", labA), " vs ", gsub("/","|", labB), ".pdf"), width = 6, height = 5)
  
}

rm(plt, ann, cs, j, labA, lavB, pvs, p)





