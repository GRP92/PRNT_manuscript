library("readxl") # read excel
library(GGally) # ggpairs() 
library(tidyr)
library(ggplot2)
#install.packages('devtools')
#devtools::install_github('VPetukhov/ggrastr')
library(ggbeeswarm) # geom_quasirandom()

# pip install mofapy2 (on terminal)
# devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
# https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/getting_started.md
# https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/downstream_analysis.html
# (Optional) set up reticulate connection with Python
# reticulate::use_python("/Users/ricard/anaconda3/envs/base_new/bin/python", required = T)
## Warning: replacing previous import 'DelayedArray::pnorm' by 'stats::pnorm'
## when loading 'MOFA2'
library(MOFA2)



mycols <- c('PRNT+' = 'coral2', 'PRNT-' = 'aquamarine4')
# mycols <- c('IGG+' = 'coral2', 'IGG-' = 'aquamarine4')



# ================================================== META MATRIX ==================================================

meta <- read.table("DataIN/20200902 metadati PRNT.txt", row.names = 1, sep = "\t", header = T, check.names=F, stringsAsFactors = F)
meta <- meta[,c(2, 3, 9)]


# ================================================== PROTEOMICS MATRIX ==================================================

PRT <- read.table("DataIN/Olink_KD.csv", sep = ",", header = T, check.names=F, stringsAsFactors = F)
PRT$ID <- gsub(" ","", PRT$ID)

PRTcov <- PRT[grepl("CACTUS", PRT$ID2, fixed = TRUE),]
rownames(PRTcov) <- PRTcov$ID
PRTcov <- PRTcov[,-c(1:2)]
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



# ================================================== MERGE ==================================================

CD40Lpos <- merge(CD40LposUS, CD40LposST, by = "row.names")
rownames(CD40Lpos) <- CD40Lpos[,1]
CD40Lpos <- CD40Lpos[,-1]
rm(CD40LposST, CD40LposUS)

FACS <- merge(CD40Lpos, BposUS, by = "row.names")
rownames(FACS) <- FACS[,1]
FACS <- FACS[,-1]

ALL <- merge(FACS, PRTcov, by = "row.names")
rownames(ALL) <- ALL[,1]
ALL <- ALL[,-1]
ALL <- merge(meta, ALL, by = "row.names")

Group <- ALL$PRNT
Group[Group > -1 ] <- 1
Group <- gsub("-1","PRNT-", Group)
Group <- gsub("1","PRNT+", Group)

rownames(ALL) <- ALL[,1]
ALL <- ALL[,-c(1:4)]


# =================== MATRIX 

CD40L <- ALL[, c(1:ncol(CD40Lpos))]
Bcells <- ALL[, c((ncol(CD40Lpos) + 1):(ncol(CD40Lpos) + ncol(BposUS)))]
PRT <- ALL[, c((ncol(CD40Lpos) + ncol(BposUS) + 1):ncol(ALL))]

CD40L <- as.matrix(CD40L)
Bcells <- as.matrix(Bcells)
PRT <- as.matrix(PRT)

class(CD40L) <- "numeric"
class(Bcells) <- "numeric"
class(PRT) <- "numeric"

rm(CD40Lpos, BposUS, ALL, FACS, PRTcov)


# ================================================== NORMALIZATION ==================================================

CD40L <- log2(CD40L + 1)
Bcells <- log2(Bcells + 1)
# PRT is already in logarithmic scale


# ================================================== FILTERING ==================================================

if(F){
  
  limIQR <- 0.3
  
  # Filter for IQR (Q3-Q1)
  iqr <- apply(CD8, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xlab="IQR", ylab="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) CD8 <- CD8[, -ir]
  
  iqr <- apply(CD40L, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xlab="IQR", ylab="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) CD40L <- CD40L[, -ir]
  
  iqr <- apply(Bcells, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xlab="IQR", ylab="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) Bcells <- Bcells[, -ir]
  
  iqr <- apply(NK, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xlab="IQR", ylab="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) NK <- NK[, -ir]
  
  iqr <- apply(PRT, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xlab="IQR", ylab="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) PRT <- PRT[, -ir]
  
}



# ================================================== DATA LOAD ==================================================

# Multiple formats are allowed for the input data:

## -- Option 1 -- ##
# nested list of matrices, where the first index refers to the view and the second index refers to the group.
# samples are stored in the rows and features are stored in the columns.
# Missing values must be filled with NAs, including samples missing an entire view

# (...)

## -- Option 2 -- ##
# data.frame with columns ["sample","feature","view","group","value"]
# In this case there is no need to have missing values in the data.frame,
# they will be automatically filled in when creating the corresponding matrices


#ldata <- list(CD40L = t(CD40L), NK = t(NK), CD8 = t(CD8), Bcells = t(Bcells), Proteomics = t(PRT))
# ldata <- list('FACS UM' = t(cbind(CD40L, Bcells)), 'FACS CL' = t(cbind(NK, CD8)), Proteomics = t(PRT))
# ldata <- list(FACS = t(cbind(CD40L, Bcells, NK, CD8)), Proteomics = t(PRT))

ldata <- list(Tcells = t(CD40L), Bcells = t(Bcells), Proteomics = t(PRT))



# ================================================== CREATE MOFA OBJECT ==================================================

# The aim of the multi-group framework is to identify the sources of variability *within* the groups. If your aim is to find a factor that 'separates' the groups, 
# you DO NOT want to use the multi-group framework. Please see the FAQ (https://github.com/bioFAM/MOFA2#2-faq-on-the-multi-group-functionality)

# MOFAobject <- create_mofa(ldata, groups=Group)
MOFAobject <- create_mofa(ldata)

# Visualise data structure
plot_data_overview(MOFAobject)


# ================================================== DEFINE OPTIONS ==================================================

# Data options
# - scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance (default is FALSE)
# - scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance (default is FALSE)
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- T

# Model options
# - likelihoods: likelihood per view (options are "gaussian","poisson","bernoulli"). By default, they are guessed internally. We advise users to use “gaussian” whenever possible!
# - num_factors: number of factors. By default K=10
# - spikeslab_factors: use spike-slab sparsity prior in the factors? default is FALSE.
# - spikeslab_weights: use spike-slab sparsity prior in the weights? default is TRUE.
# - ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
# - ard_weights: use ARD prior in the weights? Default is TRUE if using multiple views.
# Only change the default model options if you are familiar with the underlying mathematical model!
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 7
# model_opts$likelihoods <- c("gaussian","gaussian","gaussian","gaussian","gaussian")
# model_opts$likelihoods <- c("gaussian","gaussian","gaussian")

# Training options
# - maxiter: number of iterations. Default is 1000.
# - convergence_mode: "fast", "medium", "slow". For exploration, the fast mode is good enough.
# - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# - freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
# - drop_factor_threshold: minimum variance explained criteria to drop factors while training. Default is -1 (no dropping of factors)
# - gpu_mode: use GPU mode? This needs cupy installed and a functional GPU, see https://cupy.chainer.org/
# - verbose: verbose mode?
# - seed: random seed
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$maxiter <- 4000

# (Optional) Set stochastic inference options
# Only recommended with very large sample size (>1e6) and when having access to GPUs
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.1 to 0.5)
# - start_stochastic: first iteration to apply stochastic inference (recommended > 5)
# stochastic_opts <- get_default_stochastic_options(MOFAobject)


# ================================================== PREPARE MOFA OBJECT ==================================================

MOFAobject <- prepare_mofa(MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts) # stochastic_options = stochastic_opts


# ================================================== TRAIN THE MODEL ==================================================

# R requires the package reticulate to communicate with Python, and this is the source of most of the installation problems. 
# https://rstudio.github.io/reticulate/


library(reticulate)
# ho creato il mio ENV chiamato myR
# source /anaconda3/etc/profile.d/conda.sh
# conda activate myR
# pip install mofapy2

use_python("/anaconda3/envs/myR/bin/python3.6", required = T)
conda_list()
py_config()

# Using a conda enviroment called myR
use_condaenv("myR", required = TRUE)

outfile <- paste0(getwd(),"/MOFA/model.hdf5")
MOFAmodel <- run_mofa(MOFAobject, outfile)



# ================================================== LOAD TRAINED MODEL AND ADD METADATA ==================================================

model <- load_model(paste0(getwd(),"/MOFA/model.hdf5"))

sample_metadata <- model@samples_metadata
# head(sample_metadata, n=3)
sample_metadata$condition <- Group
# sample_metadata$age <- Age
# sample_metadata$sex <- Sex
# sample_metadata

samples_metadata(model) <- sample_metadata
# head(model@samples_metadata, n=3)


# ================================================== VARIANCE DECOMPOSITION ==================================================

# Total variance explained per view and group
# head(model@cache$variance_explained$r2_total[[1]]) # group 1

# Variance explained for every factor in per view and group
# head(model@cache$variance_explained$r2_per_factor[[1]]) # group 1

plot_variance_explained(model, x="view", y="factor") +
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust= 1, size = 15), axis.text.y = element_text(size = 15), legend.title = element_text(size = 13))
ggsave(("MOFA/Variance Decomposition.pdf"), width = 4, height = 7)


plot_variance_explained(model, x="view", y="factor", plot_total = T)[[2]] +
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust= 1, size = 19), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 15))
ggsave(("MOFA/Variance Decomposition Tot.pdf"), width = 3, height = 7)



# ================================================== VISUALIZATION OF SINGLE FACTORS ==================================================


# Get factor values
Z <- get_factors(model, factors = 1:7 , groups = "all", as.data.frame=TRUE)
Z$factor <- as.factor(Z$factor)
df <- merge(Z,  cbind(sample = rownames(Bcells), condition = Group), by="sample")

df$condition <- factor(df$condition, levels = c("PRNT+", "PRNT-"))


# Beeswarm Plots v1
ggplot(df, aes(x=factor, y=value, col=condition)) + facet_wrap(~factor, nrow=1, scales="free_x") +
  labs(x="", y="Factor value") + geom_quasirandom(size = 3, dodge.width = 0, alpha = 0.8, method = "quasirandom", width = 0.4) + 
  scale_color_manual(values=mycols) + theme_classic() + geom_hline(yintercept=0, linetype="dashed", size=0.4, alpha=0.5) +
  scale_shape_manual(values=c(20, 18)) +
  theme(text = element_text(size=18), panel.border = element_rect(color="black", size=0.2, fill=NA), 
        strip.background = element_rect(colour = "black", size=0.3), panel.spacing = unit(0,"lines"),
        axis.line = element_line(colour = 'grey', size = 0.2), axis.text.x = element_text(color="white"),
        axis.ticks.x = element_line(color = "white"))
ggsave(("MOFA/Beeswarm plots.pdf"), width = 14, height = 6)



if (F){
  # Violin Plots
  plot_factor(model, factor = 1:5, color_by = "condition", dot_size = 1.5, dodge = T, legend = T, add_violin = T, violin_alpha = 0.3) + 
    scale_fill_manual(values=mycols) + theme(text = element_text(size=18)) 
  ggsave(("MOFA/Violin plots.pdf"), width = 14, height = 6)
}


# Violin Plots
ggplot(df[df$factor %in% c("Factor1", "Factor2", "Factor3", "Factor4", "Factor5"),], aes(x=factor, y=value, fill=condition, col=condition)) + facet_wrap(~factor, nrow=1, scales="free_x") +
  geom_quasirandom(size = 1.8, dodge.width = 3, alpha = 1, method = "quasirandom", width = 0.4, show.legend = F) + 
  geom_violin(alpha= 0.2, trim = T, scale = "width", position = position_dodge(width=3), show.legend = T, col = "black", width = 2.8) +
  scale_shape_manual(values=c(20, 18)) + scale_color_manual(values=mycols) + scale_fill_manual(values=mycols) + 
  # geom_hline(yintercept=0, linetype="dashed", size=0.4, alpha=0.5) +
  labs(x="", y="Factor value") + theme_classic() + 
  theme(text = element_text(size=18), panel.border = element_rect(color="black", size=0.2, fill=NA), 
        strip.background = element_rect(colour = "black", size=0.3), panel.spacing = unit(0,"lines"),
        axis.line = element_line(colour = 'grey', size = 0.2), axis.text.x = element_text(color="white"),
        axis.ticks.x = element_line(color = "white"))
ggsave(("MOFA/Violin plots.pdf"), width = 14, height = 6)




# ================================================== VISUALIZATION OF COMBINATION OF FACTORS ==================================================

# plot_factors(model, factor = 1:5, color_by = "condition", shape_by = "sex") + 
#  scale_fill_manual(values=mycols) + scale_color_manual(values=mycols) + theme(text = element_text(size=18))

# spread over factors
df2 <- spread(df, key="factor", value="value")

ggpairs(df2, columns = c(4:8),
        lower = list(continuous=GGally::wrap("points", size=2, alpha = 0.6)), 
        diag = list(continuous=GGally::wrap("densityDiag", alpha = 0.5, size = 0.2)), 
        upper = list(continuous=GGally::wrap("points", size=2, alpha = 0.6)),
        mapping = ggplot2::aes(color=condition)) + 
  scale_fill_manual(values=mycols) + scale_color_manual(values=mycols) + 
  theme(text = element_text(size=16), panel.grid.major = element_blank(),
        panel.border = element_rect(color="black", size=0.3, fill=NA),
        panel.background = element_rect(fill = "white"), strip.background = element_rect(color = "black", size=0.3, fill = "gray94"))
ggsave(("MOFA/Combinations of factors.pdf"), width = 7, height = 7)


# axis.text = element_blank(), , axis.ticks = element_blank()

# ================================================== VISUALIZATION OF FEATURE WEIGHTS ==================================================

#plot_weights(model, view = "FACS", factor = 1, nfeatures = 7,  scale = T, text_size = 4)

var <- model@cache$variance_explained$r2_per_factor[[1]]
var <- cbind(Factor = 1:nrow(var), var)

for(z in 1:length(ldata)){
  
  plot_top_weights(object = model, view = names(ldata)[z], factor = var[which.max(var[,z + 1]), 1], nfeatures = 10) + 
    theme(text = element_text(size=15), strip.background = element_rect(color = "black", fill = "gray94"), plot.margin = unit(c(0.5, 0.5, 0.5, wtd)))
  ggsave(paste0("MOFA/Feature weights F", var[which.max(var[,z + 1]), 1]," (", names(ldata)[z], ").pdf"), width = 10, height = 4)
}


z <- 2
plot_top_weights(object = model, view = names(ldata)[z], factor = 3, nfeatures = 10) + 
  theme(text = element_text(size=15), strip.background = element_rect(color = "black", fill = "gray94"))
ggsave(paste0("MOFA/Feature weights F", 3," (", names(ldata)[z], ").pdf"), width = 10, height = 4)

# ================================================== VISUALIZATION OF PATTERNS IN THE DATA ==================================================

# df3 <- merge(df2, cbind(sample = rownames(freq), freq), by = "sample")
# m <- df3[, c(1:3,10,17)]
# names(m)[5] <- "CD3_Freq_of_Parent"
# ggplot(m, aes(Factor5, CD3_Freq_of_Parent, color= condition)) + geom_point(alpha = 0.4, size = 3) + geom_smooth(formula=y~x, method=lm) + scale_color_manual(values=mycols)

library(ggpubr)

plot_data_scatter(model, view = "NK", factor = 2,
                  features = 6,           # number of features to plot (they are selected by loading)
                  add_lm = TRUE,          # add linear regression
                  color_by = "condition", lm_per_group = F) + scale_fill_manual(values=mycols) 
ggsave(("MOFA/Scatter plots NK_values vs Factor1_values.pdf"), width = 10, height = 6)


plot_data_scatter(model, view = "Proteomics", factor = 1,
                  features = 6,           # number of features to plot (they are selected by loading)
                  add_lm = TRUE,          # add linear regression
                  color_by = "condition", lm_per_group = F) + scale_fill_manual(values=mycols) 
ggsave(("MOFA/Scatter plots Proteomics_values vs Factor1_values.pdf"), width = 10, height = 6)




library(uwot)

set.seed(42)
model <- run_umap(model)




model <- load_model(paste0(getwd(),"/MOFA/model.hdf5"))

sample_metadata <- model@samples_metadata
head(sample_metadata, n=3)
sample_metadata$condition <- Group
sample_metadata

samples_metadata(model) <- sample_metadata
head(model@samples_metadata, n=3)

model <- run_tsne(model, dims = 2, perplexity = 7, verbose=F, max_iter = 5000, theta = 0.0, pca_scale = F, normalize = T)

Z <- model@dim_red[["TSNE"]]

ggplot(Z, aes(TSNE1, TSNE2, col=Group)) +  geom_point(size = 6, alpha = 0.8)  + scale_color_manual(values=mycols) + 
  theme_bw(base_size = 17) + theme(legend.title = element_blank(), panel.grid.major = element_blank()) 
ggsave(("MOFA/TSE.pdf"), width = 9, height = 8)


