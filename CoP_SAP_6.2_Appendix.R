## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(kableExtra)
library(ggplot2)
library(RColorBrewer)
library(boot)
library(dplyr)
library(readxl)
library(GGally)
library(cowplot)
library(stringr)
library(ggpubr)


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------
# At the request of a reviewer, we have created a 6.2 Appendix that redoes the
# correlation Tables 6.2.1 and 6.2.2 except on only the males (6.2.1m, 6.2.2m)
# or the females (6.2.1f, 6.2.2f). We keep most of the code the same as in 
# CoP_SAP_6.2_Final.Rmd and make as few changes as needed.
#
# changes are marked in sections labeled "Appendix Changes"


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
#Functions 

#Function to get bootstrap measures of Spearman's correlation between 2 measures
boot.spearman_ci <- function(df, i) {
  df2 <- df[i, ]
  return(cor(df2$logx, df2$logy, use = "complete.obs"))
}

#Function to get Spearman's correlation and 95% CI by bootstrap
cor_return_boot <- function(logx, logy, R=5000) {
  #Get Spearman's correlation
  df1 <- data.frame(logx=logx, logy=logy)
  df1 <- df1[complete.cases(df1), ]
  df <- data.frame(cbind(rank(df1$logx), rank(df1$logy)))
  names(df) <- c("logx", "logy")
  cor_num <- cor.test(df$logx, df$logy)$estimate
  #Do 5000 bootstrap samples to calculate Spearman's CI, and get 2.5 %ile and 97.5 %ile as
  #95% bootstrap confidence bands
  boot_samp <- boot(df, boot.spearman_ci, R=R)$t
  lower <- quantile(boot_samp, 0.025)
  upper <- quantile(boot_samp, 0.975)
  
  cor_return <- data.frame(result = paste0(round(cor_num, 2), " (", round(lower, 2), ", ", round(upper, 2), ")"), row.names = NULL)
  
  return(cor_return)
}

#Function to get Spearman's correlation and 95% CI
cor_return <- function(logx, logy) {
  #Get Spearman's correlation
  #Determine ranks for X and Y
  
  df1 <- data.frame(logx=logx, logy=logy)
  df1 <- df1[complete.cases(df1), ]
  df <- data.frame(cbind(rank(df1$logx), rank(df1$logy)))
  names(df) <- c("logx", "logy")
  #Get upper and lower confidence bounds using Pearson's correlation on the ranks
  out <- cor.test(df$logx, df$logy, use="complete.obs")
  lower <- out$conf.int[1]
  upper <- out$conf.int[2]
  cor_num <- out$estimate
  #Return Spearman's cor and lower & upper CI bounds as data frame
  cor_return <- data.frame(result = paste0(round(cor_num, 2), " (", round(lower, 2), ", ", round(upper, 2), ")"), row.names = NULL)
  return(cor_return)
}

#This is a function for the pairs plot to remove the middle blank diagonal
gpairs_lower <- function(g){
  g$plots <- g$plots[-(1:g$nrow)]
  g$yAxisLabels <- g$yAxisLabels[-1]
  g$nrow <- g$nrow -1

  g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
  g$xAxisLabels <- g$xAxisLabels[-g$ncol]
  g$ncol <- g$ncol - 1

  g
}


#Function to make paired plots of immune markers
#Input dataframe and columns of immune markers on regular scale

immune_markers_pairs <- function(dat, cols, caption.num) {
  
  #For each manufacturer, plot a paired scatter plot for the immune correlates, with different colors for dose groups, log transform X and Y axes
  
  plots <- (gpairs_lower(ggpairs(dat, columns = cols, ggplot2::aes(col = Subgroup), upper  = list(continuous = "blank"), diag  = list(continuous = "blankDiag"), lower=list(continuous = wrap("points",size=0.5)), columnLabels = gsub('\\.|-|_', ' ', colnames(dat)[cols]), labeller = label_wrap_gen(18), legend = c(2, 1)) + scale_color_manual(values = c("red", "blue", "forestgreen", "tan1", "orchid1", "black", "purple3"), name = "Dose Group", drop=FALSE)  + theme_bw() + theme(strip.background=element_rect(fill="white"), strip.text = element_text(size=5), axis.text = element_text(size=5), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.position = "bottom") + ggtitle(dat[1, ]$manuf)))  
  
  plots[1,1]<-plots[1,1]+scale_y_continuous(trans="log10",limits=c(10,70000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[2,1]<-plots[2,1]+scale_y_continuous(trans="log10",limits=c(100,4000000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[3,1]<-plots[3,1]+scale_y_continuous(trans="log10",limits=c(100,800000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[4,1]<-plots[4,1]+scale_y_continuous(trans="log10",limits=c(5,90000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[5,1]<-plots[5,1]+scale_y_continuous(trans="log10",limits=c(5,70000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[6,1]<-plots[6,1]+scale_y_continuous(trans="log10",limits=c(5,35000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[7,1]<-plots[7,1]+scale_y_continuous(trans="log10",limits=c(5,15000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  
  
  plots[2,2]<-plots[2,2]+scale_y_continuous(trans="log10",limits=c(100,4000000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[3,2]<-plots[3,2]+scale_y_continuous(trans="log10",limits=c(100,800000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[4,2]<-plots[4,2]+scale_y_continuous(trans="log10",limits=c(5,90000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[5,2]<-plots[5,2]+scale_y_continuous(trans="log10",limits=c(5,70000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[6,2]<-plots[6,2]+scale_y_continuous(trans="log10",limits=c(5,35000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  plots[7,2]<-plots[7,2]+scale_y_continuous(trans="log10",limits=c(5,15000))+scale_x_continuous(trans="log10",limits=c(10,70000))
  
   plots[3,3]<-plots[3,3]+scale_y_continuous(trans="log10",limits=c(100,800000))+scale_x_continuous(trans="log10",limits=c(100,4000000))
  plots[4,3]<-plots[4,3]+scale_y_continuous(trans="log10",limits=c(5,90000))+scale_x_continuous(trans="log10",limits=c(100,4000000))
  plots[5,3]<-plots[5,3]+scale_y_continuous(trans="log10",limits=c(5,70000))+scale_x_continuous(trans="log10",limits=c(100,4000000))
  plots[6,3]<-plots[6,3]+scale_y_continuous(trans="log10",limits=c(5,35000))+scale_x_continuous(trans="log10",limits=c(100,4000000))
  plots[7,3]<-plots[7,3]+scale_y_continuous(trans="log10",limits=c(5,15000))+scale_x_continuous(trans="log10",limits=c(100,4000000))
  
  plots[4,4]<-plots[4,4]+scale_y_continuous(trans="log10",limits=c(5,90000))+scale_x_continuous(trans="log10",limits=c(100,800000))
  plots[5,4]<-plots[5,4]+scale_y_continuous(trans="log10",limits=c(5,70000))+scale_x_continuous(trans="log10",limits=c(100,800000))
  plots[6,4]<-plots[6,4]+scale_y_continuous(trans="log10",limits=c(5,35000))+scale_x_continuous(trans="log10",limits=c(100,800000))
  plots[7,4]<-plots[7,4]+scale_y_continuous(trans="log10",limits=c(5,15000))+scale_x_continuous(trans="log10",limits=c(100,800000))
  
  
  plots[5,5]<-plots[5,5]+scale_y_continuous(trans="log10",limits=c(5,70000))+scale_x_continuous(trans="log10",limits=c(5,90000))
  plots[6,5]<-plots[6,5]+scale_y_continuous(trans="log10",limits=c(5,35000))+scale_x_continuous(trans="log10",limits=c(5,90000))
  plots[7,5]<-plots[7,5]+scale_y_continuous(trans="log10",limits=c(5,15000))+scale_x_continuous(trans="log10",limits=c(5,90000))
  
  plots[6,6]<-plots[6,6]+scale_y_continuous(trans="log10",limits=c(5,35000))+scale_x_continuous(trans="log10",limits=c(5,70000))
  plots[7,6]<-plots[7,6]+scale_y_continuous(trans="log10",limits=c(5,15000))+scale_x_continuous(trans="log10",limits=c(5,70000))
  
  plots[7,7]<-plots[7,7]+scale_y_continuous(trans="log10",limits=c(5,15000))+scale_x_continuous(trans="log10",limits=c(5,35000))
  
  plots <- ggmatrix_gtable(plots) 
  
  p <- ggarrange(plots)
  annotate_figure(p, fig.lab = paste0("Figure 6.2.", caption.num), fig.lab.pos = "bottom.right")
}

#Function to make paired plots of immune markers with outcomes
#Input dataframe,  columns of immune markers on regular scale, and columns of outcomes on regular scale
markers_outcomes_pairs <- function(dat, Xcols, Ycols, caption.num) {
  #For each manufacturer, plot a paired scatter plot for the immune correlates x outcomes, with different colors for dose groups, log transform X and Y axes

  p <- gpairs_lower(ggpairs(dat, columns = c(Xcols, Ycols), ggplot2::aes(col = Subgroup), upper  = list(continuous = "blank"), diag  = list(continuous = "blankDiag"), lower=list(continuous = wrap("points",size=0.5)), columnLabels = gsub('\\.|-|_', ' ', colnames(dat)[c(Xcols, Ycols)]), labeller = label_wrap_gen(8), legend = c(2, 1)) + theme_bw() + theme(strip.background=element_rect(fill="white"), strip.text = element_text(size=5), axis.text = element_text(size=5), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.position = "bottom") + ggtitle(dat[1, ]$manuf) + scale_color_manual(values = c("red", "blue", "forestgreen", "tan1", "orchid1", "black", "purple3"), name = "Dose Group", drop=FALSE))
  
    
    #p produces a plot with all pairs of immune correlates and outcomes - however, we are not interested in pairs of immune correlates together or outcomes only
    #Get only pairs of immune correlates with outcomes
  p_y <- vector(mode = "list", length = length(Ycols))
  for (k in 1:length(Ycols)) {
    p_y[[k]] <- lapply(1:length(Xcols), function(j) getPlot(p, i = length(Xcols) + k-1, j = j))
  }
  
  new <- vector(mode = "list", length = length(Ycols)*length(Xcols))
  for (i in 1:length(Ycols)) {
    new[c(c(1:8) + (8*(i-1)))] <- p_y[[i]][c(1:8)]
  }
  
  xmin <- rep(c(10,10,100,100,5,5,5,5), 12)
  xmax <- rep(c(70000,70000,4000000,800000,90000,70000,35000,15000),12)
  ymin <- rep(c(200,100,200,100,200,100,15,15,25,25,25,25),each=8)
  ymax <- rep(c(4e7,5e8,2500000,80000000,1.5e7,1.5e8,60,50,90,60,90,80),each=8)
  
  for (i in 1:length(new)) {
    if (i < 49) {
      new[[i]] <- new[[i]] + scale_x_continuous(trans = "log10", limits = c(xmin[i], xmax[i])) + scale_y_continuous(trans = "log10", limits = c(ymin[i], ymax[i]))
    } else {
      new[[i]] <- new[[i]] + scale_x_continuous(trans = "log10", limits = c(xmin[i], xmax[i])) + scale_y_continuous(limits = c(ymin[i], ymax[i]))
    }
    
  }

  plots <- ggmatrix_gtable(ggmatrix(new, nrow=length(Ycols), ncol=length(Xcols), xAxisLabels = p$xAxisLabels[c(1:length(Xcols))], yAxisLabels = p$yAxisLabels[c(length(Xcols):(length(Xcols)+length(Ycols)-1))], labeller = label_wrap_gen(8), legend = c(2, 1), title = p$title) + theme(strip.background=element_rect(fill="white", color = "black"), legend.position = "bottom"))

  p <- ggarrange(plots)
  annotate_figure(p, fig.lab = paste0("Figure 6.2.", caption.num), fig.lab.pos = "bottom.right")
  
}

#Function to calculate geometric mean using log10 input - calculate mean of log10 values, then exponentiate back
geo_mean <- function(logX) {
  #n <- length(logX)
  geo_mean <- 10^mean(logX, na.rm=T)
  #geo_mean <- prod(x)^(1/n)
  return(geo_mean)
}

#Determine geometric mean ratio - function inputs are x and y on the log scale and regular scale, and the limit of detection for y
gmr <- function(logX, binY) {
  d <- data.frame(logX = logX, binY=binY)
  #Split the data based on whether y is above the LOD, or at or below the LOD
  #Get the geometric mean of the x (correlate) values of each of the 2 groups, and
  #determine the ratio
  gmr <- tryCatch(geo_mean(d[d$binY==1, "logX"])/geo_mean(d[d$binY==0, "logX"]), error = function(e) NA)
  gmr <- tryCatch(if (gmr > 10000) {
    gmr <- formatC(gmr, format = "e", digits = 2)
  } else {
    gmr <- round(gmr, 2)
  }, error=function(e) NA)
  #T-test on the log10-correlate value between the 2 groups (y above LOD, or y at/below LOD)
  #Do same t-test to get confidence interaval bounds, then back transform by raising the upper and lower
  #CI limits to the 10th power
  tst <- tryCatch(t.test(d[d$binY==1, "logX"], d[d$binY==0, "logX"]), error=function(e) NA)
  lower <- tryCatch(10^(tst)$conf.int[1], error=function(e) NA)
  upper <- tryCatch(10^(tst)$conf.int[2], error=function(e) NA)
  upper <- tryCatch(if (upper > 10000) {
    upper <- formatC(upper, format = "e", digits = 2)
  } else {
    upper <- round(upper, 2)
  }, error=function(e) NA)
  lower <- tryCatch(if (lower > 10000) {
    lower <- formatC(lower, format = "e", digits = 2)
  } else {
    lower <- round(lower, 2)
  }, error=function(e) NA)
  df <- data.frame(result = ifelse(!is.na(gmr) & !is.na(lower) & !is.na(upper), yes = paste0(gmr, " (", lower, ", ", upper, ")"), no = ifelse((is.na(lower) | is.na(upper)) & !is.na(gmr), yes = paste0(gmr, " (NA)"), no=NA)), row.names = NULL)
  return(df)
}


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
D <- as.data.frame(readxl::read_xlsx("./data/allDat.xlsx"))
D$manuf <- factor(D$`Vaccine Candidate`, levels = c("Janssen Single Dose", "Janssen Two Dose", "Moderna", "Novavax", "Sanofi", "Placebo"), labels = c("Janssen", "Janssen", "Moderna", "Novavax", "Sanofi", "Placebo"))
D$`Challenge Day` <- factor(D$`Challenge Day`)


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------
# Appendix Changes
D.m<- subset(D,Sex=="Male")
D.f<- subset(D,Sex=="Female")


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
cont.outcomes.log <- grep("log_VL|AUC", names(D))
cont.outcomes.reg <- grep("VL Day 2|AUC", names(D))
bin.outcomes <- grep("bin", names(D))

immune.markers.log <- grep("log_MN|log_PsVNA|log_MSD", names(D))
immune.markers.reg <- grep("(?<!log_)MN|(?<!log_)PsVNA|(?<!log_)MSD", names(D), perl=T)


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
# Appendix Changes
cont.outcomes.log.m <- grep("log_VL|AUC", names(D.m))
cont.outcomes.reg.m <- grep("VL Day 2|AUC", names(D.m))
bin.outcomes.m <- grep("bin", names(D.m))

immune.markers.log.m <- grep("log_MN|log_PsVNA|log_MSD", names(D.m))
immune.markers.reg.m <- grep("(?<!log_)MN|(?<!log_)PsVNA|(?<!log_)MSD", names(D.m), perl=T)

cont.outcomes.log.f <- grep("log_VL|AUC", names(D.f))
cont.outcomes.reg.f <- grep("VL Day 2|AUC", names(D.f))
bin.outcomes.f <- grep("bin", names(D.f))

immune.markers.log.f <- grep("log_MN|log_PsVNA|log_MSD", names(D.f))
immune.markers.reg.f <- grep("(?<!log_)MN|(?<!log_)PsVNA|(?<!log_)MSD", names(D.f), perl=T)


## ----echo = FALSE, warning=F, message=F,   fig.width=8, fig.height=10---------------------------------------------------------------------
# Appendix Changes

#for (i in 1:length(unique(D$manuf))) {
#  print(immune_markers_pairs(D[D$manuf == (unique(D$manuf)[i]), ], cols=immune.markers.reg, i))
#  
#}


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
# Appendix Changes

# Male 
combs.m <- combn(immune.markers.log.m, 2)
outcomes.m <- matrix(NA, nrow=length(immune.markers.log.m), ncol = length(immune.markers.log.m))

for (i in 1:ncol(combs.m)) {
  outcomes.m[which(immune.markers.log.m == combs.m[2, i]), which(immune.markers.log.m == combs.m[1, i])] <- cor_return(D.m[ ,combs.m[1, i]], D.m[ ,combs.m[2, i]])[1, 1]
}

outcomes.m <- as.data.frame(outcomes.m)
names(outcomes.m) <- gsub("_", " ", names(D.m)[immune.markers.log.m])
names(outcomes.m) <- gsub(".", " ", names(outcomes.m), fixed=T)

row.names(outcomes.m) <- gsub("_", " ", names(D.m)[immune.markers.log.m])
row.names(outcomes.m) <- gsub(".", " ", row.names(outcomes.m), fixed=T)

outcomes.m <- apply(outcomes.m, c(1, 2), function(x){ifelse(is.na(x), "-", x)})

kable(outcomes.m, row.names = TRUE, caption = "Table 6.2.1m: (Male Only) Spearman Correlations of Virus-Specific Antibodies", align = "c", escape = F) %>% kable_styling(latex_options=c(position = "HOLD_position")) %>% column_spec(1, width = "5em") %>% column_spec(c(2:(length(immune.markers.log.m)+1)), width = "3em")


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
# Appendix Changes
# Female 

combs.f <- combn(immune.markers.log.f, 2)
outcomes.f <- matrix(NA, nrow=length(immune.markers.log.f), ncol = length(immune.markers.log.f))

for (i in 1:ncol(combs.f)) {
  outcomes.f[which(immune.markers.log.f == combs.f[2, i]), which(immune.markers.log.f == combs.f[1, i])] <- cor_return(D.f[ ,combs.f[1, i]], D.f[ ,combs.f[2, i]])[1, 1]
}

outcomes.f <- as.data.frame(outcomes.f)
names(outcomes.f) <- gsub("_", " ", names(D.f)[immune.markers.log.f])
names(outcomes.f) <- gsub(".", " ", names(outcomes.f), fixed=T)

row.names(outcomes.f) <- gsub("_", " ", names(D.f)[immune.markers.log.f])
row.names(outcomes.f) <- gsub(".", " ", row.names(outcomes.f), fixed=T)

outcomes.f <- apply(outcomes.f, c(1, 2), function(x){ifelse(is.na(x), "-", x)})

kable(outcomes.f, row.names = TRUE, caption = "6.2.1f: (Female Only) Spearman Correlations of Virus-Specific Antibodies", align = "c", escape = F) %>% kable_styling(latex_options=c(position = "HOLD_position")) %>% column_spec(1, width = "5em") %>% column_spec(c(2:(length(immune.markers.log.f)+1)), width = "3em")




## ----echo = FALSE,  warning = FALSE, message=FALSE, fig.width=8, fig.height=12------------------------------------------------------------
# Appendix Changes

#for (i in 1:length(unique(D$manuf))) {
#  print(markers_outcomes_pairs(D[D$manuf == (unique(D$manuf)[i]), ], Xcols=immune.markers.reg, Ycols=cont.outcomes.reg, #i+length(unique(D$manuf))))
#  
#}


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
# Appendix Changes 

# Male
combs.m <- combn(cont.outcomes.log.m, 2)
outcomes.m <- matrix(NA, nrow=length(cont.outcomes.log.m), ncol = length(cont.outcomes.log.m))

for (i in 1:ncol(combs.m)) {
  outcomes.m[which(cont.outcomes.log.m == combs.m[2, i]), which(cont.outcomes.log.m == combs.m[1, i])] <- cor_return(D.m[ ,combs.m[1, i]], D.m[ ,combs.m[2, i]])[1, 1]
}

outcomes.m <- as.data.frame(outcomes.m)
names(outcomes.m) <- gsub("_", " ", names(D.m)[cont.outcomes.log.m])
names(outcomes.m) <- gsub(".", " ", names(outcomes.m), fixed=T)
#names(outcomes) <- str_wrap(names(outcomes), 20)
names(outcomes.m) <- gsub("PostChallenge", "Post Challenge", names(outcomes.m))

row.names(outcomes.m) <- gsub("_", " ", names(D.m)[cont.outcomes.log.m])
row.names(outcomes.m) <- gsub(".", " ", row.names(outcomes.m), fixed=T)

outcomes.m <- apply(outcomes.m, c(1, 2), function(x){ifelse(is.na(x), "-", x)})

kable(outcomes.m, row.names = TRUE, caption = "6.2.2m: (Male Only) Spearman Correlations of Viral Load Outcomes", align = "c", escape = F) %>% kable_styling(latex_options=c(position = "HOLD_position")) %>% column_spec(1, width = "6em") %>% column_spec(c(2:(length(cont.outcomes.log.m)+1)), width = "2.5em")


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
# Appendix Changes
# Female
combs.f <- combn(cont.outcomes.log.f, 2)
outcomes.f <- matrix(NA, nrow=length(cont.outcomes.log.f), ncol = length(cont.outcomes.log.f))

for (i in 1:ncol(combs.f)) {
  outcomes.f[which(cont.outcomes.log.f == combs.f[2, i]), which(cont.outcomes.log.f == combs.f[1, i])] <- cor_return(D.f[ ,combs.f[1, i]], D.f[ ,combs.f[2, i]])[1, 1]
}

outcomes.f <- as.data.frame(outcomes.f)
names(outcomes.f) <- gsub("_", " ", names(D.f)[cont.outcomes.log.f])
names(outcomes.f) <- gsub(".", " ", names(outcomes.f), fixed=T)
#names(outcomes) <- str_wrap(names(outcomes), 20)
names(outcomes.f) <- gsub("PostChallenge", "Post Challenge", names(outcomes.f))

row.names(outcomes.f) <- gsub("_", " ", names(D.f)[cont.outcomes.log.f])
row.names(outcomes.f) <- gsub(".", " ", row.names(outcomes.f), fixed=T)

outcomes.f <- apply(outcomes.f, c(1, 2), function(x){ifelse(is.na(x), "-", x)})

kable(outcomes.f, row.names = TRUE, caption = "6.2.2f: (Female Only) Spearman Correlations of Viral Load Outcomes", align = "c", escape = F) %>% kable_styling(latex_options=c(position = "HOLD_position")) %>% column_spec(1, width = "6em") %>% column_spec(c(2:(length(cont.outcomes.log.f)+1)), width = "2.5em") 




## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
# Appendix Changes

# Male
combs2.m <- expand.grid(cont.outcomes.log.m, immune.markers.log.m)

outcomes2.m <- matrix(NA, ncol=length(immune.markers.log.m), nrow = length(cont.outcomes.log.m))


for (i in 1:nrow(combs2.m)) {
  outcomes2.m[which(cont.outcomes.log.m == combs2.m[i, 1]), which(immune.markers.log.m == combs2.m[i, 2])] <- cor_return(D.m[ ,combs2.m[i, 1]], D.m[ ,combs2.m[i, 2]])[1, 1]
}

outcomes2.m <- as.data.frame(outcomes2.m)
names(outcomes2.m) <- gsub("_", " ", names(D.m)[immune.markers.log.m])
names(outcomes2.m) <- gsub(".", " ", names(outcomes2.m), fixed=T)

row.names(outcomes2.m) <- gsub("_", " ", names(D)[cont.outcomes.log.m])
row.names(outcomes2.m) <- gsub(".", " ", row.names(outcomes2.m), fixed=T)

kable(outcomes2.m, row.names = TRUE, caption = "6.2.3m: (Male Only) Spearman Correlations of Virus-Specific Antibodies with Viral Load Measurements", align = "c", escape = F) %>% kable_styling(latex_options=c(position = "HOLD_position")) %>% column_spec(1, width = "6em") %>% column_spec(c(2:(length(immune.markers.log.m)+1)), width = "3em")


## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
# Appendix Changes
# Female
combs2.f <- expand.grid(cont.outcomes.log.f, immune.markers.log.f)

outcomes2.f <- matrix(NA, ncol=length(immune.markers.log.f), nrow = length(cont.outcomes.log.f))


for (i in 1:nrow(combs2.f)) {
  outcomes2.f[which(cont.outcomes.log.f == combs2.f[i, 1]), which(immune.markers.log.f == combs2.f[i, 2])] <- cor_return(D.f[ ,combs2.f[i, 1]], D.f[ ,combs2.f[i, 2]])[1, 1]
}

outcomes2.f <- as.data.frame(outcomes2.f)
names(outcomes2.f) <- gsub("_", " ", names(D.f)[immune.markers.log.f])
names(outcomes2.f) <- gsub(".", " ", names(outcomes2.f), fixed=T)

row.names(outcomes2.f) <- gsub("_", " ", names(D)[cont.outcomes.log.f])
row.names(outcomes2.f) <- gsub(".", " ", row.names(outcomes2.f), fixed=T)

kable(outcomes2.f, row.names = TRUE, caption = "6.2.3f: (Female Only) Spearman Correlations of Virus-Specific Antibodies with Viral Load Measurements", align = "c", escape = F) %>% kable_styling(latex_options=c(position = "HOLD_position")) %>% column_spec(1, width = "6em") %>% column_spec(c(2:(length(immune.markers.log.f)+1)), width = "3em")



## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------
# Appendix Changes
#combs3 <- expand.grid(immune.markers.log, bin.outcomes)

#outcomes3 <- matrix(NA, nrow=length(immune.markers.log), ncol = length(bin.outcomes))


#for (i in 1:nrow(combs3)) {
#  outcomes3[which(immune.markers.log == combs3[i, 1]), which(bin.outcomes == combs3[i, 2])] <- gmr(D[ ,combs3[i, 1]], D[ ,combs3[i, #2]])[1, 1]
#}

#outcomes3 <- as.data.frame(outcomes3)
#names(outcomes3) <- gsub("_", " ", names(D)[bin.outcomes])
#names(outcomes3) <- gsub(".", " ", names(outcomes3), fixed=T)

#row.names(outcomes3) <- gsub("_", " ", names(D)[immune.markers.log])
#row.names(outcomes3) <- gsub(".", " ", row.names(outcomes3), fixed=T)

#kable(outcomes3, row.names = TRUE, caption = "Geometric Mean Ratios of Virus-Specific Antibodies with Binary Viral Load Measurements", #align = "c", escape = F) %>% kable_styling(latex_options=c(position = "HOLD_position")) %>% column_spec(1, width = "5em") %>% #column_spec(c(2:(length(bin.outcomes)+1)), width = "5em") %>% footnote(general = "NA is indicated when there is only 1 negative signal #and so a 95% CI can not be calculated.")

