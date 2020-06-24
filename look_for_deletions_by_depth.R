##Script to look at large deletion in A. fum. Chr 6 in the Clade 2 isolates.


#load modules
require(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)



#set dir
setwd("~/Desktop/IHE/depth_files")
options(stringsAsFactors = TRUE)

#load files 

#too big to load locally
#depth_df <-read.csv("mosdepth.1bp.gg_nonorm.tab.gz",sep="\t",header=F)
#depth_df <-read.csv("mosdepth.1bp.gg.tab.gz",sep="\t",header=F)

#read in read coverage file
depth_df <-read.csv("mosdepth.100bp.gg_nonorm.tab.gz",sep="\t",header=F)
#depth_df <-read.csv("mosdepth.100bp.gg.tab.gz",sep="\t",header=F)

colnames(depth_df) <- c("Chrom","START","END","DEPTH","TYPE","STRAIN")
#rownames(depth_df) <- depth_df$STRAIN

#get mean depth for each isolate
mean_depth_ea_strain<- aggregate(depth_df$DEPTH ~depth_df$STRAIN, FUN=mean)
colnames(mean_depth_ea_strain)<- c("strain", "depth")

#get depth variance across isolates
var_depth_ea_strain<- aggregate(depth_df$DEPTH ~ depth_df$STRAIN, FUN=var)
colnames(var_depth_ea_strain)<- c("strain", "var")

#combine into single df
mean_depth_ea_strain$var<- var_depth_ea_strain$var

#sort the df by mean, then by variance
mean_depth_ea_strain_sorted<- mean_depth_ea_strain[order(mean_depth_ea_strain$depth, mean_depth_ea_strain$var),]


#plot this 
#pdf(<X>_quality_cutoffs.pdf)
depth_plot<- barplot(mean_depth_ea_strain_sorted$depth, 
                     names.arg = mean_depth_ea_strain_sorted$strain, 
                     las = 2, 
                     cex.names = 0.25,
                     cex.axis = .5,
                     ylab = "mean read depth",
                     ylim = c(0,200),
                     main = "mean read depth over 5kb windows")
abline(h=15, col = "red")


min(mean_depth_ea_strain_sorted$depth)
max(mean_depth_ea_strain_sorted$depth)

##do the same for variance 
#variance graph
variance_plot<- barplot(mean_depth_ea_strain_sorted$var,
                        names.arg = mean_depth_ea_strain_sorted$strain,
                        las = 2,
                        cex.names = 0.25,
                        cex.axis = .5,
                        ylab = "depth variance",
                        main= "mean variance in depth over 5kb windows")
abline(h = 5000, col = "red")


###look at depth across all Chrs for each strain 

#isolate strain 
C1_TP1<- depth_df[depth_df$STRAIN == "DMC2_AF100-1_3",]
C1_TP1 <- C1_TP1[C1_TP1$DEPTH > 0,]

#assign color pallet
chr_colors<- c("#D53E4F",	
               "#F46D43",
               "#FDAE61",
               "#FEE08B",
               "#E6F598",
               "#66C2A5",
               "#3288BD",
               "#3E5F83",
               "#683D84",
               "#A03A7E")[C1_TP1$Chrom]


#range of concern 2337367:2355185, make highlight here
plot(C1_TP1$START, C1_TP1$DEPTH, 
     col=chr_colors, 
     xlab="", 
     ylab="", 
     xlim = c(min(C1_TP1$START),max(C1_TP1$END)), 
     pch=16,
     cex=.3)
legend(4500000, 6000, legend=c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr 7", "Chr 9", "Mito"),
       col=chr_colors, fill = rev(unique(chr_colors)), cex=0.4)



###look at Chr 6 specifically. 
depth_df<- depth_df[depth_df$Chrom == "Chr6_A_fumigatus_Af293",]
#zoom in to 2337367:2355185 (+/- 10,000)
depth_df_test<- depth_df[depth_df$START > 2324367,]
depth_df<- depth_df_test[depth_df_test$START < 2365185,]


#split out the strains and remove depth of zero counts to make graphing cleaner
C1_TP1<- depth_df[depth_df$STRAIN == "DMC2_AF100-1_3",]
C1_TP1 <- C1_TP1[C1_TP1$DEPTH > 0,]
C1_TP3<- depth_df[depth_df$STRAIN == "DMC2_AF100-3B",]
C1_TP3 <- C1_TP3[C1_TP3$DEPTH > 0,]
C1_TP4<- depth_df[depth_df$STRAIN == "DMC2_AF100-4B",]
C1_TP4 <- C1_TP4[C1_TP4$DEPTH > 0,]
C1_TP5<- depth_df[depth_df$STRAIN == "DMC2_AF100-5B",]
C1_TP5 <- C1_TP5[C1_TP5$DEPTH > 0,]
C1_TP6<- depth_df[depth_df$STRAIN == "DMC2_AF100-6B",]
C1_TP6 <- C1_TP6[C1_TP6$DEPTH > 0,]
C1_TP10<- depth_df[depth_df$STRAIN == "AF100-10B",]
C1_TP10 <- C1_TP10[C1_TP10$DEPTH > 0,]
C1_TP12<- depth_df[depth_df$STRAIN == "AF100-12_7",]
C1_TP12 <- C1_TP12[C1_TP12$DEPTH > 0,]

C2_TP8<- depth_df[depth_df$STRAIN == "DMC2_AF100-8B",]
C2_TP8 <- C2_TP8[C2_TP8$DEPTH > 0,]
C2_TP9<- depth_df[depth_df$STRAIN == "DMC2_AF100-9B",]
C2_TP9 <- C2_TP9[C2_TP9$DEPTH > 0,]
C2_TP10<- depth_df[depth_df$STRAIN == "AF100-10_5",]
C2_TP10 <- C2_TP10[C2_TP10$DEPTH > 0,]
C2_TP11<- depth_df[depth_df$STRAIN == "AF100-11_3",]
C2_TP11 <- C2_TP11[C2_TP11$DEPTH > 0,]
C2_TP12<- depth_df[depth_df$STRAIN == "DMC2_AF100-12_9",]
C2_TP12 <- C2_TP12[C2_TP12$DEPTH > 0,]

space<- depth_df[depth_df$STRAIN == "1F1SW_F4",]
space<- space[space$DEPTH >0,]
#graph and split chrs by factor (8 chrs + mito)


#range of concern 2337367:2355185, make highlight here

#view just the first one to adjust plotting
#plot(C1_TP1$START, C1_TP1$DEPTH, 
#     col=chr_colors, 
#     xlab="", 
#     ylab="", 
#     xlim = c(min(C1_TP1$START),max(C1_TP1$END)), 
#     pch=16,
#     cex=.3,
#     rect(2337367, 0, 2355185, max(C1_TP1$DEPTH), col = "#FEE08B99", border=F))

#write funtion for plotting
plot_fun <- function(data){
  plot(data$START, data$DEPTH, 
       col=chr_colors, 
       xlab="", 
       ylab="", 
       pch=16,
       cex=2, #set to .3 for viewing 
       main = deparse(substitute(data)),
       rect(2337367, 0, 2355185, max(data$DEPTH), col = "#FEE08B99", border=F))}

#plot
par(mfrow=c(4,2), mar=c(2,2,2,2)) 
plot_fun(C1_TP1)
plot_fun(C1_TP3)
plot_fun(C1_TP4)
plot_fun(C1_TP5)
plot_fun(C1_TP6)
plot_fun(C1_TP10)
plot_fun(C1_TP12)

plot_fun(C2_TP8)
plot_fun(C2_TP9)
plot_fun(C2_TP10)
plot_fun(C2_TP11)
plot_fun(C2_TP12)
plot_fun(space)

dev.off()
#there does not appear to be any large deletions in chr 6. . . . 

#genes within this region (from fungi DB)
#Afu6g09600 - zinc metallopeptadase
#Afu6g09610 - NRPS9
#Afu6g09620 - hyp protein
#Afu6g09630 - TF (C6 finger domain transcription factor gliZ)
#Afu6g09640 - Aminotransferase gliI, putative
#Afu6g09650 - Dipeptidase gliJ
#Afu6g09660 - Nonribosomal peptide synthetase gliP

HOG_genes<- c("Afu2g04680", 
              "Afu2g05740",
              "Afu5g08420",
              "Afu5g06420",
              "Afu1g15950",
              "Afu1g12940",
              "Afu2g00660",
              "Afu4g01020",
              "Afu6g10240",
              "Afu4g10280",
              "Afu5g08390",
              "Afu6g12522",
              "Afu1g10940",
              "Afu5g09100")


genes_over_deletion<- c("Afu6g09590",
                        "Afu6g09600", 
                        "Afu6g09610", 
                        "Afu6g09620", 
                        "Afu6g09630", 
                        "Afu6g09640", 
                        "Afu6g09650",
                        "Afu6g09660") 


intersect(genes_over_deletion,HOG_genes) #- none. 






