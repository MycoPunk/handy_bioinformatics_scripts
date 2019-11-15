##this script takes snpEff data (or a smilarly formatted tab delimited file from BCF tools)
##and a depth sequencing profile from mosdepth 
##and computes the relationship between the number of no-calls and the sequencing depth 
##with the intention of determining which strains to target for resequencing for optimized SNP calling

#set wd
setwd("<>")
options(stringsAsFactors = TRUE)

#load libraries
library(dplyr)

#load files 
#read in read coverage file
depth_df <-read.csv("<>.tab.gz",sep="\t",header=F)
colnames(depth_df) <- c("Chrom","START","END","DEPTH","TYPE","STRAIN")
rownames(depth_df) <- depth_df$STRAIN


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
                     ylim = c(0,550),
                     main = "Candida: mean read depth over 5kb windows")
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

#dev.off()    

#names of isolates below the cut off 
isolates_below_detpth_cutoff<- mean_depth_ea_strain_sorted[mean_depth_ea_strain_sorted$depth <15,]
isolates_below_detpth_cutoff

#names of isolates above variance cut off 
isolates_above_variance_cutoff<- mean_depth_ea_strain_sorted[mean_depth_ea_strain_sorted$var > 5000,]
isolates_above_variance_cutoff 

#non-overlapping names of isolates that need to be resequenced
combined_fail<- as.data.frame(rbind(isolates_below_detpth_cutoff, isolates_above_variance_cutoff))
isolates_that_dont_pass_QC<- as.data.frame(unique(combined_fail$strain))
isolates_that_dont_pass_QC

#get a full table of data for the isolates that pass QC and don't need resequencing 
isolates_that_do_pass_QC_full<- mean_depth_ea_strain_sorted%>% filter(!strain %in% isolates_that_dont_pass_QC)

#print the isoaltes that didn't pass QC to a file
write.table(isolates_that_dont_pass_QC, "isolates_that_do_not_pass_QC.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)



###get the number of missed calls from the snpEFF file
#load main dat files from snpEff 
snpEff_data<-read.delim("ClusNanoL1B_201902.snpEff.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

##clean the data for easy processing##

#NOTE - this is only if the input is from the snpEFF matrix file - otherwise skip 
#remove the brackets form ALT 
#snpEff_data$ALT<-gsub("\\[|\\]", "", snpEff_clade_1$ALT)
#remove the "/" from all the the SNP calls 
#snpEff_data[] <- lapply(snpEff_data, gsub, pattern='/', replacement="")

#subset to get only the columns of interest - note you may have to change the starting col position depnding on how the data was generated
subset_snpEff_data<- snpEff_data[,5:(ncol(snpEff_data) -1)]

#get number of no-calls per strain
no_calls_ea_strain <- as.data.frame(colSums(subset_snpEff_data[, ] == "."))
no_calls_ea_strain$strain<- row.names(no_calls_ea_strain)
colnames(no_calls_ea_strain) <- c("n_no_calls", "strain")

#check that the strains are in the same order and that the naming conventions are the same between strains 
#note - depending on how the data was generated, snpEFF attaches col. numbers, adds GT designators etc. 

#test<- cbind(mean_depth_ea_strain, no_calls_ea_strain$strain)
#View(test) NOPE, not in the same order. Need to rename

#annoying .. need to reformat- change no_calls_ea strain, to match mean depth naming scheme
#remove GT from the naming scheme
#no_calls_ea_strain$strain_reformat<- lapply(no_calls_ea_strain$strain, function(x) {gsub(".GT", "", x)})
#git rid of "X."
#no_calls_ea_strain$strain_reformat2<- lapply(no_calls_ea_strain$strain_reformat, function(x) {gsub("X.", "", x)})
#git rid of the <numeric>.
#no_calls_ea_strain$strain_reformat3<- lapply(no_calls_ea_strain$strain_reformat2, function(x) {gsub("^[^.]*.", "", x)})
##reorder the no_calls_ea_strain df based on the new naming scheme:
#no_calls_ea_strain$strain_reformat3<-as.character(no_calls_ea_strain$strain_reformat3)
#no_calls_ea_strain_orderd<- no_calls_ea_strain[order(no_calls_ea_strain$strain_reformat3),]

##reorder mean_depth_ea_strain by strain name
#mean_depth_ea_strain_ordered<- mean_depth_ea_strain[order(mean_depth_ea_strain$strain),]
##in this example - there was a problem, an extra entry in the depth file that wasn't in the snpEFF file 
##remember to call dim and check before proceeing
#grep("ATCC42720", names(snpEff_data), value = TRUE)
##the strain ATCC42720 is in the depth file, but not in the snpEff file foe some reason - remove this isolate to compare the two sets. 
#mean_depth_ea_strain_ordered_2<- mean_depth_ea_strain_ordered[mean_depth_ea_strain_ordered$strain != "ATCC42720",]
#class(no_calls_ea_strain$strain_reformat3)

#check that they're the same now
test<- cbind(as.data.frame(no_calls_ea_strain_orderd$strain_reformat3), as.data.frame(mean_depth_ea_strain_ordered_2$strain)) 
#View(test) #OK

#bind the resultant datasets, once the naming conventions are the same, and the datasets are in the same order
depth_vs_missing_calls_df<- as.data.frame(cbind(n_no_calls= no_calls_ea_strain_orderd$n_no_calls, depth=mean_depth_ea_strain_ordered_2$depth))
row.names(depth_vs_missing_calls_df)<- no_calls_ea_strain_orderd$strain_reformat3

#plot the relationship between depth and no-calls in the dataset  
plot(depth_vs_missing_calls_df$depth, depth_vs_missing_calls_df$n_no_calls,
     main = "read depth by no-calls", 
     xlab = "read depth", 
     ylab = "n missing calls")

lmt<- lm(depth_vs_missing_calls_df$n_no_calls ~ depth_vs_missing_calls_df$depth)
sum1<- summary(lmt)
abline(lmt)
rsq <- bquote(italic(R)^2 == .(format(summary(lm)$adj.r.squared, digits=2)))
text(x = 350, y = 20000, labels = rsq)
pval<- round(sum1$coefficients[2,4], digits = 3)
pval <- bquote(italic(P) == .(pval))
text(x = 350, y = 25000, labels = pval)
