setwd("/Volumes/My_life/Plankton/Plankton/plankton2.0")

# DOWNLOADING (installing) packages ----
# YOU ONLY HAVE TO execute just once when first using a script:
# source("https://bioconductor.org/biocLite.R")
# biocLite("labdsv")
# biocLite("dada2")
# biocLite('phyloseq')
# biocLite('ShortRead')

# Load the libraries you need (DO THIS EVERY TIME YOU OPEN A SCRIPT)
library(labdsv)  # this loads MASS, which conflicts with the "select" function in dplyr(tidyverse). Load tidyverse last.
library(dada2)
library(ShortRead)
library(phyloseq)
library(tidyverse) # for data wrangling and ggplot2

# Making sample information table -----

sample_info <- read.delim("orig_raw_data_total.mapping.txt")
head(sample_info)
summary(sample_info$Site)

# Make a vector called `inshore_sites` that lists all of the inshore sites
inshore_sites <- c("PuntaDonato", "STRIPoint", "Cristobal", "PuntaLaurel")

# Make a new column called `siteType` that (as a factor) enters the text "inshore" if the site name is contained in the vector `inshore_sites` and enters the text "offshore" if it isn't
sample_info$siteType <- as.factor(ifelse(sample_info$Site %in% inshore_sites, "inshore","offshore"))
summary(sample_info)

# Rename the column called `Number` to `tech_rep` and make it a factor
names(sample_info)
colnames(sample_info)[8] <- "techRep"

# Get rid of columns you don't need. Only keep SampleID, Site, Time, techRep, and siteType
names(sample_info)
sam_info <- sample_info %>% 
            dplyr::select(SampleID, Site, Time, techRep, siteType)
head(sam_info)

# Set path to unzipped, renamed fastq files (sequencing samples)
path <- "/Volumes/My_life/Plankton/Plankton/DADA2/Plankton_data/"
fns <- list.files(path)
fns

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files

# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq; OTHERWISE MODIFY
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) #the last number will select the field for renaming
head(sample.names)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Visualize raw data 

# First, lets look at quality profile of R1 reads. Plot the first and last 4 samples.
# This function plots a visual summary of the distribution of quality scores as a function of sequence position for the input fastq file.
plotQualityProfile(fnFs[c(1:4)])
plotQualityProfile(fnFs[c(74:77)])
# Where do the base call qualities get lower than ~30? -----> ~250 bp in forward reads

# Then look at quality profile of R2 reads
plotQualityProfile(fnRs[c(1,2,3,4)])
plotQualityProfile(fnRs[c(74:77)])
# Where do the base call qualities get lower than ~30? -----> ~200 bp in forward reads


# The reverse reads are significantly worse quality, especially at the end, common in Illumina sequencing.
# This isn’t too worrisome, DADA2 incorporates quality information into its error model which makes the algorithm more robust, 
# but trimming as the average qualities crash is still a good idea as long as our reads will still overlap. 

# Recommend trimming where quality profile crashes

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter and Trim (this takes awhile) 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
              truncLen=c(250,200),
              maxN=0, # DADA does not allow Ns
              maxEE=c(1,1), # allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
              truncQ=2, # truncate reads at the first instance of a quality score less than or equal to 2
              trimLeft=c(24,19), #N nucleotides to remove from the start of each read to remove sequencing primers
              rm.phix=TRUE, # remove reads matching phiX genome
              matchIDs=TRUE, # enforce matching between id-line sequence identifiers of F and R reads
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)
class(out)

# How many reads did we lose?
summary(out)

out_stats <- as.data.frame(out) %>% mutate(perc_reads_remaining = reads.out/reads.in*100)
mean(out_stats$perc_reads_remaining) # we only lost 20% of the reads
sum(out_stats)

# Save the out file 
#save(sam_info, out, filtFs, filtRs, sample.names, file="outData.RData")

# For DADA2 alaysis-------
load("outData.RData")

# A word on Expected Errors vs a blanket quality threshold
# Take a simple example: a read of length two with quality scores Q3 and Q40, corresponding to error 
# probabilities P=0.5 and P=0.0001. The base with Q3 is much more likely to have an error than the base with Q40 
# (0.5/0.0001 = 5,000 times more likely), so we can ignore the Q40 base to a good approximation. Consider a large 
# sample of reads with (Q3, Q40), then approximately half of them will have an error (because of the P=0.5 from 
# the Q2 base). We express this by saying that the expected number of errors in a read with quality scores (Q3, Q40) is 0.5.
# As this example shows, low Q scores (high error probabilities) dominate expected errors, but this information 
# is lost by averaging if low Qs appear in a read with mostly high Q scores. This explains why expected errors 
# is a much better indicator of read accuracy than average Q.

# Learn Error Rates
# DADA2 learns its error model from the data itself by alternating estimation of the error rates 
# and the composition of the sample until they converge on a jointly consistent solution (this is similar to the E-M algorithm)
# As in many optimization problems, the algorithm must begin with an initial guess, for which the maximum possible 
# error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).


setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
# Forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
# Reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate reads
# Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding 
# “abundance”: the number of reads with that unique sequence. 
# Dereplication substantially reduces computation time by eliminating redundant comparisons.
# DADA2 retains a summary of the quality information associated with each unique sequence. The consensus 
# quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. 
# These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer Sequence Variants 
#takes a long time

setDadaOpt(BAND_SIZE=32)
# DADA analysis on forward and reverse reads 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# now, look at the dada class objects by sample
# will tell how many 'real' variants in unique input seqs
# By default, the dada function processes each sample independently, but pooled processing is available with 
# pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. 
# See our discussion about pooling samples for sample inference. 
dadaFs[[70]]
dadaRs[[70]]


# Merge paired reads. Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample

head(mergers[[70]])
summary((mergers[[70]]))

# We now have a data.frame for each sample with the merged $sequence, its $abundance, and the 
# indices of the merged $forward and $reverse denoised sequences. Paired reads that did not 
# exactly overlap were removed by mergePairs.

# Construct sequence table
# a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab))), xlab="BP", ylab="abundance", main="Histogram of sequence lengths")

# 365-386 based on histogram of sequence lengths and where there were more (at a faily conservative length)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(365 ,386)] 
table(nchar(getSequences(seqtab2)))
dim(seqtab2)

 
# Remove chimeras 

# The core dada method removes substitution and indel errors, but chimeras remain. 
# Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
# than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
# a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

# The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
# but can be substantial. Here chimeras make up about 36% of the inferred sequence variants (138-89 = 49 => 49/138), 
# BUT those variants account for only about 0.5% of the total sequence reads
# Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)

# Track Read Stats -----

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track) 

write.csv(track,file="final_sequence.csv",row.names=TRUE,quote=FALSE)
# Checking how much was lost over all through out the analysis
tracklost <- as.data.frame(track) %>% 
  mutate(remaining_afterfilter = nonchim/input*100)
mean(tracklost$remaining_afterfilter) 
# lost 61% of reads

# How much was lost at each step all looks good 
tracklostind <- as.data.frame(track) %>% 
 mutate(remaining_filtered = filtered/input*100) %>% # lost 20% of reads
 mutate(remaining_denoised = denoised/filtered*100) %>% # lost 0% of reads
 mutate(remaining_merged = merged/denoised*100) %>% # lost 49% of reads
 mutate(remaining_tabled = tabled/merged*100) %>% # lost <1% of reads
 mutate(remaining_nonchim = nonchim/tabled*100) # lost <1% of reads
mean(tracklostind$remaining_nonchim) 

head(tracklost)

# Assign Taxonomy
# It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to classify sequence variants taxonomically. 
# DADA2 provides a native implementation of the RDP's naive Bayesian classifier. 
# The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, and outputs 
# the taxonomic assignments with at least minBoot bootstrap confidence.
# Here, I have supplied a modified version of the GeoSymbio ITS2 database (Franklin et al. 2012)

# Silva downloaded from database and then matched 
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa", 
                       minBoot = 5,
                       multithread = TRUE,
                       tryRC = TRUE,
                       outputBootstraps = FALSE)
summary(taxa)

# All Eukaryotes!!!
# Most class = Maxillopoda

# Rownames of the sample variable table (sam_info) and the OTU table (seqtab.nochim) don't match
rownames(seqtab.nochim)
rownames(sam_info) <- sam_info$SampleID
rownames(sam_info)

# Make them match
rownames(seqtab.nochim) <- sub("-",".",rownames(seqtab.nochim))
rownames(seqtab.nochim)

identical(sort(rownames(seqtab.nochim)),sort(rownames(sam_info)))
# they match now!

# START PHYLOSEQ (load) ------
# save(sam_info, seqtab.nochim, taxa, file = "dada2_output.Rdata")

load("dada2_output.Rdata") 
# CLEANED TILL HERE----
# Make OTU - sequence - taxa table for later
colnames(seqtab.nochim)[c(1:3)]
colnames(taxa)[c(1:10)]
rownames(taxa)[c(1)]

seqtab.trans <- as.data.frame(t(seqtab.nochim))
head(seqtab.trans)

# create short OTU IDs
seqtab.trans$ids <- paste0("OTU", seq(1, length(colnames(seqtab.nochim))))
head(seqtab.trans)

# merge with taxa
otu_taxa_seq <- merge(seqtab.trans, taxa, by = 0)
otu_taxa <- select(otu_taxa_seq, "Order", "ids")

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sam_info), 
               tax_table(taxa))
ps


# Let's look at what we have for phyloseq
otu_table(ps)[1:5, 1:5]
# yikes, nasty column names. Change that.


colnames(seqtab.nochim) <- ids
head(seqtab.nochim)[2,]

# Replace taxa names in the phyloseq object
taxa_names(ps) <- ids

# Try again...
otu_table(ps)[1:5, 1:5]
# Much better! Rows = samples. Columns = OTUs. Abundances (counts) fill the cells.

# What are the sample varaibles in our 'ps' object?
sample_variables(ps)
# "SampleID" "Site"     "Time"     "techRep"  "siteType"

# taxonomic ranks are "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"  
rank_names(ps)

# How many samples do we have?
nsamples(ps)
# 77, correct

# What does the taxonomic table look like?
tax_table(ps)[1:5, 1:6]
# Each OTU is associated with six levels of taxonomy (KPCOFG --- no species)

# Start stats ----

# Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset**
# total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)

# Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).

# Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index.

# plot Shannon and Simpson Diveristy by Site, color by Site Type (inshore vs. offshore)
plot_richness(ps, 
              x="Site", 
              measures=c("Shannon", "Simpson"), 
              color="siteType") +
              geom_jitter()+
              scale_color_manual(values=c("salmon","royalblue1"))+
              theme_bw()

# plot Shannon and Simpson Diveristy by Site Type
plot_richness(ps, 
              x="siteType", 
              measures=c("Shannon", "Simpson"), 
              color="siteType") + 
              geom_jitter()+
              scale_color_manual(values=c("salmon","royalblue1"))+
              theme_bw()



# Ordinate Samples
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray",k=20)

# *** No convergence -- monoMDS stopping criteria:
# ^^^^^^^^ NMDS bray doesn't converge...

# Bar-plots

# The top 30 most abundant OTUs

top30 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:30]
ps.top30 <- transform_sample_counts(ps, function(x) x/sum(x))
ps.top30 <- prune_taxa(top30, ps.top30)

plot_bar(ps.top30, x="Site", fill="Class") + 
  facet_wrap(~siteType, scales="free_y") + 
  theme_bw()

btm30 <- names(sort(taxa_sums(ps), decreasing=FALSE))[1:30]
ps.btm30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.btm30 <- prune_taxa(btm30, ps.btm30)

plot_bar(ps.btm30, x="Site", fill="Phylum") + 
  facet_wrap(~siteType, scales="free_x") +
  theme_bw()

# Save
# save(sam_info, seqtab.nochim, taxa, ps, otu_taxa_seq, goods.log, goods, file="startHere4PCoA.Rdata")

# START HERE FOR Principal coordinate analysis ----
library(vegan)
library(MCMC.OTU)
library(ggfortify)
library(cluster)
library(labdsv)
library(tidyverse)

# Load in data
load("startHere4PCoA.Rdata")
# DIFFERENCE BY SITE/STRI TIME

# Read in data 
alldat <- as.data.frame(seqtab.nochim)
summary(alldat)[,1:4]

# names are OTUs
names(alldat)
str(alldat)

alldat$sample <- as.factor(row.names(alldat))
summary(alldat)[,(ncol(alldat)-4):ncol(alldat)]

# subset alldat for samples in "mid" only----
summary(sam_info)

Mid_samples <- sam_info %>%
  filter(Time=="Mid")
summary(Mid_samples)

alldat.Mid <- alldat %>%
  filter(sample %in% Mid_samples$SampleID)

# purging under-sequenced samples; 
# and OTUs represented in less than 3% of all samples

goods <- purgeOutliers(alldat.Mid,
                       count.columns = c(1:ncol(alldat)-1),
                       otu.cut = 0.00001,
                       zero.cut = 0.01)

summary(goods)[,1:6]
summary(goods)[,495:499]

# creating a log-transfromed normalized dataset for PCoA:
goods.log <- logLin(data = goods,
                    count.columns = 2:length(names(goods)))
summary(goods.log)[,1:6]
# DIFFERENCE BY SITE (REMOVE STRI TIMEPOINTS) -------

# Read in data 
alldat <- as.data.frame(seqtab.nochim)
summary(alldat)[,1:4]

# names are OTUs
names(alldat)
str(alldat)

alldat$sample <- as.factor(row.names(alldat))
summary(alldat)[,(ncol(alldat)-4):ncol(alldat)]

# subset alldat for samples in "mid" only (no early or late day collections because they are all from STRI) 
summary(sam_info)

Mid_samples <- sam_info %>%
  filter(Time=="Mid")
summary(Mid_samples)

alldat.Mid <- alldat %>%
  filter(sample %in% Mid_samples$SampleID)

# purging under-sequenced samples; 
# and OTUs represented in less than 3% of all samples

goods <- purgeOutliers(alldat.Mid,
                        count.columns = c(1:ncol(alldat)-1),
                        otu.cut = 0.00001,
                        zero.cut = 0.01)

summary(goods)[,1:6]
summary(goods)[,495:499]

# creating a log-transfromed normalized dataset for PCoA:
goods.log <- logLin(data = goods,
                    count.columns = 2:length(names(goods)))
summary(goods.log)[,1:6]

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
goods.dist <- vegdist(goods.log, method = "manhattan")
goods.pcoa <- pcoa(goods.dist)

# make conditions
table(sam_info$SampleID %in% goods$cdat)

conditions <- sam_info %>%
  filter(SampleID %in% goods$cdat)

table(conditions$SampleID %in% goods$cdat)
head(conditions)
summary(conditions)

#make conditions with mean max time for sites we have----
tempconditions <- merge(conditions, meansbysite, by="Site", all=T)
summary(tempconditions$Site)

# plotting by type
scores <- goods.pcoa$vectors
margin <- 0.01


# play around with these numbers
xaxis <- 1
yaxis <- 2

plot(scores[,xaxis], scores[,2],type="n",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("Axis", xaxis,"(", round(goods.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("Axis", yaxis,"(", round(goods.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     main="PCoA by STRI Time Points") +
#STRIPoint by time of day 
 # points(scores[conditions$Time=="Early",xaxis],scores[conditions$Time=="Early",yaxis], col="red", pch=19) +
 # points(scores[conditions$Time=="Mid",xaxis],scores[conditions$Time=="Mid",yaxis], col="purple", pch=17) +
 # points(scores[conditions$Time=="Late",xaxis],scores[conditions$Time=="Late",yaxis], col="green", pch=15) 
  
  
# inshore sites
points(scores[conditions$Site=="PuntaDonato",xaxis],scores[conditions$Site=="PuntaDonato",yaxis], col="salmon", pch=19) +
points(scores[conditions$Site=="STRIPoint",xaxis],scores[conditions$Site=="STRIPoint",yaxis], col="salmon", pch=17) +  
points(scores[conditions$Site=="Cristobal",xaxis],scores[conditions$Site=="Cristobal",yaxis], col="salmon", pch=15) +
points(scores[conditions$Site=="PuntaLaurel",xaxis],scores[conditions$Site=="PuntaLaurel",yaxis], col="salmon", pch=18) 
# offshore sites
points(scores[conditions$Site=="DragoMar",xaxis],scores[conditions$Site=="DragoMar",yaxis], col="royalblue4", pch=1) +
points(scores[conditions$Site=="BastimentosN",xaxis],scores[conditions$Site=="BastimentosN",yaxis], col="royalblue4", pch=2) +
points(scores[conditions$Site=="BastimentosS",xaxis],scores[conditions$Site=="BastimentosS",yaxis], col="royalblue4", pch=0) +
points(scores[conditions$Site=="PopaIsland",xaxis],scores[conditions$Site=="PopaIsland",yaxis], col="royalblue4", pch=5)
  ordihull(scores,conditions$Time,label=T, draw = "polygon", col = c("royalblue4", "salmon", alpha = 255))
  legend("bottomright", c("PuntaDonato","DragoMar","STRIPoint","BastimentosN","Cristobal","BastimentosS","PuntaLaurel","PopaIsland"), pch=c(19,1,17,2,15,0,18,5), col=c("salmon","royalblue4"), cex=0.5, bty = "n")
  
# MCMC OTU analysis on site type midday only -----
  
# reformat data for mcmc.otu
goods.mcmc <- merge(goods,conditions, by = 1)
names(goods.mcmc)[c(1:3)]
names(goods.mcmc)[c(499:503)]
names(goods.mcmc)[1] <- "sample"
goods.mcmc_taxa <- merge(goods.mcmc,otu_taxa, by= "ids")
# stacking the data table
gs <- otuStack(goods.mcmc,
    count.columns=c(2:(ncol(goods.mcmc)-4)),
    condition.columns = c(1,500:503))
head(gs)

# fitting the model
mm <- mcmc.otu(
  fixed = "siteType",
  data = gs
)

summary(mm)

# selecting the OTUs that were modeled reliably
acpass <- otuByAutocorr(mm,gs)
head(acpass)

# calculating effect sizes and p-values:
ss <- OTUsummary(mm,gs,summ.plot=FALSE)

# correcting for mutliple comparisons (FDR)
ss <- padjustOTU(ss)

# getting significatly changing OTUs (FDR<0.05)
sigs <- signifOTU(ss, p.cutoff = 0.015)
head(sigs)
sigs
# plotting them
quartz()
ss2 <- OTUsummary(mm,gs_order,otus=sigs)
#wont let me edit the plot----
ss2 + labs(x= "sitetype", y="propotion") 
ss2
# bar-whiskers graph of relative changes:
gs_order <- merge(gs, otu_taxa, by= 2)
ssr=OTUsummary(mm,gs,otus=signifOTU(ss),relative=TRUE)
# displaying effect sizes and p-values for significant OTUs
ss$otuWise[sigs]


#capscale ----
# save(sam_info, goods, goods.log, sigs, ss, file="mcmcanalysis.Rdata")
load("mcmcanalysis.Rdata")
# specify groups
head(sam_info)
table(sam_info$SampleID %in% goods$cdat)

conditions <- sam_info %>%
  filter(SampleID %in% goods$cdat)
row.names(conditions) <- conditions$SampleID
head(conditions)

siteType <- as.vector(conditions$siteType)
site <- as.vector(conditions$Site)
time <- as.vector(conditions$Time)

siteShape <- ifelse(site=="Cristobal", 15, ifelse(site=="DragoMar", 1, ifelse(site=="PopaIsland", 17, ifelse(site=="PuntaDonato", 19, 5))))
siteTypeColor <- ifelse(siteType=="inshore","salmon", "royalblue4")

metaData <- data.frame(cbind(siteType, site, time))
head(metaData)

# constrained by a model:
cmd <- capscale(goods.log~siteType+site,metaData,distance="manhattan")
cmd
head(cmd)
summary(cmd)
# I'm not sure this is the correct way to calculate this...
eig <- eigenvals(cmd)[c(1:7)]
eig
perc.var <- eig/sum(eig)

cmd.scores <- scores(cmd)
cmd.scores

axes2plot <- c(1,2)  # try 2,3 too

head(cmd)
plot(cmd, choices = axes2plot, display = "sites",
     xlab = paste("CAP1 ", round(perc.var[1],2)*100, "% variance explained",sep=""), 
     ylab = paste("CAP2 ", round(perc.var[2],2)*100, "% variance explained",sep=""),
     type="points")
points(cmd, display="sites",pch=siteShape, col=siteTypeColor)
ordispider(cmd)

# inshore sites
  points(cmd.scores[metaData$Site=="PuntaDonato"],cmd.scores[metaData$Site=="PuntaDonato"], col="salmon", pch=19) 
  points(scores[conditions$Site=="STRIPoint"],scores[conditions$Site=="STRIPoint",yaxis], col="salmon", pch=17) +  
  points(scores[conditions$Site=="Cristobal"],scores[conditions$Site=="Cristobal",yaxis], col="salmon", pch=15) +
  points(scores[conditions$Site=="PuntaLaurel"],scores[conditions$Site=="PuntaLaurel",yaxis], col="salmon", pch=18) 
# offshore sites
  points(scores[conditions$Site=="DragoMar",xaxis],scores[conditions$Site=="DragoMar",yaxis], col="royalblue4", pch=1) +
  points(scores[conditions$Site=="BastimentosN",xaxis],scores[conditions$Site=="BastimentosN",yaxis], col="royalblue4", pch=2) +
  points(scores[conditions$Site=="BastimentosS",xaxis],scores[conditions$Site=="BastimentosS",yaxis], col="royalblue4", pch=0) +
  points(scores[conditions$Site=="PopaIsland",xaxis],scores[conditions$Site=="PopaIsland",yaxis], col="royalblue4", pch=5)


ordispider(cmd,choices = axes2plot, groups = metaData$site)
ordiellipse(cmd, choices = axes2plot, groups = metaData$siteType, draw="polygon", col=c("salmon","royalblue4"),label=F)


anova(cmd)
step(cmd)

# about ____% of variation is due to constraints (i.e. model factors)

adonis(goods.log~siteType+site,metaData,distance="manhattan")

#save(conditions, goods.log, otu_taxa, file="heatmap.RData")
load("heatmap.Rdata")
#heat map atempt----
# install.packages("pheatmap")
library("pheatmap")

#top 22 significant: "OTU54"  "OTU133" "OTU134" "OTU141" "OTU219" "OTU237" "OTU249" "OTU370" 
# "OTU409" "OTU413" "OTU428" "OTU439" "OTU446" "OTU449" "OTU454" "OTU466" "OTU499" "OTU508" 
#"OTU516" "OTU540" "OTU562" "OTU616"
goodssignam <- goods.log %>% select("OTU54", "OTU133", "OTU134", "OTU141", "OTU219", "OTU237", "OTU249", "OTU370", 
                                  "OTU409", "OTU413", "OTU428", "OTU439", "OTU446", "OTU449", "OTU454", "OTU466", 
                                  "OTU499", "OTU508", "OTU516", "OTU540", "OTU562", "OTU616")
row.names(otu_taxa) <- otu_taxa$ids
goodssignam <- data.frame(t(goodssignam[-1]))
named <- merge(goodssignam, otu_taxa, by= 0)
row.names(named) <- named$Row.names
named <- named %>% mutate(id = paste(Order, ids, sep="_"))
row.names(named) <- named$id
drops <- c("Order","id", "ids", "Row.names")
named <- named[, !(names(named) %in% drops)]

heat.colors = colorRampPalette(rev(c("salmon","white","royalblue4")),bias=1)(100)
pheatmap(as.matrix(named),
         scale= "column", cex=0.9,
         color = heat.colors, border_color=NA,
         clustering_distance_rows="correlation", cluster_cols=T)
#add sample names at the bottom back in 
# DESeq for Stats-------
# load deseq

library("DESeq")
library(genefilter)

# load in OTU table
otu<-read.csv("Sep21_OutputDADA_AllOTUs_FocusYesOnly.csv")
head(otu)
rownames(otu)<-otu$X #make sample IDs the rownames
length(names(otu))

# must swap rows/columns to match DESeq format

counts<-data.frame(t(otu[,2:90])) # only the columns with count data

# must remove sample with zero counts -> 2015M11r11
counts<-counts[,c(1:3,5:110)]

# Creating table of conditions for your experiment, 

Year=Type=Family=c(1:length(names(counts)))
Year[grep("2015",names(counts))]="yr2015"
Year[grep("2016",names(counts))]="yr2016" 
Type[grep("r",names(counts))]="larvae"
Type[grep("r",names(counts),invert=TRUE)]="adult"
m7<-c("M7",".7")
m9<-c("M9",".9")
m11<-c("M11",".11")
Family[grep(paste(m7,collapse="|"),names(counts))]="M7"
Family[grep(paste(m9,collapse="|"),names(counts))]="M9"
Family[grep(paste(m11,collapse="|"),names(counts))]="M11"
Family[grep("24",names(counts))]="M24"

conditions=data.frame(cbind(Year,Type,Family))
head(conditions)

real=newCountDataSet(counts,conditions) 
real=estimateSizeFactors(real)

sizeFactors(real)

# # ####all the data you ever wanted about quality control - how to find outlier samples - not sure that they're worth removing; it's the two #7 parents...
library(arrayQualityMetrics)

cds=estimateDispersions(real,method="blind")
vsdBlind=varianceStabilizingTransformation(cds)

getwd()

 v="/Users/drcarl/Dropbox/AIMSpostdoc/KateMontiSpawnExpt/deseqAQM1"

arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("Year"),force=TRUE) #check .html output file in new folder

# ###############Remove outliers as detected; repeat arrayQualityMetrics above after regenerating newCountDataSet
# ###############to confirm that all outliers have been removed 

#build the DESeq object
real=estimateDispersions(real,sharingMode="gene-est-only")  #using a pooled dispersion estimate, parametric fit
#can use gene-est-only for OTU data, should have >7 reps per factor level

goods=t(counts(real,normalized=TRUE))

# what is the proportion of samples with data for these OTUs?
withData=apply(goods,2,function(x){sum(x>0)/length(x)})

hist(withData) #Most OTUs not well represented across samples; there are <10 with counts in more than 60% of samples

# what percentage of total counts does each OTU represent?
props=apply(goods,2,function(x){sum(x)/sum(goods)})

barplot(props)

props: sq1=66%; sq2=16.8%; sq3=5.1%; sq4=4.9%; sq5=1.3%; rest<1% 

# ######################Determining quality filtering cutoffs

fit0=fitNbinomGLMs(real, count ~ 1,glmControl=list(maxit=200)) # null model: expression does not depend on anything
fit1=fitNbinomGLMs(real, count ~ Family,glmControl=list(maxit=200))
fit2=fitNbinomGLMs(real, count ~ Type,glmControl=list(maxit=200))
fit3=fitNbinomGLMs(real, count ~ Year,glmControl=list(maxit=200))


pvals.f<-nbinomGLMTest(fit1,fit0) #expression just due to family
pvals.t<-nbinomGLMTest(fit2,fit0) #expression just due to type
pvals.y<-nbinomGLMTest(fit3,fit0) #expression just due to year

par(mfrow=c(1,3))

pvalue=pvals.f; title="pvals.f" #change to each type of pval: f, t and y
theta=seq(from=0,to=0.8,by=0.02)

filterChoices=data.frame(`mean`=rowMeans(counts(real)),`median`=apply((counts(real)),1,median),`min`=rowMin(counts(real)),`max`=rowMax(counts(real)),`sd`=rowSds(counts(real)))
rejChoices=sapply(filterChoices,function(f) filtered_R(alpha=0.1,filter=f,test=pvalue,theta=theta,method="BH"))
library("RColorBrewer")
myColours=brewer.pal(ncol(filterChoices),"Set1")

matplot(theta,rejChoices,type="l",lty=1,col=myColours,lwd=2,xlab=expression(theta),ylab="number of rejections",main=title)
legend("bottomleft",legend=colnames(filterChoices),fill=myColours)

# #look for peak in graph - corresponds to correct theta and best-fit line for which metric to use 
#looks like 0.35 will be fine for all

# #######################Quality Filtering Data based on theta - get rid of genes with low variance

# #Toss OTUs where variance is too low to matter
rs=rowMeans(counts(real)) #using mean as quality filtering metric based on analyses above
theta=0.35 
use=(rs>quantile(rs,probs=theta)) ###
table(use) 
# use
# FALSE  TRUE 
   # 31    58 

realFilt=real[use,]
vsd=getVarianceStabilizedData(realFilt) #variance stabilized counts => good for 

normal=counts(realFilt,normalized=TRUE) #just extracting normalized counts


######################## Now for the real Model Testing - account for family; but interested in type*year effects; allow 30 iterations

fit0=fitNbinomGLMs(realFilt, count ~ Family, glmControl=list(maxit=30)) 
fit1=fitNbinomGLMs(realFilt, count ~ Family+Type, glmControl=list(maxit=30))
fit2=fitNbinomGLMs(realFilt, count ~ Family+Year, glmControl=list(maxit=30))

fit3=fitNbinomGLMs(realFilt, count ~ Family+Type+Year, glmControl=list(maxit=30))
fit4=fitNbinomGLMs(realFilt, count ~ Family+Type*Year, glmControl=list(maxit=30))


# testing section

pvals.y<-nbinomGLMTest(fit2,fit0) #testing significance of year term
pvals.t<-nbinomGLMTest(fit1,fit0) #testing significance of type term
pvals.i<-nbinomGLMTest(fit4,fit3) #testing significance of interaction term


#adjusting for multiple testing AND
#making non-convergent model p-values NA's -NOTE: both models must have converged 
adjp.t<-p.adjust(pvals.t,method="BH")
adjp.t=data.frame(adjp.t)
pvals.t=data.frame(pvals.t)
converged<-as.data.frame(cbind(fit0$converged,fit1$converged))
converged$test<-converged$V1+converged$V2
rownames(fit1)->rownames(pvals.t); rownames(fit1)->rownames(converged);rownames(fit1)->rownames(adjp.t);
converged$test<-apply(converged,1,all)
badmods<-subset(converged,(!converged$test))
for ( i in rownames(badmods)){pvals.t[i,1]<-NA}
for ( i in rownames(badmods)){adjp.t[i,1]<-NA}

adjp.y<-p.adjust(pvals.y,method="BH")
adjp.y=data.frame(adjp.y)
pvals.y=data.frame(pvals.y)
converged<-as.data.frame(cbind(fit0$converged,fit2$converged))
converged$test<-converged$V1+converged$V2
rownames(fit2)->rownames(pvals.y); rownames(fit2)->rownames(converged);rownames(fit2)->rownames(adjp.y);
converged$test<-apply(converged,1,all)
badmods<-subset(converged,(!converged$test))
for ( i in rownames(badmods)){pvals.y[i,1]<-NA}
for ( i in rownames(badmods)){adjp.y[i,1]<-NA}


adjp.i<-p.adjust(pvals.i,method="BH")
adjp.i=data.frame(adjp.i)
pvals.i=data.frame(pvals.i)
converged<-as.data.frame(cbind(fit4$converged,fit3$converged))
rownames(fit4)->rownames(pvals.i); rownames(fit4)->rownames(converged);rownames(fit4)->rownames(adjp.i);
converged$test<-apply(converged,1,all)
badmods<-subset(converged,(!converged$test))
for ( i in rownames(badmods)){pvals.i[i,1]<-NA}
for ( i in rownames(badmods)){adjp.i[i,1]<-NA}

summary(adjp.t)
summary(adjp.y)
summary(adjp.i)


#creating table of all multiple test corrected p-values with variance stabilized count data 
PPV_VSD<-cbind(vsd, "adjp.t" = adjp.t$adjp.t, "adjp.y" = adjp.y$adjp.y,"adjp.i" = adjp.i$adjp.i,"pval.t" = pvals.t$pvals.t, "pval.y" = pvals.y$pvals.y, "pval.i" = pvals.i$pvals.i)  

PPV_Norm<-cbind(normal, "adjp.t" = adjp.t$adjp.t, "adjp.y" = adjp.y$adjp.y,"adjp.i" = adjp.i$adjp.i,"pval.t" = pvals.t$pvals.t, "pval.y" = pvals.y$pvals.y, "pval.i" = pvals.i$pvals.i)  


write.csv(PPV_VSD, file="VSDandPVALS_FocusOnly_MeansTheta035_sep26.csv", quote=F) #writing an output file of vsd plus p-values
write.csv(PPV_Norm, file="NormCtsandPVALS_FocusOnly_MeansTheta035_sep26.csv", quote=F) #

########################## counting, venn diagram:

p<-data.frame(PPV_VSD) #OR
p<-read.csv("VSDandPVALS_FocusOnly_MeansTheta035_sep26.csv"); rownames(p)<-p$X

inter=row.names(p[p$adjp.i<=0.05 & !is.na(p$adjp.i),])#all rows where adjusted pvalue is less than or equal to 0.05 and is not an NA
year=row.names(p[p$adjp.y<=0.05 & !is.na(p$adjp.y),])
type=row.names(p[p$adjp.t<=0.05 & !is.na(p$adjp.t),])


candidates=list("Age"=type,"Year"=year,"Age*Year"=inter)
library(gplots)
quartz()
venn(candidates)
inter
year
type


#


###############################################
##### plotting  #######
###############################################


#How related are sig OTUs look like?








goods<-t(data.frame(PPV_Norm[,1:109]))




###Normalize reads for plotting
#first sum each row in OTU columns and obtain mean sum, then multiply each value in matrix by mean sum

head(goods)
names(goods)
goods$sum=apply(goods[,c(6:length(goods[1,]))], 1, function(x) {sum(x, na.rm = T)})
goods$msum = mean(goods$sum)/goods$sum

goodsNorm<-cbind(goods[,c(1:5)],((goods[,c(6:28)]*goods$msum)))

goodsNorm$NewSum=apply(goodsNorm[,c(6:28)], 1, function(x) {sum(x, na.rm = T)})
head(goodsNorm)

#NOW convert normalized reads to percents (note, all samples rescaled to have equal total reads - 18751.11)
goodsNormPerc<-cbind(goods[,c(1:5)],((goods[,c(6:28)]*goods$msum)/18751.11))
apply(goodsNormPerc[,c(6:28)], 1, function(x) {sum(x, na.rm = T)}) #test that all are now 100% = 1

######Now plot normalized read proportions
gss=otuStack(goodsNormPerc,count.columns=c(6:length(goodsNormPerc[1,])),condition.columns=c(1:5))[1:3795,] #must remove 'summ' 

gss$otu=factor(gss$otu,levels=
c("sq1","sq2","sq3","sq4","sq5","sq6","sq7","sq8","sq9","sq10","sq11","sq12","sq13","sq14","sq15","sq16","sq17","sq19","sq20","sq21","sq22","sq23", #C15
               "sq18"))#, #C15.6
               #"sq16")) #D1
               
focus<-subset(gss,Focus=="yes")

#get gradient of colors for the C15 OTUs
colfunc <- colorRampPalette(c("lightblue", "navy"))
colfunc(22)
plot(rep(1,22),col=colfunc(22),pch=19,cex=3)


p<-ggplot(focus, aes(x=sample, y=count, fill = otu))+geom_bar(position = "stack",stat="identity")+theme_bw()+scale_fill_manual(values=c(colfunc(22), rep("cyan3",1))) +facet_wrap(~ColonyID+Year,scales="free_x")
p

ggplot_build(p)$data

#Now for some stats
# stacking the data; adjust count.columns and condition.columns values for your data

gss=otuStack(goods,count.columns=c(5:length(goods[1,])),condition.columns=c(1:4))

mm=mcmc.otu(fixed="Type+Colony",
            data=gss,
            nitt=55000,thin=50,burnin=5000) #long chain to improve modeling of rare variants)


            
#selecting hte OTUs that were modeled reliably (OTUs too rare for confident parameter estimates are discarded)
acpass=otuByAutocorr(mm,gss,ac.cut=0.1)
acpass

#verifying weird OTUs
plot(mm)

#calculating differences and p-values between all pairs of factor combinations
smm0=OTUsummary(mm,gss,otus=acpass,summ.plot=FALSE)

#adjusting p-values for multiple comparisons
smmA=padjustOTU(smm0)

#significant OTUs at FDR<0.05
sigs=signifOTU(smm0,p.cutoff=0.1)
sigs

#plotting the significant ones
smm1=OTUsummary(mm,gss,otus=sigs)