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
library(dada2) # for dada2
library(ShortRead) # for dada2
library(phyloseq) # for phyloseq
library(tidyverse) # for data wrangling and ggplot2
library(vegan) # for adonis
library(MCMC.OTU) # for MCMC.oTU
library(cluster)
library(labdsv)
library(pheatmap) # for pretty heatmaps

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

# Make OTU - sequence - taxa table for later
colnames(seqtab.nochim)[c(1:3)]
colnames(taxa)[c(1:10)]
rownames(taxa)[c(1)]

seqtab.trans <- as.data.frame(t(seqtab.nochim))
head(seqtab.trans)

# create short OTU IDs
seqtab.trans$ids <- paste0("OTU", seq(1, length(colnames(seqtab.nochim))))
ids <-paste0("OTU", seq(1, length(colnames(seqtab.nochim))))
head(seqtab.trans)

# merge with taxa
otu_taxa_seq <- merge(seqtab.trans, taxa, by = 0)
# coorelate otu with taxa (removing seq info)
otu_taxa <- select(otu_taxa_seq, "ids", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# make sure it worked 
head(otu_taxa)

# relevel sam_info so that the order is Early -> Mid -> Late, instead of alphabetical (Early, Late, Mid)
head(sam_info)
levels(sam_info$Time)
sam_info$Time <- factor(sam_info$Time, levels = c("Early", "Mid", "Late"))
levels(sam_info$Time)

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

# Separate data and conds objects into two groups: STRI only (AM+midday+PM) and All sites (Midday only)
head(sam_info)

# Subset the `ps` object into STRI only (all timepoints) and all site (midday only)
psSTRI <- subset_samples(ps, Site=="STRIPoint")
psMid <- subset_samples(ps, Time=="Mid")

# Save
save(sam_info, seqtab.nochim, taxa, ps, psSTRI, psMid, otu_taxa, otu_taxa_seq, file="startHere4vegan.Rdata")

# START HERE FOR VEGAN ---------


# Load in data
load("startHere4vegan.Rdata")

# Sum technical replicates --- added 5 August 2018 ------
# Tidy up "sam_info
sam_info <- sam_info %>% 
  mutate(site = ifelse(Site=="STRIPoint", "STRI_in", 
                       ifelse(Site=="BastimentosN", "BasN_off", 
                              ifelse(Site=="BastimentosS", "BasS_off", 
                                     ifelse(Site=="Cristobal", "Cris_in", 
                                            ifelse(Site=="DragoMar", "Drag_off", 
                                                   ifelse(Site=="PopaIsland","Popa_off",
                                                          ifelse(Site=="PuntaLaurel","Laur_in", "Dona_in"))))))), 
         tow = substr(SampleID,3,3),
         samID = as.factor(paste(Time,tow,site,sep="_")))

head(sam_info)
summary(sam_info)

old2new <- sam_info %>% 
  select(SampleID, samID) %>%
  column_to_rownames("SampleID")
head(old2new)

# sum tech reps
alldat0 <- as.data.frame(seqtab.nochim)
head(alldat0)[c(1)]
alldat1 <- merge(alldat0,old2new,by=0)
head(alldat1[c(1,1528)])

alldat <- alldat1 %>% 
  group_by(samID) %>% # group by sample ID (will group technical replicates)
  mutate_all(function(x) as.numeric(as.character(x))) %>% # convert to numeric so "sum" will work
  summarise_all(funs(sum)) %>% # add up the counts
  select(-Row.names) # it makes a column called "Row.names", remove it
head(alldat)[1:10]

# manual spot check to make sure that worked
sum(alldat1[alldat1$samID=="Mid_1_STRI_in",][2])==alldat[alldat$samID=="Mid_1_STRI_in",][2]
sum(alldat1[alldat1$samID=="Early_6_STRI_in",][12])==alldat[alldat$samID=="Early_6_STRI_in",][12]

# subset alldat for samples in "mid" only(run to make goods only Mid time pt)----
summary(sam_info)

Mid_samples <- sam_info %>%
  filter(Time=="Mid")
summary(Mid_samples)
head(Mid_samples)

alldat.Mid <- alldat %>%
  filter(samID %in% Mid_samples$samID) %>%
  rename(sample = samID) # MCMC.OTU needs this column to be named "sample"
head(alldat.Mid[c(1:3)])

# subset alldat for samples in "STRI" only(run to make goods only STRI)----
STRI_samples <- sam_info %>%
  filter(Site=="STRIPoint")
summary(STRI_samples)

alldat.STRI <- alldat %>%
  filter(samID %in% STRI_samples$samID) %>%
  rename(sample = samID) # MCMC.OTU needs this column to be named "sample"
head(alldat.STRI[c(1:3)])
table(sapply(alldat.STRI,is.numeric))

# purging under-sequenced samples; 
# change what is being purged in alldat. depending on site or time 

goods.STRI <- purgeOutliers(alldat.STRI,
                       count.columns = c(2:ncol(alldat.STRI)),
                       sampleZcut = (-2.5),
                       otu.cut = 0.00001,
                       zero.cut = 0.01)
summary(goods.STRI)[,1:6]

goods.Mid <- purgeOutliers(alldat.Mid,
                            count.columns = c(2:ncol(alldat.Mid)),
                            otu.cut = 0.00001,
                            zero.cut = 0.01)
summary(goods.Mid)[,1:6]

# creating a log-transfromed normalized dataset for PCoA:
goods.log.STRI <- logLin(data = goods.STRI,
                    count.columns = 2:length(names(goods.STRI)))
summary(goods.log.STRI)[,1:6]

goods.log.Mid <- logLin(data = goods.Mid,
                         count.columns = 2:length(names(goods.Mid)))
summary(goods.log.Mid)[,1:6]

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
goods.dist.STRI <- vegdist(goods.log.STRI, method = "manhattan")
goods.pcoa.STRI <- pcoa(goods.dist.STRI)

goods.dist.Mid <- vegdist(goods.log.Mid, method = "manhattan")
goods.pcoa.Mid <- pcoa(goods.dist.Mid)

# make conditions
table(sam_info$samID %in% goods.STRI$sample)

conditions.STRI <- sam_info %>%
  select(Site,Time,siteType,site,samID) %>%
  filter(samID %in% goods.STRI$sample) %>%
  distinct()
head(conditions.STRI)

conditions.Mid <- sam_info %>%
  select(Site,Time,siteType,site,samID) %>%
  filter(samID %in% goods.Mid$sample) %>%
  distinct()
head(conditions.Mid)

table(conditions.STRI$samID %in% goods.STRI$sample)
table(conditions.Mid$samID %in% goods.Mid$sample)

#make conditions with mean max time for sites we have----
# skipping for now...
# tempconditions <- merge(conditions, meansbysite, by="Site", all=T)
# summary(tempconditions$Site)

# plotting by type
scores.STRI <- goods.pcoa.STRI$vectors
scores.Mid <- goods.pcoa.Mid$vectors
margin <- 0.01

# play around with these numbers
xaxis <- 1
yaxis <- 2

plot(scores.Mid[,xaxis], scores.Mid[,2],type="n",
     xlim=c(min(scores.Mid[,xaxis])-margin,max(scores.Mid[,xaxis])+margin),
     ylim=c(min(scores.Mid[,2])-margin,max(scores.Mid[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("Axis", xaxis,"(", round(goods.pcoa.Mid$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("Axis", yaxis,"(", round(goods.pcoa.Mid$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     main="PCoA by Site Ttype (Midday Only)") +
  ordihull(scores.Mid,conditions.Mid$siteType,label=F, draw = "polygon", col = c("salmon", "royalblue4", alpha = 255))

  # inshore sites
  points(scores.Mid[conditions.Mid$Site=="PuntaDonato",xaxis],scores.Mid[conditions.Mid$Site=="PuntaDonato",yaxis], col="salmon", pch=19) +
  points(scores.Mid[conditions.Mid$Site=="STRIPoint",xaxis],scores.Mid[conditions.Mid$Site=="STRIPoint",yaxis], col="salmon", pch=17) +  
  points(scores.Mid[conditions.Mid$Site=="Cristobal",xaxis],scores.Mid[conditions.Mid$Site=="Cristobal",yaxis], col="salmon", pch=15) +
  points(scores.Mid[conditions.Mid$Site=="PuntaLaurel",xaxis],scores.Mid[conditions.Mid$Site=="PuntaLaurel",yaxis], col="salmon", pch=18) 

  # offshore sites
  points(scores.Mid[conditions.Mid$Site=="DragoMar",xaxis],scores.Mid[conditions.Mid$Site=="DragoMar",yaxis], col="royalblue4", pch=1) +
  points(scores.Mid[conditions.Mid$Site=="BastimentosN",xaxis],scores.Mid[conditions.Mid$Site=="BastimentosN",yaxis], col="royalblue4", pch=2) +
  points(scores.Mid[conditions.Mid$Site=="BastimentosS",xaxis],scores.Mid[conditions.Mid$Site=="BastimentosS",yaxis], col="royalblue4", pch=0) +
  points(scores.Mid[conditions.Mid$Site=="PopaIsland",xaxis],scores.Mid[conditions.Mid$Site=="PopaIsland",yaxis], col="royalblue4", pch=5)
  
  legend(100,50, c("PuntaDonato","DragoMar","STRIPoint","BastimentosN","Cristobal","BastimentosS","PuntaLaurel","PopaIsland"), pch=c(19,1,17,2,15,0,18,5), col=c("salmon","royalblue4"), cex=0.5, bty = "n")
  
  
plot(scores.STRI[,xaxis], scores.STRI[,2],type="n",
       xlim=c(min(scores.STRI[,xaxis])-margin,max(scores.STRI[,xaxis])+margin),
       ylim=c(min(scores.STRI[,2])-margin,max(scores.STRI[,2])+margin),
       mgp=c(2.3,1,0),
       xlab=paste("Axis", xaxis,"(", round(goods.pcoa.STRI$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
       ylab=paste("Axis", yaxis,"(", round(goods.pcoa.STRI$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
       main="PCoA by Time (STRI Only)") +
      ordihull(scores.STRI,conditions.STRI$Time,label=F, 
               draw = "polygon", col = c("gold", "orange","blue", alpha = 255))
      # ordispider(scores.STRI,conditions.STRI$Time,label=F,lwd=2,col = c("gold", "orange", "blue"))
#STRIPoint by time of day 
  points(scores.STRI[conditions.STRI$Time=="Early",xaxis],scores.STRI[conditions.STRI$Time=="Early",yaxis], col="gold", pch=19) +
  points(scores.STRI[conditions.STRI$Time=="Mid",xaxis],scores.STRI[conditions.STRI$Time=="Mid",yaxis], col="orange", pch=17) +
  points(scores.STRI[conditions.STRI$Time=="Late",xaxis],scores.STRI[conditions.STRI$Time=="Late",yaxis], col="blue", pch=15) 
  legend("bottomright", c("Early", "Mid", "Late"), lwd=3, col=c("gold","orange","blue"), bty="n")


# Shannon diversity -----
sample <- as.vector(goods.STRI$sample)
time <- as.vector(conditions.STRI$Time)
shannonSTRI <- diversity(goods.log.STRI, "shannon")
df.shannon.STRI <- data.frame(sample,time,shannonSTRI)
df.shannon.STRI$time <- factor(df.shannon.STRI$time, levels = c("Early", "Mid", "Late"))
boxplot(shannonSTRI~time,data=df.shannon.STRI)

sample <- as.vector(goods.Mid$sample)
siteType <- as.vector(conditions.Mid$siteType)
shannonMid <- diversity(goods.log.Mid, "shannon")
df.shannon.Mid <- data.frame(sample,siteType,shannonMid)
boxplot(shannonMid~siteType,data=df.shannon.Mid)

# Simpson -------
sample <- as.vector(goods.STRI$sample)
time <- as.vector(conditions.STRI$Time)
simpsonSTRI <- diversity(goods.log.STRI, "simpson")
df.simpson.STRI <- data.frame(sample,time,simpsonSTRI,3)
df.simpson.STRI$time <- factor(df.simpson.STRI$time, levels = c("Early", "Mid", "Late"))
boxplot(simpsonSTRI~time,data=df.simpson.STRI)

sample <- as.vector(goods.Mid$sample)
siteType <- as.vector(conditions.Mid$siteType)
simpsonMid <- diversity(goods.log.Mid, "simpson")
df.simpson.Mid <- data.frame(sample,siteType,simpsonMid)
boxplot(simpsonMid~siteType,data=df.simpson.Mid)

# MCMC OTU analysis on site type midday only -----
  
# reformat data for mcmc.otu
goods.mcmc.STRI <- merge(goods.STRI,conditions.STRI, by.x = "sample", by.y="samID")
names(goods.mcmc.STRI)[c(1:3)]
names(goods.mcmc.STRI)[c(528:532)]

goods.mcmc.Mid <- merge(goods.Mid,conditions.Mid, by.x = "sample", by.y="samID")
names(goods.mcmc.Mid)[c(1:3)]
names(goods.mcmc.Mid)[c(499:506)]

# can't figure this out right now... skipping it 2August2018
# head(goods.mcmc.STRI[c(1:2)])
# head(otu_taxa)
# goods.mcmc_taxa.STRI <- merge(goods.mcmc.STRI,otu_taxa, by = "ids")

# stacking the data table
gs.STRI <- otuStack(goods.mcmc.STRI,
    count.columns=c(2:(ncol(goods.mcmc.STRI)-4)),
    condition.columns = c(1,529:532))
head(gs.STRI)

gs.Mid <- otuStack(goods.mcmc.Mid,
                    count.columns=c(2:(ncol(goods.mcmc.Mid)-4)),
                    condition.columns = c(1,503:506))
head(gs.Mid)

# fitting the model
mm.STRI <- mcmc.otu(
  fixed = "Time",
  data = gs.STRI
)

summary(mm.STRI)

mm.Mid <- mcmc.otu(
  fixed = "siteType",
  data = gs.Mid
)

summary(mm.Mid)

# selecting the OTUs that were modeled reliably
acpass.STRI <- otuByAutocorr(mm.STRI,gs.STRI)
head(acpass.STRI)

acpass.Mid <- otuByAutocorr(mm.Mid,gs.Mid)
head(acpass.Mid)

# calculating effect sizes and p-values:
ss.STRI <- OTUsummary(mm.STRI,gs.STRI,summ.plot=FALSE)
ss.Mid <- OTUsummary(mm.Mid,gs.Mid,summ.plot=FALSE)

# correcting for mutliple comparisons (FDR)
ss.STRI <- padjustOTU(ss.STRI)
ss.Mid <- padjustOTU(ss.Mid)

# getting significatly changing OTUs (FDR<0.05)
sigs.STRI <- signifOTU(ss.STRI, p.cutoff = 0.01)
length(sigs.STRI) #25 when cutoff is 0.01
sigs.STRI

sigs.Mid <- signifOTU(ss.Mid, p.cutoff = 0.05)
length(sigs.Mid) #20 when cutoff is 0.05
sigs.Mid

# Filter stack for significant OTUs
sig25_STRI <- gs.STRI %>%
  filter(otu %in% sigs.STRI)

sig20_Mid <- gs.Mid %>%
  filter(otu %in% sigs.Mid)

# Add taxonomy
gs_order_STRI <- merge(sig25_STRI, otu_taxa, by.x=2, by.y=1)
head(gs_order_STRI)
gs_order_Mid <- merge(sig20_Mid, otu_taxa, by.x=2, by.y=1)
head(gs_order_Mid)

# SAVE/LOAD MCMC -----------------
save(sam_info, goods.STRI, goods.Mid, goods.log.STRI, goods.log.Mid, 
     mm.STRI, mm.Mid, 
     sigs.STRI, sigs.Mid, ss.STRI, ss.Mid,
     gs.STRI, gs.Mid,
     gs_order_STRI, gs_order_Mid,
     sig20_Mid, sig25_STRI,
     file="mcmcanalysis.Rdata")

load("mcmcanalysis.Rdata")

# CAPSCALE ANALYSIS -------------
# specify groups
conditions.STRI <- sam_info %>%
  select(Site,Time,siteType,site,samID) %>%
  filter(samID %in% goods.STRI$sample) %>%
  distinct()
head(conditions.STRI)

conditions.Mid <- sam_info %>%
  select(Site,Time,siteType,site,samID) %>%
  filter(samID %in% goods.Mid$sample) %>%
  distinct()
head(conditions.Mid)


siteType <- as.vector(conditions.Mid$siteType)
site <- as.vector(conditions.Mid$Site)
time <- as.factor(conditions.STRI$Time)

siteShape <- ifelse(site=="Cristobal", 15, ifelse(site=="DragoMar", 1, ifelse(site=="PopaIsland", 17, ifelse(site=="PuntaDonato", 19, 5))))
siteTypeColor <- ifelse(siteType=="inshore","salmon", "royalblue4")
timeColor <- ifelse(time=="Early","gold",ifelse(time=="Mid","orange","blue"))

metaData.Mid <- data.frame(cbind(siteType, site))
head(metaData.Mid)

# constrained by a model:
cmd.STRI <- capscale(goods.log.STRI~time,time,distance="manhattan")
cmd.STRI
head(cmd.STRI)
summary(cmd.STRI)

cmd.Mid <- capscale(goods.log.Mid~site+siteType,metaData.Mid,distance="manhattan")
cmd.Mid
head(cmd.Mid)
summary(cmd.Mid)

# Calculate percentage of variance explained by CAP
eigenvals(cmd.STRI) # CAP columns 1 and 2
eig.STRI <- eigenvals(cmd.STRI)[c(1:2)]
eig.STRI
perc.var.STRI <- eig.STRI/sum(eig.STRI)

cmd.scores.STRI <- scores(cmd.STRI)
cmd.scores.STRI

eigenvals(cmd.Mid) # CAP columns 1-7
eig.Mid <- eigenvals(cmd.Mid)[c(1:7)]
eig.Mid
perc.var.Mid <- eig.Mid/sum(eig.Mid)

cmd.scores.Mid <- scores(cmd.Mid)
cmd.scores.Mid

# CAPSCALE STATS
anova(cmd.STRI)
step(cmd.STRI)
# about ____% of variation is due to constraints (i.e. model factors)
time <- as.data.frame(time)
set.seed(1)
adonis.STRI <- adonis(goods.log.STRI~time,time,distance="manhattan")
adonis.STRI$aov.tab$`Pr(>F)`[1]
# p = 0.228
set.seed(1)
adonis.Mid <- adonis(goods.log.Mid~site+siteType,metaData.Mid,distance="manhattan")
adonis.Mid$aov.tab$`Pr(>F)`[1]
# p = 0.626

# Plot CAPSCALE
axes2plot <- c(1,2)  # try 2,3 too

plot(cmd.STRI, choices = axes2plot, display = "sites",
     xlab = paste("CAP1 ", round(perc.var.STRI[1],2)*100, "% variance explained",sep=""), 
     ylab = paste("CAP2 ", round(perc.var.STRI[2],2)*100, "% variance explained",sep=""),
     type="n")
# ordihull(cmd.STRI, groups=time, draw="polygon", label=F,col=c("gold","orange","blue"))
points(cmd.STRI, display="sites", pch = 16, col=timeColor)
legend("bottomleft", inset=.02, c("Early","Mid","Late"), 
       fill=c("gold","orange","blue"), cex=0.8, bty='n')
legend("bottomright", inset=.02, paste("p = ",adonis.STRI$aov.tab$`Pr(>F)`[1], sep=" "), 
       cex=0.8, bty='n')
ordispider(cmd.STRI, col=timeColor)

plot(cmd.Mid, choices = axes2plot, display = "sites",
     xlab = paste("CAP1 ", round(perc.var.Mid[1],2)*100, "% variance explained",sep=""), 
     ylab = paste("CAP2 ", round(perc.var.Mid[2],2)*100, "% variance explained",sep=""),
     type="n")
ordihull(cmd.Mid, groups=siteType, draw="polygon", label=F,col=c("salmon","royalblue4"))
points(cmd.Mid, display="sites", pch=siteShape, col=siteTypeColor)
legend("bottomleft", inset=.02, c("Inshore","Offshore"), 
       fill=c("salmon","royalblue4"), cex=0.8, bty='n')
legend("bottomright", inset=.02, paste("p = ",adonis.Mid$aov.tab$`Pr(>F)`[1], sep=" "), 
       cex=0.8, bty='n')
# ordiellipse(cmd.Mid, groups=siteType, col=c("salmon","royalblue4"))
# ordispider(cmd.Mid, col=siteTypeColor)

# Make heatmap ------------------------------
library(pheatmap)

# Midday only 
# reformat for making the heatmap
counts_Mid <- gs_order_Mid %>%
  mutate(taxa = paste(Kingdom,Phylum,Class,Order,Family,Genus,sep="_")) %>% # make full taxa name
  mutate(taxa_otu = paste(otu,taxa,sep="_")) %>% # make unique taxa names
  select(otu,count,sample,taxa_otu) %>% # only keep columns you need
  spread(sample,count) %>% # format for heatmap
  filter(otu %in% sig20_Mid$otu) %>%# filter for significant
  column_to_rownames("taxa_otu") %>% # make the unique taxa the rownames
  select(-"otu") # get rid of otu column

head(counts_Mid)
nrow(counts_Mid)

# change order for plotting
names(counts_Mid)
order <- c(grep("in",names(counts_Mid)),
           grep("off",names(counts_Mid)))
counts_Mid <- counts_Mid[c(order)]

# Make color scale
heat.colors <- colorRampPalette(rev(c("red","white","gold")),bias=1)(100)

# Plot heatmap
pheatmap(as.matrix(counts_Mid),
         scale = "row", cex=0.9,
         color = heat.colors, border_color=NA,
         clustering_distance_rows="correlation", cluster_cols=F)

# STRI only 
# reformat for making the heatmap
counts_STRI <- gs_order_STRI %>%
  mutate(taxa = paste(Kingdom,Phylum,Class,Order,Family,Genus,sep="_")) %>% # make full taxa name
  mutate(taxa_otu = paste(otu,taxa,sep="_")) %>% # make unique taxa names
  select(otu,count,sample,taxa_otu) %>% # only keep columns you need
  spread(sample,count) %>% # format for heatmap
  filter(otu %in% sig25_STRI$otu) %>%# filter for significant
  column_to_rownames("taxa_otu") %>% # make the unique taxa the rownames
  select(-"otu") # get rid of otu column

head(counts_STRI)
nrow(counts_STRI)

# change order for plotting
names(counts_STRI)
order <- c(grep("Early",names(counts_STRI)),
           grep("Mid",names(counts_STRI)),
           grep("Late",names(counts_STRI)))
counts_STRI <- counts_STRI[c(order)]

# Make color scale
heat.colors <- colorRampPalette(rev(c("red","white","gold")),bias=1)(100)

# Plot heatmap
pheatmap(as.matrix(counts_STRI),
         scale = "row", cex=0.9,
         color = heat.colors, border_color=NA,
         clustering_distance_rows="correlation", cluster_cols=F)



