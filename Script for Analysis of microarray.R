### Annotation and normalization from Adam using custom CDF files from ProbeMatchDB: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/ensg.asp
biocLite(c("affy", "limma", "gcrma"))
biocLite("org.Mm.eg.db");biocLite("annotate")

library(annotate);library(org.Mm.eg.db)
library(affy); library(limma); library(gcrma)
# targets <- readTargets("target.txt") target file is for analysis of oligodendrocyte microarraydata
expr <- ReadAffy(filenames = targets$FileName, cdfname = "Mouse4302_Mm_ENTREZG", verbose = TRUE) #The “cdfname” should be the “Custom CDF Name” that’s listed on the BrainArray site.
expr.norm <- rma(expr)

geneIDs <- gsub("_at", "", rownames(expr.norm)) 
geneSym <- unlist(lookUp(geneIDs, "org.Mm.eg.db", "SYMBOL"))
geneNames <- unlist(lookUp(geneIDs, "org.Mm.eg.db", "GENENAME"))


# create annotation matrix, make sure that the individual objects are in the right order!
geneAnnotation <-cbind(geneIDs, as.data.frame(geneSym[match(geneIDs, names(geneSym))]),   as.data.frame(geneNames[match(geneIDs, names(geneNames))]))
colnames(geneAnnotation) <- c("EntrezID", "Gene.Symbol","GeneDesc")

write.table(expr.norm, file="expr.norm.xls", sep="\t", col.names = NA) 
featureNames(expr.norm) <- gsub("_at", "", featureNames(expr.norm))

df<- data.frame(t(as.data.frame(expr.norm)))
df1 <- cbind(EntrezID = row.names(df), df)
rownames(df1) <- NULL
EntrezID <- gsub("X", "", df1$EntrezID)
newdf<-cbind(EntrezID,df)
rownames(newdf) <- NULL
rownames(geneAnnotation)<- NULL
annotdata<-merge(geneAnnotation,newdf, by="EntrezID")
write.table(annotdata, file="annotdata.xls", sep="\t", col.names = NA)

groups<-annotdata[,4:15]

# number of columns per group (1-3, 4-6)
n <- 3
# number of groups
n_grp <- 4 ## n_grp <- ncol(dat) / n
# column indices (one vector per group)
idx_grp <- split(seq(groups), rep(seq(n_grp), each = n))

# calculate the row means for all groups
res <- lapply(idx_grp, function(i) {
  # subset of the data frame
  tmp <- groups[i]
  # calculate row means
  rowMeans(tmp, na.rm = TRUE)
})
dat <- as.data.frame(res)
colnames(dat)<-c("ZT3", "ZT9","ZT15","ZT21") 

DATA<-cbind(annotdata[,1:3], dat)
colnames(DATA)<-c("EntrezID" ,"Gene.symbol","GeneDesc","ZT3", "ZT9","ZT15","ZT21")

#combine with previous data for cell types
final.data <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/final.data.txt")
rownames(final.data)<- NULL
mergedSW<-merge(DATA,final.data, by = "Gene.symbol", all.y=TRUE)
write.table(mergedSW, file="mergedSW.xls", sep="\t", col.names = NA)

# HGNC nomenclature - ion channels
Ion.channels.HGNC <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/Nomenclature/Ion channels HGNC.txt")
HGNC.MGI <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/Nomenclature/HGNC.to.MGI")
colnames(HGNC.MGI)<-c("HGNC.ID","mgi_id")
ion.MGI<-merge(HGNC.MGI,Ion.channels.HGNC, by="HGNC.ID", all.y=TRUE)
ion.mgi<-merge(ion.MGI,geneid, by="mgi_id")
ion.symbol<-ion.mgi[,c("Gene.symbol.y","Entrez.Gene.Name")]
colnames(ion.symbol)<-c("Gene.symbol","Entrez.Gene.Name")

# NCBI nomenclature ion
Ion.channels.NCBI <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/Nomenclature/Ion channels NCBI.txt")
colnames(Ion.channels.NCBI)<-c("EntrezID","Gene.symbol")

# Ion cell types 
astro.ion1<-unique(merge(Ion.channels.NCBI,astro_sw, by="EntrezID"))
astro.ion2<-unique(merge(ion.symbol,astro_sw, by="Gene.symbol"))

endo.ion1<-unique(merge(Ion.channels.NCBI,endo_sw, by="EntrezID"))
endo.ion2<-unique(merge(ion.symbol,endo_sw, by="Gene.symbol"))

neuron.ion1<-unique(merge(Ion.channels.NCBI,neuron_sw, by="EntrezID"))
neuron.ion2<-unique(merge(ion.symbol,neuron_sw, by="Gene.symbol"))

oligo.ion1<-unique(merge(Ion.channels.NCBI,oligo_sw, by="EntrezID"))
oligo.ion2<-unique(merge(ion.symbol,oligo_sw, by="Gene.symbol"))

# HGNC nomenclature - slc transporters
genenames <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/Nomenclature/genenames.scltransport")
slc.MGI<-merge(HGNC.MGI,genenames, by="HGNC.ID", all.y=TRUE)
slc.mgi<-merge(slc.MGI,geneid, by="mgi_id")
slc.symbol<-slc.mgi[,c("Gene.symbol","Entrez.Gene.Name")]

astro.slc<-unique(merge(slc.symbol,astro_sw, by="Gene.symbol"))

endo.slc<-unique(merge(slc.symbol,endo_sw, by="Gene.symbol"))

neuron.slc<-unique(merge(slc.symbol,neuron_sw, by="Gene.symbol"))

oligo.slc<-unique(merge(slc.symbol,oligo_sw, by="Gene.symbol"))




# Combine with cell type specific data
setwd(dir = "~/Desktop/Nadia/Sleep journal club/Spreadsheet organization//Cell types & SW")

all_astro <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/all_astro.txt")
head(all_astro)
head(DATA)
astro_sw<-merge(DATA,all_astro, by="Gene.symbol")
nrow(all_astro) #263
nrow(astro_sw) #242
write.table(astro_sw, file="astro_sw.xls", sep="\t", col.names = NA)

all_neuron <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/all_neuron.txt")
neuron_sw<-merge(DATA,all_neuron, by="Gene.symbol")
nrow(all_neuron) #392
nrow(neuron_sw) #357
write.table(neuron_sw, file="neuron_sw.xls", sep="\t", col.names = NA)

all_endo <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/all_endo.txt")
endo_sw<-merge(DATA,all_endo, by="Gene.symbol")
nrow(all_endo) #595
nrow(endo_sw) #534
write.table(endo_sw, file="endo_sw.xls", sep="\t", col.names = NA)

all_oligo <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/all_oligo.txt")
oligo_sw<-merge(DATA,all_oligo, by="Gene.symbol")
nrow(all_oligo) #184
nrow(oligo_sw) #161
write.table(oligo_sw, file="oligo_sw.xls", sep="\t", col.names = NA)

## (B1) DEG Analysis for RMA Data (All timepoints)
design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4,4))) # Creates appropriate design matrix.
colnames(design) <- c("ZT3", "ZT9", "ZT15", "ZT21") # Assigns nicer column names.
contrast.matrix <- makeContrasts(ZT3-ZT9, ZT3-ZT15, ZT3-ZT21, ZT9-ZT15, ZT9-ZT21, ZT15-ZT21, levels=design) # Creates appropriate contrast matrix for pairwise comparisons.
fit <- lmFit(expr.norm, design) # Fits a linear model for each gene based on the given series of arrays.
fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
rma_deg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
rma_deg_result <- rma_deg_result[rma_deg_result$adj.P.Val<=0.05,]
write.table(rma_deg_result, "rma_deg_result.xls", col.names = NA, quote=FALSE, sep="\t")

## (B2) DEG Analysis for RMA Data (Sleep, ZT3 and ZT9 vs wake ZT15 and ZT21)
design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4,4))) # Creates appropriate design matrix.
colnames(design) <- c("ZT3", "ZT9", "ZT15", "ZT21") # Assigns nicer column names.
contrast.matrix <- makeContrasts(ZT3-ZT15, ZT9-ZT21, ZT9-ZT15, ZT9-ZT21, levels=design) # Creates appropriate contrast matrix for pairwise comparisons.
fit <- lmFit(expr.norm, design) # Fits a linear model for each gene based on the given series of arrays.
fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
#fdr=5%
rma_swdeg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
rma_swdeg_result <- rma_swdeg_result[rma_swdeg_result$adj.P.Val<=0.05,]
write.table(rma_swdeg_result, "rma_sw_deg_result.xls", col.names = NA, quote=FALSE, sep="\t")

rma_venn <- decideTests(fit2, p.value=0.05)
pdf(file="rma_venn.pdf"); vennDiagram(rma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.

#fdr=10%
rma_swdeg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
rma_swdeg_result <- rma_swdeg_result[rma_swdeg_result$adj.P.Val<=0.10,]
write.table(rma_swdeg_result, "rma_sw_deg10_result.xls", col.names = NA, quote=FALSE, sep="\t")

rma_venn <- decideTests(fit2, p.value=0.10)
pdf(file="rma_venn.pdf"); vennDiagram(rma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.

#fdr=20%, nrow=332
rma_swdeg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
rma_swdeg_result <- rma_swdeg_result[rma_swdeg_result$adj.P.Val<=0.20,]
write.table(rma_swdeg_result, "rma_sw_deg20_result.xls", col.names = NA, quote=FALSE, sep="\t")

rma_venn <- decideTests(fit2, p.value=0.20)
pdf(file="rma_venn.pdf"); vennDiagram(rma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.


## (B2) DEG Analysis for RMA Data (Sleep, ZT3 vs wake ZT15)
design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4,4))) # Creates appropriate design matrix.
colnames(design) <- c("ZT3", "ZT9", "ZT15", "ZT21") # Assigns nicer column names.
contrast.matrix <- makeContrasts(ZT3-ZT15, levels=design) # Creates appropriate contrast matrix for pairwise comparisons.
fit <- lmFit(expr.norm, design) # Fits a linear model for each gene based on the given series of arrays.
fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
#fdr=5% = 41 genes
rma_swdeg5_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
rma_swdeg5_result <- rma_swdeg5_result[rma_swdeg5_result$adj.P.Val<=0.05,]
write.table(rma_swdeg5_result, "rma_s3w15_deg_result.xls", col.names = NA, quote=FALSE, sep="\t")

rma_venn <- decideTests(fit2, p.value=0.05)
pdf(file="rma_venn.pdf"); vennDiagram(rma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.

#fdr=10% = 124 genes
rma_swdeg10_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
rma_swdeg10_result <- rma_swdeg10_result[rma_swdeg10_result$adj.P.Val<=0.10,]
write.table(rma_swdeg10_result, "rma_s3w15_deg10_result.xls", col.names = NA, quote=FALSE, sep="\t")

rma_venn <- decideTests(fit2, p.value=0.10)
pdf(file="rma_venn.pdf"); vennDiagram(rma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.

#fdr=20% = 332 genes
rma_swdeg20_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
rma_swdeg20_result <- rma_swdeg20_result[rma_swdeg20_result$adj.P.Val<=0.20,]
write.table(rma_swdeg20_result, "rma_s3w15_deg20_result.xls", col.names = NA, quote=FALSE, sep="\t")

rma_venn <- decideTests(fit2, p.value=0.20)
pdf(file="rma_venn.pdf"); vennDiagram(rma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.


#Annotation of sw DEG RMA data 
df<- data.frame(rma_swdeg20_result)
df1 <- cbind(EntrezID = row.names(df), df)
rownames(df1) <- NULL
swdeg_annotdata<-merge(geneAnnotation,df1, by="EntrezID")

names(swdeg_annotdata)[names(swdeg_annotdata) == 'GeneSymbol'] <- 'Gene.symbol'
write.table(swdeg_annotdata, "swdeg_annotdata.xls", col.names = NA, quote=FALSE, sep="\t")
# MERGE WITH BARRES and SW DATA, genes =264 after merge
SW_DEG_AllDATA<-merge(swdeg_annotdata, mergedSW, by="EntrezID", all.y=TRUE)
write.table(SW_DEG_AllDATA, "S3W15_DEG20_AllDATA.xls", col.names = NA, quote=FALSE, sep="\t")

# NOT MERGED in SW DEG DATASET = 68 = 20%
present<-swdeg_annotdata$EntrezID %in% mergedSW$EntrezID
sum(present=="TRUE")
absent<-swdeg_annotdata$Gene.symbol[present=="FALSE"]
length(absent) 

# MERGE WITH CELL SPECIFIC DATA
#all_astro <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/all_astro.txt")
head(all_astro)
head(DATA)
astro_deg<-merge(swdeg_annotdata,all_astro, by="Gene.symbol",all.y=TRUE)
astro_s3w15<-merge(DATA,astro_deg, by="Gene.symbol",all.y=TRUE)
nrow(all_astro) #263
nrow(astro_s3w15) #242
write.table(astro_s3w15, file="astro_s3w15_deg20.xls", sep="\t", col.names = NA)

# TEST FOR # of common genes!
present<-swdeg_annotdata$Gene.symbol %in% all_astro$Gene.symbol
sum(present=="TRUE")
genes.present<-swdeg_annotdata$Gene.symbol[present=="TRUE"]
absent<-swdeg_annotdata$Gene.symbol[present=="FALSE"]

#all_neuron <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/all_neuron.txt")
neuron_deg<-merge(swdeg_annotdata,all_neuron, by="Gene.symbol",all.y=TRUE)
neuron_s3w15<-merge(DATA,neuron_deg, by="Gene.symbol",all.y=TRUE)
nrow(all_neuron) #392
nrow(neuron_s3w15) #357
write.table(neuron_s3w15, file="neuron_s3w15_deg20.xls", sep="\t", col.names = NA)

#all_endo <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/all_endo.txt")
endo_deg<-merge(swdeg_annotdata,all_endo, by="Gene.symbol",all.y=TRUE)
endo_s3w15<-merge(DATA,endo_deg, by="Gene.symbol",all.y=TRUE)
nrow(all_endo) #595
nrow(endo_s3w15) #534
write.table(endo_s3w15, file="endo_s3w15_deg20.xls", sep="\t", col.names = NA)

#all_oligo <- read.delim("~/Desktop/Nadia/Sleep journal club/Spreadsheet organization/all_oligo.txt")
oligo_deg<-merge(swdeg_annotdata,all_oligo, by="Gene.symbol",all.y=TRUE)
oligo_s3w15<-merge(DATA,oligo_deg, by="Gene.symbol",all.y=TRUE)
nrow(all_oligo) #184
nrow(oligo_s3w15) #161
write.table(oligo_s3w15, file="oligo_s3w15_deg20.xls", sep="\t", col.names = NA)


#######################
## Homework R Script ##
#######################
## How to run this script:
#   source("homework_script.R")

## Dependencies:
# 	R working directory needs to contain
# 	all *.cel files and targets.txt file
#       R packages:
source("http://bioconductor.org/biocLite.R")
biocLite(c("affy", "limma", "gcrma"))
         
         ## RMA Commands
         ## (A) Normalization: RMA
setwd("~/Desktop/Nadia/Sleep journal club/SW_microarray/Yang_microarray_analysis/")         
library(affy); library(limma); library(gcrma) # Loads required libraries.
         targets <- readTargets("target.txt") # Import targets information.
         data <- ReadAffy(filenames=targets$FileName) # Import expression raw data and stores them as AffyBatch object.
         probedata <- data.frame(ID=probeNames(data), PM=pm(data), MM=mm(data)) # Access proble level data.         
         eset_rma <- rma(data) # Normalizes the data with 'rma' function and assigns them to exprSet object.
         #exprs(eset) <- log2(exprs(eset)) # Only MAS5 stores absolute intensities. GCRMA and RMA methods store log2 intensities.
         pData(eset_rma) # Lists the analyzed file names.
         write.exprs(eset_rma, file="rma_all.xls") # Exports all affy expression values to tab delimited text file. 
         ## Create Box Plots for Raw Data and Normalized Data
         pdf(file="raw_boxplot.pdf"); boxplot(data, col="red", main="Raw Data"); dev.off() # Generates boxplot for un-normalized log intensity values.
         pdf(file="rma_boxplot.pdf"); boxplot(data.frame(exprs(eset_rma)), col="blue", main="RMA Normalized Data"); dev.off() # Generates boxplot for RMA normalized log intensity values.
         
         ## (B) DEG Analysis for RMA Data (All timepoints)
         design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4,4))) # Creates appropriate design matrix.
         colnames(design) <- c("ZT3", "ZT9", "ZT15", "ZT21") # Assigns nicer column names.
         contrast.matrix <- makeContrasts(ZT3-ZT15, levels=design) # Creates appropriate contrast matrix for pairwise comparisons.
         fit <- lmFit(eset_rma, design) # Fits a linear model for each gene based on the given series of arrays.
         fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
         fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
         rma_deg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
         rma_deg_result <- rma_deg_result[rma_deg_result$adj.P.Val<=0.10,]
         write.table(rma_deg_result, "rma_deg10_result.xls", col.names = NA, quote=FALSE, sep="\t")
         

         ## (C) Create Venn Diagram for RMA Data
         rma_venn <- decideTests(fit2, p.value=0.05)
         pdf(file="rma_venn.pdf"); vennDiagram(rma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.
         
         ## (A) Normalization: GCRMA
         eset_gcrma <- gcrma(data) # Normalizes the data with 'gcrma' function and assigns them to exprSet object.
          exprs(eset) <- log2(exprs(eset)) # Only MAS5 stores absolute intensities. GCRMA and RMA methods store log2 intensities.
          pData(eset_gcrma) # Lists the analyzed file names.
          write.exprs(eset_gcrma, file="gcrma_all.xls") # Exports all affy expression values to tab delimited text file. 
          ## Create Box Plots for Raw Data and Normalized Data
          pdf(file="raw_boxplot.pdf"); boxplot(data, col="red", main="Raw Data"); dev.off() # Generates boxplot for un-normalized log intensity values.
          pdf(file="gcrma_boxplot.pdf"); boxplot(data.frame(exprs(eset_gcrma)), col="blue", main="gcRMA Normalized Data"); dev.off() # Generates boxplot for RMA normalized log intensity values.

### Annotation from Girke
biocLite("mouse4302.db")
library("mouse4302.db") # Loads required annotation package.
Annot <- data.frame(ACCNUM=sapply(contents(mouse4302ACCNUM), paste, collapse=", "), 
                    SYMBOL=sapply(contents(mouse4302SYMBOL), paste, collapse=", "), 
                    DESC=sapply(contents(mouse4302GENENAME), paste, collapse=", ")) # Constructs a data frame containing the gene IDs, gene symbols and descriptions for all probe sets on the chip.

all <- merge(Annot, eset_gcrma, by.x=0, by.y=0, all=T) # Merges everything with above expression data.
write.table(all, file="my_annot_file.xls", sep="\t", col.names = NA) # Exports data to text file that can be imported into Excel.     


## (B) DEG Analysis for GCRMA Data
         fit <- lmFit(eset_gcrma, design) # Fits a linear model for each gene based on the given series of arrays.
         fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
         fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
         gcrma_deg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
         gcrma_deg_result <- gcrma_deg_result[gcrma_deg_result$adj.P.Val<=0.05,]
         write.table(gcrma_deg_result, "gcrma_deg_result.xls", quote=FALSE, col.names = NA, sep="\t")
         
         ## (C) Create Venn Diagram for GCRMA Data
         gcrma_venn <- decideTests(fit2, p.value=0.05)
         pdf(file="gcrma_venn.pdf"); vennDiagram(gcrma_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.
         
         ## (A) Normalization: MAS5
         eset_mas5 <- mas5(data) # Normalizes the data with 'mas5' function and assigns them to exprSet object.
         exprs(eset_mas5) <- log2(exprs(eset_mas5)) # Only MAS5 stores absolute intensities. Other methods store log2 intensities.
         
         ## (B) DEG Analysis for MAS5 Data
         fit <- lmFit(eset_mas5, design) # Fits a linear model for each gene based on the given series of arrays.
         fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
         fit2 <- eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
         mas5_deg_result <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000)
         mas5_deg_result <- mas5_deg_result[mas5_deg_result$adj.P.Val<=0.05,]
         write.table(mas5_deg_result, "mas5_deg_result.xls", quote=FALSE, col.names = NA, sep="\t")
         
         ## (C) Create Venn Diagram for MAS5 Data
         mas5_venn <- decideTests(fit2, p.value=0.05)
         pdf(file="mas5_venn.pdf"); vennDiagram(mas5_venn); dev.off() # Creates venn diagram of all changed genes with p-value equal or less than 0.05.
         
         ## (D) Identifiy the Overlap Between the Three Methods
         overlap2 <-(merge(rma_deg_result, gcrma_deg_result, by.x = "ID", by.y = "ID", all = FALSE))
         overlap3 <-(merge(overlap2, mas5_deg_result, by.x = "ID", by.y = "ID", all = FALSE))
         write.table(overlap3, "overlap.xls", quote=FALSE, col.names = NA, sep="\t")
         
         