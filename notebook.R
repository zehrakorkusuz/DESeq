if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))

## Upload all the packages
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

# Read the data into R
seqdata <- read.delim("data/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("data/SampleInfo.txt", stringsAsFactors = TRUE)


head(seqdata)
dim(seqdata)

sampleinfo

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
colnames(countdata)

# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)

head(countdata)

table(colnames(countdata)==sampleinfo$SampleName)

### Convert counts to DGEList object ###
y <- DGEList(countdata)
# have a look at y
y

# See what slots are stored in y
names(y)

y$samples


#groups

group <- paste(sampleinfo$CellType, sampleinfo$Status, sep=".")
group

group <- factor(group)
group

y$samples$group <- group
y$samples

### Adding Annotation

columns(org.Mm.eg.db)
ann <- select(org.Mm.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))
head(ann)

# double check that entrez id column matches our y$counts row name
table(ann$ENTREZID==rownames(y$counts))

y$genes <- ann


### filtering lowly expressed genes

# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

"""Quality control
Now that we have got rid of the lowly expressed genes and have 
our counts stored in a DGEList object, we can look at a few different plots 
to check that the data is good quality, and that the samples are as we would expect."""
y <- y[keep, keep.lib.sizes=FALSE]
y$samples$lib.size

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# we can also adjust the labelling if we want
barplot(y$samples$lib.size/1e06, names=colnames(y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")


plotMDS(y)
# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType)
## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(sampleinfo$Status)

col.status <- c("blue","red","black")[sampleinfo$Status]
col.status

plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

#We need to update sample info

sampleinfo

# sampleinfo object with the corrected sample info
sampleinfo <- read.delim("data/SampleInfo.txt", stringsAsFactors = TRUE)
sampleinfo

# We need to correct the info for the groups
group <- factor(paste(sampleinfo$CellType,sampleinfo$Status,sep="."))
y$samples$group <- group

# Redo the MDSplot with corrected information
par(mfrow=c(1,2))
col.cell <- c("purple","orange")[sampleinfo$CellType]
col.status <- c("blue","red","black")[sampleinfo$Status]
plotMDS(y,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

