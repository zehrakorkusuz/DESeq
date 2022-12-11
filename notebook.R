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