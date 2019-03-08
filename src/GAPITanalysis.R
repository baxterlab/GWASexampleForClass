###Run once to install multtest
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("multtest")
rm(list=ls())
library(reshape2)
library("multtest") 
library("gplots") 
#install.packages("LDheatmap")
library("LDheatmap") 
library("genetics") 
#install.packages("ape")
library("ape")
#install.packages("EMMREML")
library("EMMREML")
library(compiler)
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt") 
source("http://zzlab.net/GAPIT/emma.txt") 

####Choose phenotype file to run####
#myY <- read.table("../data/standardizedMaizeDataFL06.commonNames.csv", head = TRUE, sep=",") 
myY <- read.table("../data/standardizedMaizeDataPU09.commonNames.csv", head = TRUE, sep=",") 

load("../data/Maize282.0minor2major.maf05missing.2.missingImptoMajor.rda") #loads 'allGeno' object

#make snp info table
snp_info = cbind(data.frame(Name=colnames(allGeno)),colsplit(colnames(allGeno), "_", c("Chromosome", "Position","other"))[,1:2])
#remove characters from chromsome names
snp_info$Chromosome <- as.numeric(gsub("\\D{1,2}",'', snp_info$Chromosome, perl=TRUE))    

#add taxa column to beginning of allGeno
allGeno <- cbind(data.frame(Taxa=row.names(allGeno),stringsAsFactors = FALSE),allGeno)

#get pca file from previous gapit run
#myQ <- read.table("~ekesinger/GAPIT Genotype FL06/GAPIT.PCA2.csv",head=TRUE,sep=",")
#myQ <- myQ[which(myQ$taxa %in% myY$Taxa),1:4]

###Narrow down genotype file to just what intersects between the files###
allGeno <- allGeno[which(rownames(allGeno) %in% intersect(rownames(allGeno), myY$Taxa)),]

###Do same for phenotype file###
myY <- myY[which(myY$Taxa %in% intersect(rownames(allGeno), myY$Taxa)),]

#myQ <- myQ[which(myQ$taxa %in% intersect(myQ$taxa, rownames(allGeno))),1:4]
###Phenotype and genotype files should now have same number of rows###
stopifnot(nrow(allGeno)==nrow(myY))


#make sure everything is in the same order...
allGeno <- allGeno[order(row.names(allGeno)),]
myY <- myY[order(myY$Taxa),]
#myQ <- myQ[order(myQ$taxa),]

##Choose trait(s) to analyze here##
colnames(myY)
myYtraitToAnalyze <- myY[,c("Taxa","Mo98","Mn55","Cd111","Cu65")]

###This will subsamples genotype file so it runs quickly (but gets less results)
# keepCols <- sample(2:ncol(allGeno),30000,replace = FALSE)
# allGeno <- allGeno[,c(1,keepCols)]
# snp_info <- snp_info[keepCols,]

####Gapit will spit everything into current working directory, 
######so change to desired output directory
setwd("../results/")
myGAPIT <- GAPIT( 
  Y=myYtraitToAnalyze, 
  GD=allGeno,
  GM=snp_info,
  PCA.total=3, #calculate and use 3 PCs 
  SNP.fraction=0.5, #only use half of the SNPs to calculate PCs and Kinship
  #  file.fragment = 1000,
  SNP.MAF=0.05
)

#parse results#
####Just keep SNPs with small pval####
#####Not necessarily significant, just narrowing it down
pvalThresh <- 1e-4
traitsAnalyzed <- colnames(myYtraitToAnalyze)[2:ncol(myYtraitToAnalyze)]
allResults <- data.frame(stringsAsFactors = FALSE)
for(i in traitsAnalyzed){
  thisResult <- read.table(paste0("GAPIT.MLM.",i,".GWAS.Results.csv"),sep=",",stringsAsFactors = FALSE, header=TRUE)
  thisResult <- thisResult[thisResult$P.value <= pvalThresh,]
  thisResult$Trait <- i
  allResults <- rbind(allResults,thisResult)
}

write.table(allResults, file="allGAPITresults.csv", sep=",", row.names=FALSE, col.names=TRUE)

