
##################################################################
# Full Script to Run WGCNA - with user input required throughout #
##################################################################

# ensure ssh forwarding is allowed


### User-defined details ###
############################

# Ensure module installation
library("WGCNA")
library("gtools")
library("RColorBrewer")
library("png")

# number of allowed threads
threads <- 10

## Minimum module size
# minimum number of modules to be included in any one module
minModuleSize <- 30

# how distinct do you want module colours to be from one another?
# (Max 100, min 0. May have to reduce if error appears, as not enough colours for modules)
distinction <- 100

# Set working directory
setwd("/net/isi-scratch/harryc/purkinje/6_WGCNA/final_runs")


## Required Files for WGCNA (all with headers)

# 1. A .txt numeric (except column 1) table of information about each of the samples (for this script, where the first column is sample names). Equivalent to ClinicalTraits.csv from the WGCNA sample files.
SampleInfo_file <- "sample_traits.txt"

# 2. A .csv table of information about each of the genes (for this script, where the first column is Ensembl.Gene.ID). Equivalent to GeneAnnotation.csv from the WGCNA sample files.
GeneInfo_file <- "mart_export.csv"

# 3. A .txt table of normalized expression data, where columns are sample names and rows are gene names. Equivalent to LiverFemale3600.csv from the WGCNA sample files.
ExpressionData_file <- "counts_noAd_normalized.txt"


## SELECTED POWER

# This can be selected after running part 1 of the main script
# and looking at the soft thresholding plots (scale independence and mean connectivity).

# The power selected should be the lowest power for which the scale-free
# topology fit index curve flattens out upon reaching a high value.

# if selected power is defined as NA, script can auto-pick the power as the first
# above a user-specified index (e.g. 0.85)

selected_power <- 3
power_threshold <- 0.85 # if required


## SELECTED MODULE MERGE HEIGHT

# This can be selected after running part 2 of the main script
# and looking at the module eigengene cluster plot
# (3_ClusteringAndCorrelation/PremergedModuleEigengeneClustering.png)

# merge height selected should be based on causing only those modules closest to merge

selected_MMH <- 0.5



##############################
##############################
##############################


### Main Script ### PART 1 ###
# run script in two halves


# reads in files
options(stringsAsFactors = FALSE)
SampleInfo <- read.table(SampleInfo_file,header=TRUE,row.names=1)
GeneInfo <- read.csv(GeneInfo_file,header=TRUE)
raw_ExpressionData <- read.table(ExpressionData_file,header=TRUE)
preQC_ExpressionData <- t(raw_ExpressionData)


## QC

dir.create("1_QC",showWarnings=FALSE)

# performs auto QC based on missing entries and zero-variance genes
gsg <- goodSamplesGenes(preQC_ExpressionData, verbose = 3)

# writes resulting gene lists
write.table(colnames(preQC_ExpressionData)[gsg$goodGenes],file="1_QC/genes_retained.txt",quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
write.table(colnames(preQC_ExpressionData)[!gsg$goodGenes],file="1_QC/genes_removed.txt",quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
write.table(rownames(preQC_ExpressionData)[gsg$goodSamples],file="1_QC/samples_retained.txt",quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
write.table(rownames(preQC_ExpressionData)[!gsg$goodSamples],file="1_QC/samples_removed.txt",quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)

# Remove the offending genes and samples from the data:
ExpressionData <- preQC_ExpressionData[gsg$goodSamples,gsg$goodGenes]

# clustering of samples to check for outliers
sampleTree <- flashClust(dist(ExpressionData),method="average")
# Plot the sample tree
png(file="1_QC/SampleClustering_RemovalCheck.png")
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)
dev.off()

## Use of Clinical Trait Data

# Clustering using Clinical Trait Data
sampleTree2 <- flashClust(dist(ExpressionData),method="average")
# Convert traits to a colour representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(SampleInfo,signed=FALSE)
# Plot the sample dendrogram and the colours underneath
png(file="1_QC/SampleClustering_RemovalCheck_ColourByTraitData.png")
plotDendroAndColors(sampleTree2,traitColors,groupLabels=names(SampleInfo),main="Sample dendrogram and trait heatmap")
dev.off()


## Network Construction

dir.create("2_PowerSelection",showWarnings=FALSE)
enableWGCNAThreads(nThreads=threads)


# Constructing a weighted gene network entails the choice of the soft thresholding power Beta to which co-expression similarity is raised to calculate adjacency [1].The authors of [1] have proposed to choose the soft thresholding power based on the criterion of approximate scale-free topology. We refer the reader to that work for more details; here we illustrate the use of the function pickSoftThreshold that performs the analysis of network topology and aids the user in choosing a proper soft-thresholding power. The user chooses a set of candidate powers (the function provides suitable default values), and the function returns a set of network indices that should be inspected.

# [1] B. Zhang and S. Horvath. A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology, 4(1):Article 17, 2005

# scale-free topology fit index = a measure of scale-free topology, determined from a Log-log plot of whole-network connectivity distribution. If the x-axis is the logarithm of whole network connectivity [log10(k)], and the y-axis is the logarithm of the corresponding frequency distribution [log10(p(k))], a scale-free topology would follow a straight line
# See Figure 2A - Langfelder and Horvath (2008) "WGCNA: an R package for weighted correlation network analysis" BMC Bioinformatics

# soft thresholding power = the value to which the original coexpression values are raised, so that high correlations are emphasized at the expense of low correlations (this does so in a "soft" manner, compared to "hard" thresholding, which would binarily separate correlations into coexpression and non-coexpression)
# you want to select the lowest power for which the scale-free topology fit index reaches a high value (e.g. 0.9)

# Choose a set of soft-thresholding powers
powers <- c(c(1:10),seq(from=12,to=20,by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(ExpressionData,powerVector=powers,verbose=5)

# Plot the results
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
png(filename="2_PowerSelection/SoftThresholdingPlots.png",width=960)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
        type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        labels=powers,cex=cex1,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1],sft$fitIndices[,5],
        xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="n",
        main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red")
dev.off()


##############################
##############################
##############################


### Main Script ### PART 2 ###
# can use once power has been selected


# if power not selected, auto-picks one
if(is.na(selected_power)){
    abline_flag <- TRUE
    selected_power <- powers[min(which(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] > power_threshold))]
}else{
    abline_flag <- FALSE
}

# Replots with selected power highlighted
cols <- rep("red",length(powers))
cols[which(powers==selected_power)] <- "green"

png(filename="2_PowerSelection/SoftThresholdingPlots.png",width=960)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
        type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        labels=powers,cex=cex1,col=cols)
if(abline_flag){abline(h=power_threshold,col="green")}
plot(sft$fitIndices[,1],sft$fitIndices[,5],
        xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="n",
        main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col=cols)
dev.off()


# calculates adjacency matrix
dir.create("3_ClusteringAndCorrelation",showWarnings=FALSE)
adjacency <- adjacency(ExpressionData,power=selected_power,type="signed")

# turn adjacency into topological overlap
# "to minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity"

# Topological overlap matrix = Matrix of Topological Overlap Measure. While the adjacency matrix considers each pair of genes in isolation, topological overlap considers each pair of genes in relation to all other genes in the network. More specifically, genes will have a high topological overlap if they are connected to roughly the same group of genes in the network (i.e. they share the same neighbourhood).
# See: Yip and Horvath (2007) "Gene network interconnectedness and the generalized topological overlap measure" BMC Bioinformatics


TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# hierarchical clustering of all genes, with dissimilarity based on topological overlap
# (flashclust used as it is fastest)
geneTree <- flashClust(as.dist(dissTOM),method="average")
# plot the resulting tree (dendrogram)
png(filename="3_ClusteringAndCorrelation/TOMClusteringDendrogram.png",width=960)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity",labels=FALSE, hang=0.04)
dev.off()

# The clustering dendrogram plotted by the last command is shown in Figure 2. In the clustering tree (dendrogram), each leaf, that is a short vertical line, corresponds to a gene. Branches of the dendrogram group together densely interconnected, highly co-expressed genes. Module identification amounts to the identification of individual branches ("cutting the branches of the dendrogram"). There are several methods for branch cutting; our standard method is the Dynamic Tree Cut from the package dynamicTreeCut. The next snippet of code illustrates its use.

# module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro=geneTree,distM=dissTOM,deepSplit=2,pamRespectsDendro=FALSE,minClusterSize=minModuleSize)

# number of genes in each module (module zero is reserved for unassigned genes)
table(dynamicMods)

## added own line here, to reduce colours to only those with highly differential rgb 
distinct_cols <- colors(distinct=T)[ apply( t(col2rgb(colors(distinct=T))) ,1, function(x) (max(x)-min(x))>distinction ) ]
distinct_cols <- distinct_cols[grep("^[[:alpha:]]*$", distinct_cols)]

# convert numeric labels into colors
dynamicColors=labels2colors(dynamicMods,colorSeq=distinct_cols)
table(dynamicColors)

# plot the dendrogram and colors underneath
png(filename="3_ClusteringAndCorrelation/DynamicTOMClusteringDendrogram.png",width=960)
plotDendroAndColors(geneTree,dynamicColors,"Dynamic Tree Cut",dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05,main="Gene dendrogram and module colors")
dev.off()

# calculate eigengenes
MEList <- moduleEigengenes(ExpressionData,colors=dynamicColors)
MEs <- MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- flashClust(as.dist(MEDiss),method="average")
png(filename="3_ClusteringAndCorrelation/PremergedModuleEigengeneClustering.png")
plot(METree, main="Clustering of module eigengenes",xlab="",sub="")
dev.off()


##############################
##############################
##############################


### Main Script ### PART 3 ###
# can use once cut-off for module merging has been decided

MEDissThres <- selected_MMH

# Call an automatic merging function
merge = mergeCloseModules(ExpressionData,dynamicColors,cutHeight=MEDissThres,verbose=3)
# the merged module colors
mergedColors=merge$colors
# Eigengenes of the new merged modules
mergedMEs <- merge$newMEs

# calculate dissimilarity of new eigengenes
mergedMEDiss <- 1-cor(mergedMEs)
# cluster new eigengenes
mergedMETree <- flashClust(as.dist(mergedMEDiss),method="average")
png(filename="3_ClusteringAndCorrelation/PostmergedModuleEigengeneClustering.png")
plot(mergedMETree, main="Merged clustering of module eigengenes",xlab="",sub="")
dev.off()

# plot new gene dendrogram with merged module colors
png(filename="3_ClusteringAndCorrelation/MergedDynamicTOMClusteringDendrogram.png",width=960)
plotDendroAndColors(geneTree,cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut","Merged Dynamic"), dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05,main="Gene dendrogram and module colors")
dev.off()

# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey",standardColors(50))
moduleLables <- match(moduleColors,colorOrder)-1
MEs <- mergedMEs

# table of number of genes in each module
moduleCounts <- table(moduleColors)
write.table(moduleCounts,file="3_ClusteringAndCorrelation/ClusteringGeneCounts.txt",quote=FALSE,row.names=FALSE)

# all genes and their assigned module
geneAssignment <- cbind(colnames(ExpressionData),moduleColors)
write.table(geneAssignment,file="3_ClusteringAndCorrelation/ClusteringGeneAssignment.txt",quote=FALSE,row.names=FALSE)


## Gene Identification

nGenes <- ncol(ExpressionData)
nSamples <- nrow(ExpressionData)
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(ExpressionData,moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs,SampleInfo,use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
# obtain variances explained by eigengenes (for plotting later)
var_exp <- moduleEigengenes(ExpressionData,moduleColors)$varExplained
names(var_exp) <- names(MEs0)

## Correlation heatmap plot
png(filename="3_ClusteringAndCorrelation/CorrelationHeatmap.png",height=960,width=600)
# Will display correlations and their p-values
orderedModuleCounts <- moduleCounts[gsub("ME","",rownames(moduleTraitCor))]
textMatrix <- paste(signif(moduleTraitCor,2),"  (p-value ",
                signif(moduleTraitPvalue,1),")\n",as.vector(orderedModuleCounts)," genes",sep="")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix=moduleTraitCor, xLabels=names(SampleInfo),
                yLabels=names(MEs), ySymbols=names(MEs),
                colorLabels=FALSE, colors=blueWhiteRed(50),
                textMatrix=textMatrix, setStdMargins=FALSE,
                cex.text=1, zlim=c(-1,1),
                main=paste("Module-trait relationships"))
dev.off()



## Module-trait correlatons

dir.create("4_CorrelationResults",showWarnings=FALSE)

# calculates correlation data
geneModuleMembership <- as.data.frame(cor(ExpressionData,MEs,use="p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))

modNames <- substring(names(MEs),3)
names(geneModuleMembership) <- paste("MM",modNames,sep="")
names(MMPvalue) <- paste("p.MM",modNames,sep="")

# loops through traits
for(trait in colnames(SampleInfo)){
    
    tmp_trait <- as.data.frame(SampleInfo[trait])
    dir.create(paste("4_CorrelationResults",trait,sep="/"),showWarnings=FALSE)
    
    # average MEs (for plotting later)
    MEs_mean <- c()
    for(i in unique(unlist(tmp_trait))){
        MEs_mean <- rbind(MEs_mean, apply( MEs[which(unlist(tmp_trait)==i),], 2, mean) )
    }
    MEs_mean <- as.data.frame(MEs_mean)
    
    # calculates trait significance in all genes
    geneTraitSignificance <- as.data.frame(cor(ExpressionData,tmp_trait,use="p"))
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
    names(geneTraitSignificance) <- paste("GS.",names(tmp_trait),sep="")
    names(GSPvalue) <- paste("p.GS.",names(tmp_trait),sep="")
    
    write.table(abs(geneTraitSignificance),
                        file=paste("4_CorrelationResults/",trait,"/",trait,"_GeneTraitSignificance.txt",sep=""),
                        quote=FALSE)
    
    # loops through all modules
    tmp_cor <- rownames(moduleTraitCor)
    for(module in tmp_cor){
        
        tmp_module <- gsub("ME","",module)
        dir.create(paste("4_CorrelationResults",trait,tmp_module,sep="/"),showWarnings=FALSE)
        
        column <- match(tmp_module,modNames)
        moduleGenes <- moduleColors==tmp_module
        
        # plots correlation of genes within this module
        png(filename=paste("4_CorrelationResults/",trait,"/",tmp_module,"/",trait,"-",tmp_module,"_AllGenesCorrelation.png",sep=""))
        par(mfrow=c(1,1))
        verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                        abs(geneTraitSignificance[moduleGenes,1]),
                        xlab=paste("Module Membership in",tmp_module,"module"),
                        ylab=paste("Gene significance for",trait),
                        main=paste("Module membership vs. gene significance\n"),
                        cex.main=1.2, cex.lab=1.2, cex.axis=1.2, col=tmp_module)
        dev.off()
        
        png(filename=paste("4_CorrelationResults/",trait,"/",tmp_module,"/",trait,"-",tmp_module,"_AllGenesCorrelation(not_absolutes).png",sep=""))
        par(mfrow=c(1,1))
        verboseScatterplot(geneModuleMembership[moduleGenes,column],
                        geneTraitSignificance[moduleGenes,1],
                        xlab=paste("Module Membership in",tmp_module,"module"),
                        ylab=paste("Gene significance for",trait),
                        main=paste("Module membership vs. gene significance\n(not using absolute values)\n"),
                        cex.main=1.2, cex.lab=1.2, cex.axis=1.2, col=tmp_module)
        dev.off()
        
        # writes table of genes within module, and membership to module
        write.table(abs(geneModuleMembership[moduleGenes,column,drop=FALSE]),file=paste("4_CorrelationResults/",trait,"/",tmp_module,"/",trait,"-",tmp_module,"_ModuleGenesMembership.txt",sep=""),quote=FALSE,col.names=FALSE)
        
        # plots raw PC1 vs. trait for this module
        png(filename=paste("4_CorrelationResults/",trait,"/",tmp_module,"/",trait,"-",tmp_module,"_PC1.png",sep=""))
        par(mar=c(5.1,6.6,4.1,2.1))
        plot(unlist(tmp_trait),unlist(MEs[module]),xlab=trait,col=tmp_module,
                        ylab=paste("Module Eigengenes\n(i.e. PC1 of module)\nExplaining ",round(var_exp[module]*100,2),"% of the Variance",sep=""),
                        main=paste("Module (",tmp_module,") Eigengenes Vs. ",trait,sep=""))
        dev.off()
        
        # plots mean PC1 vs. trait for this module
        png(filename=paste("4_CorrelationResults/",trait,"/",tmp_module,"/",trait,"-",tmp_module,"_mean_PC1.png",sep=""))
        par(mar=c(5.1,6.6,4.1,2.1))
        plot(sort(unique(unlist(tmp_trait))),unlist(MEs_mean[module])[order(unique(unlist(tmp_trait)))],
            ylab=paste("Module Eigengenes\n(i.e. PC1 of module)\nExplaining ",round(var_exp[module]*100,2),"% of the Variance",sep=""),
            type="b",xlab=trait,col=tmp_module,main=paste("Module (",tmp_module,") Eigengenes Vs. ",trait,sep=""))
        dev.off()
        
        ## plots all genes in a module (as you see with MFuzz plots)
        
        # obtains list of genes hard clustered to module
        subset_genes <- geneAssignment[which(geneAssignment[,2]==tmp_module),1]
        
        # sorts expression data by decreasing Module Membership p-value (so heat colour increases as p-value decreases)
        subset_MMPvalue <- MMPvalue[subset_genes,]
        decP_genes <- rownames(subset_MMPvalue)[order(subset_MMPvalue[[paste("p.MM",tmp_module,sep="")]],decreasing=T)]
        subset_EData <- raw_ExpressionData[decP_genes,]
        
        # takes mean for each stage
        subset_EData_mean <- c()
        for(i in unique(unlist(tmp_trait))){
            subset_EData_mean <- cbind(subset_EData_mean,
                apply( subset_EData[,rownames(tmp_trait)[which(tmp_trait==i)]], 1, mean) )
        }
        colnames(subset_EData_mean) <- unique(unlist(tmp_trait))
        subset_EData_mean <- subset_EData_mean[,mixedsort(colnames(subset_EData_mean))]
        
        # standardizes for plotting
        final_EData <- t(apply(subset_EData_mean,1,function(p) (p-mean(p))/sd(p)))
        
        ## plots
        
        ylims <- c( floor(min(final_EData)) , ceiling(max(final_EData)) )
        
        # build palette based on module colour
        
        tinted_cols <- do.call("cbind",(lapply(
                            seq(0.8,0,-0.025), function(x) col2rgb(tmp_module)+(x*(255-col2rgb(tmp_module)))
                        )))
        tinted_cols_hex <- apply(tinted_cols,2,function(x) rgb(x[1],x[2],x[3],max=255))
        
        shaded_cols <- do.call("cbind",(lapply(
                            seq(1,0.8,-0.01), function(x) col2rgb(tmp_module)*x
                        )))
        shaded_cols_hex <- apply(shaded_cols,2,function(x) rgb(x[1],x[2],x[3],max=255))
        
        cols_hex <- c(tinted_cols_hex,shaded_cols_hex)
        pal_cols <- colorRampPalette(cols_hex)(nrow(final_EData))
        
        png(filename=paste("4_CorrelationResults/",trait,"/",tmp_module,"/",trait,"-",tmp_module,"_FULLPLOT.png",sep=""))
        plot( final_EData[1,],type="l",ylim=ylims, col=pal_cols[1],
                main=paste(tmp_module,"gene expression\ncoloured by p-value"),
                xaxt="n", xlab=trait,
                ylab="Standardized Mean Normalized Expression Level"
            )
        axis(1,at=axTicks(1),labels=colnames(final_EData))
        for(j in 2:nrow(final_EData)){ lines( final_EData[j,],type="l", col=pal_cols[j] )}
        dev.off()
        
    }
    
    
    ## Summary output of network analysis results
    
    probes <- colnames(ExpressionData)
    probes2annot <- match(probes,GeneInfo$Ensembl.Gene.ID)
    geneInfo0 <- GeneInfo[probes2annot,]
    geneInfo0$Ensembl.Gene.ID <- probes
    geneInfo0 <- cbind(geneInfo0,moduleColors,geneTraitSignificance,GSPvalue)
    # Order modules by their significance for trait
    modOrder <- order(-abs(cor(MEs,tmp_trait,use="p")))
    # Add module membership information in the chosen order
    for(mod in 1:ncol(geneModuleMembership)){
        oldNames <- names(geneInfo0)
        geneInfo0 <- data.frame(geneInfo0,geneModuleMembership[,modOrder[mod]],
                                MMPvalue[,modOrder[mod]])
        names(geneInfo0) <- c(oldNames,paste("MM.",modNames[modOrder[mod]],sep=""),
                                paste("p.MM.",modNames[modOrder[mod]],sep=""))
    }
    # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
    geneOrder <- order(geneInfo0$moduleColors,-abs(geneInfo0[paste("GS",trait,sep=".")]))
    geneInfo <- geneInfo0[geneOrder,]
    
    write.table(geneInfo,file=paste("4_CorrelationResults/",trait,"/",trait,"_NetworkSummary.txt",sep=""),quote=FALSE)
    
    
}


## Visualization of networks

dir.create("5_NetworkVisualization",showWarnings=FALSE)

# transform topological overview with a power to make moderately strong connections more visible in heatmap
plotTOM <- dissTOM^selected_power
diag(plotTOM) <- NA
#TOMplot(plotTOM,geneTree,moduleColors,main="Network heatmap plot, all genes")

MET <- orderMEs(cbind(MEs,tmp_trait))

png(filename="5_NetworkVisualization/EigengeneDendrogram.png")
plotEigengeneNetworks(MET,"Eigengene dendrogram",marDendro=c(1,5,3,1),plotHeatmaps=FALSE)
dev.off()

png(filename="5_NetworkVisualization/EigengeneHeatmap.png")
plotEigengeneNetworks(MET,"Eigengene adjacency heatmap",marHeatmap=c(6,6,4,4),plotDendrograms=FALSE,xLabelsAngle=90)
dev.off()


## Additional functionality - replots fullplots in one

trait <- colnames(SampleInfo)[1]

mfrow_list <- c( "","1,2","1,3","1,4","2,3","2,3","2,4","2,4","3,3","3,4","3,4","3,4","3,5","3,5","3,5","4,4","3,6","3,6","3,7","3,7","3,7","4,6","4,6","4,6","5,5" )
tmp_mfrow <- as.numeric(unlist(strsplit( mfrow_list[length(rownames(moduleTraitCor))] ,",")))

png(filename="AllModules_FullPlots.png",width=500*tmp_mfrow[2],height=500*tmp_mfrow[1])
par(mfrow=tmp_mfrow,oma=rep(0,4),mar=rep(0,4))
for(module in rownames(moduleTraitCor)){
    tmp_module <- gsub("ME","",module)
    ima <- readPNG(paste("4_CorrelationResults/",trait,"/",tmp_module,"/",trait,"-",tmp_module,"_FULLPLOT.png",sep=""))
    plot.new() ; rasterImage(ima,0,0,1,1)
}
dev.off()











