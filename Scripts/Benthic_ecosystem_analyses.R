library(lubridate)
library(birk)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(microbiome)
library(WGCNA)
library(mixOmics)
library(igraph)
library(reshape2)
library(HiveR)
library(FactoMineR)
library(ggrepel)
library(Hmisc)
library(corrplot)
library(treemap)
library(gplots)
library(ggcorrplot)
library(RColorBrewer)

##################################
# Import Environmental variables #
##################################
# Load environmental data
all_benthic_data<-read.csv("Benthic_environmental_data.csv", sep=";", h=T, dec=',') # all corrected
# Set the row names of 'all_benthic_data' to the first column of the dataset
rownames(all_benthic_data)<-all_benthic_data[,c(1)]
# Remove the first column from 'all_benthic_data' as it is now used for row names.
all_benthic_data<-all_benthic_data[,-c(1)]
  
#######################
### Import Omics data #
#######################
# Load sample metadata 
samples_df<-read.csv("Benthic/Benthic_samples.csv", sep=";", h=T, dec='.')

# Load raw omics data
ASV_Benthic<-read.csv("Benthic/18S_Benthic_data.csv", sep=";", h=T, dec='.')

# Extract the first column of 'ASV_Benthic', which contains ASV (Amplicon Sequence Variant) identifiers.
ASVs<-ASV_Benthic[,c(1)]

# Extract the columns representing abundance data for each sample
ASV_mat<-ASV_Benthic[,c(3:147)

Set the column names of 'ASV_mat' to match the sample names
colnames(ASV_mat)<-samples_df$Sample
ASV_mat<-cbind(ASVs,ASV_mat)

# Extract taxonomy information
tax_mat<-ASV_Benthic[,c(148:155)]
tax_mat<-cbind(ASVs,tax_mat)

#################################################################################
############################## Build ASV data base for analysis #################
#################################################################################
row.names(ASV_mat) <- ASV_mat$ASVs
ASV_mat <- ASV_mat[,-c(1)] 
row.names(tax_mat) <- tax_mat$ASVs
tax_mat <- tax_mat[,-c(1)] 
row.names(samples_df) <- samples_df$Sample

# Keep sample 
keep_<-match(rownames(all_benthic_data),as.character(samples_df$Trap_name))
keep_<-keep_[!is.na(keep_)]

samples<-sample_data(samples_df[c(keep_),])
ASV_mat<-ASV_mat[,c(keep_[!is.na(keep_)])]
tax_mat <- as.matrix(tax_mat)

ASV = otu_table(ASV_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
Omicsbenthic <- phyloseq(ASV, TAX, samples)

#################################################################################
############################## Summarize Data Base ##############################
#################################################################################

# Check if any ASV are not present in any samples. If True removed
any(taxa_sums(Omicsbenthic)==0)
Omicsbenthic_a <- prune_taxa(taxa_sums(Omicsbenthic) > 0, Omicsbenthic)
any(taxa_sums(Omicsbenthic_a)==0)

# Original abundance ASV versus removed
ntaxa(Omicsbenthic)-ntaxa(Omicsbenthic_a)

#################################################################################
############################## Remove Phylum unclassified #######################
#################################################################################

### Composition
Omicsbenthic_a.com<-subset_taxa(Omicsbenthic_a, !is.na(Supergroup))

taxic <- as.data.frame(Omicsbenthic_a.com@tax_table) 
taxic$ID <- paste(as.matrix(taxic)[cbind(1:nrow(taxic), max.col(!is.na(taxic),"last")-1)],"",as.matrix(taxic)[cbind(1:nrow(taxic), max.col(!is.na(taxic),"last"))])

# Add the ASV ids from ASV table into the taxa table at the end.
taxic$ASV <- rownames(taxic) 

# You can see that we now have extra taxonomy levels.
colnames(taxic)

# convert it into a matrix.
taxmat <- as.matrix(taxic)

# convert into phyloseq compaitble file.
new.tax <- tax_table(taxmat)  

# incroporate into phyloseq Object
tax_table(Omicsbenthic_a.com) <- new.tax 

##############################################################
############################## Remove Singleton  #############
##############################################################

summarize_phyloseq(Omicsbenthic_a.com) # We have 95 singleton

Benthic_matrice <- prune_samples(sample_sums(Omicsbenthic_a.com) > 4500, Omicsbenthic_a.com)

# Filter ASVs that have 50 sequences in at least 3 samples
subset_3samples = genefilter_sample(Benthic_matrice, filterfun_sample(function(x) x >= 50), A=3)

################## Remove biais
subset_3samples[subset_3samples="ASV_45"]=FALSE			# Capra_hircus
subset_3samples[subset_3samples="ASV_B_1379"]=FALSE		# Capra_hircus
subset_3samples[subset_3samples="ASV_B_11847"]=FALSE	# Capra_hircus

subset_3samples[subset_3samples="ASV_537"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_578"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_603"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_659"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_1179"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_1718"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_1765"]=FALSE		# Homo_sapiens
###############################################

Benthic_matrice2 = prune_taxa(subset_3samples, Benthic_matrice)
ntaxa(Benthic_matrice2)
nsamples(Benthic_matrice2)

z<-plot_richness(Benthic_matrice2)
z<-z$data
index<-z[which((z$variable=="Shannon")|(z$variable=="Simpson")|(z$variable=="Chao1")),c(1,7:8)]
index<-dcast(index,Sample~variable)
rownames(index)<-index[,c(1)]
index<-index[,-c(1)]

# Combine indeces in the environmental table
all_benthic_data<-merge(all_benthic_data,index, by=0)
rownames(all_benthic_data)<-all_benthic_data[,c(1)]
all_benthic_data<-all_benthic_data[,-c(1)]

# Print the filtered matrix
asvs <- otu_table(Benthic_matrice2)
asvs <- as.data.frame(asvs)

#################################################################################
######################### CLR Transformation ####################################
#################################################################################
# Counting_Table
Benthic_matrice2_clr<-microbiome::transform(Benthic_matrice2,'clr')
datASVClr<-t(otu_table(Benthic_matrice2_clr))

datASVClr_2<-t(otu_table(Benthic_matrice2))
#################################################################################
######################### Running WGCNA.     ####################################
#################################################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datASVClr, powerVector = powers, verbose = 5)

par(mfrow=c(1,2))    # set the plotting area into a 1*2 array
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red")
abline(h=0.90, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
#################################################################################

# Turn adjacency into topological overlap matrix (TOM)
adjacency <- adjacency(datASVClr, power = 6, type="signed")
adjacency[adjacency < 0] = 0
adjacency[adjacency > 1] = 1
TOMadj = TOMsimilarity(adjacency, TOMType="signed")
dissTOMadj <- 1- TOMadj

# Clustering using TOM
# Call the hierarchical clustering function 
hclustGeneTree <- hclust(as.dist(dissTOMadj), method = "average")

# Plot the resulting clustering tree (dendogram)
sizeGrWindow(12, 9)
plot(hclustGeneTree, xlab = "", sub = "", 
     main = "Gene Clustering on TOM-based disssimilarity", 
     labels = FALSE, hang = 0.04)

# Make the modules larger, so set the minimum higher
minModuleSize <- 20

# Module ID using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = hclustGeneTree, 
                             distM = dissTOMadj,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)

# Change some color
dynamicColors[dynamicColors=="yellow"]="orange"
dynamicColors[dynamicColors=="green"]="grey"

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(hclustGeneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
dynamic_MEList <- moduleEigengenes(datASVClr, colors = dynamicColors)
dynamic_MEs <- dynamic_MEList$eigengenes

# Calculate dissimilarity of module eigengenes
dynamic_MEDiss <- 1-cor(dynamic_MEs)
dynamic_METree <- hclust(as.dist(dynamic_MEDiss), method="average")

# Plot the hclust
sizeGrWindow(7,6)
plot(dynamic_METree, main = "Dynamic Clustering of module eigengenes",
     xlab = "", sub = "")
######################## MERGE SIMILAR MODULES
dynamic_MEDissThres <- 0.25

# Plot the cut line
abline(h = dynamic_MEDissThres, col = "red")

# Call an automatic merging function
merge_dynamic_MEDs <- mergeCloseModules(datASVClr, dynamicColors, cutHeight = dynamic_MEDissThres, verbose = 3)

# define number of samples
nSamples <- nrow(datASVClr)

moduleColors = merge_dynamic_MEDs $colors;

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datASVClr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Select Environmental that macth with our omics data
datTraits<-all_benthic_data
datTraits<-datTraits[c(match(as.character(rownames(sample_data(Benthic_matrice2))),rownames(datTraits))),]
datTraits<-datTraits[,c(19,18,17,16,14,1:12,21:27,29:31,47:50)]

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# PLOT
p_val_module<-signif(moduleTraitPvalue, 1)
p_val_module1<-p_val_module
cor_module<-signif(moduleTraitCor, 2)
cor_module[p_val_module1>=0.05]<-" "

p_val_module[0.0001>=p_val_module1]<-"****"
p_val_module[(0.001>=p_val_module1) & (p_val_module1>0.0001)]<-"***"
p_val_module[(0.01>=p_val_module1) & (p_val_module1>0.001)]<-"**"
p_val_module[(0.05>p_val_module1) & (p_val_module1>0.01)]<-"*"
p_val_module[p_val_module1>=0.05]<-" "

sizeGrWindow(10,6)
textMatrix <- p_val_module
dim(textMatrix) <- dim(moduleTraitCor)

# Display the correlation values within a heatmap
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = colnames(datTraits),
               yLabels = names(MEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1), xLabelsAngle=90, 
			   plotLegend=0,horizontalSeparator.lty=2)

ggcorrplot(round(moduleTraitCor,3),
		ggtheme = ggplot2::theme_bw,
		p.mat = as.matrix(moduleTraitPvalue), 
		sig.level = 0.005, 
		insig = "blank", 
		colors = c("#6D9EC1", "white", "#E46726"), 
		method = "circle",
		outline.color="black",
		tl.cex=15)
	   
# Define variable weight containing the weight column of datTrait
CORG<- as.data.frame(datTraits$CORG)
names(CORG) <- "CORG" # rename

# Calculate the correlations between modules
geneModuleMembership <- as.data.frame(WGCNA::cor(datASVClr, MEs, use = "p"))
# p-values
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# Calculate the correlation for the trait
geneTraitSignificance <- as.data.frame(WGCNA::cor(datASVClr, CORG, use = "p"))
# p-values
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nSamples))

names(geneTraitSignificance) <- paste("GS.", names(CORG), sep = "")
names(GSPvalue) <- paste("p.GS.", names(CORG), sep = "")

# Remove the first two letters in string
modNames <- substring(names(geneModuleMembership),3)

par(mfrow = c(2,4))  
# Initialize for loop and visualize Gene Significance vs module membership of the modules
for (i in names(geneModuleMembership)) {
  
  # Pull out the module we're working on
  module <- substring(i,3)
  print(module)   
  
  # Find the index in the column of the dataframe 
  column <- match(module, modNames)
  
  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors == module
  vecASVnames = rownames(geneTraitSignificance)
  print(paste("There are ", length(vecASVnames[moduleGenes]), " ASVs in the ", module, " module.", sep = ""))
  print(vecASVnames[moduleGenes])
  
  # Make the plot
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                 abs(geneTraitSignificance[moduleGenes, 1]),
                 xlab = paste("Module Membership in", module, "module"),
                 ylab = "Gene significance for Corg",
                 main = paste("Module membership vs. gene significnace \n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}    

select_C=names(geneModuleMembership)[c(2)] # Positively correlated to carbon content in sediment cores
module <- select_C
print(module)
# Select ASVs in the module selected
moduleGenes = is.finite(match(moduleColors, substring(module,3)))
vecASVnames = rownames(geneTraitSignificance)
species_module<-vecASVnames[moduleGenes]
color_<-moduleColors[moduleGenes]
taxon<-tax_table(Benthic_matrice2_clr)[match(species_module,rownames(tax_table(Benthic_matrice2_clr)))]

# Replot Gene Significance vs module membership for modules interests (here carbon organic content)
par(mfrow = c(1,2))
for (i in select_C) {
  # Pull out the module we're working on
  module <- substring(i,3)
  print(module)   
  # Find the index in the column of the dataframe 
  column <- match(module, modNames)
  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors == module
  vecASVnames = rownames(geneTraitSignificance)
  print(paste("There are ", length(vecASVnames[moduleGenes]), " ASVs in the ", module, " module.", sep = ""))
  print(vecASVnames[moduleGenes])
  # Make the plot
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                 abs(geneTraitSignificance[moduleGenes, 1]),
                 xlab = paste("Module Membership in", module, "module"),
                 ylab = "Gene significance for Corg",
                 main = paste("Module membership vs. gene significnace \n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}

#######################################################################################################
############################ Modules relative abundance across the samples ############################
#######################################################################################################

# Transformation en abundance relative
Benthic_matrice2.rel<- microbiome::transform(Benthic_matrice2, "compositional")
taxmat <- tax_table(Benthic_matrice2.rel)[,c(1:8)]
datExpr<-t(otu_table(Benthic_matrice2.rel))

annot = data.frame(taxmat)
dim(annot)
probes = colnames(datExpr)
probes2annot = match(probes, rownames(annot))

moduleColors = merge_dynamic_MEDs $colors;

ASVInfo = data.frame(ASVs = probes,
                       Kingdom = annot$Kingdom[probes2annot],
                       Supergroup = annot$Supergroup[probes2annot],
                       Division = annot$Division[probes2annot],
                       Class = annot$Class[probes2annot],
                       Order = annot$Order[probes2annot],
					   Family = annot$Family[probes2annot],
					   Genus = annot$Genus[probes2annot],
					   Species = annot$Species[probes2annot],
                       moduleColor = moduleColors)
table_module <- ASVInfo[,c(1:10)]

table_perc <- t(datExpr)
table_perc <- as.data.frame(table_perc)
table_perc$tot <- rowSums(table_perc)
table_perc<- cbind(ASVs = rownames(table_perc), table_perc)
abond <- table_perc[,c(1,141)]
abond <- as.data.frame(abond)

table_perc <- table_perc[,-c(1,141)]
table_perc <- t(table_perc)

perc_module <- merge(table_module, abond, by = "ASVs")
com_imp<-unique(moduleColors)

# Extract ASVs by modules
module_com_1 <- perc_module[perc_module$moduleColor == com_imp[1],]
module_com_2 <- perc_module[perc_module$moduleColor == com_imp[2],]
module_com_3 <- perc_module[perc_module$moduleColor == com_imp[3],]
module_com_4 <- perc_module[perc_module$moduleColor == com_imp[4],]
module_com_5 <- perc_module[perc_module$moduleColor == com_imp[5],]
module_com_6 <- perc_module[perc_module$moduleColor == com_imp[6],]


Contr_ASVs_modl <- data.frame(com_1 = apply(table_perc[,colnames(table_perc) %in% module_com_1$ASVs], 1, sum),
                       com_2 = apply(table_perc[,colnames(table_perc) %in% module_com_2$ASVs], 1, sum),
                       com_3 = apply(table_perc[,colnames(table_perc) %in% module_com_3$ASVs], 1, sum),
                       com_4 = apply(table_perc[,colnames(table_perc) %in% module_com_4$ASVs], 1, sum),
                       com_5 = apply(table_perc[,colnames(table_perc) %in% module_com_5$ASVs], 1, sum),
                       com_6 = apply(table_perc[,colnames(table_perc) %in% module_com_6$ASVs], 1, sum))	

colnames(Contr_ASVs_modl)<-com_imp

info_sample<-sample_data(Benthic_matrice2.rel)

df <- cbind(samples = rownames(Contr_ASVs_modl), Contr_ASVs_modl)
df <- cbind(Year = info_sample$Year, df)
df <- cbind(Position = info_sample$Position, df)

df_m <- melt(df, id.vars=c("samples", "Year","Position"))

df_m$Position <- factor(df_m$Position, levels = unique(df_m$Position))
df_m$Year <- factor(df_m$Year, levels = sort(unique(df_m$Year)))

mcast <- dcast(df_m, Position + Year ~ variable, mean)
df_m2 <- melt(mcast, id.vars=c("Position","Year"))

# Change the default color
df_m2$color<-col2hex(df_m2$variable)
new_color<-as.character(df_m2$variable)
new_color[which(new_color=="orange")] <- "#66c2a5"
new_color[which(new_color=="blue")] <- "#fc8d62"
new_color[which(new_color=="turquoise")] <- "#e78ac3"
new_color[which(new_color=="grey")] <- "#a6d854"
new_color[which(new_color=="brown")] <- "#e5c494"
new_color[which(new_color=="red")] <- "#b3b3b3"

# Order for visualisation
df_m2$Position<- factor(df_m2$Position, levels=c('EG1','EG2','EG3','EG4','HG9','HG8','HG7','HG6','HG5','HG4','HG3','HG2','HG1','N5','N4','N3','N2','N1','S3'))
df_m2$Gradient<-"Longitudinal"
df_m2[which((df_m2$Position=="N5")|(df_m2$Position=="N4")|(df_m2$Position=="N3")|(df_m2$Position=="N2")|(df_m2$Position=="N1")|(df_m2$Position=="S3")),c(6)]<-"Latitudinal"
df_m2$Gradient<- factor(df_m2$Gradient, levels=c('Longitudinal','Latitudinal'))

# Plot
df_m2$new_color<-new_color
ggplot(df_m2, aes(x=Year, y = value, fill=variable)) +
    geom_bar(stat = "identity", position = "fill", show.legend = FALSE) +
    facet_grid(.~Position, scales="free", space="free_x")+
	scale_fill_manual(values=unique(df_m2$new_color))+
    theme_bw()+
	coord_flip() +
    theme(strip.text.y = element_text(size = 15, hjust = 0), strip.background = element_blank(), axis.text.x = element_text(size = 10 , angle = 40, hjust=1), 
	axis.text.y = element_text(size = 10),legend.position="none")+
	ylab("Subnetwork contribution")

#############################################################
### Relative abundance of the main families across module ###
#############################################################
table <- perc_module[,c(6,10,11)]
cast <- dcast(table, Order ~ moduleColor, sum) 

##   invalid factor level, NA generated
levels(cast$Order) <- c(levels(cast$Order), "Unknown") #include a new category
cast$Order[is.na(cast$Order)] <- "Unknown"
cast$Order

rownames(cast) <- cast$Order
cast <- cast[,-1]

genre_table <- as.matrix(cast)

s_num<-dim(genre_table)[2]
asv.perc<-matrix(rep("NA",times=(dim(genre_table)[1]*s_num)),nrow=dim(genre_table)[1],ncol=s_num)
rownames(asv.perc)<-rownames(genre_table)
colnames(asv.perc)=colnames(genre_table) 
totals<-colSums(genre_table)

for(s in c(1:s_num)){
  vec<-(genre_table[,s]/totals[s])*100
  asv.perc[,s]<-vec
  rm(vec)
}

df_m <- melt(asv.perc)
colnames(df_m) <- c("Order","moduleColor","value")

# Keep the 30 most abundant families
cast$tot <- rowSums(cast)
cast_Order<- cast[order(cast$tot, decreasing=TRUE),]
cast <- cast_Order[c(1:30),]
df_m2 <- df_m[df_m$Order %in% row.names(cast),]

df_m2[df_m2 == 0] <- NA
df_m2 <- na.omit(df_m2)

df_m2$moduleColor <- factor(df_m2$moduleColor)
df_m2$Order <- factor(df_m2$Order)
df_m2$value <- as.numeric(as.character(df_m2$value))

df_m2$color<-col2hex(df_m2$moduleColor)

treemap(df_m2, index=c("moduleColor","Order"), vSize="value", vColor="color" , type="color",
 
      fontsize.labels=c(16,13),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
      fontcolor.labels=c("white","white"),    # Color of labels
      fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
      bg.labels=c("transparent"),              # Background color of labels
      align.labels=list(
          c("left", "top"), 
          c("right", "bottom")),               # Where to place labels in the rectangle?
      overlap.labels=0.5,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
      inflate.labels=F,                        # If true, labels are bigger when rectangle is bigger.
 	  border.col=c("black","white"),
  )

#####################################
############## sPLS ASV  ############
#####################################
ASVs_CLR<-t(as.data.frame(otu_table(Benthic_matrice2_clr)))
tax<-tax_table(Benthic_matrice2_clr)
taxonomy<-as.data.frame(tax)
taxonomy$subnet<-dynamicColors

new_color<-as.character(taxonomy$subnet)
new_color[which(new_color=="orange")] <- "#66c2a5"
new_color[which(new_color=="blue")] <- "#fc8d62"
new_color[which(new_color=="turquoise")] <- "#e78ac3"
new_color[which(new_color=="grey")] <- "#a6d854"
new_color[which(new_color=="brown")] <- "#e5c494"
new_color[which(new_color=="red")] <- "#b3b3b3"


my.col <- c("Subnetwork_1" = "#66c2a5",
			"Subnetwork_2" = "#a6d854",
			"Subnetwork_3" = "#b3b3b3", 
			"Subnetwork_4" = "#e5c494",
			"Subnetwork_5" = "#e78ac3", 
			"Subnetwork_6" = "#fc8d62")
taxonomy$subnet<-new_color

datTraits2<-datTraits
env<-datTraits2[,c(1,7:17,5,4,3,2,29:31)]
rownames(env)<- c(1:dim(env)[1])
rownames(ASVs_CLR)<-c(1:dim(ASVs_CLR)[1])
all_bis<-merge(env, ASVs_CLR, by=0)
all_bis$Row.names<-NULL
all_bis<-na.omit(all_bis)

Y<-all_bis[,c(1:length(env))]
X<-all_bis[,-c(1:length(env))]

pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.gene)
pca.env <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.env)

#PLS
liver.pls <- pls(X, Y, ncomp=4, mode = "regression")
pls_vip <- mixOmics::vip(liver.pls)
########################################################################################
# sPLS
spls_ASV <- mixOmics::spls(X, Y, ncomp=4, keepX = c(20,620,570,270), keepY= c(3,3,3,17), mode = "regression")

##### Figure ######
cim.spls_ASV <- mixOmics::cim(spls_ASV, 
	comp = 1:4, 
	xlab = "Environmental factors", 
	ylab = "Species", 
	margins = c(7,20), 
	dist.method = c("correlation", "correlation"), 
	clust.method=c("complete", "complete"), 
	mapping = "XY",
	row.names= taxonomy$Family,
	row.sideColors = new_color, 
	legend = list(legend=names(my.col), col=my.col),
	color=brewer.pal(n = 9, name = "Greys"))
#######################################################################################