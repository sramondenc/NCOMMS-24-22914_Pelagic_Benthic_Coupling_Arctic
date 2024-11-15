# Package
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
library(gephi)
library(treemap)
library(gplots)
library("stringr")
library(ggcorrplot)
library(RColorBrewer)

#Function
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
##################################
# Import Environmental variables #
##################################

all<-read.csv("Pelagic_environmental_data.csv", sep=";", dec=',') # all corrected
rownames(all)<-all[,c(1)]
all<-all[,-c(1)] 
 
#######################
### Import Omics data #
#######################

# Sample name table
samples_df<-read.csv("Pelagic_samples.csv", sep=";", h=T, dec='.' )
samples_df$Date_start<- as.Date(samples_df$Date_start,format = '%m/%d/%Y')
samples_df$julian<- as.numeric(format(samples_df$Date_start,"%j"))

# Omic raw data
ASV_Planktonic<-read.csv("18S_Pelagic_data.csv", sep=";", h=T, dec='.')
ASVs<-ASV_Planktonic[,c(1)]
ASV_mat<-ASV_Planktonic[,c(20:102)]
colnames(ASV_mat)<-samples_df$Traps
ASV_mat<-cbind(ASVs,ASV_mat)

#Taxonomy table
tax_mat<-ASV_Planktonic[,c(103:110)]
tax_mat<-cbind(ASVs,tax_mat)

#################################################################################
############################## Build OTU data base for analysis #################
#################################################################################
row.names(ASV_mat) <- ASV_mat$ASVs
ASV_mat <- ASV_mat[,-c(1)] 
row.names(tax_mat) <- tax_mat$ASVs
tax_mat <- tax_mat[,-c(1)] 
row.names(samples_df) <- samples_df$Traps

samples = sample_data(samples_df) 
ASV_mat <- as.matrix(ASV_mat)
tax_mat <- as.matrix(tax_mat)

ASV = otu_table(ASV_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
Omics_pelagic <- phyloseq(ASV, TAX, samples)

#################################################################################
############################## Summarize Data Base ##############################
#################################################################################

# Check if any ASV are not present in any samples. If True removed
any(taxa_sums(Omics_pelagic)==0)
Omics_pelagic_a <- prune_taxa(taxa_sums(Omics_pelagic) > 0, Omics_pelagic)
any(taxa_sums(Omics_pelagic_a)==0)

# Original abundance ASV versus removed
ntaxa(Omics_pelagic)-ntaxa(Omics_pelagic_a)

#################################################################################
############################## Remove Phylum unclassified #######################
#################################################################################

### Composition
Omics_pelagic_a.com<-subset_taxa(Omics_pelagic_a, !is.na(Supergroup))

taxic <- as.data.frame(Omics_pelagic_a.com@tax_table) 
taxic$ID <- paste(as.matrix(taxic)[cbind(1:nrow(taxic), max.col(!is.na(taxic),"last")-1)],"",as.matrix(taxic)[cbind(1:nrow(taxic), max.col(!is.na(taxic),"last"))])

# Add the ASV ids from OTU table into the taxa table at the end.
taxic$ASV <- rownames(taxic) 

# You can see that we now have extra taxonomy levels.
colnames(taxic)

# convert it into a matrix.
taxmat <- as.matrix(taxic)

# convert into phyloseq compaitble file.
new.tax <- tax_table(taxmat)  

# incroporate into phyloseq Object
tax_table(Omics_pelagic_a.com) <- new.tax 

#############################################################
############################## Remove Singleton #############
#############################################################
summarize_phyloseq(Omics_pelagic_a.com) # We have 134 singleton
Pelagic_matrice <- prune_samples(sample_sums(Omics_pelagic_a.com) > 4500, Omics_pelagic_a.com)

# Filter ASVs that have 50 sequences in at least 3 samples
subset_3samples = genefilter_sample(Pelagic_matrice, filterfun_sample(function(x) x >= 50), A=3)

################## Remove biais
subset_3samples[subset_3samples="ASV_45"]=FALSE		# Capra_hircus
subset_3samples[subset_3samples="ASV_P_520"]=FALSE		# Capra_hircus
subset_3samples[subset_3samples="ASV_P_854"]=FALSE		# Capra_hircus
subset_3samples[subset_3samples="ASV_P_1086"]=FALSE		# Capra_hircus
subset_3samples[subset_3samples="ASV_P_1456"]=FALSE		# Capra_hircus
subset_3samples[subset_3samples="ASV_P_1591"]=FALSE		# Capra_hircus
subset_3samples[subset_3samples="ASV_P_2871"]=FALSE		# Capra_hircus

subset_3samples[subset_3samples="ASV_537"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_578"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_603"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_659"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_1179"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_1718"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_1765"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_489"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_1798"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_1895"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_2534"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_2583"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_2966"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_3729"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_5117"]=FALSE		# Homo_sapiens
subset_3samples[subset_3samples="ASV_P_5973"]=FALSE		# Homo_sapiens

subset_3samples[subset_3samples="ASV_P_315"]=FALSE		# Musa_basjoo

subset_3samples[subset_3samples="ASV_P_486"]=FALSE		# Pinus_taeda
subset_3samples[subset_3samples="ASV_P_1474"]=FALSE		# Pinus_taeda
subset_3samples[subset_3samples="ASV_P_4907"]=FALSE		# Pinus_taeda

subset_3samples[subset_3samples="ASV_P_534"]=FALSE		# Phaseolus_acutifolius
subset_3samples[subset_3samples="ASV_P_4531"]=FALSE		# Phaseolus_acutifolius
##################

Pelagic_matrice2 = prune_taxa(subset_3samples, Pelagic_matrice)
ntaxa(Pelagic_matrice2)
nsamples(Pelagic_matrice2)

# Print the filtered matrix
asvs <- otu_table(Pelagic_matrice2)
asvs <- as.data.frame(asvs)

#################################################################################
######################### CLR Transformation ####################################
#################################################################################
# Counting_Table
Pelagic_matrice2_clr<-microbiome::transform(Pelagic_matrice2,'clr')
datASVClr<-t(otu_table(Pelagic_matrice2_clr))

datASVClr_2<-t(otu_table(Pelagic_matrice2))
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

# define numbers of samples
nSamples <- nrow(datASVClr)

moduleColors = merge_dynamic_MEDs $colors;

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datASVClr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Select Environmental that macth with our omics data
datTraits<-all[,c(3:7,9:32)]
datTraits<-datTraits[c(match(as.character(sample_data(Pelagic_matrice2)$Date_start),rownames(datTraits))),]

moduleTraitCor <- cor(MEs, datTraits, use = "pairwise.complete.obs")
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

par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = c("Copepod","Amphipod","Ostracod","Pteropod","Cheatognaths","ice_concentration","temp_air","v_winds","u_winds","v_ice","u_ice","Press_atm","SST","SSS","u_currents","v_currents","ice_distance","chla_sat","AMO","A0","NAO","Seasonal_var","total.flux","CaCO3","DSi.PSi","POC","PON","d15N","d13C"),
               yLabels = names(MEs), 
               ySymbols = names(MEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               cex.text = 0.6,
			   cex.lab.x=1,
			   cex.lab.y = 0.001,
               zlim = c(-1,1), xLabelsAngle=90, 
			   plotLegend=0,horizontalSeparator.lty=1)

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
POC_export<- as.data.frame(datTraits$POC)
names(POC_export) <- "POC_export"

# Calculate the correlations between modules
geneModuleMembership <- as.data.frame(WGCNA::cor(datASVClr, MEs, use = "p"))
#p-values
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# Calculate the correlation for the trait
geneTraitSignificance <- as.data.frame(WGCNA::cor(datASVClr, POC_export, use = "p"))
# p-values
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nSamples))

names(geneTraitSignificance) <- paste("GS.", names(POC_export), sep = "")
names(GSPvalue) <- paste("p.GS.", names(POC_export), sep = "")

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
  #print(column)
  
  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors == module
  vecASVnames = rownames(geneTraitSignificance)
  print(paste("There are ", length(vecASVnames[moduleGenes]), " ASVs in the ", module, " module.", sep = ""))
  print(vecASVnames[moduleGenes])
  
  # Make the plot
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                 abs(geneTraitSignificance[moduleGenes, 1]),
                 xlab = paste("Module Membership in", module, "module"),
                 ylab = "Gene significance for POC flux",
                 main = paste("Module membership vs. gene significnace \n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}    

select_POC=names(geneModuleMembership)[c(4,9:10)] # Positively & Negatively correlated to POC export
module <- select_POC
print(module)
moduleGenes = is.finite(match(moduleColors, substring(module,3)))
vecASVnames = rownames(geneTraitSignificance)
species_module<-vecASVnames[moduleGenes]
color_<-moduleColors[moduleGenes]
taxon<-tax_table(Pelagic_matrice2_clr)[match(species_module,rownames(tax_table(Pelagic_matrice2_clr)))]

# Replot Gene Significance vs module membership for modules interests (here POC export)
par(mfrow = c(1,3))  
for (i in select_POC) { 
  # Pull out the module we're working on
  module <- substring(i,3)
  if(module=='turquoise'){
	  new_color_dot<-"#fdae6b"
  } else if(module=='black'){
      new_color_dot<-"#3182bd"
  } else if(module=='orange'){
	  new_color_dot<-"#756bb1"}
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
                 ylab = "Gene significance for POC flux",
                 main = paste("Module membership vs. gene significnace \n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = new_color_dot)
}  


###################################################
## Modules relative abundance across the samples ##
###################################################

Pelagic_matrice2.rel<- microbiome::transform(Pelagic_matrice2, "compositional")
taxmat <- tax_table(Pelagic_matrice2.rel)[,c(1:8)]
datExpr<-t(otu_table(Pelagic_matrice2.rel))

annot = data.frame(taxmat)
dim(annot)
probes = colnames(datExpr)
probes2annot = match(probes, rownames(annot))

moduleColors = merge_dynamic_MEDs $colors;

ASVInfo = data.frame(ASVs = probes,
                       Phylum = annot$Division[probes2annot],
                       Class = annot$Class[probes2annot],
                       Order = annot$Order[probes2annot],
                       Family = annot$Family[probes2annot],
                       Genus = annot$Genus[probes2annot],
                       moduleColor = moduleColors)
table_module <- ASVInfo[,c(1:7)]

table_perc <- t(datExpr)
table_perc <- as.data.frame(table_perc)
table_perc$tot <- rowSums(table_perc)
table_perc<- cbind(ASVs = rownames(table_perc), table_perc)
abond <- table_perc[,c(1,83)]
abond <- as.data.frame(abond)

table_perc <- table_perc[,-c(1,83)]
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
module_com_7 <- perc_module[perc_module$moduleColor == com_imp[7],]
module_com_8 <- perc_module[perc_module$moduleColor == com_imp[8],]
module_com_9 <- perc_module[perc_module$moduleColor == com_imp[9],]
module_com_10 <- perc_module[perc_module$moduleColor == com_imp[10],]
module_com_11 <- perc_module[perc_module$moduleColor == com_imp[11],]



Contr_ASVs_modl <- data.frame(com_1 = apply(table_perc[,colnames(table_perc) %in% module_com_1$ASVs], 1, sum),
                       com_2 = apply(table_perc[,colnames(table_perc) %in% module_com_2$ASVs], 1, sum),
                       com_3 = apply(table_perc[,colnames(table_perc) %in% module_com_3$ASVs], 1, sum),
                       com_4 = apply(table_perc[,colnames(table_perc) %in% module_com_4$ASVs], 1, sum),
                       com_5 = apply(table_perc[,colnames(table_perc) %in% module_com_5$ASVs], 1, sum),
                       com_6 = apply(table_perc[,colnames(table_perc) %in% module_com_6$ASVs], 1, sum),
                       com_7 = apply(table_perc[,colnames(table_perc) %in% module_com_7$ASVs], 1, sum),
                       com_8 = apply(table_perc[,colnames(table_perc) %in% module_com_8$ASVs], 1, sum),
                       com_9 = apply(table_perc[,colnames(table_perc) %in% module_com_9$ASVs], 1, sum),
                       com_10 = apply(table_perc[,colnames(table_perc) %in% module_com_10$ASVs], 1, sum),
                       com_11 = apply(table_perc[,colnames(table_perc) %in% module_com_11$ASVs], 1, sum))

colnames(Contr_ASVs_modl)<-com_imp

info_sample<-sample_data(Pelagic_matrice2.rel)

df <- cbind(samples = rownames(Contr_ASVs_modl), Contr_ASVs_modl)
df <- cbind(Year = info_sample$Year, df)
df <- cbind(Season = info_sample$Season, df)

df_m <- melt(df, id.vars=c("samples", "Year","Season"))

df_m$Season <- factor(df_m$Season, levels = c("Winter","Spring","Summer","Faller"))
df_m$Year <- factor(df_m$Year, levels = c("2000","2001","2003","2004","2005","2006","2007","2008","2009","2011","2012"))

df_m <- na.omit(df_m)

mcast <- dcast(df_m, Season + Year ~ variable, mean)
df_m2 <- melt(mcast, id.vars=c("Season","Year"))

df_m2 <- na.omit(df_m2)

# Change the default color
df_m2$color<-col2hex(df_m2$variable)
new_color<-as.character(df_m2$variable)
new_color[which(new_color=="black")] <- "#3182bd"
new_color[which(new_color=="blue")] <- "#bcbddc"
new_color[which(new_color=="brown")] <- "#e6550d"
new_color[which(new_color=="greenyellow")] <- "#d9d9d9"
new_color[which(new_color=="grey")] <- "#31a354"
new_color[which(new_color=="magenta")] <- "#a1d99b"
new_color[which(new_color=="orange")] <- "#756bb1"
new_color[which(new_color=="pink")] <- "#bdbdbd"
new_color[which(new_color=="purple")] <- "#636363"
new_color[which(new_color=="red")] <- "#9ecae1"
new_color[which(new_color=="turquoise")] <- "#fdae6b"

df_m2$new_color<-new_color
ggplot(df_m2, aes(x=Year, y = value, fill=variable)) +
    geom_bar(stat = "identity", position = "fill", show.legend = FALSE) +
    facet_grid(.~Season, space="free_x") +
	scale_fill_manual(values=unique(df_m2$new_color))+
	theme_bw()+
    theme(strip.text = element_text(size = 15, hjust=0), strip.background = element_blank(), axis.text.x = element_text(size = 10 , angle = 40, hjust=1), 
	axis.text.y = element_text(size = 10),
	legend.position="none")+
	ylab("Subnetwork contribution")

#############################################################
### Relative abundance of the main families across module ###
#############################################################

table <- perc_module[,c(4,7,8)]
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
cast_order <- cast[order(cast$tot, decreasing=TRUE),]
cast <- cast_order[c(1:30),]
df_m2 <- df_m[df_m$Order %in% row.names(cast),]

df_m2[df_m2 == 0] <- NA
df_m2 <- na.omit(df_m2)

df_m2$moduleColor <- factor(df_m2$moduleColor)
df_m2$Order <- factor(df_m2$Order)
df_m2$value <- as.numeric(as.character(df_m2$value))

ggplot(df_m2, aes(x=moduleColor , y = Order), size=as.numeric(value), color=as.numeric(value)) +
  geom_point(aes(size=as.numeric(value)),color="#757272", position = "identity")+ 
  scale_size(limits = c(0,100), range=c(1,10), breaks=c(10,20,30,40,50,60,70,80,90)) +
  scale_color_gradient2(limits=c(0,100),low="#00ffff", mid="#383838",high="#f26846", midpoint = 45, na.value=NA)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df_m2$color<-col2hex(df_m2$moduleColor)

treemap(df_m2, index=c("moduleColor","Order"), vSize="value", vColor="color" , type="color",
 
      fontsize.labels=c(15,11),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
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
ASVs_CLR<-t(as.data.frame(otu_table(Pelagic_matrice2_clr)))
tax<-tax_table(Pelagic_matrice2_clr)
taxonomy<-as.data.frame(tax)
taxonomy$subnet<-dynamicColors

new_color<-as.character(taxonomy$subnet)
new_color[which(new_color=="black")] <- "#3182bd"
new_color[which(new_color=="blue")] <- "#bcbddc"
new_color[which(new_color=="brown")] <- "#e6550d"
new_color[which(new_color=="greenyellow")] <- "#d9d9d9"
new_color[which(new_color=="grey")] <- "#31a354"
new_color[which(new_color=="magenta")] <- "#a1d99b"
new_color[which(new_color=="orange")] <- "#756bb1"
new_color[which(new_color=="pink")] <- "#bdbdbd"
new_color[which(new_color=="purple")] <- "#636363"
new_color[which(new_color=="red")] <- "#9ecae1"
new_color[which(new_color=="turquoise")] <- "#fdae6b"

my.col <- c("Subnetwork_1" = "#31a354",
			"Subnetwork_2" = "#636363",
			"Subnetwork_3" = "#756bb1", 
			"Subnetwork_4" = "#9ecae1",
			"Subnetwork_5" = "#a1d99b", 
			"Subnetwork_6" = "#bcbddc",
			"Subnetwork_7" = "#bdbdbd",
			"Subnetwork_8" = "#d9d9d9",
			"Subnetwork_9" = "#e6550d",
			"Subnetwork_10" = "#fdae6b",
			"Subnetwork_11" = "#3182bd")
taxonomy$subnet<-new_color

datTraits2<-datTraits[,-c(18)]
env<-datTraits2
colnames(env)<-c("Copepod","Amphipod","Ostracod","Pteropod","Cheatognaths","[ice]","temp_air","v_winds","u_winds","v_ice","u_ice","pressure_atm","SST", "SSS","u_currents","v_currents","ice_distance","AMO","A0","NAO","Seasonal_var","total flux","CaCO3","Si","POC","PON","d15N","d13C")
rownames(env)<- c(1:dim(env)[1])
rownames(ASVs_CLR)<-c(1:dim(ASVs_CLR)[1])
all_bis <-merge(env, ASVs_CLR, by=0)
all_bis$Row.names<-NULL
all_bis2<-na.omit(all_bis)

Y<-all_bis2[,c(1:length(env))]
X<-all_bis2[,-c(1:length(env))]

pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.gene)
pca.env <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.env)

liver.pls <- pls(X, Y, ncomp=3, mode = "regression")
pls_vip <- mixOmics::vip(liver.pls)

#######################################################################################
# sPLS
spls_ASV <- mixOmics::spls(X, Y, ncomp=3, keepX= c(510,475,370), keepY= c(3,3,10), mode = "regression")

##### Figure 1 ######
cim.spls_ASV <- mixOmics::cim(spls_ASV, 
	comp = 1:3, 
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