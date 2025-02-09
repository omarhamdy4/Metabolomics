# Metabolomics Analysis Using R
***
> #### This notebook demonstrates a step-by-step workflow for analyzing metabolomics data using limma. The analysis includes preprocessing the expression matrix, handling metadata, Biomarker discovery for differentially expressed metabolites (DEMs), and Plotting important figures such as volcano plots, Heatmaps, PCA, OPLS-DA, PCA, and Roc curve.
## The Workflow goes as follows:
![Metabolomics Analysis Pipeline script workflow template](https://github.com/user-attachments/assets/fbb64407-a86e-49b9-ad24-61659cf1a2e9)
***
## 1- Library loading

> ------------------------------------------------------------------------

```{r setup, include=FALSE}
library(magrittr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(limma)
library(vegan)
library(cluster)
library(factoextra)
library(gridExtra)
library(PerformanceAnalytics)
library(corrplot)
library(Hmisc)
library(RColorBrewer)
library(impute)
library(pathview)
library(glmnet)
library(ggfortify)
library(ComplexHeatmap)
library(pROC)
```

## 2- Reading Data

> ------------------------------------------------------------------------

```{r}
setwd("E:/1.Fresh Grad/02_EgComBio2023/MODA Metabolomics")
negative_hilic_metabolite <- read_csv("m_MTBLS8920_LC-MS_negative_hilic_metabolite_profiling_v2_maf.csv")
negative_hilic_metabolite$charge <- "NEG"
positive_hilic_metabolite<- read_csv("m_MTBLS8920_LC-MS_positive_hilic_metabolite_profiling_v2_maf.csv")
positive_hilic_metabolite$charge <- "POS"

Metadata <- read_csv("metadata.csv")
Metadata <- Metadata[1:57,]
```

## 3- Pre-processing

> ------------------------------------------------------------------------

```{r Pre-processing intensity matrix}
All_hilic_metabolite <- rbind(positive_hilic_metabolite,negative_hilic_metabolite)
All_metabolites <- All_hilic_metabolite %>% select(c(metabolite_identification,22:78))
metabolite.name <- unlist(All_metabolites$metabolite_identification)
All_metabolites[All_metabolites==0.0|0.00|0.000] <- NA

All_metabolites <- All_metabolites[,-1]
All_metabolites_mat <- apply(All_metabolites, 2, as.numeric)
row.names(All_metabolites_mat) <-metabolite.name
class(All_metabolites_mat)

```

## 4- Quality Control

> ------------------------------------------------------------------------

```{r QC}
#remove metaboloites with zero variance
varCol=apply(All_metabolites_mat, 1, var, na.rm = T)
constCol <- (varCol == 0 | is.na(varCol))
sum(constCol)
#All_metabolites_mat <- All_metabolites_mat[,!constCol]

#Compute missinigness rate in your metaboloites
round(sum(is.na(All_metabolites_mat))/(nrow(All_metabolites_mat) * ncol(All_metabolites_mat))*100,1)
h=hist((apply(is.na(All_metabolites_mat), 1, sum)/nrow(All_metabolites_mat) ) *100,breaks=10,main="",xlab="Percentage of missingness")
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))
plot(density( (apply(is.na(All_metabolites_mat), 1, sum)/nrow(All_metabolites_mat) )*100),xlab="percentage of missingness",main='')

#Impute missing values
Imputed_All_metabolites_mat=impute.knn(All_metabolites_mat,k=10)$data

#Normalization(Log transformation)
Logged_Imputed_All_metabolites_mat <- log2(Imputed_All_metabolites_mat + 1)
head(Logged_Imputed_All_metabolites_mat)

# Scaling the data
Scaled_Logged_Imputed_All_metabolites_mat=scale(Logged_Imputed_All_metabolites_mat,center = TRUE, scale = TRUE)

#save the normalized data
write.csv(Scaled_Logged_Imputed_All_metabolites_mat,"Normalized data.csv")
```
#### Visualize QC Results
```{r Visualize QC Results}

#Visualize QC Results [Log transformation/Scaling] (Subset of Metabolites)
options(repr.plot.width=10,repr.plot.height=8)
par(mar = c(8,5,2,2),mfrow=c(1,3))
boxplot(Imputed_All_metabolites_mat[,1:20], main="Before log2" ,horizontal=T, names=rownames(Imputed_All_metabolites_mat)[1:20],las=2,
       col = "lightgreen")
boxplot(Logged_Imputed_All_metabolites_mat[,1:20], main="After log2" ,horizontal=T, names=rownames(Logged_Imputed_All_metabolites_mat)[1:20],
        las=2,col = "lightgreen")
boxplot(Scaled_Logged_Imputed_All_metabolites_mat[,1:20], main="After log2 + Scaled " ,horizontal=T,
        names=rownames(Scaled_Logged_Imputed_All_metabolites_mat)[1:20],
        las=2,col = "lightgreen")
#Visualize QC Results [Log transformation/Scaling] (Subset of Samples)
options(repr.plot.width=10,repr.plot.height=8)
par(mar = c(8,5,2,2),mfrow=c(1,3))
boxplot(Imputed_All_metabolites_mat, main="Before log2" ,horizontal=T,las=2,col = "lightgreen")    
boxplot(Logged_Imputed_All_metabolites_mat ,main="After log2" ,horizontal=T,las=2,col = "lightgreen")
boxplot(Scaled_Logged_Imputed_All_metabolites_mat, main="After log2 + Scaled" ,horizontal=T,las=2,col = "lightgreen")

```

## 5- Multivariate Analysis

> ------------------------------------------------------------------------
#### Principal Component Analysis (PCA)
PCA is an unsupervised dimensionality reduction technique used to simplify complex datasets by transforming them into a lower-dimensional space. It identifies the directions (principal components) that maximize variance in the data, allowing for visualization, noise reduction, and feature extraction. PCA is widely used in exploratory data analysis, pattern recognition, and data compression.

#### Orthogonal Projections to Latent Structures Discriminant Analysis (OPLS-DA)
OPLS-DA is a supervised multivariate analysis method used to separate classes in a dataset. It extends PLS-DA by removing variation orthogonal to the class separation, improving interpretability. OPLS-DA is commonly used in metabolomics, genomics, and other fields for biomarker discovery and classification tasks. It provides clearer insights into class-discriminatory features compared to traditional PLS-DA.

### 5- (a) PCA

```{r PCA}
options(repr.plot.width=8,repr.plot.height=6)

df_pca <- prcomp(t(Scaled_Logged_Imputed_All_metabolites_mat))
df_out <- as.data.frame(df_pca$x)
ggplot(df_out,aes(x=PC1,y=PC2,color=Metadata$`Factor Value[Disease]` ))+
geom_point()+ggtitle("")+labs(color='')+
geom_point(size=8,alpha=0.5)+ #Size and alpha just for fun
stat_ellipse(geom="polygon", level=0.95, alpha=0.2)+ #adding a stat
#ggplot2::geom_polygon(data=df_out,aes(x=PC1,y=PC2,fill = Metadata$`Factor Value[Disease]`), alpha = 0.3)+
geom_text(label=Metadata$`Source Name`,nudge_x=0.35, nudge_y=0.1,check_overlap=T)+
theme(  plot.title = element_text(hjust = 0.5,size=15,face = "bold"),
        axis.text.x = element_text( size = 15, angle = 45, hjust = .5, vjust = 0.5, face = "plain"),
        axis.text.y = element_text( size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text( size = 15, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text( size = 15, angle = 90, hjust = .5, vjust = .5, face = "bold"),
        #legend.title=element_text(size=20),
        legend.title=element_blank(), # remove legend title name
        legend.text = element_text(size=15,face="plain"),
        strip.text = element_text(size = 15,face="plain") ,
        legend.position="right",

        # for transparent background
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.line = element_line(colour = "black") # adding a black line for x and y axis
)+xlab(paste0("PC 1 (", round(df_pca$sdev[1],1),"%)"))+
    ylab(paste0("PC 2 (", round(df_pca$sdev[2],1),"%)"))


```

### 5- (b) OPLS-DA
```{r OPLS-DA}
#OPLS-DA relies on a projection of X data as PCA does, but here we are rephrasing the question from one of corresponding to a maximum variance model to become one of corresponding to a maximum separation model. This separation is then possible because the projection of the X data will be guided by known class information. We will supply class information to the OPLS algorithm and thereby we go from an unsupervised modeling approach (PCA) to supervised modeling approach (OPLS).
# OPLS-DA is an excellent tool to find “What’s the difference” between two groups (such as Good and Bad product).  The OPLS-DA model will indicate which are the driving forces among the variables and we can then make score plots to visualize the differences if they exist.
#The Variable Importance in Projection (VIP), which reflects both the loading weights for each component and the variability of the response explained by this component (Pinto, Trygg, and Gottfries 2012; Mehmood et al. 2012), can be used for feature selection

#R2Y represents the goodness of fit of the model. It quantifies how well the model explains the variance in the dependent variable(s) (Y) based on the independent variable(s) (X).
#R2Y ranges from 0 to 1, where: 0 means the model explains none of the variance in Y. 1 means the model explains all the variance in Y. A higher R2Y indicates a better fit of the model to the data.
#Q2Y represents the predictive ability of the model. It is calculated using cross-validation and measures how well the model can predict new data. Q2Y also ranges from 0 to 1, where: 0 means the model has no predictive ability. 1 means the model has perfect predictive ability. A higher Q2Y indicates better predictive performance.
#The VIP score quantifies the contribution of each independent variable (X) to the model's ability to discriminate between classes (Y).Range: VIP scores are typically non-negative, with higher values indicating greater importance.A common threshold is VIP > 1.0, which is often used to identify significant variables.Variables with VIP scores < 1.0 are considered less important.

library(ropls)  #https://www.bioconductor.org/packages/release/bioc/vignettes/ropls/inst/doc/ropls-vignette.html#1_The_ropls_package
plot.oplsda <- opls(t(Scaled_Logged_Imputed_All_metabolites_mat), Metadata$`Factor Value[Disease]`,predI = 1, orthoI = NA, subset = "odd")
VIPscores <- data.frame(
    gene_symbol = names(plot.oplsda@vipVn),  # Extract metabolite names
    VIP_Score = plot.oplsda@vipVn,          # Extract VIP scores
    row.names = NULL                 # Ensure no row names are assigned
)


```
## 6- Biomarker Discovery (*Differential Expression Analysis [Limma]*)

> ------------------------------------------------------------------------

```{r Diff. Exp. Analysis}
table(Metadata$`Factor Value[Disease]`)
#type = as.character(Metadata$`Factor Value[Disease]`)

groups = factor(Metadata$`Factor Value[Disease]`, levels = c("Breast_Cancer", "Healthy_Control"))
design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
contrast<-makeContrasts(Breast_Cancer-Healthy_Control,levels=design) # pos logfc means in cancer is more than normal
fit <- lmFit(as.matrix(Scaled_Logged_Imputed_All_metabolites_mat), design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
DEGs=topTable(fit2, adjust.method='fdr', number=999999999,p.value=1,coef = 1)
dim(DEGs)

topres <- data.frame(DEGs)
topres$gene_symbol <- topres$ID
topres$diffexpressed <- "NA"
topres$diffexpressed[topres$logFC > 0.5 & topres$P.Value < 0.05] <- "UP"
topres$diffexpressed[topres$logFC < -0.5 & topres$P.Value < 0.05] <- "DOWN"
topres_VIP <- inner_join(topres,VIPscores,by="gene_symbol")

#These are the top selected panel of metabolites depending on (logFC | p-value | VIP score)
DEMs <- topres_VIP[topres_VIP$logFC >= 0.5 & topres_VIP$P.Value < 0.05 |topres_VIP$logFC <= -0.5 & topres_VIP$P.Value < 0.05 &topres_VIP$VIP_Score >1,]
names_DEMs<- unlist(DEMs$ID)
write.csv(names_DEMs,"names_DEMs.csv")

```

## 7- Downstream Plots

> ------------------------------------------------------------------------

### 7- (a) K-means Clustering

```{r K-means Clustering}
#removing duplicated rownames
dupliTEMPO <-duplicated(rownames(Scaled_Logged_Imputed_All_metabolites_mat))
tempo <- Scaled_Logged_Imputed_All_metabolites_mat[!dupliTEMPO,]
#clustering
kmeans2 <- kmeans((tempo), centers = 2, nstart = 25)
factoextra::fviz_cluster(kmeans2, data = tempo, ellipse = T,labelsize = 1,ggtheme = theme_minimal(), main ="K-means Clustering",  outlier.color = "black",shape = NULL,outlier.shape = 19 )

```

### 7- (b) Ward Hierarchical Clustering

```{r ward Hierarchical Clustering}
# NICE resource for f=dendograms https://www.gastonsanchez.com/visually-enforced/how-to/2012/10/03/Dendrograms/
# Ward Hierarchical Clustering
d <- dist(t(Scaled_Logged_Imputed_All_metabolites_mat), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
library("ggdendro")
ggdendrogram(fit, rotate = TRUE, size = 4, theme_dendro = FALSE,ggtheme = theme_minimal())

```

### 7- (c) Heatmap

```{r Heatmap}
deg_100 <- DEGs[order(DEGs$P.Value, -abs(DEGs$logFC)),]   #sort to get top 100
deg_100 <- DEGs[1:100,]                                         #select only them
deg_100_exp <- as.matrix(Logged_Imputed_All_metabolites_mat[deg_100$ID,])                     #get their normalized expression values

# Create column annotations for the heatmap
sam_condition <- Metadata$`Factor Value[Disease]`
column_annot <- HeatmapAnnotation(
  Condition = sam_condition,
  col = list(Condition = c("Breast_Cancer" = "purple", "Healthy_Control" = "violet")))

# Generate the heatmap
Heatmap(
  matrix = deg_100_exp,
  top_annotation = column_annot,
  row_title = "Top 100 DEGs",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 4),
  column_names_rot = 45,
  column_names_centered = TRUE)
# The heatmap resulted shows no significant signal between disease and control

```

### 7- (d) Volcano plot

```{r}
res_df <- data.frame(DEGs)
res_df$gene_symbol <- res_df$ID

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)))
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (logfc respectively positive or negative)
res_df$diffexpressed <- "NA"
# if logFoldchange > 1 and pvalue < 0.05, set as "UP"
res_df$diffexpressed[res_df$logFC > 0.5 & res_df$P.Value < 0.05] <- "UP"
# if logFoldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_df$diffexpressed[res_df$logFC < -0.5 & res_df$P.Value < 0.05] <- "DOWN"
# Create a new column "delabel" to de, that will contain the name of the top differentially expressed genes (NA in case they are not) 
#*fOR dOWN GENES*
res_df$delabel <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$logFC), "gene_symbol"], 6), res_df$gene_symbol, NA)
#*fOR UPP GENES*
res_df$delabel2 <- res_df$delabel
res_df$delabel2 <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$logFC,decreasing = T), "gene_symbol"], 4), res_df$gene_symbol, NA)

# Now plot it
ggplot(res_df, aes(x = logFC, y = -log(P.Value), col = res_df$diffexpressed)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  labs(
    title = "Volcano Plot",
    x = "log Fold change",
    y = "-log10 adjusted Pvalue"
  ) +
geom_vline(xintercept = c(0, -0), col="black",linetype = "dashed") +
geom_hline(yintercept = -log10(0.05), col="black",linetype = "dashed")+
geom_point(size = 2)+
scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
coord_cartesian(ylim = c(0,20), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
labs(color = 'Regulation', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
scale_x_continuous(breaks = seq(-5, 5, 1))+ # to customise the breaks in the x axis
ggtitle('Breast Cancer vs Control')+  # Plot title 
geom_text_repel(label=res_df$delabel2,max.overlaps = Inf)+
geom_text_repel(label=res_df$delabel,max.overlaps = Inf) # To show all labels 

```

### 7- (e) ROC curve

```{r ROC curve}

#ROC Curve
x <- select(Metadata, "Factor Value[Disease]")
colnames(x) <- "disease_state"
x <- mutate(x,
            value=case_when(
              disease_state =="Breast_Cancer" ~ "1",
              disease_state =="Healthy_Control" ~ "0",
              TRUE ~ "X" #this replaces the non-mentioned other values to be x
            ))
x$value <- as.numeric(x$value)
Nx2 <- t(Scaled_Logged_Imputed_All_metabolites_mat)
CallMe <- function(q) {
  x$ex <- Nx2[,q]
  glm.fit=glm(value~ex,data=x, family=binomial)
  library(pROC)
  par(pty="s")
  roc(x$value, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE,
      xlab="100-Specificity",ylab="Sensitivity", col="#377eb8",lwd=4, print.auc=TRUE, print.auc.x=45)    
  legend("bottomright" ,legend=c(q), col=c("#377eb8"), lwd=4)
}
#ROC curves for all the Differentially expressed Metabolites
for (i in names_DEMs)    CallMe(i)

```



## 8-Pathways enrichment analysis [over-representation analysis]
> ------------------------------------------------------------------------
#### Before the next block of code
- After saving the top metabolites we had in the file 'names_DEMs.csv'
- IDs were converted to KEGG IDs using MetaboAnalyst ID converter
- The IDs were uploaded to consensus path db for the over-representation analysis
- Data was downloaded and then uploaded to R to have a better visualization


```{r ORA DEMs}
 library(readr)
overlapping2 <- read_delim("ORA_results.tab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
overlapping_filtered2 <- overlapping2 %>% filter(`p-value`<0.01)
zz=str_wrap(overlapping_filtered2$pathway,width = 50)

p6 <- ggplot() +

  geom_point( data=overlapping_filtered2, mapping=aes(x = overlapping_filtered2$size,
                                                     y =-log(overlapping_filtered2$`q-value`),
                                                     color=overlapping_filtered2$source,
                                                     size=overlapping_filtered2$size)) +
  scale_size(range = c(10, 30),guide = 'none')+

  labs(x = "Size of genes pathway", y = "-Log(q-value)",color="Pathway source",size="# of overlaped metabolites") +
  ggtitle("Pathways analysis")+
  geom_label_repel(aes(x = overlapping_filtered2$size,
                       y =-log(overlapping_filtered2$`q-value`), color=overlapping_filtered2$source,
                       label =str_wrap(overlapping_filtered2$pathway,width=20)) ,
                   min.segment.length = unit(2, 'lines'),
                   size = 3.5,force=1, arrow = arrow(length = unit(0.02, "npc")),segment.color = 'red',
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),show_guide = F) +
theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"))+
theme_minimal()

print(p6)
```
