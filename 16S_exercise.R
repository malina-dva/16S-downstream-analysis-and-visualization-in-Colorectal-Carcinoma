#' ---
#' title: "16S downstream analysis in colorectal carcinoma"
#' author: "Malina Doynova"
#' date: "5/02/2024"
#' output:
#'   html_notebook:
#'     theme: simplex # Change the theme to 'darkly'
#'     highlight: tango # Change the highlight style
#'     toc: yes # Include table of content
#'     toc_float: false # If to be floating or not 
#' ---
#' 
## ----setup, include=FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' # 1.  Install and/or load the required packages
#' 
## -----------------------------------------------------------------------------------------------------

# install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "vegan", "microbiome", "Maaslin2"))
# BiocManager::install(c("microbiomeMarker"))
# BiocManager::install(c("ggplot2", "tidyverse", "dplyr", "data.table"))

library("phyloseq")
library("vegan")
library("ggplot2")
library("microbiomeMarker")
library("tidyverse")
library("dplyr")
library("data.table")
library("microbiome")
library("Maaslin2")
library("igraph")
library("RColorBrewer")
library("nlme")


#' 
#' # 2. Set the working directory
#' 
## -----------------------------------------------------------------------------------------------------
setwd("C://Users//Malina//Desktop//Microbial compositionality")

#' 
#' # 3. Load the example dataset
#' 
## ----echo = TRUE--------------------------------------------------------------------------------------
data(kostic_crc)
kostic_crc
## Explore the dataset
dim(sample_data(kostic_crc))
head(sample_data(kostic_crc))

table(sample_data(kostic_crc)$DIAGNOSIS)

head(tax_table(kostic_crc))

head(otu_table(kostic_crc)[,1:20])

#' 
#' # 4. Function to create a phyloseq object from data files
#' 
## -----------------------------------------------------------------------------------------------------

create_phyloseq <- function(otu_file, sample_file, tax_file) {
  # Load OTU table
  otu_data <- read.csv(otu_file, row.names = 1, stringsAsFactors = FALSE)
  otu_matrix <- as.matrix(otu_data)
  
  # Load sample data
  sample_data <- read.csv(sample_file, row.names = 1)
  
  # Load taxonomic table
  tax_data <- read.csv(tax_file, row.names = 1, stringsAsFactors = FALSE)
  tax_table <- tax_data[, 1:7]  # Assuming taxonomic ranks are in columns 1 to 7
  
  # Create phyloseq object
  physeq_obj <- phyloseq(otu_table(otu_matrix, taxa_are_rows = TRUE),
                         sample_data,
                         tax_table)
  
  return(physeq_obj)
}

# this is how the function is run

# # Provide file paths
# otu_file <- "path/to/otu_table.csv"
# sample_file <- "path/to/sample_data.csv"
# tax_file <- "path/to/taxonomic_table.csv"
# 
# # Create phyloseq object
# physeq_obj <- create_phyloseq(otu_file, sample_file, tax_file)



#' # 5. Normalization of the sample counts
#' 
## ----echo=TRUE----------------------------------------------------------------------------------------
## check how large the sample sizes (number of counts) are 

max_difference = max(sample_sums(kostic_crc))/min(sample_sums(kostic_crc))
max_difference

max(sample_sums(kostic_crc))
min(sample_sums(kostic_crc))

## norm is only needed if sample sizes are more than 10x different (and they are)

## distribution of the sample sizes
hist(sample_sums(kostic_crc), breaks = 177)


## get relative abundance with the microbiome package
kostic_crc.compositional <- microbiome::transform(kostic_crc, "compositional")

head(otu_table(kostic_crc.compositional)[,1:20])

## get relative abundance with the microbiomeMarker package
kostic_crc_TSS = microbiomeMarker::normalize(kostic_crc, method="TSS")
head(otu_table(kostic_crc_TSS)[,1:20])

## output of both from above are relative abundances

## css norm (originally from metagenomeSeq packege) ## 
# Cumulative Sum Scaling (CSS) is a median-like quantile normalization which corrects differences in sampling depth (library size).
# While standard relative abundance (fraction/percentage) normalization re-scales all samples to the same total sum (100%), CSS keeps a variation in total counts between samples. 
# CSS re-scales the samples based on a subset (quartile) of lower abundant taxa 
#(relatively constant and independent), 
# thereby excluding the impact of (study dominating) high abundant taxa.

kostic_crc_CSS = microbiomeMarker::normalize(kostic_crc, method="CSS")
head(otu_table(kostic_crc_CSS)[,1:20])

## cpm norm originally lefse package ###
kostic_crc_CPM = microbiomeMarker::normalize(kostic_crc, method="CPM")

head(otu_table(kostic_crc_CPM)[,1:20])
#tail(otu_table(kostic_crc_CSS))


## rarefy norm ###
## from microbiomeMarker package ##
# 
# kostic_crc_rare = microbiomeMarker::normalize(kostic_crc, method="rarefy")
# 
# head(otu_table(kostic_crc_rare)[,1:20])

## from phyloseq package 

ps.rarefied = phyloseq::rarefy_even_depth(kostic_crc, rngseed=123, sample.size=524, replace=F)



#' 
#' 
#' # 6. Plot stack bar plots of the top 15 taxa based on abund
#' 
## ----echo=TRUE----------------------------------------------------------------------------------------
##  extract the top 15 taxa on this taxa level for CSS 
top_phy = tax_glom(kostic_crc_CSS, "Genus")
top15 = names(sort(taxa_sums(top_phy), decreasing=TRUE)[1:15])
top15_css_ps = prune_taxa(top15, top_phy)


phyloseq::plot_bar(top15_css_ps, "DIAGNOSIS", fill="Genus")+
  theme( axis.ticks.x=element_blank(), panel.background=element_rect(fill=NA), panel.grid.major=element_line(colour="#ebebeb")) + 
  labs(x=NULL)  

## get rid of the black squares and get better colors
plot_bar(top15_css_ps, "DIAGNOSIS", fill="Genus")+
  theme( axis.ticks.x=element_blank(), panel.background=element_rect(fill=NA), panel.grid.major=element_line(colour="#ebebeb")) + 
  labs(x=NULL) + geom_bar(stat="identity") 


## Generate 15 distinct colors from the palette

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
sampled_colors <- sample(color, 15)

plot_bar(top15_css_ps, "DIAGNOSIS", fill = "Genus") +
  theme(axis.ticks.x = element_blank(), panel.background = element_rect(fill = NA), panel.grid.major = element_line(colour = "#ebebeb")) +
  labs(x = NULL) +
  geom_bar(stat = "identity", aes(fill = Genus)) +  # Specify fill aesthetic here
  scale_fill_manual(values = sampled_colors )  # Use scale_fill_manual for fill colors

## plot each genus in a separate facet

# plot_bar(top15_css_ps, "DIAGNOSIS", fill="DIAGNOSIS", facet_grid=~Genus)  + # # #geom_bar(stat="identity")

## plot top 5 ##

top_phy5 = tax_glom(kostic_crc_CSS, "Genus")
top5 = names(sort(taxa_sums(top_phy5), decreasing=TRUE)[1:5])
top5_css_ps = prune_taxa(top5, top_phy5)

tax_table(top5_css_ps)

## plot each genus in a separate facet
plot_bar(top5_css_ps, "DIAGNOSIS", fill="DIAGNOSIS", facet_grid=~Genus) + geom_bar(stat="identity")

## changing the appearance 
## https://stackoverflow.com/questions/63133884/removing-the-original-color-outline-in-r-when-using-a-new-pallette-in-a-barplot
## remove black boxes 
## https://github.com/joey711/phyloseq/issues/721

##  extract the top 15 taxa on this taxa level for TSS
top_phy = tax_glom(kostic_crc_TSS, "Genus")
top15 = names(sort(taxa_sums(top_phy), decreasing=TRUE)[1:15])
top15_tss_ps = prune_taxa(top15, top_phy)


# prepare to create a nice plot for saving
top15_tss_G  =plot_bar(top15_tss_ps, "DIAGNOSIS", fill="Genus")+
  theme( axis.ticks.x=element_blank(), panel.background=element_rect(fill=NA), panel.grid.major=element_line(colour="#ebebeb")) + 
  labs(x=NULL)+
  geom_bar(stat = "identity", aes(fill = Genus)) +  
  scale_fill_manual(values = sampled_colors ) +
  theme_bw() +
  theme(
    text = element_text(size = 12, face = "bold"),  # Make all text bold
    axis.title = element_text(size = 14, face = "bold"),  # Make axis titles bold
    axis.text = element_text(size = 12),  # Adjust axis text size
    legend.title = element_text(size = 14, face = "bold"),  # Make legend title bold
    legend.text = element_text(size = 12),  # Adjust legend text size
    plot.title = element_text(size = 16, face = "bold")  # Make plot title bold
  )

  
top15_tss_G

#save the image
#ggsave("top15_Genera_plot_TSS.png", top15_tss_G, width = 6, height = 8, dpi = 300)

# plot rarefaction curve 
# Extract OTU abundance data and transpose
# otu_table_transposed<- as.data.frame(t(otu_table(kostic_crc)))
# head(otu_table_transposed)
# rarecurve(as.matrix(otu_table_transposed), step=500, cex=0.1)


#' # 7. Plot alpha diversity between groups
#' 
## ----echo=TRUE----------------------------------------------------------------------------------------

## estimation of alpha div only works with integers and not rel abund values (float numbers)
## what each of the alpha div means 
## Observed - counts the number of unique species (OTUs or ASVs) present in a sample
## Chao1 -  It takes into account the number of rare species (singletons and doubletons) observed in the sample and extrapolates the total richness based on the frequency of these rare species
## ACE -  Takes into account both the number of rare species observed and their relative abundances. inflates the number of rare taxa and inflates again the number of taxa with abundance 1.
## Shannon's diversity index (also known as Shannon-Wiener index) considers both species richness and evenness in a community. It takes into account the number of different species present as well as the relative abundance of each species.
## Simpson's diversity index measures the probability that two randomly selected sequences are of the same species


## CSS ##
alpha_div = plot_richness(kostic_crc_CSS, x="DIAGNOSIS", color="DIAGNOSIS", title="CSS", measures=c("Observed", "Chao1", "ACE", "Simpson", "Shannon")) + 
  geom_boxplot() +
  theme_bw() +
  theme(
    text = element_text(size = 12, face = "bold"),  # Make all text bold
    axis.title = element_text(size = 14, face = "bold"),  # Make axis titles bold
    axis.text = element_text(size = 12),  # Adjust axis text size
    legend.title = element_text(size = 14, face = "bold"),  # Make legend title bold
    legend.text = element_text(size = 12),  # Adjust legend text size
    plot.title = element_text(size = 16, face = "bold")  # Make plot title bold
  )

alpha_div

#save the image
#ggsave("alpha_diversity_plot_CSS.png", alpha_div, width = 10, height = 6, dpi = 300)

## rarefies ##

plot_richness( ps.rarefied, x="DIAGNOSIS", color="DIAGNOSIS", title="rarefied", measures=c("Observed", "Chao1", "ACE", "Simpson", "Shannon")) + 
  geom_boxplot() + 
  theme_bw() 

## stat sign diffrences in alpha?
##https://microbiome.github.io/course_2021_radboud/alpha-diversity.html

## https://rpubs.com/lconteville/713954
richness <- estimate_richness(kostic_crc_CSS)
head(richness)

## for normally distributed shannon index values 
hist(richness$Shannon, breaks =  25)
shapiro_test <- shapiro.test(richness$Shannon)
shapiro_test
## almost normally distributed but because of the outliers it is not 
anova.sh = aov(richness$Shannon ~ sample_data(kostic_crc_CSS)$DIAGNOSIS)
summary(anova.sh)

kruskal.test(richness$Shannon ~ sample_data(kostic_crc_CSS)$DIAGNOSIS)

## because the alpha div is one number per sample we can append this to the metadata of 
## our phyloseq object 

## accounting for the multiple sampling from the same individial
# Load the nlme package

tail(sample_data(kostic_crc_CSS)$ANONYMIZED_NAME)

richness_exp <- cbind(richness, DIAGNOSIS = sample_data(kostic_crc_CSS)$DIAGNOSIS, ANONYMIZED_NAME = as.factor(sample_data(kostic_crc_CSS)$ANONYMIZED_NAME))

mixed_model <- lme(Shannon ~ DIAGNOSIS, data = richness_exp, random = ~ 1 | ANONYMIZED_NAME)

summary(mixed_model)


## if you have experimental design where the group you check has more than two levels
## you will have to perform posthoc test in addition to anova and KW

## here is an online tutorial that can help with this
## https://scienceparkstudygroup.github.io/microbiome-lesson/04-alpha-diversity/index.html


## Faith's diversity or PD (Phylogenetic diversity) are based on phylogenetic information for each sample (which we don't have for this data)

# library(biomeUtils)
# data("FuentesIliGutData")
# # reduce size for example
# ps1 <- subset_samples(FuentesIliGutData, ILI == "C")
# ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
# 
# meta_tib <- calculatePD(ps1, justDF = TRUE)
# # check
# meta_tib[c(1, 2, 3), c("PD", "SR")]
# #>                 PD  SR
# #> sample_1  967.3003 259
# #> sample_2 1035.1305 291
# #> sample_3 1189.0248 336

## calculate Pielou’s evenness
data_evenness <- vegan::diversity(data.frame(otu_table(kostic_crc_CSS))) / log(specnumber(data.frame(otu_table(kostic_crc_CSS))))
head(data_evenness)


#' 
#' # 8. Estimate and plot the overall dissimilarity (beta diversity) in the microbial composition between the groups
#' 
## ----echo=TRUE, message = FALSE-----------------------------------------------------------------------

## beta diversity ##
dist_methods = unlist(distanceMethodList)
dist_methods

## Jaccard and Bray-Curtis dissimilarities are non-phylogenetic measures that consider the presence/absence and abundance of features in samples, respectively, while UniFrac dissimilarities are phylogenetic measures that 
## incorporate information about the evolutionary relatedness of microbial taxa

physeq2.ord <- phyloseq::ordinate(kostic_crc_CSS, "NMDS", "bray")

## NMDS is non-linear and focuses on preserving rank order of dissimilarities, PCA is linear and focuses on capturing maximum variance, 
## and PCoA can capture both linear and non-linear relationships and different types of distances measures

#shape is also a parameter we can assign metadata variable to 
p = plot_ordination(kostic_crc_CSS, physeq2.ord, type="samples", color="DIAGNOSIS",
                    title="OTUs")

# Add ellipses
beta_div <- p + stat_ellipse(level = 0.95, type = "norm", geom = "polygon", alpha = 0, aes(color = DIAGNOSIS), size = 1) +
  geom_point(size = 4) + # Adjust the size of points+
  theme_bw() +
  guides(size = "none") + # 
  theme(
    text = element_text(size = 12, face = "bold"),  # Make all text bold
    axis.title = element_text(size = 14, face = "bold"),  # Make axis titles bold
    axis.text = element_text(size = 12),  # Adjust axis text size
    legend.title = element_text(size = 14, face = "bold"),  # Make legend title bold
    legend.text = element_text(size = 12),  # Adjust legend text size
    plot.title = element_text(size = 16, face = "bold")  # Make plot title bold
  )

beta_div

#save image 
#ggsave("beta_diversity_plot_CSS.png", beta_div , width = 10, height = 6, dpi = 300)


## stat of the beta diversity 

metadata = data.frame(sample_data(kostic_crc_CSS))
phy_beta = phyloseq::distance(kostic_crc_CSS, method = "bray")
adonis2(phy_beta ~ DIAGNOSIS, data = metadata, perm=999)

## The adonis2 function fits a multivariate linear model to the distance matrix using the grouping variable (DIAGNOSIS) as a predictor. 
## It then assesses the significance of the relationship between the grouping variable and the structure of the distance matrix through permutation testing.
## it assumes homogeneity of dispersion among groups

## A small p-value suggests that the grouping variable significantly explains the variation in the distance matrix,
## indicating differences in community composition or structure between the groups defined by the DIAGNOSIS variable

## betadisper function checks if the within groups dispersion are homogeneous(compositions vary similarly within the group)

bd = betadisper(phy_beta, metadata$'DIAGNOSIS')
anova(bd)


## perform the analyses correcting for multiple sampling from the same individual

## change the sample id to be a factor 

# Extract sample data from the phyloseq object
sample_data <- sample_data(kostic_crc_CSS)
# write.csv(data.frame(sample_data), "sample_data.csv")

# Modify the variable of interest to be a factor
sample_data$ANONYMIZED_NAME <- as.factor(sample_data$ANONYMIZED_NAME)
levels(sample_data$ANONYMIZED_NAME)
# Replace the modified sample data in the phyloseq object
sample_data(kostic_crc_CSS) <- sample_data


## check relationship between individuals 
plot_ordination(kostic_crc_CSS, physeq2.ord, type="samples", color="ANONYMIZED_NAME",
                title="OTUs")+
  guides(color = FALSE) # no legend 


# Extract the NMDS results from physeq2.ord
nmds_results <- scores(physeq2.ord)

# Combine the NMDS results with the sample IDs
plot_data <- cbind(sample_data(kostic_crc_CSS), nmds_results$sites)

# Plot the ordination
p <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = X.SampleID)) +
  geom_point() +
  labs(title = "OTUs") +
  guides(color = FALSE)

p + geom_text(aes(label = ANONYMIZED_NAME), size =2.5, vjust = -0.5)

adonis2(phy_beta ~ metadata$DIAGNOSIS, data = metadata, perm=999, strata = metadata$ANONYMIZED_NAME)

## try weighted unifrac for diagnosis 

# phy_tree(kostic_crc_CSS)
# physeq2.ord1<- ordinate(kostic_crc_CSS, "NMDS", "wunifrac") #or unifrac
# we don't have the phy_tree 
# but if we did, we could have also explored this beta_div measure 


#' # 9. Identify DA taxa between groups
## ----eval= FALSE--------------------------------------------------------------------------------------
## 
## ## microbiomeMarker is a wrapper for a lot of R packages specifically designed to
## ## profile DA taxa from microbial compositional data
## 
## ## lefse ##
## 
## ## KW detect stat sign features
## ## Wolcoxon and multigrp_strat test for consistency of stat sign across groups of replicates
## ## LDA analyses are finally performed to identify the exact effect size of the difference
## 
## ## CPM is the default norm method for lefse ##
## 
## df_kostic_cpm =normalize(kostic_crc, method = "CPM")%>%microbiomeMarker::run_lefse(
##                                        wilcoxon_cutoff = 0.05,
##                                        norm = "none",
##                                        group = "DIAGNOSIS",
##                                        kw_cutoff = 0.05,
##                                        multigrp_strat = TRUE,
##                                        lda_cutoff = 1.5
## )
## 
## df_kostic_cpm_table  = (marker_table(df_kostic_cpm))
## dim(df_kostic_cpm_table)
## # 158 DA taxa across taxa levels
## head(df_kostic_cpm_table)
## 
## #write the table
## #write.csv(data.frame(df_kostic_cpm_table), "df_kostic_cpm_table_DA.csv")
## 
## ## results are obtained not only for species but for all taxa levels
## ## We can't easily correct for co-founders or repeated measures with these analyses
## 
## 
## ### prepare the data for Maaslin2
## ## MaAsLin2 relies on general linear models to accommodate most modern epidemiological study designs,
## ## including cross-sectional and longitudinal, and offers a variety of data exploration, normalization, and transformation methods
## ## it allows correction for co-founders and repeated measures (correcting for intraindividual variability)
## 
## 
## kostic_crc_CPM = microbiomeMarker::normalize(kostic_crc, method = "CPM")
## 
## crc_CPM_OTU = data.frame(otu_table(kostic_crc_CPM))
## dim(crc_CPM_OTU)
## head(crc_CPM_OTU)
## crc_CPM_TAXA = data.frame((tax_table(kostic_crc_CPM)))
## dim(tax_table(kostic_crc_CPM))
## head(crc_CPM_TAXA)
## crc_CPM_META = data.frame(sample_data(kostic_crc_CPM))
## head(crc_CPM_META)
## # These are how they should be
## crc_CPM_TAXA = crc_CPM_TAXA %>% dplyr::mutate(Species1 = paste(rownames(crc_CPM_TAXA), Genus, Species, sep = "_"))
## crc_CPM_TAXA = crc_CPM_TAXA %>% dplyr::select(Species1)
## 
## head(crc_CPM_TAXA)
## head(crc_CPM_OTU[, 1:20])
## 
## # Add row names as a column in both data frames
## crc_CPM_OTU <- rownames_to_column(crc_CPM_OTU, "RowNames")
## crc_CPM_TAXA <- rownames_to_column(crc_CPM_TAXA, "RowNames")
## 
## # Merge the two data frames by the new column "RowNames"
## crc_CPM_OTU_S <- merge(crc_CPM_OTU, crc_CPM_TAXA, by.x = "RowNames", by.y = "RowNames")
## head(crc_CPM_OTU_S)
## # Optionally, remove the redundant "RowNames" column after merging
## crc_CPM_OTU_S <- crc_CPM_OTU_S %>% dplyr::select(-RowNames)
## crc_CPM_OTU_S = column_to_rownames(crc_CPM_OTU_S, "Species1")
## 
## crc_CPM_OTU_S = data.frame(t(crc_CPM_OTU_S))
## head(crc_CPM_OTU_S)
## colnames(crc_CPM_OTU_S) <- sub("X", "", colnames(crc_CPM_OTU_S))
## 
## crc_CPM_META$ANONYMIZED_NAME = as.factor(crc_CPM_META$ANONYMIZED_NAME)
## 
## 
## head(crc_CPM_OTU_S[1:20, 1:20])
## head(crc_CPM_META[1:20, 1:20])
## 
## ## We will need crc_CPM_OTU_S and crc_CPM_META
## 
## fit_data <- Maaslin2(
##   crc_CPM_OTU_S, crc_CPM_META, 'C://Users//Malina//Desktop//Microbial compositionality/maaslin_results_random_effect/',
##   fixed_effects = c('DIAGNOSIS'),
##   normalization = "NONE",
##   reference = c("DIAGNOSIS,Healthy"),
##   random_effects = c("ANONYMIZED_NAME")
##   )
## 
## 
## ## DA analyses of the same dataset with deseq2
## ## https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
## 

#' # 10. Get co-abundance interaction network for either Tumor or Healthy samples  
## ----echo= TRUE---------------------------------------------------------------------------------------

### will subset the initial CPM norm phyloseq obj to get only the Tumor samples and will
### generate graph from it 
## will subset species that are present at least in 50% of the samples (filtering on prevalence)
## threshold for correlation will be 0.4 (spearman)

kostic_crc_CPM

subset_kostic_crc <- subset_samples(kostic_crc_CPM, DIAGNOSIS == "Tumor")
# Define the prevalence threshold (e.g., taxa present in at least 50% of samples)
prevalence_threshold <- 0.5

# Filter taxa based on prevalence
subset_kostic_crc <- filter_taxa(subset_kostic_crc , function(x) sum(x > 0) / length(x) >= prevalence_threshold, TRUE)

## create OTU input table with species or genus (when NA is present for species) names
crc_CPM_OTU = data.frame(otu_table(subset_kostic_crc))
dim(crc_CPM_OTU)
crc_CPM_TAXA = data.frame((tax_table(subset_kostic_crc)))
dim(tax_table(subset_kostic_crc))

crc_CPM_META = data.frame(sample_data(subset_kostic_crc))
head(crc_CPM_META)
# These are how they should be 
crc_CPM_TAXA = crc_CPM_TAXA %>% dplyr::mutate(Species1= paste(rownames(crc_CPM_TAXA), Genus, Species, sep = "_"))
crc_CPM_TAXA = crc_CPM_TAXA %>% dplyr::select(Species1)

head(crc_CPM_TAXA)
head(crc_CPM_OTU[, 1:20])

# Add row names as a column in both data frames
crc_CPM_OTU <- rownames_to_column(crc_CPM_OTU, "RowNames")
crc_CPM_TAXA <- rownames_to_column(crc_CPM_TAXA, "RowNames")

# Merge the two data frames by the new column "RowNames"
crc_CPM_OTU_S <- merge(crc_CPM_OTU, crc_CPM_TAXA, by.x = "RowNames", by.y = "RowNames")
head(crc_CPM_OTU_S)
# Optionally, remove the redundant "RowNames" column after merging
crc_CPM_OTU_S <- crc_CPM_OTU_S %>% dplyr::select(-RowNames)
crc_CPM_OTU_S = column_to_rownames(crc_CPM_OTU_S, "Species1")

rownames(crc_CPM_OTU_S) <- sub("X", "", rownames(crc_CPM_OTU_S))
head(crc_CPM_OTU_S[,1:20])

dim(crc_CPM_OTU_S)

### build unweighted graph 
correlation_matrix <- cor(data.frame(t(crc_CPM_OTU_S)),method = "spearman")
head(correlation_matrix[,1:5])
threshold <- 0.4

# Create an adjacency matrix based on the correlation threshold
adjacency_matrix <- ifelse(abs(correlation_matrix) > threshold, 1, 0)
head(adjacency_matrix[,1:5])

# Remove the diagonal values 
diag(adjacency_matrix) <- 0
head(adjacency_matrix[,1:5])

graph <- graph_from_adjacency_matrix(adjacency_matrix , 
                                     mode = "undirected")
graph
#plot(graph)

#Cluster based on edge betweenness
ceb = cluster_edge_betweenness(graph)

#Extract community memberships from the clustering
community_membership <- membership(ceb)
head(community_membership)
length(ceb)
# Count the number of nodes in each cluster
cluster_sizes <- table(community_membership)
cluster_sizes

# Add cluster memberships as a vertex attribute
igraph::V(graph)$cluster_membership <- community_membership

# Map cluster memberships to colors using a color palette function
cluster_colors <- rainbow(max(community_membership))
head(cluster_colors)
# Create a color variable based on cluster memberships
V(graph)$color_variable <- cluster_colors[community_membership]

# vetrex degree
vertex_degrees <- igraph::degree(graph)
head(vertex_degrees)
# Add vertex degree as a vertex attribute
igraph::V(graph)$degree <- vertex_degrees


vertex.attr = list(
  cluster_membership = V(graph)$cluster_membership,
  name = V(graph)$name,
  degrees = V(graph)$degree,
  color = V(graph)$color_variable)

# Specify the file name for saving the GraphML file
graphml_file <- "T1graph_kostic_CPM_with_metadata_unweighted_undirected.graphml"

# Write the igraph object to GraphML format with metadata
write_graph(
  graph,
  graphml_file,
  format  = "graphml")

# Combine vertex attributes into a data frame
vertex_data <- data.frame(vertex.attr)

head(vertex_data)
write.csv(vertex_data, "T1vertex_data_kostic_CPM_with_metadata_unweighted_undirected.csv", row.names = FALSE)

## create a function for building co-abundance network

get_interaction_network = function(norm_phyloseq_obj, group, graphml_file, vertex_data_file){
  
  subset_kostic_crc <- subset_samples(norm_phyloseq_obj, DIAGNOSIS == paste(group))
  # Define the prevalence threshold (e.g., taxa present in at least 50% of samples)
  prevalence_threshold <- 0.5
  
  # Filter taxa based on prevalence
  subset_kostic_crc <- filter_taxa(subset_kostic_crc , function(x) sum(x > 0) / length(x) >= prevalence_threshold, TRUE)
  
  ## create OTU input table with species or genus (when NA is present for species) names
  crc_CPM_OTU = data.frame(otu_table(subset_kostic_crc))

  crc_CPM_TAXA = data.frame((tax_table(subset_kostic_crc)))
 
  
  crc_CPM_META = data.frame(sample_data(subset_kostic_crc))
  
  # These are how they should be 
  crc_CPM_TAXA = crc_CPM_TAXA %>% dplyr::mutate(Species1= paste(rownames(crc_CPM_TAXA), Genus, Species, sep = "_"))
  crc_CPM_TAXA = crc_CPM_TAXA %>% dplyr::select(Species1)
  
  
  # Add row names as a column in both data frames
  crc_CPM_OTU <- rownames_to_column(crc_CPM_OTU, "RowNames")
  crc_CPM_TAXA <- rownames_to_column(crc_CPM_TAXA, "RowNames")
  
  # Merge the two data frames by the new column "RowNames"
  crc_CPM_OTU_S <- merge(crc_CPM_OTU, crc_CPM_TAXA, by.x = "RowNames", by.y = "RowNames")

  # Optionally, remove the redundant "RowNames" column after merging
  crc_CPM_OTU_S <- crc_CPM_OTU_S %>% dplyr::select(-RowNames)
  crc_CPM_OTU_S = column_to_rownames(crc_CPM_OTU_S, "Species1")
  
  rownames(crc_CPM_OTU_S) <- sub("X", "", rownames(crc_CPM_OTU_S))
  
  ### build unweighted graph 
  correlation_matrix <- cor(data.frame(t(crc_CPM_OTU_S)),method = "spearman")
 
  threshold <- 0.4
  
  # Create an adjacency matrix based on the correlation threshold
  adjacency_matrix <- ifelse(abs(correlation_matrix) > threshold, 1, 0)
 
  
  # Remove the diagonal values 
  diag(adjacency_matrix) <- 0
  
  
  graph <- graph_from_adjacency_matrix(adjacency_matrix , 
                                       mode = "undirected")
  
  
  #Cluster based on edge betweenness
  ceb = cluster_edge_betweenness(graph)
  
  #Extract community memberships from the clustering
  community_membership <- membership(ceb)
  
  # Count the number of nodes in each cluster
  cluster_sizes <- table(community_membership)
  
  # Add cluster memberships as a vertex attribute
  igraph::V(graph)$cluster_membership <- community_membership
  
  # Map cluster memberships to colors using a color palette function
  cluster_colors <- rainbow(max(community_membership))
  
  # Create a color variable based on cluster memberships
  V(graph)$color_variable <- cluster_colors[community_membership]
  
  # vetrex degree
  vertex_degrees <- igraph::degree(graph)
 
  # Add vertex degree as a vertex attribute
  igraph::V(graph)$degree <- vertex_degrees
  
  vertex.attr = list(
    cluster_membership = V(graph)$cluster_membership,
    name = V(graph)$name,
    degrees = V(graph)$degree,
    color = V(graph)$color_variable)
  
  # Specify the file name for saving the GraphML file
  graphml_file <- graphml_file
  
  # Write the igraph object to GraphML format with metadata
  write_graph(
    graph,
    graphml_file,
    format  = "graphml")
  
  # Combine vertex attributes into a data frame
  vertex_data <- data.frame(vertex.attr)
  
  write.csv(vertex_data, vertex_data_file, row.names = FALSE)
  

}

## run the function to obtain the graph file and vertex metadata for the Healthy samples

get_interaction_network(
  norm_phyloseq_obj = kostic_crc_CPM,
  group = "Healthy",
  graphml_file = "H1graph_kostic_CPM_with_metadata_unweighted_undirected.graphml",
  vertex_data_file = "H1vertex_data_kostic_CPM_with_metadata_unweighted_undirected.csv")


