# SMG8-DGE-Analysis

*Citation: Recessive deleterious variants in SMG8 Expand the Role of Nonsense-Mediated Decay in Developmental Disorders in Humans*

Scripts for quality control and DGE analysis

## Dependencies
- perl (v5.26.2)

  - Statistics::ANOVA 0.14

  - Statistics::Descriptive

  - Statistics::Multtest

- R (3.6.3)

  - pheatmap

  - ggfortify

  - ggrepel
  
  - gridExtra
  
  - ggpubr

## Use 

### PCA
We performed sample-level quality control of the normalized count by applying principal components analysis (PCA) and analyzing the distribution of the samples on the PC1 and PC2 coordinates.  
```
perl pca.skin.pl  ###Prepare PCA Inpput
./pcr.r
```


### Quality Assessment
To assess the quality of the scaling of the TPM normalization on our RNA-seq data, we measured the association between the noise level and the stability level of highly expressed genes in both fibroblasts and LCL datasets. For this analysis, we selected genes with the constraints: mean TPM > 50 and minimum TPM > 0 (LCL, n=2,465; fibroblasts, n=2,596), and we computed the internal control gene-stability measure from geNorm (M) [Vandesompele et al., Genome Biol, 2002] of unnormalized counts and the coefficient of variation (CV) of the TPM-normalized counts

```
perl select_gene.skin.pl  ####Select genes and compute the M and CV
./ggscatter.r
```

### DGE Analysis
Prior to our DGE analysis, we computed the z-score of normalized counts for each gene in the control set.  We treated samples that have the absolute z-score > 2.0 as outliers and removed them from the control set for each gene.  To perform this DGE analysis, we developed custom perl scripts.  In the custom scripts, we used the anova method in perl module Statistics::ANOVA with parameters: independent =1; parametric= 1; and ordinal= 0 to perform ANOVA on the datasets. To obtain FDR-adjusted p-values from the ANOVA output, we used perl module Statistics::Multtest.  The log fold-change was defined to be the log2-transformed ratio of the mean of the normalized count from the case samples to the mean of the normalized count from the control samples, and it was also implemented in perl.  
```
perl anova.skin.pl     ##### removed outliers and calculate Anova p-value
perl add_info.pl skin  ####Add q-value, gene name and NMD related information
./pheatmap.all.r         ####visualize expression profile
```
