# SMG8-DGE-Analysis

**Citation: Recessive deleterious variants in SMG8 Expand the Role of Nonsense-Mediated Decay in Developmental Disorders in Humans**

Scripts for DGE analysis

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
```
perl pca.blood.pl 
```

### Quality Assessment

```
perl select_gene.skin.pl
./ggscatter.r
```

### DGE Analysis

```
perl anova.skin.pl     #####Calculate Anova p-value
perl add_info.pl skin  ####Add q-value, gene name and NMD related information
./pheatmap.all.r         ####visualize expression profile
```
