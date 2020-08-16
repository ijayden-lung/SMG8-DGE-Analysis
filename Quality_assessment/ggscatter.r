#!/usr/bin/env Rscript

library(ggpubr)

pdf("blood.tpm50.MiVsCV.pdf")

df = read.table("blood.tpm50.MiVsCV.txt",sep="\t",header=TRUE,row.names = 1)
sp <- ggscatter(df, x = "M_i", y = "CV",
				shape = 21,
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Internal control gene-stability measure",
				ylab = "Coefficient of variation with normalized count"          
				)
# Add correlation coefficient
sp + stat_cor(method = "spearman", cor.coef.name = "rho",label.x = 1.2, label.y = 3)

#> `geom_smooth()` using formula 'y ~ x'
dev.off()


pdf("skin.tpm50.MiVsCV.pdf")
df = read.table("skin.tpm50.MiVsCV.txt",sep="\t",header=TRUE,row.names = 1)
sp <- ggscatter(df, x = "M_i", y = "CV",
				shape = 21,
				add = "reg.line",  # Add regressin line
				add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
				conf.int = TRUE, # Add confidence interval
				xlab = "Internal control gene-stability measure",
				ylab = "Coefficient of variation with normalized count"          
				)
# Add correlation coefficient
sp + stat_cor(method = "spearman", cor.coef.name = "rho",label.x = 1.2, label.y = 3)

#> `geom_smooth()` using formula 'y ~ x'
dev.off()
