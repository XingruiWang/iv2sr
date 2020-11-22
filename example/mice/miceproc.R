mice <- read.csv("mice/BXH_ClinicalTraits_361mice_forNewBXH.csv")
x.imp <- t(read.delim("mice/expression_imputed_liver.txt"))
z.imp <- t(read.delim("mice/genotype_imputed.txt"))

mice.id <- intersect(mice$MiceID[!is.na(mice$WeightG)], intersect(rownames(x.imp), rownames(z.imp)))
y <- mice$WeightG[match(mice.id, mice$MiceID)]
x <- x.imp[mice.id, ]
z <- z.imp[mice.id, ]

sex <- mice$sex[match(mice.id, mice$MiceID)] - 1
coef.sex <- lm(y ~ sex)$coef[2]
y <- y - coef.sex*sex

gene.map <- read.delim("mice/gene_map.txt")
x <- x[, colnames(x) %in% gene.map$ID & apply(x, 2, sd) > 0.1]
save(x, y, z, file="mice/mice_processed_liver.rda")
