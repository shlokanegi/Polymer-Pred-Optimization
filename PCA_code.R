install.packages("factoextra")
install.packages("FactoMineR")
install.packages("readr")
library(readr)
init_dat <- read_csv("temp_PCA.csv", col_types = cols(polymer = col_skip()))
dat=init_dat[,c(-1,-2,-12)]
View(dat)

library(factoextra)
res.pca <- prcomp(dat, scale = TRUE)
fviz_eig(res.pca) ## SCREE PLOT
summary(res.pca)

## SCORE PLOT (Plot of Individuals)
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, axes = c(3,4))

## LOADING PLOT (Plot of Variables)
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, axes = c(3,4))

## BI-PLOT
fviz_pca_biplot(res.pca, repel = TRUE, col.var = "#2E9FDF", col.ind = "#696969", axes = c(1,2))
fviz_pca_biplot(res.pca, repel = TRUE, col.var = "#2E9FDF", col.ind = "#696969", axes = c(3,4))

#### EIGENVALUES
eig.val <- get_eigenvalue(res.pca)
eig.val

# Contributions of variables to PC1, PC2, PC3 and PC4
fviz_contrib(res.pca, choice = "var", axes = 1, top = 9)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 9)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 9)
fviz_contrib(res.pca, choice = "var", axes = 4, top = 9)

library(FactoMineR)
# Visualize individual cos2 on axes 1,2,3 and 4 together
fviz_cos2(res.pca, choice ="ind", axes = 1:4)

# Visualize variable categories cos2 on axes 1,2,3 and 4 together
fviz_cos2(res.pca, choice ="var", axes = 1:4)


########### REMOVING OUTLIERS AND USING TOP 5 PREDICTORS ONLY #############
library(readr)
dat <- read_csv("Work_drug_R2.csv", col_types = cols(delta_Tg = col_skip()))
ndat = dat[c(-52,-57,-97,-45,-70,-75,-29,-89,-76,-38,-15),c(1,4,6,7,8,10)]
df = ndat[,-6]
res.pca <- prcomp(df, scale = TRUE)
summary(res.pca)

## SCORE PLOT (Plot of Individuals)
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

## LOADING PLOT (Plot of Variables)
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

## BI-PLOT
fviz_pca_biplot(res.pca, repel = TRUE, col.var = "#2E9FDF", col.ind = "#696969", axes = c(1,2))


######### GROUP BY Tg
dat <- read_csv("Work_drug_R3(Tg).csv")
ndat = dat[c(-52,-57,-97,-45,-70,-75,-29,-89,-76,-38,-15),c(1,4,6,7,8,10)]
df1 = ndat[,-6]
res.pca1 <- prcomp(df1, scale = TRUE)
fviz_pca_ind(res.pca1,
             #label = "none", # hide individual labels
             habillage = ndat$delta_Tg, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             repel = TRUE
)





