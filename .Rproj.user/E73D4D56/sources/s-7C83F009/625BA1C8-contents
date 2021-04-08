library(dplyr)

### read data
overfeed.data = read.csv('Re__Obesity_data_set_2/overfeeding_Cytokine_Chemokine_luminexWithSSPG.csv', row.names = 1,
                         header = F, stringsAsFactors = FALSE, sep=',', dec='.')


pheno.data = read.csv('cytof/phenotypes.csv', 
                      header = T, stringsAsFactors = FALSE, sep=',', dec='.')
pheno.data
pheno.data$sspg.didff = pheno.data$SSPG2 - pheno.data$SSPG1
pheno.data$wt.diff = pheno.data$Wtwk4 - pheno.data$Base.Wt
pheno.data$lab. = sub('-', '_', pheno.data$lab.)
pheno.data$lab. = sub('/', '_', pheno.data$lab.)

pheno.data = rename(pheno.data, subject.id = 'lab.')

# overfeed.data = data.frame(t(overfeed.data[c(1,3,4,5,6,  grep('_FCH', rownames(overfeed.data))),]))
overfeed.data = data.frame(t(overfeed.data[c(1,3,4,5,6, grep('_B', rownames(overfeed.data)), grep('_S', rownames(overfeed.data))),]))
str(overfeed.data)

### change all columns that should be numbers to numeric
# first find exclude all columns that should be factors
num.col = seq(1,ncol(overfeed.data))[!seq(1, ncol(overfeed.data)) %in% c(1)] 
num.col
# overfeed.data[1:10,1:10]
# change all remaining to numeric
overfeed.data[,num.col] = sapply(num.col,
                                 function(x) as.numeric(as.character(overfeed.data[,x])))



# overfeed.data[,num.col] = as.numeric(overfeed.data[,num.col])

# str(overfeed.data)
# # overfeed.data = rename(overfeed.data, subject.id = 'SUBJECT.ID')
# 
# str(overfeed.data )
# overfeed.data[1:10,1:10]
# is.na(overfeed.data[1:10,1:10])
# overfeed.data[1:10,1:10]
# 
# ## change time point to numbers
# # overfeed.data$TIME.POINT = rep(c(0,2,4,8), nrow(overfeed.data)/4)
# 
# # library(lme4)
# # get the metabolite columns
# # metabolite_cols = colnames(overfeed.data)[7:ncol(overfeed.data)]
# 
# library(dplyr)
# rownames(overfeed.data) = paste0(overfeed.data$subject.id, '_', overfeed.data$TIME.POINT)
# rownames(overfeed.data)
# pca.data = overfeed.data %>% 
#   filter(TIME.POINT == 0)
# # select(7:ncol(overfeed.data)) 
# 
# all.base.data = inner_join(pheno.data, pca.data, by = 'subject.id')
# all.base.data[1:10, 1:20]
# 
# # pca.met.data = pca.data %>% select(7:ncol(overfeed.data)) 
# 
# 
# pca.met.data = all.base.data[complete.cases(all.base.data[,20:ncol(all.base.data)]),]
# sum(colSums(is.na(pca.met.data)) > 0)
# pca.met.data[,20:ncol(pca.met.data)] = Filter(function(x)(length(unique(x))>1), pca.met.data[,20:ncol(pca.met.data)])
# pca.met.data[1:10, 1:20]
# pca.met.data = data.frame(pca.met.data) %>% remove_rownames() %>% column_to_rownames(var = 'subject.id')
# pca.met.data[1:10, 1:20]

rownames(overfeed.data) = overfeed.data$nmID
data1 = overfeed.data[c(1:5, grep('_B', colnames(overfeed.data)))]
colnames(data1) = sub('_B', '', colnames(data1))
data1$type = rep('base', nrow(data1))
data2 = overfeed.data[c(1:5, grep('_S', colnames(overfeed.data)))]
colnames(data2) = sub('_S', '', colnames(data2))
data2$type = rep('peak', nrow(data2))
pca.data = bind_rows(data1, data2)
str(pca.data)

pca.res = prcomp(pca.data[,6:(ncol(pca.data)-1)], scale = T)
library("factoextra")
fviz_eig(pca.res)
# fviz_pca_ind(pca.res, axes = c(3,4),
fviz_pca_ind(pca.res,
             # col.ind = "cos2", # Color by the quality of representation
             # col.ind = as.factor(pca.met.data$GENDER), # Color by the quality of representation
             # col.ind = (pca.met.data$BMI.x), # Color by the quality of representation
             # col.ind = (pca.met.data$GROUP.NAME), # Color by the quality of representation
             col.ind = pca.data$type,
             # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # gradient.cols = c('red', 'white', 'blue'),
             repel = TRUE     # Avoid text overlapping
)


### Run pheatmap using the metadata data frame for the annotation
library(pheatmap)
library(RColorBrewer)
heat_colors <- rev(brewer.pal(7, "RdYlBu"))
pheatmap(t(scale(pca.data[,6:(ncol(pca.data)-1)])),
         annotation_col = data.frame(row.names = rownames(pca.data), 
                                     groups = pca.data$type,
                                     Gender = pca.data$Sexismale),
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         height = 20)
