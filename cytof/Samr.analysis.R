rm(list=ls())

library(dplyr)
require(aod)
cytof.data = read.csv('NewRAW.csv', row.names = 1,
                      header = T, stringsAsFactors = F, sep=',', dec='.')

groups = read.csv('cytofs_w_phenotypes.csv', row.names = 1,
                  header = F, stringsAsFactors = F, sep=',', dec='.')
groups = data.frame(t(groups))
# groups[1:10,1:3]
groups = groups[,c(1,3)]
groups

controls = read.csv('controls.csv', row.names = 1,
                    header = T, stringsAsFactors = F, sep=',', dec='.')

## change control columns to numeric
controls[, sapply(controls, class) == 'character'] = as.numeric(as.matrix(controls[, sapply(controls, class) == 'character']))

str(controls)
dim(controls)
ctrl_correct = rbind(controls[2,]/controls[1,],
                     controls[3,]/controls[1,],
                     controls[4,]/controls[1,],
                     controls[5,]/controls[1,])

## change cytof columns to numeric
cytof.data[, sapply(cytof.data, class) == 'character'] = as.numeric(as.matrix(cytof.data[, sapply(cytof.data, class) == 'character']))

## make annotation matrix
annotations = data.frame(batches = as.factor(sapply(rownames(cytof.data), function(x) strsplit(x, '_')[[1]][2])),
                         timepoints = as.factor(cytof.data$Time))
cytof.data = cytof.data %>% select(-c('Time'))
dim(cytof.data)
annotations$ids = sapply(rownames(annotations), function(x) {
  ls.x = strsplit(x, '_')[[1]]
  # print(ls.x[length(ls.x)])
  ls.x = strsplit(ls.x[length(ls.x)], '\\.')[[1]][1]
  sub('^0', '', ls.x)
}
)
head(groups)
head(annotations)
# length(intersect(unique(annotations$ids), unique(groups$V1)))
annotations$groups = groups$insulin.status[match(annotations$ids, groups$V1)]

## delete columns with one values
non.unique.col = which(apply(cytof.data, 2, function(x) length(unique(x))) <= 1)
cytof.data = cytof.data[,-non.unique.col]
ctrl_correct = ctrl_correct[,-non.unique.col]

# delete na columns
col_na = which(colSums(is.na(cytof.data)) > 0)
cytof.data = cytof.data[,-col_na]
dim(cytof.data)
ctrl_correct = ctrl_correct[,-col_na]
ctrl_correct %>% dim()
cytof.data %>% dim()

### original pca
pca.cytof.res = prcomp(cytof.data,  scale = T)
library("factoextra")
fviz_eig(pca.cytof.res)
pca.plot = fviz_pca_ind(pca.cytof.res,
                        geom="point",
                        col.ind = annotations$batches, # Color by the quality of representation
                        repel = TRUE     # Avoid text overlapping
)

print(pca.plot)

### combat

### all the negative values are really small, so i am turning them into 0
# neg.ind = data.frame(which(cytof.data < 0, arr.ind = T))
# # neg.ind[1,]
# neg.res = sapply(1:nrow(neg.ind), function (row){
#   (cytof.data[neg.ind[row, 1], neg.ind[row, 2]])
# })

cytof.data
cytof.data[cytof.data < 0] <- 0

### all the negative values are really small, so i am turning them into 0

library(sva)
# cytof.combat = t(ComBat(t(log(cytof.data+0.00001)), batch = annotations$batch))
cytof.combat = t(ComBat(t(cytof.data), batch = annotations$batch))

pca.cytof.combat.res = prcomp(cytof.combat,  scale = T)
library("factoextra")
fviz_eig(pca.cytof.combat.res)
pca.plot = fviz_pca_ind(pca.cytof.combat.res,
                        geom="point",
                        col.ind = annotations$batches, # Color by the quality of representation
                        repel = TRUE     # Avoid text overlapping
)
print(pca.plot)


### logfold changes
pipe = function(timepoint, in.group, datatype, pval) {
  
  print(c(timepoint, in.group, datatype, pval))
  # browser()
  # data1 = data.frame(cytof.data) %>% slice(c(which(annotations$timepoints == 1)))
  if (nchar(in.group) > 0){
    data1 = data.frame(cytof.combat) %>% slice(c(which(annotations$timepoints == 1 & annotations$groups == in.group)))
    data2 = data.frame(cytof.combat) %>% slice(c(which(annotations$timepoints == timepoint & annotations$groups == in.group)))
  }else{
    data1 = data.frame(cytof.combat) %>% slice(c(which(annotations$timepoints == 1)))
    data2 = data.frame(cytof.combat) %>% slice(c(which(annotations$timepoints == timepoint)))
  }
  
  data1 = data1 %>% select(grep(datatype, colnames(data1)))
  data2 = data2 %>% select(grep(datatype, colnames(data2)))
  
  ## diagnostics ##
  ## need normal idstribution for logistic analysis
  # library(tidyr)
  # scale(data1)
  # data_long <- gather(data.frame(scale(data1)), gene, measurement, factor_key=TRUE)
  # head(data_long)
  # ggplot(data_long, aes(gene, measurement)) +
  #   geom_point()
  ## diagnostics ##
  
  # find name of data1
  data1.id = sapply(rownames(data1), function(x) {
    ls.x = strsplit(x, '_')[[1]]
    ls.x[length(ls.x)]
  }
  )
  rownames(data1) = data1.id
  # find name of data2
  data2.id = sapply(rownames(data2), function(x) {
    ls.x = strsplit(x, '_')[[1]]
    ls.x[length(ls.x)]
  }
  )
  rownames(data2) = data2.id
  
  # find intersect of two data
  inter_row = intersect(data1.id, data2.id)
  inter_row
  
  #subset data1 and 2
  data1 = data1[inter_row,]
  data2 = data2[inter_row,]
  
  if (datatype == 'Median'){
    library(samr)
    samr.data = bind_rows(data1, data2)
    samr.data = scale(samr.data)
    samr.data = t(samr.data)
    samr.data
    # y = c(seq(1, nrow(data1)), seq(-1, -nrow(data2)))
    y = c(seq(-1, -nrow(data1)), seq(1, nrow(data2)))
    y
    d <- list(x=samr.data, y=y, geneid=row.names(samr.data), genenames=paste("g", as.character(1:nrow(samr.data)),sep=""), logged2=TRUE)
    samr.obj<-samr(d,  resp.type="Two class paired", nperms=1000)
    
    ### set delta, plot and get delta table 
    delta=.4
    # samr.plot(samr.obj,delta)
    delta.table <- samr.compute.delta.table(samr.obj)
    delta.table
    
    sam.pval = samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar) 
    print(sam.pval[sam.pval < 0.05])
    
    
    ### create significant genes table
    siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
    siggenes.table
    siggenes.table.up = siggenes.table$genes.up %>% data.frame()
    siggenes.table.lo = siggenes.table$genes.lo %>% data.frame()
    
    # siggenes.table.lo$pval = sam.pval[which(names(sam.pval) %in% siggenes.table.lo$Gene.Name)] 
    siggenes.table.lo$pval = sam.pval[match(siggenes.table.lo$Gene.Name, names(sam.pval))] %>% na.omit()
    siggenes.table.up$pval = sam.pval[match(siggenes.table.up$Gene.Name, names(sam.pval))] %>% na.omit()
    # head(siggenes.table.lo)
    # sam.pval[grep('Non.Basophils.Leukocytes.Lymphocytes.CD3..CD4.T.cels.HLA.DR.CD38.._.Median.pCREB', names(sam.pval))]
    
    head(siggenes.table$genes.lo)
    head(siggenes.table$genes.up)
    write.csv(siggenes.table.up, paste('SAM_res/', timepoint, '_UP_', in.group, '_', datatype, '.csv', sep = ''))
    write.csv(siggenes.table.lo, paste('SAM_res/', timepoint, '_LO_', in.group, '_', datatype, '.csv', sep = ''))
    # (siggenes.table$genes.up) %>% View()
    
  }
  
}

for (time in c(2,3,4)){
  # pipe(time, 'IS',  'Freq', 'padj')
  # pipe(time, 'IS', 'Median', 'padj')
  # 
  # pipe(time, 'IR',  'Freq', 'padj')
  # pipe(time, 'IR', 'Median', 'padj')
  
  # pipe(time, '',  'Freq', 'padj')
  # pipe(time, '', 'Median', 'padj')
  
  # pipe(time, 'IS',  'Freq', 'p.val')
  # pipe(time, 'IR',  'Freq', 'p.val')
  # pipe(time, '',  'Freq', 'p.val')
  # pipe(time, 'IS', 'Median', 'p.val')
  # pipe(time, 'IR', 'Median', 'p.val')
  pipe(time, '', 'Median', 'p.val')
}

