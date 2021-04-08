### old version
rm(list=ls())

library(dplyr)
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

library(sva)
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
  
  res = sapply(colnames(data1), function(x){
    test = t.test(data1[,x], data2[,x], paired = TRUE, alternative = "two.sided")
    logFC = log(mean(data2[,x])/mean(data1[,x]), 2)
    c(LogFC = logFC, p.val = test$p.value)
  }
  )
  # res
  # dim(res)
  # str(data.frame(res))
  res = data.frame(t(res))
  # res[1:10,1:2]
  # str(res)
  
  res$padj = p.adjust(res$p.val, method = 'BH')
  res[which(res$padj < 1),]
  head(res)

  write.csv(res[res$p.val < 0.05,], paste('volcano/combat_correct/', timepoint, pval, '_', in.group, '_', datatype, '.csv', sep = ''))
  
  
  library(EnhancedVolcano)
  jpeg(paste('volcano/combat_correct/', datatype, '/', timepoint, pval, '_', in.group, '_', datatype, '_ctrl_corrected.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  volcano.plot = EnhancedVolcano(res,
                                 lab = rownames(res),
                                 pCutoff = 0.05,
                                 FCcutoff = 0,
                                 x = 'LogFC',
                                 y = pval,
                                 xlim = c(min(res$LogFC, na.rm = T) - 5, max(res$LogFC, na.rm = T) + 5),
                                 drawConnectors = TRUE,
                                 legendPosition = 'bottom')
  print(volcano.plot)
  dev.off()
}

for (time in c(2,3,4)){
  pipe(time, 'IS',  'Freq', 'padj')
  pipe(time, 'IS', 'Median', 'padj')
  
  pipe(time, 'IR',  'Freq', 'padj')
  pipe(time, 'IR', 'Median', 'padj')
  
  pipe(time, '',  'Freq', 'padj')
  pipe(time, '', 'Median', 'padj')
  
  pipe(time, 'IS',  'Freq', 'p.val')
  pipe(time, 'IS', 'Median', 'p.val')
  pipe(time, 'IR',  'Freq', 'p.val')
  pipe(time, 'IR', 'Median', 'p.val')
  pipe(time, '',  'Freq', 'p.val')
  pipe(time, '', 'Median', 'p.val')
}
