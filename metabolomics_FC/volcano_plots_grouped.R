library(dplyr)
library(tibble)


### read data
overfeed.data = read.csv('../Re__Obesity_data_set/overfeeding_metabolomics.csv', row.names = 1,
                         header = F, stringsAsFactors = FALSE, sep=',', dec='.')

pheno.data = read.csv('../cytof/phenotypes.csv', 
                      header = T, stringsAsFactors = FALSE, sep=',', dec='.')
head(pheno.data)
pheno.data$sspg.didff = pheno.data$SSPG2 - pheno.data$SSPG1
pheno.data$wt.diff = pheno.data$Wtwk4 - pheno.data$Base.Wt
pheno.data$lab. = sub('-', '_', pheno.data$lab.)
pheno.data$lab. = sub('/', '_', pheno.data$lab.)

pheno.data = rename(pheno.data, subject.id = 'lab.')

### get rid of last column which contains 1 number
overfeed.data = data.frame(t(overfeed.data[,-ncol(overfeed.data)]))

### change all columns that should be numbers to numeric
# first find exclude all columns that should be factors
num.col = seq(1,ncol(overfeed.data))[!seq(1, ncol(overfeed.data)) %in% c(1,4,5,6)] 
# overfeed.data[1:10,1:10]
# change all remaining to numeric
overfeed.data[,num.col] = sapply(num.col,
                                 function(x) as.numeric(as.character(overfeed.data[,x])))

overfeed.data = rename(overfeed.data, subject.id = 'SUBJECT.ID')

str(overfeed.data )
overfeed.data[1:10,1:10]
is.na(overfeed.data[1:10,1:10])
overfeed.data[1:10,1:10]
overfeed.data[1:3,1:10]

## change time point to numbers
overfeed.data$TIME.POINT = rep(c(0,2,4,8), nrow(overfeed.data)/4)

library(lme4)
# get the metabolite columns
metabolite_cols = colnames(overfeed.data)[7:ncol(overfeed.data)]

library(dplyr)
rownames(overfeed.data) = paste0(overfeed.data$subject.id, '_', overfeed.data$TIME.POINT)
rownames(overfeed.data)

pipe = function(timepoint, group, pval) {
  
  data1 = overfeed.data %>% filter(TIME.POINT == 0) %>% filter(GROUP.NAME == group)
  data2 = overfeed.data %>% filter(TIME.POINT == timepoint) %>% filter(GROUP.NAME == group)
  
  delete.row = c(which(rowSums(is.na(data1))>20), which(rowSums(is.na(data2))>20))
  
  delete.col = c(which(apply(data2 %>% select(-c(1:6)), 2, function(x) length(unique(x))) <= 1),
                 which(apply(data1 %>% select(-c(1:6)), 2, function(x) length(unique(x))) <= 1))
  delete.col = delete.col + 6
  unique(delete.col)
  colnames(data1)[delete.col[1]]
  
  data1 = data1 %>% slice(-c(delete.row))
  data2 = data2 %>% slice(-c(delete.row))
  data1 = data1 %>% select(-c(delete.col))
  data2 = data2 %>% select(-c(delete.col))
  
  res = sapply(colnames(data1)[7:ncol(data1)], function(x){
    test = t.test(data1[,x], data2[,x], paired = TRUE, alternative = "two.sided")
    logFC = log(mean(data2[,x])/mean(data1[,x]), 2)
    c(LogFC = logFC, p.val = test$p.value)
  }
  )
  # dim(res)
  # str(data.frame(res))
  res = data.frame(t(res))
  res[1:10,1:2]
  str(res)
  
  res$padj = p.adjust(res$p.val, method = 'BH')
  res[which(res$padj < 1),]
  
  library(EnhancedVolcano)
  jpeg(paste('volcano/', timepoint, pval, group, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  volcano.plot = EnhancedVolcano(res,
                                 lab = rownames(res),
                                 pCutoff = 0.05,
                                 FCcutoff = 0,
                                 x = 'LogFC',
                                 y = pval,
                                 legendPosition = 'bottom')
  print(volcano.plot)
  dev.off()
}

pipe(2, 'Insulin Sensitive', 'padj')
pipe(4, 'Insulin Sensitive', 'padj')
pipe(8, 'Insulin Sensitive', 'padj')
pipe(2, 'Insulin Sensitive', 'p.val')
pipe(4, 'Insulin Sensitive', 'p.val')
pipe(8, 'Insulin Sensitive', 'p.val')

pipe(2, 'Insulin Resistance', 'padj')
pipe(4, 'Insulin Resistance', 'padj')
pipe(8, 'Insulin Resistance', 'padj')
pipe(2, 'Insulin Resistance', 'p.val')
pipe(4, 'Insulin Resistance', 'p.val')
pipe(8, 'Insulin Resistance', 'p.val')

pipe(2, 'Intermediate', 'padj')
pipe(4, 'Intermediate', 'padj')
pipe(8, 'Intermediate', 'padj')
pipe(2, 'Intermediate', 'p.val')
pipe(4, 'Intermediate', 'p.val')
pipe(8, 'Intermediate', 'p.val')