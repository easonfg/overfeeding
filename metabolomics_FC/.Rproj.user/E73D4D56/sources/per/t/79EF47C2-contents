library(dplyr)
library(tibble)


### read data
overfeed.data = read.csv('../Re__Obesity_data_set/overfeeding_metabolomics.csv', row.names = 1,
                         header = F, stringsAsFactors = FALSE, sep=',', dec='.')

# pheno.data = read.csv('../cytof/phenotypes.csv', 
#                       header = T, stringsAsFactors = FALSE, sep=',', dec='.')
# pheno.data
# pheno.data$sspg.didff = pheno.data$SSPG2 - pheno.data$SSPG1
# pheno.data$wt.diff = pheno.data$Wtwk4 - pheno.data$Base.Wt
# pheno.data$lab. = sub('-', '_', pheno.data$lab.)
# pheno.data$lab. = sub('/', '_', pheno.data$lab.)
# 
# pheno.data = rename(pheno.data, subject.id = 'lab.')

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

## change time point to numbers
overfeed.data$TIME.POINT = rep(c(0,2,4,8), nrow(overfeed.data)/4)

library(lme4)
# get the metabolite columns
metabolite_cols = colnames(overfeed.data)[7:ncol(overfeed.data)]

library(dplyr)
rownames(overfeed.data) = paste0(overfeed.data$subject.id, '_', overfeed.data$TIME.POINT)
rownames(overfeed.data)

pipe = function(timepoint, pval) {
  
  # library(tidyr)
  # log.overfeed.data = overfeed.data
  # # column_to_rownames(log.overfeed.data, var = '')
  # log.overfeed.data.graph = gather(log.overfeed.data, metabolites, measurement,  7:ncol(log.overfeed.data), factor_key=TRUE)
  # ggplot(log.overfeed.data.graph, aes(metabolites, measurement)) +
  #   geom_point()
  # 
  # log.overfeed.data[,7:ncol(log.overfeed.data)] = log(log.overfeed.data[,7:ncol(log.overfeed.data)], 2)
  # scale(log.overfeed.data[,7:ncol(log.overfeed.data)])
  # log.overfeed.data[,7:ncol(log.overfeed.data)] = log.overfeed.data[,7:ncol(log.overfeed.data)] + 10
  # 
  
  # overfeed.data[,7:ncol(overfeed.data)] = scaling(overfeed.data[,7:ncol(overfeed.data)], type = 'pareto')
  
  data1 = overfeed.data %>% filter(TIME.POINT == 0)
  data2 = overfeed.data %>% filter(TIME.POINT == timepoint)
  
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
  
  # data1 = data1 %>% select(-c('X3..cystein.S.yl.acetaminophen.', 'allopurinol','allopurinol.riboside','desmethylnaproxen'))
  # data2 = data2 %>% select(-c('X3..cystein.S.yl.acetaminophen.', 'allopurinol','allopurinol.riboside','desmethylnaproxen'))

  comb.data = bind_rows(data1, data2)
  comb.data[comb.data$TIME.POINT != 0, 'TIME.POINT'] = 1
  comb.data$TIME.POINT = as.factor(comb.data$TIME.POINT)
  comb.data[,7:ncol(comb.data)] = scaling(comb.data[,7:ncol(comb.data)], type = 'pareto')
  infCol = c(which(is.infinite(colSums(comb.data[,7:ncol(comb.data)]))))
  infCol = unique(infCol)
  infCol
  comb.data = comb.data %>% select(-c(infCol+6))
  
  data1 = comb.data %>% filter(TIME.POINT == 0)
  data2 = comb.data %>% filter(TIME.POINT != 0)
  
  res = sapply(colnames(data1)[7:ncol(data1)], function(x){
    print(x)
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
  
  res$padj = p.adjust(res$p.val, method = 'fdr')
  res[which(res$padj < 0.1),]
  
  # library(EnhancedVolcano)
  # jpeg(paste('volcano/', timepoint, pval, '.jpeg', sep = ''),
  #      units="in", width=10, height=10, res=500)
  # volcano.plot = EnhancedVolcano(res,
  #                                lab = rownames(res),
  #                                pCutoff = 0.1,
  #                                FCcutoff = 0,
  #                                x = 'LogFC',
  #                                y = pval,
  #                                legendPosition = 'bottom')
  # print(volcano.plot)
  # dev.off()
}

pipe(2, 'padj')
pipe(4, 'padj')
pipe(8, 'padj')

pipe(2, 'p.val')
pipe(4, 'p.val')
pipe(8, 'p.val')
