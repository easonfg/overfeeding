rm(list=ls())

library(dplyr)
require(aod)
cytof.data = read.csv('NewRAW.csv', row.names = 1,
                      header = T, stringsAsFactors = F, sep=',', dec='.')

phenotype = read.csv('phenotypes.csv', row.names = 1,
                     header = T, stringsAsFactors = F, sep=',', dec='.')
head(phenotype)

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

## add age and gender and bmi
rownames(phenotype) = gsub('1792-', '', rownames(phenotype))
rownames(phenotype) = gsub('^0', '', rownames(phenotype))
rownames(phenotype)
annotations$age = phenotype$Age[match(annotations$ids, rownames(phenotype))]
annotations[annotations$ids == '443', 'age'] = 55
annotations$bmi = phenotype$BMI[match(annotations$ids, rownames(phenotype))]
annotations[annotations$ids == '443', 'bmi'] = 28.8
annotations$gender = phenotype$Gender[match(annotations$ids, rownames(phenotype))]
annotations[annotations$ids == '443', 'gender'] = 1
annotations$gender = annotations$gender %>% factor()

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
cytof.combat = cytof.data
# neg.ind = data.frame(which(cytof.combat < 0, arr.ind = T))
# # neg.ind[1,]
# neg.res = sapply(1:nrow(neg.ind), function (row){
#   (cytof.combat[neg.ind[row, 1], neg.ind[row, 2]])
# })
# neg.res

# pca.cytof.combat.res = prcomp(cytof.combat,  scale = T)
# library("factoextra")
# fviz_eig(pca.cytof.combat.res)
# pca.plot = fviz_pca_ind(pca.cytof.combat.res,
#                         geom="point",
#                         col.ind = annotations$batches, # Color by the quality of representation
#                         repel = TRUE     # Avoid text overlapping
# )
# print(pca.plot)

### diagnosistics
### diagnosistics

### logfold changes
pipe = function(timepoint1, timepoint2, in.group, datatype, pval) {
  
  print(c(timepoint1, timepoint2, in.group, datatype, pval))
  # browser()
  # data1 = data.frame(cytof.data) %>% slice(c(which(annotations$timepoints == 1)))
  if (nchar(in.group) > 0){
    data1 = data.frame(cytof.combat) %>% slice(c(which(annotations$timepoints == timepoint1 & annotations$groups == in.group)))
    data2 = data.frame(cytof.combat) %>% slice(c(which(annotations$timepoints == timepoint2 & annotations$groups == in.group)))
    sub.annotate1 = annotations %>% slice(c(which(annotations$timepoints == timepoint1 & annotations$groups == in.group)))
    sub.annotate2 = annotations %>% slice(c(which(annotations$timepoints == timepoint2 & annotations$groups == in.group)))
  }else{
    data1 = data.frame(cytof.combat) %>% slice(c(which(annotations$timepoints == timepoint1)))
    data2 = data.frame(cytof.combat) %>% slice(c(which(annotations$timepoints == timepoint2)))
    sub.annotate1 = annotations %>% slice(c(which(annotations$timepoints == timepoint1)))
    sub.annotate2 = annotations %>% slice(c(which(annotations$timepoints == timepoint2)))
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
  
  annotate1.id = sapply(rownames(sub.annotate1), function(x) {
    ls.x = strsplit(x, '_')[[1]]
    ls.x[length(ls.x)]
  }
  )
  rownames(sub.annotate1) = annotate1.id
  
  annotate2.id = sapply(rownames(sub.annotate2), function(x) {
    ls.x = strsplit(x, '_')[[1]]
    ls.x[length(ls.x)]
  }
  )
  rownames(sub.annotate2) = annotate2.id
  
  # find intersect of two data
  inter_row = intersect(data1.id, data2.id)
  inter_row
  
  #subset data1 and 2
  data1 = data1[inter_row,]
  data2 = data2[inter_row,]
  sub.annotate1 = sub.annotate1[inter_row, ]
  sub.annotate2 = sub.annotate2[inter_row, ]
  
  # if (datatype == 'Median'){
  #   library(samr)
  #   samr.data = bind_rows(data1, data2)
  #   samr.data = scale(samr.data)
  #   samr.data = t(samr.data)
  #   samr.data
  #   y = c(seq(1, nrow(data1)), seq(-1, -nrow(data2)))
  #   y
  #   d <- list(x=samr.data, y=y, geneid=row.names(samr.data), genenames=paste("g", as.character(1:nrow(samr.data)),sep=""), logged2=TRUE)
  #   samr.obj<-samr(d,  resp.type="Two class paired", nperms=1000)
  #   
  #   ### set delta, plot and get delta table 
  #   delta=.4
  #   # samr.plot(samr.obj,delta)
  #   delta.table <- samr.compute.delta.table(samr.obj)
  #   delta.table
  #   
  #   ### create significant genes table
  #   siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
  #   siggenes.table
  #   head(siggenes.table$genes.lo)
  #   head(siggenes.table$genes.up)
  #   write.csv(siggenes.table$genes.up, paste('SAM_res/', timepoint, '_UP_', in.group, '_', datatype, '.csv', sep = ''))
  #   write.csv(siggenes.table$genes.lo, paste('SAM_res/', timepoint, '_LO_', in.group, '_', datatype, '.csv', sep = ''))
  #   # (siggenes.table$genes.up) %>% View()
  #   
  # }
  
  # res = sapply(colnames(data1), function(x){
  #   test = t.test(data1[,x], data2[,x], paired = TRUE, alternative = "two.sided")
  #   logFC = log(mean(data2[,x])/mean(data1[,x]), 2)
  #   c(LogFC = logFC, p.val = test$p.value)
  # }
  # )
  
  library(tibble)
  logitstic.data = bind_rows(data1 %>% add_column(timepoint = factor(rep(0, nrow(data1)))),
                             data2 %>% add_column(timepoint = factor(rep(1, nrow(data2)))))
  # logitstic.data[,-ncol(logitstic.data)] = scale(logitstic.data[,-ncol(logitstic.data)])
  logitstic.data = Filter(function(x)(length(unique(x))>1), logitstic.data)
  
  logitstic.data = add_column(logitstic.data, age = c(sub.annotate1$age, sub.annotate2$age), .before = "timepoint")
  logitstic.data = add_column(logitstic.data, bmi = c(sub.annotate1$bmi, sub.annotate2$bmi), .before = "timepoint")
  
  logitstic.data$timepoint = as.factor(as.numeric(logitstic.data$timepoint))
  logitstic.data[,-ncol(logitstic.data)] = log(logitstic.data[,-ncol(logitstic.data)] + 0.000001, 2)
  logitstic.data$group = c(sub.annotate1$groups, sub.annotate2$groups)
  logitstic.data$batch = c(sub.annotate1$batches, sub.annotate2$batches)
  logitstic.data$gender = c(sub.annotate1$gender, sub.annotate2$gender)
  
  res = sapply(colnames(logitstic.data[,1:(ncol(logitstic.data)-6)]), function(x){
    # res = sapply(colnames(logitstic.data[,1:(ncol(logitstic.data)-3)]), function(x){
    print(x)
    
    # shapiro.res = shapiro.test(logitstic.data[,x])
    # if (shapiro.res$p.value > 0.05){
    #   
    #   pair.model = lm(as.formula(paste(x, '~timepoint + group + batch + age+bmi+gender')), data = logitstic.data)
    #   # betas = pair.model$coefficients
    #   # log2fc = log(sum(betas)/betas[1], 2)
    #   wt.res = sapply(2:length(coef(pair.model)), function(x){
    #     temp.wt.res = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(x))
    #     temp.wt.res$result$chi2[3]
    #   })
    #   names(wt.res) = coef(pair.model) %>% names() %>% .[2:length(coef(pair.model))] %>% paste0(., '.p.val')
    #   # c(Log2FC = unname(pair.model$coefficients[2]), wt.res)
    #   c(Log2FC = unname(pair.model$coefficients[2]), wt.res)
    # }else{
    #   browser()
    #   ancova.np <- sm.ancova(Prewt, Postwt, Treat, model="equal", logitstic.data)
    # }
    

    # pair.model = lm(as.formula(paste(x, '~timepoint + group + batch + age+bmi+gender')), data = logitstic.data)
    # 
    # if (shapiro.test(logitstic.data[,x])$p.val < 0.05 & ols_test_normality(pair.model)$shapiro$p.value > 0.05){
    #   browser()
    # }
    # if (shapiro.test(logitstic.data[,x])$p.val > 0.05 & ols_test_normality(pair.model)$shapiro$p.value < 0.05){
    #   browser()
    # }
    
    pair.model = lm(as.formula(paste(x, '~timepoint + group + batch + age+bmi+gender')), data = logitstic.data)
    if (ols_test_normality(pair.model)$shapiro$p.value > 0.05){
      
      wt.res = sapply(2:length(coef(pair.model)), function(x){
        temp.wt.res = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(x))
        temp.wt.res$result$chi2[3]
      })
      names(wt.res) = coef(pair.model) %>% names() %>% .[2:length(coef(pair.model))] %>% paste0(., '.p.val')
      # c(Log2FC = unname(pair.model$coefficients[2]), wt.res)
      c(Log2FC = unname(pair.model$coefficients[2]), wt.res)
    }
    
    
    # pair.model = glm(as.formula(paste(x, '~timepoint')),
    #                  data=logitstic.data, family=binomial(link="logit"))
    # pair.model = lm(as.formula(paste(x, '~timepoint + group+batch')), data = logitstic.data)
    # wt.res
    # wt.res.timepoint = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(2))
    # wt.res.ir = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(3))
    # wt.res.is = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(4))
    # wt.res.batch = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(5))
    # wt.res.age = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(6))
    # wt.res.bmi   = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(7))
    # wt.res.gender = wald.test(b=coef(pair.model), Sigma=vcov(pair.model), Terms=c(8))
    # sum.model = summary(pair.model)
    # c(Log2FC = unname(pair.model$coefficients[2]),
    #   timepoint.p.val = wt.res.timepoint$result$chi2[3],
    #   ir.p.val = wt.res.ir$result$chi2[3],
    #   is.p.val = wt.res.is$result$chi2[3],
    #   batch.p.val = wt.res.batch$result$chi2[3],
    #   age.p.val = wt.res.age$result$chi2[3],
    #   bmi.p.val = wt.res.bmi$result$chi2[3],
    #   gender.p.val = wt.res.gender$result$chi2[3]
    # )
  }
  )
  
  # res
  # dim(res)
  # str(data.frame(res))
  # res = data.frame(t(res))
  res = Filter(Negate(is.null), res) %>% data.frame() %>% t() %>% data.frame()
  res
  head(res)
  # res[1:10,1:2]
  # str(res)
  
  for( col in colnames(res)[2:length(colnames(res))]){
    res[,paste0(col, '.adj')] = p.adjust(res[,col], method = 'fdr')
  }
  
  # res$timepoint.padj = p.adjust(res$timepoint.p.val.P, method = 'BH')
  # res$ir.padj = p.adjust(res$ir.p.val.P, method = 'BH')
  # res$is.padj = p.adjust(res$is.p.val.P, method = 'BH')
  # res$batch.padj = p.adjust(res$batch.p.val.P, method = 'BH')
  # res$age.padj = p.adjust(res$age.p.val.P, method = 'BH')
  # res$bmi.padj = p.adjust(res$bmi.p.val.P, method = 'BH')
  # res$gender.padj = p.adjust(res$gender.p.val.P, method = 'BH')
  res = res[order(res[,2]),]
  # head(res) %>% View()
  
  write.csv(res[res[,2]  < 0.05,], paste('logistic2_more_ctrl_batchCorrectinLM/', timepoint1, '_', timepoint2, '_', in.group, '_', datatype, '_more_ctrl_BatchCinLM.csv', sep = ''))
  
  
  # library(EnhancedVolcano)
  # # jpeg(paste('volcano/combat_correct/', datatype, '/', timepoint, pval, '_', in.group, '_', datatype, '_ctrl_corrected.jpeg', sep = ''),
  # jpeg(paste('logistic2_controlled/', timepoint, pval, '_', in.group, '_', datatype, '.jpeg', sep = ''),
  #      units="in", width=10, height=10, res=500)
  # res.plot = res
  # res.plot$p.val = -log(res.plot$p.val)
  # res.plot = res.plot %>% rownames_to_column('name') 
  # head(res.plot)
  # # volcano.plot = EnhancedVolcano(res,
  # #                                lab = rownames(res),
  # #                                pCutoff = 0.05,
  # #                                FCcutoff = 0,
  # #                                x = 'Estimates',
  # #                                y = pval,
  # #                                drawConnectors = TRUE,
  # #                                xlim = c(-2,2),
  # #                                xlab = 'Estimates',
  # #                                legendPosition = 'bottom')
  # volcano.plot = ggplot(res.plot, aes(Estimates, p.val, label = name)) +
  #   geom_point(aes(color = ifelse(p.val > -log(0.05), 'sig', 'non-Sig'))) +
  #   geom_hline(yintercept = -log(0.05)) +
  #   geom_text_repel(data = subset(res.plot, p.val > -log(0.05))) +
  #   xlab("Estimates") + 
  #   ylab("-log10(P)") +
  #   xlim(-2,2) +
  #   labs(color='') + ## legend title
  #   theme_bw()
  # print(volcano.plot)
  # dev.off()
}

for (time in list(list(1,2), list(1,3), list(1,4), list(2,3), list(2,4), list(3,4))){
  
  time1 = time[[1]]
  time2 = time[[2]]
  # for (time2 in c(2,3,4)){
  # pipe(time, 'IS',  'Freq', 'padj')
  # pipe(time, 'IS', 'Median', 'padj')
  # 
  # pipe(time, 'IR',  'Freq', 'padj')
  # pipe(time, 'IR', 'Median', 'padj')
  
  # pipe(time, '',  'Freq', 'padj')
  # pipe(time, '', 'Median', 'padj')
  
  # pipe(time, 'IS',  'Freq', 'p.val')
  # pipe(time, 'IR',  'Freq', 'p.val')
  # pipe(time1, time2, '',  'Freq', 'p.val')
  # pipe(time, 'IS', 'Median', 'p.val')
  # pipe(time, 'IR', 'Median', 'p.val')
  pipe(time1, time2, '', 'Median', 'p.val')
  # }
}

