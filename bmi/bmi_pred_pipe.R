library(dplyr)
library(ggplot2)

###### delta sspg ##############################

# bmi.lm = lm(BMI ~., data = combined.data)
# summary(bmi.lm)


### multiple linear
# linear.reg = lm(percdeltaSSPG ~., data = data.frame(scale(sspg)))
# summary(linear.reg)
### multiple linear

pipe = function(dataset, outcome, control_var, graph_title){
  ### univariate linear
  lin.f.stat.p.val.coeff = sapply(colnames(dataset %>% select(-outcome)), function(x){
    if (length(control_var > 0)){
      uni.linear = lm(as.formula(paste0(outcome, ' ~ ', x, paste0('+', control_var, collapse = ''))), data = data.frame(scale(dataset)))
    }else{
      uni.linear = lm(as.formula(paste0(outcome, '~ ', x)), data = data.frame(scale(dataset)))
    }
    uni.sum.res = summary(uni.linear)
    ##check for  collinearity
    if (length(uni.linear$coefficients) > 2){
      if (car::vif(uni.linear)[1] > 2){
        browser()
      }
    }
    # browser()
    f.stat.p.val = pf(uni.sum.res$fstatistic[1],uni.sum.res$fstatistic[2],uni.sum.res$fstatistic[3],lower.tail=FALSE)
    if (f.stat.p.val < 0.05){
      # browser()
      return (c(f.stat.p.val = f.stat.p.val, coeff = uni.sum.res$coefficients[2,1],  coeff.p.val = uni.sum.res$coefficients[2,4]))
    }
  })
  lin.f.stat.p.val.coeff
  
  if (sum(sapply(lin.f.stat.p.val.coeff, function(x) is.null(x)) == FALSE) == 0){
  }else{
    
    if (sum(sapply(lin.f.stat.p.val.coeff, function(x) is.null(x))) == 0){
      sig.genes = data.frame(t(lin.f.stat.p.val.coeff))
    }else{
      sig.genes = lin.f.stat.p.val.coeff[sapply(lin.f.stat.p.val.coeff, function(x) !is.null(x))]
      sig.genes = data.frame(t(as.data.frame(sig.genes)))
      head(sig.genes)
    }
    
    ## adjusted pvalues
    sig.genes = sig.genes[order(sig.genes$coeff, decreasing = T), ]
    sig.genes$fdr = p.adjust(sig.genes$f.stat.p.val.value, method = 'fdr')
    sig.genes
    ## adjusted pvalues
    
    ##select significant ones
    sig.genes = sig.genes[sig.genes$coeff.p.val < 0.05, ]
    ##select significant ones
    
    ## line graphs
    sub_dataset = dataset[,c(outcome, rownames(sig.genes))]
    sub_dataset
    library(tidyr)
    sub_dataset = gather(sub_dataset, genes, measurements, rownames(sig.genes))
    ggplot(sub_dataset, aes(measurements, sub_dataset[,outcome])) +
      geom_point() +
      geom_smooth(method='lm') +
      facet_grid(. ~ genes)
    ## line graphs
    
    ## dot plots
    sig.genes$genes = factor(rownames(sig.genes), levels = c(rownames(sig.genes)))
    jpeg(paste('bmi/results/dotplot/', graph_title, '.jpeg', sep = ''),
         units="in", width=5, height=5, res=500)
    gg.res = ggplot(sig.genes, aes(coeff, genes)) + 
      geom_point(aes(size = -log10(fdr), fill = coeff), 
                 colour = 'black', shape = 21) + 
      scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
      ylab('') + xlab('') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
      labs(size = '-log10(FDR)', fill = 'Estimates') +
      geom_vline(xintercept=c(0), linetype="dotted") +
      ggtitle(graph_title)
    print(gg.res)
    dev.off()
    ## dot plots
    
    ## csv output
    write.csv(sig.genes, paste('bmi/results/univariate_csv/', graph_title, '_univariate.csv', sep = ''))
    ## csv output
    ### univariate linear
  }
  
  
  #################### glmnet#########################################
  library(glmnet)
  cv <- cv.glmnet(as.matrix(scale(dataset %>% select(-c(outcome)))),
                  dataset[,outcome], nfolds = nrow(dataset))
  
  summary(cv)
  jpeg(paste('bmi/results/glmnet_cv/', graph_title, '.jpeg', sep = ''),
       units="in", width=8, height=8, res=500)
  plot(cv)
  dev.off()
  
  return_features = function(model, lambda){
    coeff = coef(model, s = lambda)
    top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
    top_features = top_features[order(top_features$coefficient),]
    # print(top_features)
    return(top_features)
  }
  min.coeff = return_features(cv, 'lambda.min') 
  return_features(cv, 'lambda.1se')
  write.csv(min.coeff, paste('bmi/results/glmnet_csv/', graph_title, '_glmnet.csv', sep = ''))
}


sspg = read.csv('Re__Obesity_data_set_2/obesityData_luminexWithSSPG.csv', row.names = 1)
sspg = t(sspg)
head(sspg)
sspg_deltasspg = sspg %>% data.frame()%>% select(c(1,2,4,5, grep('_FCH', colnames(sspg))))
sspg_deltaweight = sspg %>% data.frame()%>% select(c(1,3,4,5, grep('_FCH', colnames(sspg))))
sspg_bmi = sspg %>% data.frame()%>% select(c(1,4,5, grep('_FCH', colnames(sspg))))
# head(sspg)
# sspg[1:10,1:10]
pipe(sspg_deltasspg, 'percdeltaSSPG', c('Sexismale', 'Age'), 'sspg_deltasspg_controlled')
pipe(sspg_deltasspg, 'percdeltaSSPG', c(), 'sspgsspg_deltasspg')
pipe(sspg_deltaweight, 'percdeltawt', c('Sexismale', 'Age'), 'sspg_deltaweight_controlled')
pipe(sspg_deltaweight, 'percdeltawt', c(), 'sspg_deltaweight')
pipe(sspg_bmi, 'BMI1', c('Sexismale', 'Age'), 'sspg_bmi')

# bmi_all = read.csv('Re__Obesity_data_set_2/bmi_all.csv', row.names = 1)
# bmi_all = t(bmi_all)
# bmi_all = bmi_all %>% data.frame()%>% select(-c(3,5,6))
# bmi_all = bmi_all %>% select(-grep('^mod_', colnames(bmi_all)))
# pipe(bmi_all, 'BMI', c('Gender', 'Age', 'CMV', 'EBV'), 'bmi_all')
#
# kip.cytokines = read.csv('Re__Obesity_data_set/full_1KIP.csv', row.names = 1)
# sub.kip.cytokines = kip.cytokines %>% select(-c(1,2,3, 8:53))
# pipe(sub.kip.cytokines, 'BMI', c('CMV', 'GENDER', 'AGE'), '1kip.cytokines')
# 
kip.cytokines = read.csv('Re__Obesity_data_set/full_1KIP.csv', row.names = 1)
sub.kip.cytokines = kip.cytokines %>% select(-c(1,2,3, 8:53))
sub.kip.cytokines = sub.kip.cytokines[grep('SLV', rownames(sub.kip.cytokines)),]
dim(sub.kip.cytokines)
# pipe(sub.kip.cytokines, 'BMI', c('CMV', 'GENDER', 'AGE'), '1kip.cytokines.SLV')
# 
kip.cytokines = read.csv('Re__Obesity_data_set/full_1KIP.csv', row.names = 1)
sub.kip.cytokines = kip.cytokines %>% select(-c(1,2,3, 8:53))
sub.kip.cytokines = sub.kip.cytokines[grep('CFS', rownames(sub.kip.cytokines)),]
dim(sub.kip.cytokines)
# pipe(sub.kip.cytokines, 'BMI', c('CMV', 'GENDER', 'AGE'), '1kip.cytokines.CFS')

# montoya = read.csv('Re__Obesity_data_set_2/montoya_GE.csv', row.names = 1)
# pipe(montoya, 'bmi_montoya...1.', c("CMV", "EBV","GENDER","AGE"  ), 'montoya')

