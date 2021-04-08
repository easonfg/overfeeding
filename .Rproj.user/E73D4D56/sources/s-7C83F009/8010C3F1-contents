library(dplyr)

###### delta sspg ##############################
# # sspg = read.csv('Re__Obesity_data_set_2/obesityData_luminexWithSSPG.csv', row.names = 1)
# sspg = read.csv('Re__Obesity_data_set_2/overfeeding_Cytokine_Chemokine_luminexWithSSPG.csv', 
#                 skip = 1, row.names = 1)
# # row.names = 1)
# sspg
# str(sspg)
# sspg = t(sspg)
# sspg[1:3, 1:10]
# str(sspg)
# # weigts
# sspg = sspg %>% data.frame()%>% select(c(1,2,4,5, grep('_FCH', colnames(sspg))))
# head(sspg)
# 
# 
# # bmi.lm = lm(BMI ~., data = combined.data)
# # summary(bmi.lm)
# 
# ### pls
# library(caret)
# library(pls)
# set.seed(123)
# model <- train(
#   percdeltaSSPG~., data = scale(sspg), method = "pls",
#   scale = TRUE,
#   trControl = trainControl("cv", number = 10),
#   tuneLength = 10
# )
# # Plot model RMSE vs different values of components
# plot(model)
# # Print the best tuning parameter ncomp that
# # minimize the cross-validation error, RMSE
# model$bestTune
# summary(model$finalModel)
# ### pls
# 
# 
# ### multiple linear
# linear.reg = lm(percdeltaSSPG ~., data = data.frame(scale(sspg)))
# summary(linear.reg)
# 
# ### multiple linear
# 
# ### univariate linear
# lin.f.stat.p.val.coeff = sapply(colnames(sspg %>% select(-'percdeltaSSPG')), function(x){
#   # uni.linear = lm(as.formula(paste0('percdeltaSSPG ~ ',x, '+BMI1')), data = data.frame(sspg))
#   uni.linear = lm(as.formula(paste0('percdeltaSSPG ~ ',x, '+Age +Sexismale')), data = data.frame(sspg))
#   # uni.linear = lm(as.formula(paste0('percdeltaSSPG ~ ', x)), data = data.frame(scale(sspg)))
#   uni.sum.res = summary(uni.linear)
#   # browser()
#   f.stat.p.val = pf(uni.sum.res$fstatistic[1],uni.sum.res$fstatistic[2],uni.sum.res$fstatistic[3],lower.tail=FALSE)
#   if (f.stat.p.val < 0.05){
#     # browser()
#     return (c(f.stat.p.val = f.stat.p.val, coeff = uni.sum.res$coefficients[2,1]))
#   }
# })
# lin.f.stat.p.val.coeff
# 
# sig.genes = lin.f.stat.p.val.coeff[sapply(lin.f.stat.p.val.coeff, function(x) !is.null(x))]
# sig.genes = data.frame(t(as.data.frame(sig.genes)))
# head(sig.genes)
# 
# ## adjusted pvalues
# sig.genes = sig.genes[order(sig.genes$coeff, decreasing = T), ]
# sig.genes$fdr = p.adjust(sig.genes$f.stat.p.val.value, method = 'fdr')
# sig.genes
# ## adjusted pvalues
# 
# ## line graphs
# sub_sspg = sspg[,c('percdeltaSSPG', rownames(sig.genes))]
# sub_sspg
# library(tidyr)
# sub_sspg = gather(sub_sspg, genes, measurements, rownames(sig.genes))
# ggplot(sub_sspg, aes(measurements, sub_sspg[,'percdeltaSSPG'])) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   facet_grid(. ~ genes)
# ## line graphs
# 
# ## dot plots
# sig.genes$genes = factor(rownames(sig.genes), levels = c(rownames(sig.genes)))
# jpeg(paste('bmi/results/', deparse(substitute(sspg)), '.jpeg', sep = ''),
#      units="in", width=5, height=5, res=500)
# gg.res = ggplot(sig.genes, aes(coeff, genes)) + 
#   geom_point(aes(size = -log10(fdr), fill = coeff), 
#              colour = 'black', shape = 21) + 
#   scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
#   ylab('') + xlab('') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
#   labs(size = '-log10(FDR)', fill = 'Estimates') +
#   geom_vline(xintercept=c(0), linetype="dotted") +
#   ggtitle(deparse(substitute(sspg)))
# print(gg.res)
# dev.off()
# ## dot plots
# ### univariate linear
# 
# #################### glmnet#########################################
# library(glmnet)
# cv <- cv.glmnet(as.matrix(scale(sspg %>% select(-percdeltaSSPG))),
#                 sspg$percdeltaSSPG, nfolds = nrow(sspg))
# 
# dim(sspg)
# 
# summary(cv)
# plot(cv)
# 
# return_features = function(model, lambda){
#   coeff = coef(model, s = lambda)
#   top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
#   top_features = top_features[order(top_features$coefficient),]
#   # print(top_features)
#   return(top_features)
# }
# return_features(cv, 'lambda.min') 
# return_features(cv, 'lambda.1se')

####################delta sspg#########################################

####################delta weight#########################################
library(dplyr)
# sspg = read.csv('Re__Obesity_data_set_2/obesityData_luminexWithSSPG.csv', row.names = 1)
sspg = read.csv('Re__Obesity_data_set_2/overfeeding_Cytokine_Chemokine_luminexWithSSPG.csv', 
                skip = 1, row.names = 1)
                 # row.names = 1)
sspg
str(sspg)
sspg = t(sspg)
sspg[1:3, 1:10]
str(sspg)
# weigts
sspg = sspg %>% data.frame()%>% select(c(1,3,4,5, grep('_FCH', colnames(sspg))))

### linear regression
linear.reg = lm(percdeltawt ~., data = data.frame(scale(sspg)))
linear.reg
summary(linear.reg)
uni.linear = lm(as.formula(paste0('percdeltawt ~ ', colnames(sspg)[grep('IL.12P70_FCH', colnames(sspg))])), data = data.frame(sspg))
uni.linear = lm(as.formula(paste0('percdeltawt ~ BMI1+', colnames(sspg)[grep('IL.12P70_FCH', colnames(sspg))])), data = data.frame(sspg))
uni.linear = lm(as.formula(paste0('percdeltawt ~ Age+', colnames(sspg)[grep('IL.12P70_FCH', colnames(sspg))])), data = data.frame(sspg))
uni.linear = lm(as.formula(paste0('percdeltawt ~ Age + Sexismale+', colnames(sspg)[grep('IL.12P70_FCH', colnames(sspg))])), data = data.frame(sspg))
uni.linear = lm(as.formula(paste0('percdeltawt ~ BMI1+', colnames(sspg)[grep('IL.12P70_FCH', colnames(sspg))])), data = data.frame(sspg))
uni.linear = lm(as.formula(('percdeltawt ~ BMI1')), data = data.frame(sspg))
# car::vif(uni.linear)
summary(uni.linear)


### univariate linear
lin.f.stat.p.val.coeff = sapply(colnames(sspg %>% select(-'percdeltawt')), function(x){
  # uni.linear = lm(as.formula(paste0('percdeltawt ~ ',x)), data = data.frame(sspg))
  # uni.linear = lm(as.formula(paste0('percdeltawt ~ ',x '+BMI1')), data = data.frame(sspg))
  uni.linear = lm(as.formula(paste0('percdeltawt ~',x, '+ Sexismale + Age')), data = data.frame(sspg))
  uni.sum.res = summary(uni.linear)
  # browser()
  f.stat.p.val = pf(uni.sum.res$fstatistic[1],uni.sum.res$fstatistic[2],uni.sum.res$fstatistic[3],lower.tail=FALSE)
  if (f.stat.p.val < 0.05){
    return (c(f.stat.p.val = f.stat.p.val, coeff = uni.sum.res$coefficients[2,1]))
  }
})

sig.genes = lin.f.stat.p.val.coeff[sapply(lin.f.stat.p.val.coeff, function(x) !is.null(x))]
sig.genes = t(as.data.frame(sig.genes))
sig.genes
sub_sspg = sspg[,c('percdeltawt', rownames(sig.genes))]
library(tidyr)
library(ggplot2)
sub_sspg = gather(sub_sspg, genes, measurements, rownames(sig.genes))
ggplot(sub_sspg, aes(measurements, percdeltawt)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(. ~ genes)

# adjusted p values
lin.f.stat.p.val = sapply(colnames(sspg %>% select(-'percdeltawt')), function(x){
  uni.linear = lm(as.formula(paste0('percdeltawt ~ ',x)), data = data.frame(sspg))
  uni.sum.res = summary(uni.linear)
  # browser()
  pf(uni.sum.res$fstatistic[1],uni.sum.res$fstatistic[2],uni.sum.res$fstatistic[3],lower.tail=FALSE)
})
lin.f.stat.p.val[which(lin.f.stat.p.val < 0.05)]
lin.f.stat.p.val
p.adjust(lin.f.stat.p.val, method = 'BH')
lin.f.stat.p.val[p.adjust(lin.f.stat.p.val, method = 'fdr') < 0.05]
### univariate linear


# bmi.lm = lm(BMI ~., data = combined.data)
# summary(bmi.lm)

### glmnet
library(glmnet)
cv <- cv.glmnet(as.matrix(scale(sspg %>% select(-percdeltawt))),
                # sspg$percdeltawt, nfolds = nrow(sspg))
                sspg$percdeltawt, nfolds = nrow(sspg), alpha = 0.5)

dim(sspg)

summary(cv)
plot(cv)

return_features = function(model, lambda){
  coeff = coef(model, s = lambda)
  top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
  top_features = top_features[order(top_features$coefficient),]
  # print(top_features)
  return(top_features)
}
return_features(cv, 'lambda.min') 
return_features(cv, 'lambda.1se')
####################delta weight#########################################