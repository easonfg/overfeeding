bmi_all = read.csv('Re__Obesity_data_set_2/bmi_all.csv', row.names = 1)
bmi_all = t(bmi_all)
bmi_all = bmi_all %>% data.frame()%>% select(-c(3,5,6))
bmi_all = bmi_all %>% select(-grep('^mod_', colnames(bmi_all)))


# bmi.lm = lm(BMI ~., data = combined.data)
# summary(bmi.lm)

library(glmnet)
cv <- cv.glmnet(as.matrix(scale(bmi_all %>% select(-BMI))),
                bmi_all$BMI,nfolds = nrow(bmi_all))

summary(cv)
plot(cv)

return_features = function(model, lambda){
  coeff = coef(model, s = lambda)
  top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
  top_features = top_features[order(top_features$coefficient),]
  # print(top_features)
  return(top_features)
}
return_features(cv, 'lambda.min') %>% View()
return_features(cv, 'lambda.1se')


############################################################
### multiple linear
bmi_all = bmi_all %>% select(-grep('^mod_', colnames(bmi_all)))
linear.reg = lm(BMI ~., data = data.frame(scale(bmi_all)))
summary(linear.reg)

### multiple linear

### univariate linear

lin.f.stat.p.val.coeff = sapply(colnames(bmi_all %>% select(-'BMI')), function(x){
  # uni.linear = lm(as.formula(paste0('BMI ~ ', x)), data = data.frame(bmi_all))
  uni.linear = lm(as.formula(paste0('BMI ~ ', x, '+Gender', '+Age', '+CMV', '+EBV')), data = data.frame(bmi_all))
  uni.sum.res = summary(uni.linear)
  # browser()
  f.stat.p.val = pf(uni.sum.res$fstatistic[1],uni.sum.res$fstatistic[2],uni.sum.res$fstatistic[3],lower.tail=FALSE)
  if (f.stat.p.val < 0.05){
    # browser()
    return (c(f.stat.p.val = f.stat.p.val, coeff = uni.sum.res$coefficients[2,1], coeff.p.val = uni.sum.res$coefficients[2,4]))
  }
})
lin.f.stat.p.val.coeff

sig.genes = lin.f.stat.p.val.coeff[sapply(lin.f.stat.p.val.coeff, function(x) !is.null(x))]
sig.genes = data.frame(t(as.data.frame(sig.genes)))
head(sig.genes)

## adjusted pvalues
sig.genes = sig.genes[order(sig.genes$coeff, decreasing = T), ]
sig.genes$fdr = p.adjust(sig.genes$f.stat.p.val.value, method = 'fdr')
sig.genes
## adjusted pvalues

## line graphs
sub_sspg = bmi_all[,c('BMI', rownames(sig.genes))]
sub_sspg
library(tidyr)
sub_sspg = gather(sub_sspg, genes, measurements, rownames(sig.genes))
ggplot(sub_sspg, aes(measurements, BMI)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(. ~ genes)
## line graphs

## dot plots
sig.genes$genes = factor(rownames(sig.genes), levels = c(rownames(sig.genes)))
jpeg(paste('bmi/results/', deparse(substitute(bmi_all)), '.jpeg', sep = ''),
     units="in", width=5, height=5, res=500)
gg.res = ggplot(sig.genes, aes(coeff, genes)) + 
  geom_point(aes(size = -log10(fdr), fill = coeff), 
             colour = 'black', shape = 21) + 
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(size = '-log10(FDR)', fill = 'Estimates') +
  geom_vline(xintercept=c(0), linetype="dotted") +
  ggtitle('bmi_all')
print(gg.res)
dev.off()
## dot plots

### univariate linear


############## vaccine response ##############################
bmi_all = read.csv('Re__Obesity_data_set_2/bmi_all.csv', row.names = 1)
bmi_all = t(bmi_all)
bmi_all = bmi_all %>% data.frame()%>% select(-c(5,6))

# bmi.lm = lm(BMI ~., data = combined.data)
# summary(bmi.lm)

library(glmnet)
cv <- cv.glmnet(as.matrix(scale(bmi_all %>% select(-resp))),
                bmi_all$resp, nfolds = nrow(bmi_all))

summary(cv)
plot(cv)

return_features = function(model, lambda){
  coeff = coef(model, s = lambda)
  top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
  top_features = top_features[order(top_features$coefficient),]
  # print(top_features)
  return(top_features)
}
return_features(cv, 'lambda.min') %>% View()
return_features(cv, 'lambda.1se')

############## response ##############################
############## response ##############################


