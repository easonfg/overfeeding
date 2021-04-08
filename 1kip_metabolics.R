kip.phenotypes = read.csv('Re__Obesity_data_set/full_1KIP.csv', row.names = 1)
kip.metabolomics = read.csv('Re__Obesity_data_set/metabolomics_1KIP.csv', row.names = 1)
kip.metabolomics = t(kip.metabolomics)

common.names = intersect(rownames(kip.phenotypes), rownames(kip.metabolomics))

combined.data = cbind(kip.phenotypes[common.names, c('CMV', 'GENDER', 'AGE', 'BMI')], kip.metabolomics[common.names,])
# kip.phenotypes[common.names,]
# View(kip.phenotypes[common.names,])

bmi.lm = lm(BMI ~., data = combined.data)
summary(bmi.lm)

library(glmnet)
cv <- cv.glmnet(as.matrix(combined.data %>% select(-BMI)),
                combined.data$BMI)

summary(cv)
plot(cv)

return_features = function(model, lambda){
  coeff = coef(model, s = lambda)
  top_features = data.frame(name = coeff@Dimnames[[1]][coeff@i + 1], coefficient = coeff@x)
  top_features = top_features[order(top_features$coefficient),]
  print(top_features)
  return(top_features)
}
return_features(cv, 'lambda.min')
return_features(cv, 'lambda.1se')
