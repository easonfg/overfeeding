sspg = read.csv('Re__Obesity_data_set_2/obesityData_luminexWithSSPG.csv', row.names = 1)
sspg = t(sspg)
sspg = sspg %>% data.frame()%>% select(c(1,4,5, grep('_FCH', colnames(sspg))))
sspg


# bmi.lm = lm(BMI ~., data = combined.data)
# summary(bmi.lm)

library(glmnet)
cv <- cv.glmnet(as.matrix(scale(sspg %>% select(-BMI1))),
                sspg$BMI1, nfolds = nrow(sspg))

single <- glmnet(as.matrix(scale(sspg %>% select(-BMI1))),
                sspg$BMI1)

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
return_features(cv, 'lambda.min') %>% View()
return_features(cv, 'lambda.1se')
