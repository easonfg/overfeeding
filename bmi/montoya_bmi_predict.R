montoya = read.csv('Re__Obesity_data_set_2/montoya_GE.csv', row.names = 1)

# bmi.lm = lm(BMI ~., data = combined.data)
# summary(bmi.lm)

library(glmnet)
cv <- cv.glmnet(as.matrix(scale(montoya %>% select(-bmi_montoya...1.))),
                  montoya$bmi_montoya...1., nfolds = nrow(montoya))

elastic.net.reg <- glmnet(as.matrix(scale(montoya %>% select(-bmi_montoya...1.))),
                  montoya$bmi_montoya...1.)

elastic.coeff = coef(elastic.net.reg)[,ncol(coef(elastic.net.reg))]
elastic.coeff = elastic.coeff[elastic.coeff > 0]
elastic.coeff[order(elastic.coeff, decreasing = T)]

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
