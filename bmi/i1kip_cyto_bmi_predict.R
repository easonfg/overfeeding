kip.cytokines = read.csv('Re__Obesity_data_set/full_1KIP.csv', row.names = 1)

# kip.phenotypes[common.names,]
# View(kip.phenotypes[common.names,])

sub.kip.cytokines = kip.cytokines %>% select(-c(1,2,3, 8:53))
sub.kip.cytokines

pca.data = sub.kip.cytokines
pca.data$rnames = sapply(rownames(pca.data), function(x) strsplit(x, '_')[[1]][1])
pca.data$rnames = sapply(pca.data$rnames, function(x) strsplit(x, '0')[[1]][1])
length(pca.data$rnames)
(pca.data$rnames)
pc.res = prcomp(pca.data[,5:(ncol(pca.data)-1)], scale = T)
fviz_pca_ind(pca.res,
             col.ind = pca.data$rnames,
             # repel = TRUE     # Avoid text overlapping
)
dim(pca.data)
(pc.res$x)


# bmi.lm = lm(BMI ~., data = combined.data)
# summary(bmi.lm)

library(glmnet)
cv <- cv.glmnet(as.matrix(scale(sub.kip.cytokines %>% select(-BMI))),
                sub.kip.cytokines$BMI, nfolds = nrow(sub.kip.cytokines))

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
