## clean data

### read data
overfeed.data = read.csv('Re__Obesity_data_set/overfeeding_metabolomics.csv', row.names = 1,
                         header = F, stringsAsFactors = FALSE, sep=',', dec='.')


### get rid of last column which contains 1 number
overfeed.data = data.frame(t(overfeed.data[,-ncol(overfeed.data)]))

### change all columns that should be numbers to numeric
# first find exclude all columns that should be factors
num.col = seq(1,ncol(overfeed.data))[!seq(1, ncol(overfeed.data)) %in% c(1,4,5,6)] 
# overfeed.data[1:10,1:10]
# change all remaining to numeric
overfeed.data[,num.col] = sapply(num.col,
                         function(x) as.numeric(as.character(overfeed.data[,x])))

str(overfeed.data )
overfeed.data[1:10,1:10]
is.na(overfeed.data[1:10,1:10])

## change time point to numbers
overfeed.data$TIME.POINT = rep(c(0,2,4,8), nrow(overfeed.data)/4)

library(lme4)
# get the metabolite columns
metabolite_cols = colnames(overfeed.data)[7:ncol(overfeed.data)]

## compare the linear mixed model to linear models and see which one gives the best AIC
compare_lme = function(i){
  # browser()
  lm1.formula = as.formula(paste(metabolite_cols[i], '~ AGE + BMI + GENDER + GROUP.NAME + TIME.POINT'))
  overfeed.lm = lm(lm1.formula, data=overfeed.data)
  
  lme.formula = as.formula(paste(metabolite_cols[i], '~ AGE + BMI + GENDER + GROUP.NAME + TIME.POINT +
        (1 + TIME.POINT|SUBJECT.ID)'))
  overfeed.lme = lmer(lme.formula, data=overfeed.data)
  
  lme2.formula = as.formula(paste(metabolite_cols[i], '~ AGE + BMI + GENDER + GROUP.NAME + TIME.POINT:GROUP.NAME+
       (1 + TIME.POINT|SUBJECT.ID)'))
  overfeed.lme2 = lmer(lme2.formula, data=overfeed.data)
  
  lme3.formula = as.formula(paste(metabolite_cols[i], '~ AGE + BMI + GENDER + TIME.POINT:GROUP.NAME  +
       (1 + TIME.POINT|SUBJECT.ID)'))
  overfeed.lme3 = lmer(lme3.formula, data=overfeed.data)
  
  library(arm)
  # display(overfeed.lm)
  # display(overfeed.lme)
  # display(overfeed.lme2)
  # display(overfeed.lme3)
  # 
  # summary(overfeed.lme)
  # summary(overfeed.lme2)
  # summary(overfeed.lme3)
  # View(overfeed.data)
  
  # print(BIC(overfeed.lm, overfeed.lme, overfeed.lme2, overfeed.lme3))
  print(AIC(overfeed.lm, overfeed.lme, overfeed.lme2, overfeed.lme3))
}

for (i in 1:length(metabolite_cols)){
  compare_lme(i)
}
## linear modesl gave the best AIC


## compare different linear models
compare_lm = function(i){
  
  # browser()
  lm1.formula = as.formula(paste(metabolite_cols[i], '~ AGE + BMI + GENDER + GROUP.NAME + TIME.POINT'))
  overfeed.lm = lm(lm1.formula, data=overfeed.data)
  
  lm2.formula = as.formula(paste(metabolite_cols[i], '~ AGE + BMI + GENDER + GROUP.NAME + TIME.POINT:GROUP.NAME'))
  overfeed.lm2 = lm(lm2.formula, data=overfeed.data)
  
  lm3.formula = as.formula(paste(metabolite_cols[i], '~ AGE + BMI + GENDER + TIME.POINT:GROUP.NAME'))
  overfeed.lm3 = lm(lm3.formula, data=overfeed.data)
  library(arm)
  # display(overfeed.lm)
  # display(overfeed.lme)
  # display(overfeed.lme2)
  # display(overfeed.lme3)
  # 
  # summary(overfeed.lme)
  # summary(overfeed.lme2)
  # summary(overfeed.lme3)
  # View(overfeed.data)
  
  # print(BIC(overfeed.lm, overfeed.lme, overfeed.lme2, overfeed.lme3))
  print(AIC(overfeed.lm, overfeed.lm2, overfeed.lm3))
}

for (i in 1:length(metabolite_cols)){
  compare_lm(i)
}
## linear models without nesting generally give better AIC

## see what variable is the most positve correlation with the expression of metabolite is
i = 500
res.variable = sapply(1:length(metabolite_cols), function(x) {
  # print('###################')
  lm1.formula = as.formula(paste(metabolite_cols[x], '~ AGE + BMI + GENDER + GROUP.NAME + TIME.POINT'))
  overfeed.lm = lm(lm1.formula, data=overfeed.data)
  # summary(overfeed.lm)
  dis.overfeed.coef = overfeed.lm$coefficients[-1]
  # dis.overfeed.coef
  # print(dis.overfeed.coef[which.max(dis.overfeed.coef)])
  return(dis.overfeed.coef[which.max(dis.overfeed.coef)])
})
res.variable
str(overfeed.data)

## see what metabolite changed most positively with time
res.pos = sapply(1:length(metabolite_cols), function(x) {
  # print('###################')
  lm1.formula = as.formula(paste(metabolite_cols[x], '~ AGE + BMI + GENDER + GROUP.NAME + TIME.POINT'))
  overfeed.lm = lm(lm1.formula, data=overfeed.data)
  return(overfeed.lm$coefficients['TIME.POINT'])
})
res.pos[which.max(res.pos)]
metabolite_cols[which.max(res.pos)]

# no.na.overfeed.data = overfeed.data[which(rowSums(is.na(overfeed.data[,7:ncol(overfeed.data)])) == 0),]
# no.na.overfeed.data[1:10,1:10]
# i = 500
# no.na.lm1.formula = as.formula(paste(metabolite_cols[i], '~ AGE + BMI + GENDER + GROUP.NAME + TIME.POINT'))
# no.na.overfeed.lm = lm(lm1.formula, data=no.na.overfeed.data)
# summary(no.na.overfeed.lm)
# display(no.na.overfeed.lm)


### finding overlaps
cytof = read.csv('cytof/NewRAW.csv', row.names = 1)
rownames(cytof)

overfeed = read.csv('Re__Obesity_data_set/overfeeding_metabolomics.csv', header = F)
overfeed

res <- str_match(rownames(cytof), "_\\s*(.*?)\\s*\\.")

ls.strings = sapply(rownames(cytof), function(x) strsplit(x, '_'))
ls.strings = sapply(ls.strings, function(x) x[length(x)])
ls.strings = sapply(ls.strings, function(x) strsplit(x, '\\.')[[1]][1])
ls.strings

strsplit(as.vector(overfeed[1,3]), '_')
feed.ls.strings = sapply(overfeed[1,-ncol(overfeed)], function(x) strsplit(as.vector(x), '_'))
feed.ls.strings = sapply(feed.ls.strings, function(x) x[2])

length(intersect(unique(feed.ls.strings), ls.strings))
setdiff(unique(feed.ls.strings), ls.strings)
setdiff(ls.strings, unique(feed.ls.strings))

