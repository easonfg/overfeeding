
## regression model
genex_s= as.data.frame(genex_s)

genex_s$BMI.30.Bin = as.factor(ifelse(genex$bmi_montoya...1. >=30, 1, 0))
genex_s$AGE = genex$AGE
genex_s$GENDER = as.factor(genex$GENDER)


glm.gene = sapply(colnames(genex_s[,1:6149]),function(x){
  print(x)
  fit = lm(as.formula(paste(x, '~BMI.30.Bin + AGE + GENDER ')),data = genex_s)
  #browser()
  summary(fit)$coefficient[2,c(1,4)]
})

# 
 fit = lm( A2LD1 ~ BMI.30.Bin +AGE * GENDER,data = genex_s)
# 
# grep('BMI.30.Bin', colnames(genex_s))

#fit = glm(as.formula(paste( 'BMI.30.Bin~', x,'+AGE + GENDER ')),data = genex_s, family= 'binomial') 
#fit = glm(BMI.30.Bin~A2LD1 +AGE +GENDER, data= genex_s, family='binomial')
# 
#fit = lm(A2LD1 ~BMI.30.Bin+AGE +GENDER, data= genex_s )

glm.gene = as.data.frame(t(glm.gene))
glm.gene = glm.gene%>% rename('p.value'= 'Pr(>|t|)')
glm.gene$p.adjust = p.adjust(glm.gene$p.value,method='fdr')
attach(glm.gene)
glm.gene = glm.gene[order(p.value),,drop =F]


library(EnhancedVolcano)

volcano.plot = EnhancedVolcano(glm.gene,
                               lab = rownames(glm.gene),
                               pCutoff = 0.05,
                               FCcutoff = 0,
                               x = 'Estimate',
                               y = 'p.adjust',
                               legendPosition = 'bottom')
print(volcano.plot)







