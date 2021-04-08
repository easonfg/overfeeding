
genex=read.csv('Gene_Expression_montoya.csv',header=T, row.names=1)

cytok = read.csv('inflamazome_BloodScreen__Gene_Modules_Metabolites.csv',header=T, row.names=1)
cytok = as.data.frame(t(cytok))

genex_s = scale(genex[,6:6154])

x = as.matrix(genex[,c(1:5)])
y= as.matrix(genex_s)
fit_reg<-regressNPermuteFast(x, y, numPerms = 1000,fullOut=TRUE)
beta_AGE <-fit_reg$r$beta[,5]
beta_Gender <-fit_reg$r$beta[,4]
beta_CMV = fit_reg$r$beta[,2]
beta_EBV = fit_reg$r$beta[,3]

x_ = t(x)
y_ = t(y)
xAGE <- as.matrix(beta_AGE) %*% x_[4,]
xGender <- as.matrix(beta_Gender) %*% x_[3,]
xCMV = as.matrix(beta_CMV) %*% x_[1,]
xEBV = as.matrix(beta_EBV) %*% x_[2,]

y_corr <- y_-as.matrix(xAGE)-as.matrix(xGender)

corr_genex = as.data.frame(y_corr)
#SAM
library(samr)
a = as.matrix(corr_genex)
b = ifelse(genex$bmi_montoya...1. >=30, 2, 1)
d = list(x=a, y=b, geneid=row.names(a), genenames=paste("g", as.character(1:nrow(x)),sep=""), logged2=FALSE)
samr.obj = samr(d, resp.type="Two class unpaired",nperms=1000)
delta=.8
samr.plot(samr.obj,delta)
delta.table <- samr.compute.delta.table(samr.obj)
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
upg = siggenes.table$genes.up
log = siggenes.table$genes.lo
all_out = rbind(upg,log)


plot(all_out[,8]~all_out[,4], xlab = "Score(d)", ylab = "FDR")
abline(h=13,col=3,lty=2)

gsam = rbind(upg,log)
ta = data.frame(apply(gsam[,2:4], 2, function(x) as.numeric(as.character(x))))
gsam = cbind(as.character(gsam[,1]),ta)
gsam = gsam %>% rename('q-value' = "q.value...")
gsam$logFC= log2(gsam$Fold.Change)





library(EnhancedVolcano)

volcano.plot = EnhancedVolcano(gsam,
                               lab = gsam$gene,
                               pCutoff = 0.1,
                               FCcutoff = 0,
                               x = 'Score.d.',
                               y = 'q-value',
                               legendPosition = 'bottom')
print(volcano.plot)








cytok_s = scale(cytok[,9:284])

a = as.matrix(cytok[,c(1,2,4)])
b= as.matrix(cytok_s)
fit_reg<-regressNPermuteFast(x, y, numPerms = 1000,fullOut=TRUE)

beta_AGE <-fit_reg$r$beta[,4]
beta_Gender <-fit_reg$r$beta[,3]

a_ = t(a)
b_ = t(b)

xAGE <- as.matrix(beta_AGE) %*% a_[3,]
xGender <- as.matrix(beta_Gender) %*% a_[2,]

cyto_corr <- b_-as.matrix(xAGE)-as.matrix(xGender)

corr_cytok = as.data.frame(cyto_corr)

###
# install.packages(c("shiny", "MASS", "preprocessCore"), dependencies = TRUE)
# 
# shiny::runGitHub("ABIS", user="giannimonaco")

rownames(y_corr) = make.names(rownames(y_corr),unique = T)

y_corr

write.table(y_corr,file = 'transpose_montoya.txt', sep='\t', quote = FALSE)



