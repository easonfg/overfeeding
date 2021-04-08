set.seed(123)
# x <- c(rpois(50,100), rpois(50,75))
x <- c(rpois(100,100))
y <- rep(factor(c("A", "B")), each=50)
df = data.frame(x, y, age = sample(1:2000, 50), group1 = c(rpois(50,90), rpois(50,65)), cov = c(rpois(50,95), rpois(50,70)))
df

(lm.mdl <- lm(x~y+age+group, data = df))
summary(lm.mdl)
betas = lm.mdl$coefficients
# (betas <- lm(x~y)$coefficients)

(means <- tapply(x,y,mean))
(means)


betas[1] + betas[2]  
betas[1]
betas[2]
betas

means[2]

log(sum(betas)/betas[1], 2) 
log(sum(betas)/betas[1], 2) 

lm.mdl2 = lm(x~y+age+group1+cov, data = df%>%mutate(x = log(x,2), age = log(age,2), group1 = log(group1,2), cov = log(cov,2)))
lm.mdl2$coefficients
summary(lm.mdl2)

wt.res = sapply(2:length(coef(lm.mdl2)), function(x){
  temp.wt.res = wald.test(b=coef(lm.mdl2), Sigma=vcov(lm.mdl2), Terms=c(x))
  temp.wt.res$result$chi2[3]
})
wt.res
