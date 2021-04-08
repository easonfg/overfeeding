data("Orthodont", package="nlme")

library(ggplot2)
theme_set(theme_bw())
pl <- ggplot(data=Orthodont) + geom_point(aes(x=age,y=distance), color="red", size=3) 
pl

lm1 <- lm(distance~age, data=Orthodont)
summary(lm1)
predict(lm1)

pl + geom_line(aes(x=age,y=predict(lm1)))
pl + geom_line(aes(x=age,y=predict(lm1))) + facet_grid(.~ Sex )

## same slope, different intercept
lm2 <- lm(distance~age+Sex, data=Orthodont)
summary(lm2)

Orthodont$pred.lm2 <- predict(lm2)
Orthodont$pred.lm2

pl + geom_line(data=Orthodont,aes(x=age,y=pred.lm2)) + facet_grid(.~ Sex )

## different slope, same intercept
lm3 <- lm(distance~age:Sex , data=Orthodont)
summary(lm3)

Orthodont$pred.lm3 <- predict(lm3)
pl + geom_line(data=Orthodont,aes(x=age,y=pred.lm3)) + facet_grid(.~ Sex )

## different slope, different intercept
lm4 <- lm(distance~age:Sex+Sex, data=Orthodont)
summary(lm4)

Orthodont$pred.lm4 <- predict(lm4)
pl + geom_line(data=Orthodont,aes(x=age,y=pred.lm4)) + facet_grid(.~ Sex )