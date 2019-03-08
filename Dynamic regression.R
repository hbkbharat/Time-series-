## Dyanamic regression example 
## And check for seasonality & serial correlation
##stationarity
library(ggplot2)
df <- read.csv("2assign.csv", header = TRUE)
df1<- read.csv("plot.csv")
head(df)

fed <- df$fedrate
ten <- df$tenyr
cpi <- df$cpi
par(mfrow=c(2,2))
plot(fed,type='l',col=4)
plot(ten,type='l',col=3)
plot(cpi,type='l',col=2)
##log convert
lfed <- log(fed)
lten <- log(ten)


##
library(ggplot2)
ggplot(df1,aes(x=time,y = value)) + geom_line()

par(mfrow=c(1,1))
plot(lsnp,type='l',col=4)
plot(lten,type='l',col=2)
plot(lwti,type='l',col=3)



## check for stationary
source(file="intord.R")
intord(fed)
intord(ten)
intord(cpi)


##convert to ts
install.packages("xts")
library(xts)
fedout = ts(fed,frequency = 12,start=c(1969,8),end=c(2018,12)) #monthly data starting august 1969
plot(fedout,type='l',col=4)
tenout = ts(ten,frequency = 12,start=c(1969,8),end=c(2018,12)) #monthly data starting august 1969
plot(tenout,type='l',col=3)
cpiout = ts(cpi,frequency = 12,start=c(1969,8),end=c(2018,12)) #monthly data starting august 1969
plot(cpiout,type='l',col=7)

## all seem to have frst lag stationarity
#convert to frst difference
dfed<-diff(fedout) 
dcpi<-diff(tenout)
dten<-diff(cpiout)


## dynamic regression :-


cc = cbind(fedout,tenout,cpiout)
summary(cc)
length(cpiout)
library(dynlm)
r2= dynlm(dten~L(dten,1:12)+L(dcpi,0:12)+L(dfed,0:12),start=c(1969,8), end=c(2018,12))
summary(r2)
AIC(r2);BIC(r2)
## kind of close to best below
r1 = dynlm(dten~L(dten,7:12)+L(dten,1:3)+L(dcpi,0:6)+L(dfed,1:12),start=c(1969,8), end=c(2018,12))
summary(r1)
AIC(r1);BIC(r1)
## lets try other combos
r3 = dynlm(dten~L(dten,7:12)+L(dten,1:3)+L(dcpi,0:1)+L(dcpi,6:6)+L(dfed,1:8)+L(dfed,11:12),start=c(1969,8), end=c(2018,12))
summary(r3)
AIC(r3);BIC(r3)
## check combo using Anova

anova(r1,r2)
##check with RMSE

RMSPE(r1)
RMSPE(r2)
RMSPE(r3)

### check for seasonality
dten_lags = embed(dten,11)
dtent = dten_lags[,1] ## this will be our y if we are using lags 
dtent_1 = dten_lags [,2]
dtent_2 = dten_lags [,3]
dtent_3 = dten_lags [,4]
dtent_4 = dten_lags [,5]
dtent_5 = dten_lags [,6]
dtent_6 = dten_lags [,7]
dtent_7 = dten_lags [,8]
dtent_8 = dten_lags [,9]
dtent_9 = dten_lags [,10]
dten_10 = dten_lags[,11]
rs0 <- lm(dtent~dtent_1+dtent_2+dtent_3+dtent_4+dtent_5+dtent_6+dtent_7+dtent_8+dtent_9+dten_10)  ## frst lag on second lag
summary(rs0)
dten_fit = r0$fitted
pred_act = cbind(dtent,dten_fit)
matplot(pred_act,type='l',col=1:2)


rs1 <- lm(dtent~dten_10) ## frst lag on 11th lag
summary(rs1)
length(dten_1)
lynx_fit_sdums = r1$fitted
pred_act = cbind(lynxt,lynx_fit)
matplot(pred_act,type='l',col=1:2)

##only with seasonal dummies in the dynlm model
n = length(dten)
source("seas.R")
s = seas(n,11)  # number of obs., frequency of data (12 = monthly)
str(s) 
rs = dynlm(dten~L(dten,7:12)+L(dten,1:3)+L(dcpi,0:6)+L(dfed,1:12)+season(dten),start=c(1969,8), end=c(2018,12))
summary(rs)
anova(rs,r1)
rs_fit <- rs$fitted.values

##only with lags 11th and 2nd 
rr<- dynlm(dtent~ dten_10 +dtent_1)
summary(rr)
rr_fit<- rr$fitted.values
pred_act = cbind(dtent,rr_fit)
matplot (pred_act,type = "l",col = 1:2)
#with both lags and dummies (dynlm doesnt work here so ols)
n = length(dtent)
s = seas(n,10)  # number of obs., frequency of data (12 = monthly)
str(s)
rrs<- lm(dtent~ dten_10 +dtent_1 +s$seas[,1:9])
summary(rrs)
length(dten_10)
length(dtent_1)
length(s$seas[,1:9])
dten_fit_dums_lags = rrs$fitted
pred_act = cbind(dtent,rr_fit,rs_fit,dten_fit_dums_lags)
matplot(pred_act,type='l',col=c("black", "red", "blue", "green", "purple"))

## compare all models
AIC(rs0);BIC(rs0)
AIC(rs);BIC(rs)
AIC(rr);BIC(rr)
AIC(rrs);BIC(rrs)

RMSPE(rs)
RMSPE(rs0)
RMSPE(rrs)

#serial correlation check
### serial correlation (start with ur best model of Q2)
r1 = dynlm(dten~L(dten,1:24)+L(dcpi,0:12)+L(dfed,1:12),start=c(1969,8), end=c(2018,12))
summary(r1)
## BG test
bgtest(r1)
bgtest(r1,order=2)
bgtest(r1,order=3)
bgtest(r1,order=4)
bgtest(r1,order=8)
bgtest(r1,order=12)

regmodel = r1
bgs = rep(0,20)
for (i in 1:20){
  b=bgtest(regmodel,order = i)
  bgs[i] = b$p.value 
}
bgs
## if the above p values are small, there is evidence of serial correlation. to avoid that, add more lags

## box test
res = r1$resid
blt  = rep(0,20)
for (i in 1:20){
  b = Box.test(res,lag=i,type ="Ljung-Box")
  blt[i]=b$p.value
}
blt
## if the above p values are small, there is evidence of serial correlation. to avoid that, add more lags

###  granger causality
##CPI check first
r1 = dynlm(dten~L(dten,7:12)+L(dten,1:3)+L(dcpi,0:0)+L(dcpi,6:6)+L(dfed,1:12),start=c(1969,8), end=c(2018,12))
summary(r1)
r_dcpi<-dynlm(dten~L(dten,7:12)+L(dten,1:3)+L(dfed,1:12),start=c(1969,8), end=c(2018,12)) ## removing cpi variable
anova(r_dcpi,r1,test="F")

## fed rate check
r1 = dynlm(dten~L(dten,7:12)+L(dten,1:3)+L(dcpi,0:0)+L(dcpi,6:6)+L(dfed,1:12),start=c(1969,8), end=c(2018,12))
summary(r1)
r_dfed<-dynlm(dten~L(dten,7:12)+L(dten,1:3)+L(dcpi,0:0)+L(dcpi,6:6),start=c(1969,8), end=c(2018,12)) ## removing fed ariable
anova(r_dfed,r1,test="F")

