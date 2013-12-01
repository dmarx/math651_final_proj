library(DMwR)
data(algae)
dat = knnImputation(algae)

###################################################

# general data epxloration

pairs(dat) # data overload
# Let's focus on numeric predictors
dim(dat) # 20018
names(dat)
pairs(dat[,4:11])
# Shared outlier in CL, NO3, NH4
# oPO4 and PO4 appear to be highly correlated
with(dat, cor(oPO4, PO4)) # 0.9149
symnum(cor(dat[,4:18])) # This is a bit nicer than "image"
par(mfrow=c(1,1))
image(t(cor(dat[,4:18])))

# Wonder if there's any correlation in the responses. Not that it matters
pairs(dat[,12:18])

cov.numeric = cov(dat[,-c(1:3)])
image(cov.numeric) # one region is EXTREMELY highly correlated.
max(cov.numeric) #  3813835
cov.numeric == max(cov.numeric) # 65
# The high covariance is the variance of NH4, which is stupid high.
var(dat$NH4) # yeah... wtf is this?

summary(dat$NH4)
plot(dat$NH4) # we've got a serious outlier in the NH4 data
which(dat$NH4==max(dat$NH4)) # 153
var(dat$NH4[-153]) # 1014321... variance reduced by a factor of 3, still pretty bad.

# Let's drop the variance terms (diagonals) of the covariance matrix to focus
# on intervariable covariances
cov.intervar = cov.numeric
for(i in 1:nrow(cov.numeric)){
  cov.intervar[i,i] = 0
}
image(abs(cov.intervar)) # ok, we've got one variable that's tight with 2 others

cov.intervar
# looks like it's PO4 covarying with NH4 and oP04

par(mfrow=c(1,1))
# I think this'll be easier to read if I build a separate plot for each column (variable)
#lapply(data.frame(abs(cov.intervar)), function(column) image(t(cbind(rep(0,length(column)),rev(column)))))
lapply(data.frame(abs(cov.intervar)), function(column) image(t(as.matrix(rev(column)))))

abs(cov.intervar)

par(mfrow=c(2,8))
vars = colnames(cov.intervar)
n=length(vars)
for(i in 1:n){
  image(t(rev(abs(cov.intervar[,i]))), axes=F, main=vars[i])
  axis(side=2, labels = vars, at=seq(1,0, length.out=n), las=1)
  text(x=rep(0, n), y=seq(1,0, length.out=n),
       labels=round(abs(cov.numeric[,i]))
       )
}
par(mfrow=c(1,1))



dat.no.out = dat[-153,]

par(mfrow=c(2,4))
for(i in 12:18){
  #plot(dat[,8], dat[,i], main=colnames(dat)[i])
  plot(dat.no.out[,8], dat.no.out[,i], main=colnames(dat)[i])
}
vars[8]

length(vars)
dim(dat)


###############################################################################

# Outlier analysis

# Let's first look for global outliers using mahalanobis distance just to get a feel.
# i.e. points that are outliers w.r.t several variables.
mal = mahalanobis(dat[,4:18], center=colMeans(dat[,4:18]), cov=cov(dat[,4:18]))
par(mfrow=c(1,1))
plot(mal)
summary(mal)
boxplot(mal)
hist(mal)
summary(mal) # 75% cutoff is 17.11
p = seq(0,1,.01)
q = quantile(mal, p)
plot(p,q)
abline(v=.9)
abline(v=.95) # even better cutoff.
# Looks like 90% quantile is a good cut off for the most typical data.
# If we dropped everything above this threshhold, how much would we lose?

q[p==.9] # 26.84
# q[p==.95] # This isn't working for some reason.
q[96] # 40.1
sum(mal>=q[p==.9]) # 20 points out of 200. Oh right, duh. 10% of our data.
# Not saying we necessarily *should* drop all 20 of these points, but we should probably
# investigate them to see why they're so strange relative to the rest of our data.

out.ix = which(mal>=q[91])  # .9 cutoff
out.ix2 = which(mal>=q[96]) # .95
# not surprisingly, our friend 153 is in here.
mal[153] # 177.7868 . Yup, this is our worst outlier.

#################################################################################

# Variable selection

# Can we do simultaneous inference in R? Why yes, yes we can!
mult.model = lm(a1+a2+a3+a4+a5+a6+a7~., dat)
summary(mult.model)
mult.step = step(mult.model)
summary(mult.step) 
# backward selection gives model:
# a1 + a2 + a3 + a4 + a5 + a6 + a7 ~ speed + mxPH + NO3 + oPO4

library(leaps)
# leaps can't do multiple inference
# we'll need to apply alg to each response variable separately
leaps(x=dat$a1, dat)
dat[,-c("a1")]


1 in c(1,2)
