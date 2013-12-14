library(DMwR) # for knnImputation and algae data
library(MASS) # for studres
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
# I think it might be interesting to color my earlier pairs plot using a color gradient
# that highlights outliers based on their mahalanobis distance.

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

# Pairs plot with white-yellow-orange-red palette highlighting our outliers.
pal = colorRampPalette(brewer.pal(9,"YlOrRd"))(max(mal))
plot(dat[,4:11], col=pal[ceiling(mal)])

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

####################################################3

# Unconstrained full first-order model with bivariate interactions
# Focusing on algae a6

names(dat)[1:11]
# Dynamically build up formula
formula1 = 'a6 ~ '
formula2 = 'a6 ~ '
for(n in names(dat)[1:11]){
  formula1 = paste(formula1, n, '+')
  formula2 = paste(formula2, n, '+') # full model, no interactions
  for(n2 in names(dat)[1:11]){
    if(n < n2){ # alphabetic deduplication
      interaction = paste(n,':',n2, sep="")
      formula1 = paste(formula1, interaction, '+')
    }
  }
}

# Get rid of trailing plus sign
formula.full.interact = substr(formula1, 1, nchar(formula1)-2)
formula.full = substr(formula2, 1, nchar(formula2)-2)

mod.full          = lm(formula.full, dat)          # adj.rsq = .355
mod.full.interact = lm(formula.full.interact, dat) # adj.rsq = .575

# Step regression to trim out weakest interactions from full model
mod.full.selectInteract = step(mod.full.interact
                              ,scope=list(end=formula.full, start=formula.full.interact)
                              ,direction='backward')

# Full model with select interactions
# --> all univariate terms kept in model to ensure that output of backwards selection 
#     on interactions model gives valid model. Otherwise, result is all interactions.
summary(mod.full)
summary(mod.full.selectInteract) # .664


# Call:
# lm(formula = a6 ~ season + size + speed + mxPH + mnO2 + Cl + 
# NO3 + NH4 + oPO4 + PO4 + Chla + season:speed + speed:mxPH + 
# mxPH:NO3 + season:mnO2 + size:mnO2 + speed:mnO2 + mnO2:NO3 + 
# mnO2:NH4 + mnO2:oPO4 + size:Cl + speed:Cl + mxPH:Cl + Cl:NH4 + 
# Cl:PO4 + NO3:oPO4 + season:NH4 + size:NH4 + speed:NH4 + NO3:NH4 + 
# NH4:oPO4 + season:oPO4 + speed:oPO4 + oPO4:PO4 + season:PO4 + 
# size:PO4 + speed:PO4 + size:Chla + mxPH:Chla + mnO2:Chla + 
# oPO4:Chla, data = dat)

# Residuals:
# Min       1Q   Median       3Q      Max 
# -12.4749  -3.1048  -0.2255   2.5452  22.7681 

# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -2.397e+01  2.646e+01  -0.906 0.366644    
# seasonspring             -6.493e+00  1.588e+01  -0.409 0.683274    
# seasonsummer             -1.798e+01  1.554e+01  -1.157 0.249269    
# seasonwinter             -2.889e+01  1.466e+01  -1.970 0.050936 .  
# sizemedium                3.471e+01  1.048e+01   3.313 0.001195 ** 
# sizesmall                 3.979e+01  1.254e+01   3.174 0.001875 ** 
# speedlow                  1.265e+02  5.122e+01   2.469 0.014847 *  
# speedmedium              -4.643e+01  2.092e+01  -2.219 0.028197 *  
# mxPH                      3.155e-01  2.425e+00   0.130 0.896667    
# mnO2                      1.899e+00  1.611e+00   1.179 0.240693    
# Cl                       -5.220e-01  3.781e-01  -1.381 0.169706    
# NO3                       8.743e+00  5.907e+00   1.480 0.141261    
# NH4                       4.098e-02  1.618e-02   2.533 0.012491 *  
# oPO4                     -2.498e-01  1.118e-01  -2.235 0.027096 *  
# PO4                       2.652e-01  5.384e-02   4.926 2.51e-06 ***
# Chla                     -1.000e+00  6.966e-01  -1.436 0.153439    
# seasonspring:speedlow    -1.990e+01  8.747e+00  -2.275 0.024547 *  
# seasonsummer:speedlow    -1.655e+01  6.489e+00  -2.550 0.011919 *  
# seasonwinter:speedlow    -6.017e+00  4.984e+00  -1.207 0.229571    
# seasonspring:speedmedium  8.800e-01  4.312e+00   0.204 0.838605    
# seasonsummer:speedmedium  3.014e+00  3.929e+00   0.767 0.444288    
# seasonwinter:speedmedium  1.335e+01  4.000e+00   3.339 0.001098 ** 
# speedlow:mxPH            -1.035e+01  5.646e+00  -1.832 0.069220 .  
# speedmedium:mxPH          3.627e+00  2.211e+00   1.640 0.103378    
# mxPH:NO3                 -1.524e+00  7.534e-01  -2.022 0.045192 *  
# seasonspring:mnO2         6.431e-01  1.390e+00   0.463 0.644290    
# seasonsummer:mnO2         1.909e+00  1.366e+00   1.397 0.164696    
# seasonwinter:mnO2         2.523e+00  1.261e+00   2.000 0.047533 *  
# sizemedium:mnO2          -3.088e+00  9.350e-01  -3.302 0.001239 ** 
# sizesmall:mnO2           -3.668e+00  1.102e+00  -3.328 0.001139 ** 
# speedlow:mnO2            -3.173e+00  1.258e+00  -2.522 0.012888 *  
# speedmedium:mnO2          1.546e+00  8.650e-01   1.787 0.076302 .  
# mnO2:NO3                  2.940e-01  1.507e-01   1.951 0.053193 .  
# mnO2:NH4                 -4.187e-03  1.153e-03  -3.630 0.000406 ***
# mnO2:oPO4                 8.314e-03  6.558e-03   1.268 0.207161    
# sizemedium:Cl             8.630e-02  5.192e-02   1.662 0.098900 .  
# sizesmall:Cl              1.216e-01  8.325e-02   1.461 0.146443    
# speedlow:Cl              -6.484e-02  1.230e-01  -0.527 0.598898    
# speedmedium:Cl            1.078e-01  7.363e-02   1.464 0.145660    
# mxPH:Cl                   6.049e-02  4.472e-02   1.353 0.178473    
# Cl:NH4                    1.926e-04  3.301e-05   5.836 4.04e-08 ***
# Cl:PO4                   -7.623e-04  1.947e-04  -3.914 0.000146 ***
# NO3:oPO4                  1.447e-02  5.789e-03   2.499 0.013700 *  
# seasonspring:NH4         -1.450e-02  5.820e-03  -2.492 0.013977 *  
# seasonsummer:NH4         -7.518e-03  3.240e-03  -2.320 0.021893 *  
# seasonwinter:NH4         -6.715e-03  3.485e-03  -1.927 0.056208 .  
# sizemedium:NH4           -1.370e-02  7.940e-03  -1.726 0.086806 .  
# sizesmall:NH4            -1.759e-03  8.281e-03  -0.212 0.832105    
# speedlow:NH4             -1.693e-02  9.846e-03  -1.719 0.088010 .  
# speedmedium:NH4          -6.666e-03  6.962e-03  -0.957 0.340141    
# NO3:NH4                   3.508e-04  1.633e-04   2.148 0.033556 *  
# NH4:oPO4                 -3.128e-05  1.309e-05  -2.390 0.018289 *  
# seasonspring:oPO4         3.389e-02  7.118e-02   0.476 0.634787    
# seasonsummer:oPO4        -2.213e-01  5.776e-02  -3.831 0.000198 ***
# seasonwinter:oPO4        -2.938e-02  5.932e-02  -0.495 0.621224    
# speedlow:oPO4             2.150e-01  6.444e-02   3.336 0.001108 ** 
# speedmedium:oPO4          8.682e-02  5.417e-02   1.603 0.111452    
# oPO4:PO4                  1.728e-04  9.094e-05   1.900 0.059618 .  
# seasonspring:PO4         -5.793e-02  4.815e-02  -1.203 0.231120    
# seasonsummer:PO4          1.411e-01  3.649e-02   3.867 0.000173 ***
# seasonwinter:PO4         -1.137e-02  4.113e-02  -0.276 0.782660    
# sizemedium:PO4           -7.630e-02  3.297e-02  -2.314 0.022212 *  
# sizesmall:PO4            -1.292e-01  3.646e-02  -3.544 0.000548 ***
# speedlow:PO4             -1.399e-01  4.894e-02  -2.859 0.004948 ** 
# speedmedium:PO4          -8.215e-02  4.186e-02  -1.963 0.051811 .  
# sizemedium:Chla          -1.291e-01  1.005e-01  -1.284 0.201412    
# sizesmall:Chla            1.851e-01  1.517e-01   1.220 0.224564    
# mxPH:Chla                 1.773e-01  7.712e-02   2.299 0.023111 *  
# mnO2:Chla                -4.371e-02  1.961e-02  -2.228 0.027572 *  
# oPO4:Chla                -2.397e-03  8.074e-04  -2.968 0.003566 ** 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

# Residual standard error: 6.724 on 130 degrees of freedom
# Multiple R-squared: 0.7828,  Adjusted R-squared: 0.6674 
# F-statistic: 6.788 on 69 and 130 DF,  p-value: < 2.2e-16 

############################################################################
# Season is not significant by itself, but has some very significant interactions with 
# PO4 and oPO4, in particular the summer factor. Similar situation with the Cl predictor.

# Let's use the full model without interactions to find outliers, and repeat the analysis
# with interactions and compare results.

####################################33

# Studentized deleted residuals
del.res.full = studres(mod.full)
del.res.full.interact = studres(mod.full.interact) # NaNs produced???

sum(is.na(del.res.full.interact)) # 1 Na
which(is.na(del.res.full.interact)) # 48
dat[48,] # a6 = 0
summary(mod.full.interact)$residual[48] # -5e-15, effectively 0
min(abs(summary(mod.full.interact)$residual)) # Looks like this is the smallest absolute residual
length(which(dat$a6==0)) # a6 is zero for 107 points. Why is THIS residual NA? Because it's close to zero, I guess...

del.res.full.interact[48]=0


# Bonferroni test
n = nrow(dat)
p = summary(mod.full.interact)$df[1]
a = 0.01 
crit.b = qt(1 - a/(2 * n), n - p - 1)  # 4.279
max(abs(del.res.full.interact))  # 4.964
which(del.res.full.interact>crit.b) # only observation 73

# X outliers from leverage
H.int = hatvalues(mod.full.interact)
h.bar.int = p/n #.58
which(H.int > 2 * h.bar.int)
plot(sort(H.int))
hist(H.int)
# H takes values from 0 up to 1, pretty unifromly distributed

crit.vals.full.int = list(cooks = 4/n
                      , dffits = 2 * sqrt(p/n)
                      , levg=2*h.bar.int
                      , bonf=crit.b
                      , dfbetas = 2/sqrt(n)
)


# repeat bonferroni test for model without interactions
n = nrow(dat)
p = summary(mod.full)$df[1]
a = 0.01 
crit.b = qt(1 - a/(2 * n), n - p - 1)  # 4.154
max(abs(del.res.full))  # 5.938
which(del.res.full>crit) # only observation 133, 170


# X outliers from leverage
H.full = hatvalues(mod.full)
h.bar.full = p/n #.58
H.full[which(H.full > 2 * h.bar.full) ] # crit = 0.16
# 20        21        35        88        89        127       128       134       153 
# 0.3309145 0.1609968 0.2072736 0.1689092 0.2429690 0.1753647 0.1974960 0.3630732 0.8931062
H.full[133]
plot(sort(H.full))
hist(H.full)
plot(mod.full, which=5)
abline(v=.22, lty=2)
abline(v=.4, lty=2)
abline(v=2*h.bar.full, lty=3) 

crit.vals.full = list(cooks = 4/n
                      , dffits = 2 * sqrt(p/n)
                      , levg=2*h.bar.full
                      , bonf=crit.b
                      , dfbetas = 2/sqrt(n)
                      )

# I think most of these aren't outliers by leverage. 
# 153 is outlier by leverage. Leverages less than .4 seem close enough to main cluster for me. 
# Maybe could go up to h=.22
##
# By cooks distance, 170, 133, and 134. 133 and 170 were outliers by studentized deleted
# residual as well, we should probably get rid of these two, as well as 153.

###########################################################################

# narrowing in on 12 most problematic observations for additional analysis

ix = c(20,21,35,73, 88,89,127,128,133,134,153,170)

df.diag.full = data.frame(index = ix
                          ,cooks = cooks.distance(mod.full)[ix]
                          ,dffits = dffits(mod.full)[ix]
                          ,leverage = H.full[ix]
                          ,abs.stud.del.res = abs(del.res.full[ix]))

df.eval.full = as.data.frame(abs(df.diag.full[,2:5]) > crit.vals.full[-5])
df.eval.full = cbind(df.diag.full$index, df.eval.full)

# Cooks and DFFITS coorborates influence of observations 89, 133, 134, 153, 170
# NB cooks distance heuristic critical threshhold is TINY! Just .02. Ramping the
# threshhold up to .2 narrows in on observations 133, 134, and 170, the points 
# suggested by the residuals vs leverage plot of the full model without interactions.

############################################################################

# Repeat analysis for model with interaction terms.


ix = c(20,21,35,73, 88,89,127,128,133,134,153,170)

df.diag.full.int = data.frame(index = ix
                          ,cooks = cooks.distance(mod.full.interact)[ix]
                          ,dffits = dffits(mod.full.interact)[ix]
                          ,leverage = H.int[ix]
                          ,abs.stud.del.res = abs(del.res.full.interact[ix]))

df.eval.full.int = as.data.frame(abs(df.diag.full.int[,2:5]) > crit.vals.full.int[-5])
df.eval.full.int = cbind(df.diag.full.int$index, df.eval.full.int)

# In the model with interaction terms, only observation 73 was highlighted as a problematic
# point based on the deleted residuals, but it looks like we actually have many influential 
# observations by cook's distance and DFFITS, including observation 73.
# This analysis corroborates the influence of 89, 133, 134, 153, 170 observed above.
#
# Oddities:
# Very high cook's distance observed in observations 153, 88, and 133. 
# Very high |DFFITS| in observations 153, 88, 20, 133
#
# In addition to the dropping terms 133, 134, and 170, I also recommend we drop observations
# 88 and 153 to mitigate influence in models including interaction terms.

############################################################################

# Repeat analysis for full model with select interaction terms. Cause why not.

p.sel = summary(mod.full.selectInteract)$df[1]
H.sel = hatvalues(mod.full.selectInteract)
h.bar.sel = p.sel/n
sum(H.sel > 2*h.bar.sel) # 21 observations... wowza. WTF?
plot(sort(H.sel))
abline(h=h.bar.sel, lty=2) # this is way too liberal of a cut off. 
abline(h=.62, lty=2) # 0.62 is a way better cutoff
sum(H.sel > .62) # 28 observations by this criterion. 
sort(H.sel, decreasing=T)[1:28]
# The usual suspects appear at the top of the list:
# 153, 20, 133, 88, (35), 89, (34), 134, (...), 170.
# Going down to obs 170 gives 15 observations to investigate by leverage.

# Crazy how much the leverages change with each model
# Clearly leverage is very, very model specific, especially in elevated dimensions.
plot(H.full)
points(H.int, col='blue')
points(H.sel, col='red')

del.res.sel = studres(mod.full.selectInteract)

# Bonferroni test
crit.b.sel = qt(1 - a/(2 * n), n - p.sel - 1)  # 4.197
max(abs(del.res.sel))  # 4.08016 , no outliers by bonferroni at a=0.01
plot(sort(del.res.sel))


crit.vals.sel = list(cooks = 4/n
                       , dffits = 2 * sqrt(p.sel/n)
                       , levg=2*h.bar.sel
                       , bonf=crit.b.sel
                       , dfbetas = 2/sqrt(n))

ix = c(20,21,35,73, 88,89,127,128,133,134,153,170, 35, 34) # added 35 and 34

df.diag.sel = data.frame(index = ix
                              ,cooks = cooks.distance(mod.full.selectInteract)[ix]
                              ,dffits = dffits(mod.full.selectInteract)[ix]
                              ,leverage = H.sel[ix]
                              ,abs.stud.del.res = abs(del.res.sel[ix]))

df.eval.sel = as.data.frame(abs(df.diag.sel[,2:5]) > crit.vals.sel[-5])
df.eval.sel = cbind(df.diag.sel$index, df.eval.sel)

# observations 35 and 34 have significant influence in this model by cooks, DFFITS and leverage,
# but not in any other models. Influence of observations 153, 170, 133, 88, 73, 35, and 20
# corroborated. Observations 21, 89, 127, and 134 had high leverage, but not significant
# influence, and observation 128 had neither.
#
# Observation 153 had highest leverage and extremely high DFFITS and cooks distance in 
# this model. Observation 133 had next highest cooks and DFFITS, and the third highest leverage.
#
# 35 and 34 don't seem like serious concerns, but should probably investigate their influence
# on the earlier models. At the very least, will keep them in mind for the next (and last) model

############################################################################3

# Use step regression to select non-categorical variables
formula.cat.only = 'a6 ~ season + size + speed '

mod.reduced = step(mod.full
                  ,scope=list(end=formula.cat.only, start=formula.full)
                  ,direction='backward')

mod.reduced$call # this isn't respecting keeping "speed" in the model

# Gives same model. Whatever, I'll just work with this one. Doesn't really matter.
mod.reduced = step(mod.full
                   #,scope=list(end=formula.cat.only, start=formula.full)
                   ,direction='backward')


# Studentized deleted residuals
del.res.reduced = studres(mod.reduced)

# Bonferroni test
n = nrow(dat)
p = summary(mod.reduced)$df[1]
a = 0.01 
crit.b.reduced = qt(1 - a/(2 * n), n - p - 1)  # 4.1527
max(abs(del.res.reduced))  # 4.964
which(del.res.reduced>crit.b.reduced) # only observations 133 and 170

# X outliers from leverage
H.red = hatvalues(mod.reduced)
h.bar.red = p/n #.58
which(H.red > 2 * h.bar.red)
plot(sort(H.red))
hist(H.red)
# usual suspects. 98, 163 and 175 are new
# 20  21  35  88  89  98 127 128 133 134 153 163 175 
ix = c(20,21,35,73, 88,89,127,128,133,134,153,170, 34, 98, 163, 175)


crit.vals.red = list(cooks = 4/n
                          , dffits = 2 * sqrt(p/n)
                          , levg=2*h.bar.red
                          , bonf=crit.b.reduced
                          , dfbetas = 2/sqrt(n)
)



df.diag.red = data.frame(index = ix
                         ,cooks = cooks.distance(mod.reduced)[ix]
                         ,dffits = dffits(mod.reduced)[ix]
                         ,leverage = H.red[ix]
                         ,abs.stud.del.res = abs(del.res.reduced[ix]))

df.eval.red = as.data.frame(abs(df.diag.red[,2:5]) > crit.vals.red[-5])
df.eval.red = cbind(df.diag.red$index, df.eval.red)

# Our new points are 98, 163, and 175. They all meet the leverage criteria
# but their cooks and DFFITS are below threshhold. By Cooks and DFFITS and
# leverage, the most influential points in the model are 89, 133, 134, and
# 153. 170 did not have high leverage, but it had high cooks and DFFTIS.


# Let's bring it all together and rebuild each analysis with all
# the considered observations from the earlier analyses:

df.diag.all = data.frame(row.names = ix #index = ix
                     ## REDUCED (BACKWARDS ELIMINATION VIA AIC) ## 
                     ,cooks.red = cooks.distance(mod.reduced)[ix]
                     ,dffits.red = dffits(mod.reduced)[ix]
                     ,leverage.red = H.red[ix]
                     ,studres.red = del.res.reduced[ix]
                     ## FULL ##
                    ,cooks.ful = cooks.distance(mod.full)[ix]
                    ,dffits.ful = dffits(mod.full)[ix]
                    ,leverage.ful = H.full[ix]
                    ,studres.ful = del.res.full[ix]
                    ## SELECTED INTERACTIONS
                    ,cooks.sel = cooks.distance(mod.full.selectInteract)[ix]
                    ,dffits.sel = dffits(mod.full.selectInteract)[ix]
                    ,leverage.sel = H.sel[ix]
                    ,studres.sel = del.res.sel[ix]
                    ## ALL INTERACTIONS ##
                    ,cooks.int = cooks.distance(mod.full.interact)[ix]
                    ,dffits.int = dffits(mod.full.interact)[ix]
                    ,leverage.int = H.int[ix]
                    ,studres.int = del.res.full.interact[ix]
                    )

dim(df.diag.all) # 16x16, that's a funny coincidence...


write.csv(df.diag.all, 'math651_final_outliers.csv')


# Conclusion: observations 







##################################################################

# testing addition of the interaction term where one of the variables is strong and the other
# is weak.

formula.seasonPO4 = paste(formula.full, " + season:PO4")
  
mod.seasonPO4 = lm(formula.seasonPO4, dat)

anova(mod.full, mod.seasonPO4) # p value = .3736

plot(mod.full$residuals, mod.seasonPO4$residuals)

# Also, try the "canary in the coal min:" Try an interaction with a completely random variable

