# Interaction Mining

data(algae)
dat = knnImputation(algae, k=14)
dat = dat[-c(170, 133, 134),] # remove outliers

# Begin by writing a function that takes in a model specification
# and returns a dataframe of p-values and f statistics for added
# variable tests of each second order interaction, tested individually against the
# naive full model. 

# I think the cox paper isn't just suggesting we add single interactions, but permute
# up to all second order interactions being included. I'm not doing that. 55 is enough.

mine.interactions = function(response, predictors=NA, data=NA, intercept = F){
  # response: a single character string defining the response variable
  # predictors: a vector of strings defining the predictors
  # data: a dataframe. Any variables that need to be treated as factors should
  #       be converted into factors in this dataframe already.
  
  # build the formula for the first order model
  # For now, let's assume response is a single character vector
  #if(length(sapply(strsplit(response,','),c))!=1){
  if(length(response)!=1){
    stop("Only single response variable currently supported")
  }
  
  left = paste(response, '~')
  right = paste(predictors, collapse = '+')
  full.formula = paste(left, right)
  if(!intercept){
    full.formula = paste(full.formula,'-1')
  }
  model.full = lm(full.formula, data)
  print(summary(model.full))
  
  # build character strings for each interaction
  n = length(predictors)
  interactions = c() # should preallocate 
  k=0
  for(i in 1:n){
    for(j in c(1:n)){      
      p1 = predictors[i]
      p2 = predictors[j]
      if(p1<p2){ # alphabetic deduplification
        k=k+1
        int.k = paste(p1, p2, sep=':')
        interactions = c(interactions, int.k)
      }
    }
  }

  # generate the models and store the relevant statistics.
  out = matrix(NA, k, 4)
  out = as.data.frame(out)
  names(out) = c("term","p-value","F", "Adjusted.R2")
  for(i in 1:k){
    int.formula = paste(full.formula, interactions[i], sep=' + ')
    model.int = lm(int.formula, data)
    A = anova(model.int, model.full)
    S = summary(model.int)
    out[i,] = c(interactions[i], A$"Pr(>F)"[2], A$"F"[2], S$adj.r.squared)
  }
  out[,2]=as.numeric(out[,2])
  out[,3]=as.numeric(out[,3])
  out
}

test.full = mine.interactions('a6', predictors=names(dat)[1:11], dat, intercept=T)
# R-squared:  0.4737;  Adjusted R-squared:  0.4301 
# In the "Outlier mitigation" section, they report "wihtout outliers" R2, should be
# adjusted R2. 
# The R2 values reported in this section are not consistent with the R2 values reported
# for the full model in the model building section. I get the same R2 they get in 
# the outlier section, and ajd-r2 as above.
# also, dropping PO4 and oPO4 from AIC model should be accompanied by F-test.

hist(test.full[,2], xlab="p-value", main="Histogram of p-values from F-tests adding\nindividual interaction terms to the full model") # under null hypothesis this should be uniform. Looks bimodal. Small sample size though.
plot(test.full[,2])
abline(h=.05)
abline(h=.01, lty=2)
nrow(test.full)
# 6 interactions significant at a=0.05
test.full[test.full[,2]<.05,]
#######################################################
term     p-value        F       Adjusted.R2
15  mnO2:NO3 0.008049715 7.181407 0.448878815805768
23   Cl:mnO2 0.036503754 4.439418 0.440685509993636
30 NO3:speed 0.015826692 4.243590 0.449777638852647
32   NO3:PO4 0.030835530 4.736276  0.44158428891519
42  oPO4:PO4 0.026928862 4.976503 0.442309498335494
#######################################################
# Although addition of terms is significant, it appears
# that the improvement to the models is marginal.
# All of these terms are present in the reduced model except speed.
mod.full = lm(a6~., dat[,c(1:11,17)])
mod.full.explicit = lm(a6~season+size+speed+mxPH+mnO2+Cl+NO3+NH4+oPO4+PO4+Chla, dat)
summary(mod.full.explicit)$adj.r.squared # .4300572
mod.explicit.aic = step(mod.full.explicit)
mod.explicit.aic$call
summary(mod.explicit.aic)$adj.r.squared
# 0.4315243
mod.explicit.aic.int = step(lm(a6~season+size+speed+mxPH+mnO2+Cl+NO3+NH4+oPO4+PO4+Chla+mnO2:NO3, dat))
summary(mod.explicit.aic.int)$adj.r.squared # 0.4531913

mod.aic = step(mod.full) 
mod.aic$call # a6 ~ season + size + speed + mnO2 + Cl + NO3 + oPO4 + PO4
# Chla and speed dropped when regressing through the origin.
summary(mod.full)$adj.r.squared # .43
summary(mod.aic)$adj.r.squared # 0.4315
pred.aic2 = c("season","size","speed", "mnO2","Cl","NO3","NH4","oPO4","PO4",'Chla')
test.aic2 = mine.interactions('a6', predictors=pred.aic2, dat, intercept=T)
test.aic2[test.aic2[,2]<.05,]
###################################################
term     p-value        F       Adjusted.R2
7   mnO2:NO3 0.009035632 6.965241  0.45037587336297
14   Cl:mnO2 0.048709975 3.938282 0.441379956080115
21 NO3:speed 0.019359959 4.032267   0.4506660573857
23   NO3:PO4 0.031165687 4.717049 0.443722415708355
33  oPO4:PO4 0.035107073 4.507426 0.443093822287251
####################################################
# These are all relative to 0.4315. No huge gains by adding single interactions.


int.plot = function(t1, t2, resp, lab1=NA, lab2=NA, bins=3, leg=T){
  fact = cut(t2, bins)  
  lev = levels(fact)
  plot(t1, resp, type = 'n', xlab=lab1, ylab='a6', main=paste("Interaction plot for",bins,"binned levels of",lab2))
  for(i in 1:bins){
    ix = fact==lev[i]
    r=resp[ix]
    t=t1[ix]
    mod.temp = lm(r~t)
    points(t, r, pch=i, col=i)
    abline(mod.temp, col=i)
  }
  if(leg){
    legend("bottomright", lev, pch=1:bins, col=1:bins, lty=c(1,1,1))
  }
}


with(dat, int.plot(NO3, PO4, a6, lab1="Cl", lab2="mn02", leg=F))
with(dat, int.plot(Cl, mnO2, a6, lab1="Cl", lab2="mn02", leg=F))

mod.full.smartInt = lm(a6~.+mnO2:NO3+NO3:speed, dat[,c(1:11,17)])

mod.full.smartInt = lm(a6~.+mnO2:NO3+NO3:speed, dat[,c(1:11,17)])
summary(mod.full.smartInt) # .4487, nothing special
mod.aic.smartInt = step(mod.full.smartInt)
mod.aic.smartInt$call
# a6 ~ season + size + speed + mnO2 + Cl + NO3 + oPO4 + PO4 + Chla + mnO2:NO3
# Valid model, thank god.
# In addition to the interaction term, this differs from the first order
# step() result by the addition of the speed predictor
summary(mod.aic.smartInt) # .4532

mod.full.smartInt2 = lm(a6~.+mnO2:NO3 + Cl:mnO2 + NO3:speed + NO3:PO4 + oPO4:PO4, dat[,c(1:11,17)])
mod.aic.smartInt2 = step(mod.full.smartInt2)
mod.aic.smartInt2$call
# a6 ~ season + size + speed + mnO2 + Cl + NO3 + oPO4 + PO4 + Chla + mnO2:Cl + speed:NO3 + oPO4:PO4
# Valid model, sweet.
summary(mod.aic.smartInt2) # 0.4621 nice

install.packages('DAAG') # for press()
library('DAAG')

press(mod.full)           # 20234.9
press(mod.aic)            # 13957.7
press(mod.full.smartInt)  # 13560.46
press(mod.aic.smartInt2)  # 13313.17 # adjr2 goes up due to overfitting.
press(mod.aic.smartInt)   # 12894.52
mod.aic.smartInt$call

# Let's visualize the interactions in the aic.smartInt model
with(dat, int.plot(NO3, mnO2, a6, lab1="NO3", lab2="mn02", bins=10, leg=T))
with(dat, int.plot(NO3, mnO2, a6, lab1="NO3", lab2="mn02", bins=3))
# Looks like there's essentially two levels of the interaction
#with(dat, int.plot(NO3, mnO2, a6, lab1="NO3", lab2="mn02", bins=2))

score.model = function(model){
  data.frame(
    r2 = summary(model)$r.squared,
    adjr2 = summary(model)$adj.r.squared,
    AIC = AIC(model),
    #MSE = mean(model$residuals^2), # This is a biased estimator
    MSE = sum(model$residuals^2)/196, # Unbiased estimator divides by n-1
    PRESS = press(model)
    )
}

score.models=function(models){
  n=length(models)
  df = matrix(NA, n, 6)
  df = data.frame(df)
  df[,1] = as.character(df[,1])
  df[,2] = as.numeric(df[,2])
  df[,3] = as.numeric(df[,3])
  df[,4] = as.numeric(df[,4])
  df[,5] = as.numeric(df[,5])
  df[,5] = as.numeric(df[,6])
  names(df) = c('model','r2', 'adjr2','AIC', 'MSE', 'PRESS')
  for(i in 1:n){
    model = models[[i]]
    print("we're here")
    df[i,1] = as.character(model$call)[[2]]
    print("we're there")
    df[i,-1] = score.model(model)
  }
  df
}


full.formula = "a6~season+size+speed+mxPH+mnO2+Cl+NO3+NH4+oPO4+PO4+Chla"
interaction.candidates = c("mnO2:NO3","NO3:speed","oPO4:PO4","NO3:PO4","Cl:mnO2")
test.formulas = c(full.formula)
for(i in 1:length(interaction.candidates)){
  inter = paste(interaction.candidates[1:i], collapse=" + ")
  cand.formula = paste(full.formula, inter, sep= " + ")
  test.formulas = c(test.formulas, cand.formula)
}

aic.models = as.list(rep(NA,6))
for(i in 1:6){
  aic.models[[i]] = step(lm(test.formulas[i], dat))
}
score.models(aic.models)
score.model(mod.full.explicit)

int1 = aic.models[[2]]

plot(cooks.distance(int1), main="Cook's Distances under Int1 model.", ylab="Cook's Distance")
ix=c(89,100)
text(ix+8, cooks.distance(int1)[ix], ix, cex=.6)
abline(h=4/197, lty=2) # 0.02030457

plot(int1, which=4)
which(cooks.distance(int1)>.1)
5/200
dat2 = dat[-c(89,100),]
int1.formula = as.character(int1$call)[[2]]
int1.out = lm(int1.formula, dat2)
summary(int1.out)$adj.r.squared # .4665495
press(int1.out) # 11383
