## ----setup, include=FALSE-------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


## ---- fig.width=4, fig.height=4, message=FALSE----------------
library(Distance)
data("PTExercise")
head(PTExercise, n=3)
conversion.factor <- convert_units("meter", NULL, "hectare")
bin.cutpoints <- seq(0, 35, by = 2.5)
# Fit half-normal detection function, no truncation, 2.5m cutpoints
PTExercise.hn <- ds(data=PTExercise, transect="point", key="hn", convert.units=conversion.factor,
                    cutpoints=bin.cutpoints)
plot(PTExercise.hn, pdf=TRUE, main="Simulated pt transect data\nHalf normal key function")


## ---- trunc20, message=FALSE----------------------------------
bin.cutpoints.20m <- seq(0, 20, by = 2.5)
# Half normal, no adjustments
PTExercise.hn.t20m <- ds(data=PTExercise, transect="point", key="hn",
                    convert.units=conversion.factor, cutpoints=bin.cutpoints.20m)
# Hazard rate, no adjustments
PTExercise.hr.t20m <- ds(data=PTExercise, transect="point", key="hr", truncation=20,
                    convert.units=conversion.factor, cutpoints=bin.cutpoints.20m)
# Uniform, cosine adjustments
PTExercise.uf.cos.t20m <- ds(data=PTExercise, transect="point", key="unif", adjustment="cos",
                             convert.units=conversion.factor, cutpoints=bin.cutpoints.20m)


## ---- echo=TRUE, eval=TRUE------------------------------------
gof_ds(PTExercise.hn.t20m)


## ---- echo=FALSE----------------------------------------------
#   Do not get excited about the code in this chunk; 
#   it is not necessary for your understanding of distance sampling.
pt.tab <- data.frame(DetectionFunction=c("Half-normal","Half-normal",
                                         "Hazard rate","Uniform"), 
                     Adjustments=c("None","None","None","Cosine"), Truncation=c(34.2,20,20,20), 
                     AIC=rep(NA,4), Density=rep(NA,4), D.CV=rep(NA,4), Lower.CI=rep(NA,4), Upper.CI=rep(NA,4))

get.results.f <- function(fit.model) {
  return(c(AIC=summary(fit.model$ddf)$aic,
         D=fit.model$dht$individuals$D$Estimate,
         D.CV=fit.model$dht$individuals$D$cv,
         lCL=fit.model$dht$individuals$D$lcl,
         uCL=fit.model$dht$individuals$D$ucl))
}
pt.tab[1,4:8] <- get.results.f(PTExercise.hn)
pt.tab[2,4:8] <- get.results.f(PTExercise.hn.t20m)
pt.tab[3,4:8] <- get.results.f(PTExercise.hr.t20m)
pt.tab[4,4:8] <- get.results.f(PTExercise.uf.cos.t20m)
knitr::kable(pt.tab, caption="Results from simulated point transect data.", digits=3)


## ---- fig.height=6--------------------------------------------
par(mfrow=c(2,2))
plot(PTExercise.hn, main="Half normal, no truncation", pdf=TRUE)
plot(PTExercise.hn.t20m, main="Half normal, truncation 20m", pdf=TRUE)
plot(PTExercise.hr.t20m, main="Hazard rate, truncation 20m", pdf=TRUE)
plot(PTExercise.uf.cos.t20m, main="Uniform with cosine, truncation 20m", pdf=TRUE)


## ---- echo=TRUE, eval=TRUE------------------------------------
data("wren_5min")
conversion.factor <- convert_units("meter", NULL, "hectare")

#Make an initial fit to look at truncation - half normal with no adjustments
bin.cutpoints <- c(0, 10, 20, 30, 40, 60, 80, 100, 120)
wren5min.hn <- ds(data=wren_5min, key="hn", adjustment = NULL,
                  transect="point", cutpoints=bin.cutpoints,
                  convert.units=conversion.factor)
plot(wren5min.hn)
plot(wren5min.hn, pdf = TRUE)


## ---- echo=TRUE, eval=TRUE------------------------------------
bin.cutpoints.100m <- bin.cutpoints <- c(0, 10, 20, 30, 40, 60, 80, 100)
#Note - code below issues a warning about non-monotonicity -- will want to
# check it out if this becomes the selected model
wren5min.uf.cos.t100 <- ds(data=wren_5min, key="unif", adjustment="cos", 
                        transect="point", cutpoints=bin.cutpoints.100m,
                        convert.units=conversion.factor)
#Note - the code below throws an error when trying the second adjustment term
# and so just returns the results for key plus one adjustment
wren5min.hn.herm.t100 <- ds(data=wren_5min, key="hn", adjustment="herm", 
                        transect="point", cutpoints=bin.cutpoints.100m,
                        convert.units=conversion.factor)
wren5min.hr.poly.t100 <- ds(data=wren_5min, key="hr", adjustment="poly", 
                        transect="point", cutpoints=bin.cutpoints.100m,
                        convert.units=conversion.factor)
AIC(wren5min.uf.cos.t100, wren5min.hn.herm.t100, wren5min.hr.poly.t100)


## ---- echo=TRUE, eval=TRUE------------------------------------
data("wren_snapshot")
wrensnap.uf.cos.t100 <- ds(data=wren_snapshot, key="unif", adjustment="cos", 
                        transect="point", cutpoints=bin.cutpoints.100m,
                        convert.units=conversion.factor)
#Note - the half-normal with Hermite adjustments often throws an error
# and will return the key function without adjustment terms, check # params
wrensnap.hn.herm.t100 <- ds(data=wren_snapshot, key="hn", adjustment="herm", 
                        transect="point", cutpoints=bin.cutpoints.100m,
                        convert.units=conversion.factor)
wrensnap.hr.poly.t100 <- ds(data=wren_snapshot, key="hr", adjustment="poly", 
                        transect="point", cutpoints=bin.cutpoints.100m,
                        convert.units=conversion.factor)
AIC(wrensnap.uf.cos.t100, wrensnap.hn.herm.t100, wrensnap.hr.poly.t100)



## ---- fig.height=5, fig.cap="Wren avoidance of observer is evident from surplus of detections 30-60m from point stations."----
# Plot detection functions
par(mfrow=c(1,2))
plot(wren5min.hr.poly.t100, main="5 minute count", pdf = TRUE)
plot(wrensnap.hr.poly.t100, main="Snapshot moment", pdf = TRUE)


## ---- fig.height=5, fig.cap="Corresponding detection function plots."----
# Plot detection functions
par(mfrow=c(1,2))
plot(wren5min.hr.poly.t100, main="5 minute count")
plot(wrensnap.hr.poly.t100, main="Snapshot moment")


## ---- fig.height=5, fig.cap="Even with evasive movement, detection function models for both the 5-minute and snapshot data pass the goodness of fit test."----
gof_ds(wren5min.hr.poly.t100)
gof_ds(wrensnap.hr.poly.t100)


## ---- echo=FALSE----------------------------------------------
# Harvest results
n <- 2
wren.tab <- data.frame(Method=c("Five minute","Snapshot"), Density=rep(NA,n), 
                       Lower.CI=rep(NA,n), Upper.CI=rep(NA,n))

get.results.f <- function(fit.model) { return(c(D=fit.model$dht$individuals$D$Estimate,
         lCL=fit.model$dht$individuals$D$lcl,
         uCL=fit.model$dht$individuals$D$ucl))
}
wren.tab[1,2:4] <- get.results.f(wren5min.hr.poly.t100)
wren.tab[2,2:4] <- get.results.f(wrensnap.hr.poly.t100)
knitr::kable(wren.tab, caption="Winter wren density estimates from 5 minute counts and snapshot moment.", digits=3)

