## ----setup, include=FALSE------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


## ---- fig.width=4, fig.height=4, message=FALSE---------------------------------------------------------------------
library(Distance)
data("PTExercise")
head(PTExercise, n=3)
conversion.factor <- convert_units("meter", NULL, "hectare")
# Fit half-normal detection function, no truncation
PTExercise.hn <- ds(data=PTExercise, transect="point", key="hn", convert.units=conversion.factor)
plot(PTExercise.hn, pdf=TRUE, main="Simulated pt transect data\nHalf normal key function")


## ---- trunc20, message=FALSE---------------------------------------------------------------------------------------
# Half normal, no adjustments
PTExercise.hn.t20m <- ds(data=PTExercise, transect="point", key="hn", truncation=20,
                    convert.units=conversion.factor)
# Hazard rate, no adjustments
PTExercise.hr.t20m <- ds(data=PTExercise, transect="point", key="hr", truncation=20,
                    convert.units=conversion.factor)
# Uniform, cosine adjustments
PTExercise.uf.cos.t20m <- ds(data=PTExercise, transect="point", key="unif", 
                        adjustment="cos", truncation=20,convert.units=conversion.factor)


## ---- echo=TRUE, eval=TRUE-----------------------------------------------------------------------------------------
gof_ds(PTExercise.hn.t20m)


## ---- echo=FALSE---------------------------------------------------------------------------------------------------
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


## ---- fig.height=6-------------------------------------------------------------------------------------------------
par(mfrow=c(2,2))
plot(PTExercise.hn, main="Half normal, no truncation", pdf=TRUE)
plot(PTExercise.hn.t20m, main="Half normal, truncation 20m", pdf=TRUE)
plot(PTExercise.hr.t20m, main="Hazard rate, truncation 20m", pdf=TRUE)
plot(PTExercise.uf.cos.t20m, main="Uniform with cosine, truncation 20m", pdf=TRUE)


## ---- echo=TRUE, eval=TRUE-----------------------------------------------------------------------------------------
data("wren_5min")
data("wren_snapshot")
conversion.factor <- convert_units("meter", NULL, "hectare")
wren5min.uf.cos.t110 <- ds(data=wren_5min, key="unif", adjustment="cos", 
                        transect="point", truncation=110, 
                        convert.units=conversion.factor)
wrensnap.hr.cos.t110 <- ds(data=wren_snapshot, key="hr", adjustment=NULL, 
                        transect="point", truncation=110, 
                        convert.units=conversion.factor)


## ---- echo=FALSE---------------------------------------------------------------------------------------------------
# Harvest results
n <- 2
wren.tab <- data.frame(Method=c("Five minute","Snapshot"), Density=rep(NA,n), 
                       Lower.CI=rep(NA,n), Upper.CI=rep(NA,n))

get.results.f <- function(fit.model) { return(c(D=fit.model$dht$individuals$D$Estimate,
         lCL=fit.model$dht$individuals$D$lcl,
         uCL=fit.model$dht$individuals$D$ucl))
}
wren.tab[1,2:4] <- get.results.f(wren5min.uf.cos.t110)
wren.tab[2,2:4] <- get.results.f(wrensnap.hr.cos.t110)
knitr::kable(wren.tab, caption="Winter wren density estimates from 5 minute counts and snapshot moment.", digits=3)


## ---- fig.height=5, fig.cap="Wren avoidance of observer is evident from surplus of detections 40-50m from point stations."----
# Plot detection functions
par(mfrow=c(1,2))
plot(wren5min.uf.cos.t110, main="5 minute count")
plot(wrensnap.hr.cos.t110, main="Snapshot moment")


## ---- fig.height=5, fig.cap="Even with evasive movement, detection function models for both the 5-minute and snapshot data pass the goodness of fit test."----
gof_ds(wren5min.uf.cos.t110)
gof_ds(wrensnap.hr.cos.t110)

