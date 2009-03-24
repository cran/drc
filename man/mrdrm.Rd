\name{mrdrm}

\alias{mrdrm}

\title{Model-robust dose-response modelling}

\description{
  Model-robust dose-response modelling is based on an optimal linear convex combination of two model fits, 
  a parametric and a non-parametric model fit. The current implementation relies on local linear regression
  (loess) for the non-parametric part.
}

\usage{
  mrdrm(object1, object2, lambda = (0:10)/10, criterion = c("gcv", "lcv"), critFct = c("ls", "ll"), 
  ls.weights = c("nonpar", "ad hoc", "none", "par", "response"), fixedEnd = FALSE, unitScale = FALSE)
}

\arguments{
  \item{object1}{object of class 'drc' (the parametric fit).}
  \item{object2}{object of class 'loess' (the non-parametric fit).}
  \item{lambda}{numeric vector of potential mixing values (between 0 and 1).}
  \item{criterion}{character string specifying the criterion to use in the PRESS* procedure.}
  \item{critFct}{character string specifying the criterion function to use in the PRESS* procedure.}  
  \item{ls.weights}{character string specifying the type weights to use in the PRESS* criterion.}
  \item{fixedEnd}{logical indicating whether or not the leave-one-out predictions of the non-parametric fit should
  be equal to the average response at the boundary dose values. If not, no such predictions are obtained at all.}
  \item{unitScale}{logical indicating if the dose values should be transformed to the unit interval in order to
  improve the local regression fit.}
}

\details{
  The PRESS* leave-one-out criterion is used to determine the optimal mixing of the parametric and non-parametric 
  model fits (Nottingham and Birch, 2000)).
}

\value{
  A list of components from the fit.
}

\references{  
  Notttingham, Q. J. and Birch, J. B. (2000) A semiparametric approach to analysing dose-response data,
  \emph{Statist. Med.}, \bold{19}, 389--404.
}

\author{Christian Ritz}

\note{
  Currently, there is no check ensuring that the parametric and non-parametric fits are based on the same data!
}

%\seealso{}

\examples{

## deguelin data from Nottingham and Birch (2000)
deguelin.m1 <- drm(r/n~dose, weights=n, data=deguelin, fct=LL.2(), type="binomial")
deguelin.m2 <- loess(r/n~dose, data=deguelin, degree=1)  # local linear regression

deguelin.mr <- mrdrm(deguelin.m1, deguelin.m2)
deguelin.mr

predict(deguelin.mr, interval = "confidence")
ED(deguelin.mr, c(10, 20, 50, 80, 90), interval = "approximate")
ED(deguelin.m1, c(10, 20, 50, 80, 90), ci = "delta")

plot(deguelin.m1, ylim=c(0,1))
plot(deguelin.mr, add = TRUE, lty = 2)

## The same results  (loess fit automatically supplied)
deguelin.mr2 <- mrdrm(deguelin.m1)
ED(deguelin.mr2, c(10, 20, 50, 80, 90), interval = "approximate")

## With fixed lambda
deguelin.mr3 <- mrdrm(deguelin.m1, deguelin.m2, lambda = 0.8)
plot(deguelin.mr3, add = TRUE, lty = 3)

## Purely non-parametric fit
deguelin.mr4 <- mrdrm(deguelin.m1, deguelin.m2, lambda = 1)
plot(deguelin.mr4, add = TRUE, lty = 4)

## On log scale (completely different results)
deguelin.m2b <- loess(r/n ~ log(dose), data = deguelin, degree = 1)
deguelin.mr2b1 <- mrdrm(deguelin.m1, deguelin.m2b, critFct = "ll")
deguelin.mr2b1

deguelin.mr2b2 <- mrdrm(deguelin.m1, deguelin.m2b, critFct = "ls")
deguelin.mr2b2

deguelin.mr2b3 <- mrdrm(deguelin.m1, deguelin.m2b, critFct = "ls", fixedEnd = TRUE)
deguelin.mr2b3


## daphnids dataset at 24 hours
daphnids1.m1<-drm(no/total~dose, weights = total, data = daphnids[1:8,], fct = LL.2(), type = "binomial")
daphnids1.m2<-loess(no/total~dose, data = daphnids[1:8,], degree = 1)

daphnids1.mr<-mrdrm(daphnids1.m1, daphnids1.m2)
daphnids1.mr

plot(daphnids1.m1)
plot(daphnids1.mr, add=TRUE, type="none", lty=2)

## daphnids dataset at 48 hours
daphnids2.m1<-drm(no/total~dose, weights = total, data = daphnids[9:16,], fct = LL.2(), type = "binomial")
daphnids2.m2<-loess(no/total~dose, data = daphnids[9:16,], degree = 1)

daphnids2.mr<-mrdrm(daphnids2.m1, daphnids2.m2)
daphnids2.mr

plot(daphnids2.m1)
plot(daphnids2.mr, add=TRUE, type="none", lty=2)


## fly dataset from Nottingham & Birch (1996)
fly<-data.frame(
conc = c(0.1,0.15,0.2,0.3,0.5,0.7,0.95),
total = c(47,53,55,52,46,54,52),
killed = c(8,14,24,32,38,50,50))

fly.m1 <- drm(killed/total~conc, weights = total, data = fly, fct = LL.2(), type = "binomial")
fly.m2 <- loess(killed/total~conc, data = fly, degree = 1)

fly.mr1 <- mrdrm(fly.m1, fly.m2)
fly.mr2 <- mrdrm(fly.m1, fly.m2, criterion="lcv")

plot(fly.m1, ylim = c(0,1))
plot(fly.mr1, add = TRUE, type = "none", lty = 3)
plot(fly.mr2, add = TRUE, type = "none", lty = 2)

fly.mr1
fly.mr2
AIC(fly.m1)

## part of bin.mat (currently in 'mrdrc')
if (FALSE)
{
bin.mat.m1<-drm(matured/total~conc, data=bin.mat[c(3,6,9,12,15),],fct=LL.2(), type = "binomial")
bin.mat.m2<-loess(matured/total~conc, data=bin.mat[c(3,6,9,12,15),], degree=1)
bin.mat.mr<-mrdrm(bin.mat.m1, bin.mat.m2)
bin.mat.mr
plot(bin.mat.mr)  # oversmoothed using GCV!

bin.mat.mr2<-mrdrm(bm1, bm.loess, criterion = "lcv")
bin.mat.mr2
plot(bin.mat.mr2)

## exp.a data (currrently in 'mrdrc')
exp.a.m1 <- drm(y ~ x, data = exp.a, fct = LL.3())
exp.a.m2 <- loess(y ~ x, data = exp.a, degree = 1, span = 0.35)
exp.a.mr <- mrdrm(exp.a.m1, exp.a.m2)
exp.a.mr

plot(exp.a.mr, type = "all", conLevel=1, broken=TRUE, xlab="Dose", ylab="Response")
predict(exp.a.mr, se.fit = TRUE)
ED(exp.a.mr, c(10, 50, 90), interval = "approximate")

## Using a given lower reference limit
ED(exp.a.mr, c(10, 50, 90), interval = "approximate", lower = 0)
ED(exp.a.mr, c(10, 50, 90), interval = "approximate", method="grid", lower = 0)


## dataset "x"
exp.x.m1 <- drm(response ~ conc, data = exp.az[exp.az$exp=="x",], fct = LL.3())
exp.x.m2 <- loess(response ~ conc, data = exp.az[exp.az$exp=="x",], degree = 1, span = 0.5)
exp.x.mr <- mrdrm(exp.x.m1, exp.x.m2)
exp.x.mr

plot(exp.x.mr, type="all")
ED(exp.x.mr, c(10, 50, 90), interval = "approximate")
ED(exp.x.mr, c(10, 50, 90), interval = "approximate", cgridsize=100)



#uniConc <- sort(unique(exp.az[exp.az$exp=="x",]$conc))
#doseLoess <- loess((0:8)/9 ~ uniConc)
#unitConc <- predict(doseLoess, data.frame(uniConc = exp.az[exp.az$exp=="x",]$conc))
#
#exp.x.m3 <- loess(exp.az[exp.az$exp=="x",]$response~unitConc, degree=1, span=0.35)
#exp.x.mr2 <- mrdrm(exp.x.m1, exp.x.m3)
#exp.x.mr2
#
#plot(exp.x.m1, type="all")
#lines(seq(0, 2, length=101), 
#predict(exp.x.m3, data.frame(unitConc = predict(doseLoess, data.frame(uniConc = seq(0, 2, length=101))))), lty = 2)



## dataset "w"
exp.w.m1 <- drm(response ~ conc, data = exp.az[exp.az$exp=="w",], fct = LL.3())
exp.w.m2 <- loess(response ~ conc, data = exp.az[exp.az$exp=="w",], degree = 1)
exp.w.mr <- mrdrm(exp.w.m1, exp.w.m2)
exp.w.mr

plot(exp.w.mr, type="all")

## dataset "u"
exp.u.m1 <- drm(response ~ conc, data = exp.az[exp.az$exp=="u",], fct = LL.3())
exp.u.m2 <- loess(response ~ conc, data = exp.az[exp.az$exp=="u",], degree = 1)
exp.u.mr <- mrdrm(exp.u.m1, exp.u.m2)
exp.u.mr

plot(exp.u.mr)

## dataset "j"
exp.j.m1 <- drm(response ~ conc, data = exp.az[exp.az$exp=="j",], fct = LL.3())
exp.j.m2 <- loess(response ~ conc, data = exp.az[exp.az$exp=="j",], degree = 1)
exp.j.mr <- mrdrm(exp.j.m1, exp.j.m2)
exp.j.mr

plot(exp.j.mr)

## dataset "k"
exp.k <- exp.az[exp.az$exp=="k",]

## Correcting doses
exp.k[43:48,1]<-0.0317826731495085
exp.k[exp.k$conc==0.1,1] <- 0.100433247152447
exp.k[exp.k$conc==0.3,1] <- 0.317369061001732
exp.k[exp.k$conc==1,1] <- 1.00288623276547
exp.k[exp.k$conc==3.2,1] <- 3.16912049553889
exp.k[exp.k$conc==10,1] <- 10.0144207659029
exp.k[exp.k$conc==31.6,1] <- 31.6455696202532 


exp.k.m1 <- drm(response ~ conc, data = exp.k, fct = LL.3())
exp.k.m2 <- loess(response ~ conc, data = exp.k, degree = 1, span = 0.5)
exp.k.mr <- mrdrm(exp.k.m1, exp.k.m2)
exp.k.mr

plot(exp.k.mr, type = "all")


## dataset "l"
exp.l.m1 <- drm(response ~ conc, data = exp.az[exp.az$exp=="l",], fct = LL.3())
exp.l.m2 <- loess(response ~ conc, data = exp.az[exp.az$exp=="l",], degree = 1)
exp.l.mr <- mrdrm(exp.l.m1, exp.l.m2)
exp.l.mr

plot(exp.l.mr)

## dataset "z"
exp.z.m1 <- drm(response ~ conc, data = exp.az[exp.az$exp=="z",], fct = LL.3())
exp.z.m2 <- loess(response ~ conc, data = exp.az[exp.az$exp=="z",], degree = 1, span = 0.35)
exp.z.mr <- mrdrm(exp.z.m1, exp.z.m2)
exp.z.mr

plot(exp.z.mr, type = "all", conLevel = 5, broken = TRUE)
plot(exp.z.m1, add = TRUE, type = "none", conLevel = 5, broken = TRUE, lty = 2)

ED(exp.z.mr, 50, interval="approximate")
ED(exp.z.m1, 50, ci="delta")

## Figure 1
par(mfrow = c(2, 2))

plot(exp.a.m1, type = "all", conLevel=1, broken=TRUE, xlab="", ylab="Response", main = "A")
plot(exp.k.m1, type = "all", broken=TRUE, xlab="", ylab="", main = "B")
plot(exp.x.m1, type = "all", broken=TRUE, xlab="Dose", ylab="Response", main = "C")
plot(exp.z.m1, type = "all", conLevel = 5, broken = TRUE, xlab="Dose", ylab="", main = "D")

par(mfrow = c(1, 1))


## Figure 3
par(mfrow = c(1, 2))

plot(exp.x.mr, type="all", broken = TRUE, ylim=c(0.05,0.84), xlab = "Dose", ylab = "Response", main = "C")
plot(exp.x.m1, add = TRUE, type = "none", broken = TRUE, lty = 2)

plot(exp.z.mr, pava = TRUE, type = "all", conLevel = 5, broken = TRUE, ylim=c(0.05,0.84), xlab = "Dose", ylab = "", main = "D")
plot(exp.z.m1, add = TRUE, type = "none", conLevel = 5, broken = TRUE, lty = 2)


par(mfrow = c(1, 1))


## Table 3
ED(exp.x.m1, c(10,20,50), ci="delta")
ED(exp.x.mr, c(10,20,50), interval="approximate")

ED(exp.z.m1, c(10,20,50), ci="delta")
ED(exp.z.mr, c(10,20,50), interval="approximate", cgridsize=100)


}

## ryegrass dataset (continuous response)
ryegrass.m1 <- drm(rootl~conc, data=ryegrass, fct=LL.4())
ryegrass.m2 <- loess(rootl~conc, data=ryegrass, degree=1)
ryegrass.mr <- mrdrm(ryegrass.m1, ryegrass.m2)
ryegrass.mr

predict(ryegrass.mr)
ED(ryegrass.mr, c(10, 50, 90), interval = "approximate")

## lettuce dataset (continuous response)
lettuce.m1 <- drm(weight~conc, data = lettuce, fct = LL.3())
lettuce.m2 <- loess(weight~conc, data = lettuce, degree = 1, span = 0.5)
lettuce.mr <- mrdrm(lettuce.m1, lettuce.m2)
lettuce.mr

plot(lettuce.mr, type = "all")

## Obtaining ED values (not working with bisection method)
ED(lettuce.mr, c(10,50), interval = "approximate", method = "grid", 
upper = predict(lettuce.mr, data.frame(conc=0)))
}

\keyword{models}
\keyword{nonlinear}

