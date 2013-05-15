\name{spsmooth-package}
\alias{spsmooth-package}
\alias{spsmooth}
\docType{package}
\title{An Extension Package for mgcv.}
\description{Produces Slepian smoother basis objects for use
in generalized additive models. Requires \link{mgcv}.}
\details{
\tabular{ll}{
Package: \tab spsmooth\cr
Type: \tab Package\cr
Version: \tab 1.1-0\cr
Date: \tab 2013-05-15\cr
License: \tab GPL (>= 2)\cr
Depends: \tab mgcv (>= 1.7.18), R (>= 2.15)
}
}
\author{
Wesley Burr
}
\references{
    Hastie T.J. & Tibshirani, R.J. (1990) Generalized Additive Models.
    Chapman and Hall.

    Wood S.N. (2006) Generalized Additive Models: An Introduction with 
    R. Chapman and Hall/CRC Press.
}
\note{
    If you are using a \code{sp} smoother on a data set which is missing
    a significant portion (i.e. 8 in 10 points) \code{\link[mgcv]{gam}}
    may have convergence issues due to lack of identifiability. This can
    typically be corrected by increasing \code{epsilon} and \code{maxit}
    slightly; empirical results indicate that changing \code{epsilon}
    from \code{1e-6} to \code{2e-6} and changing \code{maxit} from
    \code{50} to \code{150} is sufficient for most cases. This change 
    can be done in the call to \code{\link[mgcv]{gam}} by including
    \code{control=list(maxit=150, epsilon=2e-6)}. If convergence
    does fail there will be a warning.
}
\keyword{ package }\keyword{ smooth }\keyword{ model }
\examples{
    # The only function included in this package which users are intended
    # to call directly is summary().

    # This example requires the mgcv package.
    library("mgcv")
    data("chicagoNMMAPS")   # loads as object 'chic'

    # Model 1 includes a standard cubic regression fixed-df smooth function of time,
    # with k chosen to give the same EDF as the second model
    Model1 <- as.formula("death ~ pm10tmean + s(tmpd, k=6, fx=TRUE) + s(dptp, k=3, fx=TRUE) + 
              as.factor(dow) + s(time, bs='cr', k=196, fx=TRUE)")

    # Model 2 uses spsmooth and changes the smooth function of time to use 
    # the Slepian basis, 2NW-1 basis vectors.
    #
    # Choosing W: there are two approaches that can be taken. The first is 
    # to select the df first, i.e. to constrain 2NW-1=df. As N is fixed
    # for a given problem, this constrains W. The second is to choose W
    # based on model concerns, and then given N, this constrains df. In
    # the case of Model 2, we have chosen W = 7/365, i.e. 7 cycles/year,
    # and there is no correction factor for delta-t needed as our data
    # is daily.
    #
    # W needs to be in the same units as your data series. 
    #
    Model2 <- as.formula("death ~ pm10tmean + s(tmpd, k=6, fx=TRUE) + s(dptp, k=3, fx=TRUE) + 
              as.factor(dow) + s(time, bs='sp', xt=list(W=7/365, mask=chic$mask))")

    # Using the quasi(log, mu) link, and demonstrating the use of control arguments
    # that can be used in the case of poor convergence
    fit1 <- gam(Model1, data=chic, na.action=na.exclude, family=quasi(log,mu), 
                control=list(epsilon=1e-6, maxit=50))
    fit2 <- gam(Model2, data=chic, na.action=na.exclude, family=quasi(log,mu), 
                control=list(epsilon=1e-6, maxit=50))

    summary(fit1)
    summary(fit2)
}
