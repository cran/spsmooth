%     The spsmooth R package
%     An extension for mgcv.
%     Copyright (C) 2012 Wesley Burr
%
%     Written by Wesley Burr.
%
%     The data included with this package is a subselection of the 
%     NMMAPSdata / NMMAPSlite package database. It is included for 
%     demonstration purposes only.
%
%     This file is part of the spsmooth package for R.
%
%     The spsmooth package is free software: you can redistribute it and
%     or modify it under the terms of the GNU General Public License as 
%     published by the Free Software Foundation, either version 2 of the 
%     License, or any later version.
%
%     The spsmooth package is distributed in the hope that it will be 
%     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
%     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with spsmooth.  If not, see <http://www.gnu.org/licenses/>.
%
%     If you wish to report bugs please contact the author. 
%     Wesley Burr
%     wburr@mast.queensu.ca
%     239 Jeffery Hall, Queen's University, Kingston, Ontario
%     Canada, K7L 3N6

\name{chicagoNMMAPS}
\alias{chicagoNMMAPS}
\alias{chic}
\title{Chicago 1987-2000 from the NMMAPS database}
\description{
  Daily measurements of death, temperature, dewpoint, and PM10 pollution
  from Chicago, 1987 to 2000. Data taken from the NMMAPS database
  using \code{\link[NMMAPSlite]{NMMAPS}}, a package written by Roger Peng
  and available on CRAN: http://cran.r-project.org/web/packages/NMMAPSlite/index.html.
}
\usage{chicagoNMMAPS}
\format{ A data frame indicating the date, temperature, dewpoint, and PM10
    pollution from Chicago. Frame also includes the mask variable needed
    for estimating a sp basis, see \code{\link{spsmooth}}.
}
\examples{
    # This provides the code that generated this set originally
    # library("NMMAPSlite")
    # initDB()   # ** this creates a directory for a data archive
    # chic <- readCity("chic", collapseAge=TRUE)
    # chic <- chic[, c("date", "death", "time", "tmpd", "dptp", "dow", "pm10tmean")]
    # mask <- rep(FALSE, length(chic[, "pm10tmean"]))
    # mask[!is.na(chic[, "pm10tmean"])] <- TRUE
    # chic <- cbind(chic, mask)

    # save(chic, file="chicagoNMMAPS.RData")
}
\keyword{datasets}
