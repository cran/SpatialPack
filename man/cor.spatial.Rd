\name{cor.spatial}
\alias{cor.spatial}
\title{Tjostheim's Coefficient}
\description{
Computes Tjostheim's coefficient for two spatial sequences observed over
the same locations on the plane.
}
\usage{
cor.spatial(x, y, coords)
}
\arguments{
  \item{x}{an n-dimensional vector of data values.}
  \item{y}{an n-dimensional vector of data values.}
  \item{coords}{an n x 2 matrix containing coordinates
  of the n data locations in each row.}
}
\details{
The implemented  technique is a nonparametric coefficient that summarizes
the association between two spatial variables. This coefficient was first
introduced by Tjostheim (1978) and later generalized by Hubert and Golledge (1992).
The computation of the coefficient is based on the construction of ranks
associated to suitable modifications of the coordinates. Tjostheim's coefficient
is a variant of the correlation coefficient (\code{\link{cor}}) to be used
in a spatial statistics context.
}
\value{
  Tjostheim's coefficient. The variance is returned as the attribute \code{"variance"}.
}
\references{
  Tjostheim, D., (1978),
  A measure of association for spatial variables.
  \emph{Biometrika} \bold{65}, 109--114.
  
 Hubert, L., Golledge, R.G., (1982),
 Measuring association between spatially defined variables: Tjostheim's coefficient index and some extensions.
 \emph{Geographical Analysis} \bold{14}, 273--278.
}
\examples{
## This example uses the camg data set, available in the geoR Package
## (Ribeiro and Diggle, 2001).

## required to loading the camg data set from geoR package 
library(geoR)
data(camg)

## defining the variables elevation and ca020 from the camg data set
x <- camg$elevation
y <- camg$ca020

## extracting the coordinates from the camg data set
coords <- camg[,1:2]

## computing Tjostheim's coefficient
cor.spatial(x, y, coords)
}
\keyword{multivariate}
