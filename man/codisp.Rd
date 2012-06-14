\name{codisp}
\alias{codisp}
\title{Codispersion Coefficient}
\description{
Computes the codispersion coefficient between two spatial variables for
a given number of classes for the lag distance.}
\usage{
codisp(x, y, coords, nclass = 13)
}
\arguments{
  \item{x}{an n-dimensional vector of data values.}
  \item{y}{an n-dimensional vector of data values.}
  \item{coords}{an n x 2 matrix containing coordinates of the n data
  locations in each row.}
  \item{nclass}{a single number giving the number of cells for the codispersion
  coefficient. The default is 13. If this argument is NULL Sturges' formula is used.}
}
\details{
The procedure computes the codispersion coefficient for two spatial sequences
defined on general (non-rectangular) grids. First, a given number of bins are
constructed for the lag distance. Then the codispersion is computed for each bin. 
}
\value{
  A list with class \code{"codisp"} containing the following components:
  \item{coef}{a vector of size \code{nclass} containing the values of the codispersion
  coefficient.}
  \item{upper.bounds}{upper bounds of the intervals constructed to compute the
  codispersion coefficient.}
  \item{card}{number of elements in each interval generated to compute the
  codispersion coefficient.}
  
  The function \code{plot} can be used to obtain a graph of the codispersion coefficient
  versus the lag distance.
}
\references{
  Matheron, G. (1965),
  \emph{Les Variables Regionalisees et leur Estimation}.
  Masson, Paris.
  
  Rukhin, A., Vallejos, R. (2008),
  Codispersion coefficient for spatial and temporal series.
  \emph{Statistics and Probability Letters} \bold{78}, 1290--1300.
  
  Vallejos, R. (2008).
  Assessing the association between two spatial or temporal sequences.
  \emph{Journal of Applied Statistics} \bold{35}, 1323--1343.
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

## computing the codispersion coefficient 
z <- codisp(x, y, coords)
z

## plotting the codispersion coefficient vs. the lag distance
plot(z)
}
\keyword{multivariate}
