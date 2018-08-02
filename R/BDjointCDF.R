#' @importFrom stats runif
#'
#' @name BDjointCDF
#' @aliases pbdbivgeo
#'
#' @title Joint Cumulative Function for the Basu-Dhar Bivariate Geometric Distribution
#'
#' @description This function computes the joint cumulative function of the Basu-Dhar bivariate geometric distribution assuming arbitrary parameter values.
#'
#' @author Ricardo P. Oliveira \email{rpuziol.oliveira@gmail.com}
#' @author Jorge Alberto Achcar \email{achcar@fmrp.usp.br}
#'
#' @references
#'
#' Basu, A. P., & Dhar, S. K. (1995). Bivariate geometric distribution. \emph{Journal of Applied Statistical Science}, \bold{2}, 1, 33-44.
#'
#' Achcar, J. A., Davarzani, N., & Souza, R. M. (2016). Basu–Dhar bivariate geometric distribution in the presence of covariates and censored data: a Bayesian approach. \emph{Journal of Applied Statistics}, \bold{43}, 9, 1636-1648.
#'
#' de Oliveira, R. P., & Achcar, J. A. (2018). Basu-Dhar's bivariate geometric distribution in presence of censored data and covariates: some computational aspects. \emph{Electronic Journal of Applied Statistical Analysis}, \bold{11}, 1, 108-136.
#'
#' @param x matrix or vector containing the data. If x is a matrix then it is considered as x the first column and y the second column (y argument need be setted to NULL). Additional columns and y are ignored.
#'
#' @param y vector containing the data of y. It is used only if x is also a vector. Vectors x and y should be of equal length.
#'
#' @param theta vector (of length 3) containing values of the parameters \eqn{\theta_1, \theta_2} and \eqn{\theta_{12}} of the Basu-Dhar bivariate Geometric distribution. The parameters are restricted to \eqn{0 < \theta_i < 1, i = 1,2} and \eqn{0 < \theta_{12} \le 1}.
#'
#' @param lower.tail logical; If TRUE (default), probabilities are \eqn{P(X \le x, Y \le y)} otherwise \eqn{P(X > x, Y > y)}.
#'
#' @return \code{\link[BivGeo]{pbdbivgeo}} gives the values of the cumulative function.
#'
#' @return Invalid arguments will return an error message.
#'
#' @usage
#'
#' pbdbivgeo(x, y, theta, lower.tail = TRUE)
#'
#' @details
#'
#' The joint cumulative function for a random vector (\eqn{X}, \eqn{Y}) following a Basu-Dhar bivariate geometric distribution could be written as:
#' \deqn{P(X \le x, Y \le y) = 1 - (\theta_{1}\theta_{12})^{x} - (\theta_{2}\theta_{12})^{y} + \theta_{12}^{\max(x,y)}}
#' and the joint survival function is given by:
#' \deqn{P(X > x, Y > y) = \theta_{1}^{x}\theta_{2}^{y} \theta_{12}^{\max(x,y)}}
#'
#' @examples
#'
#' # If x and y are integer numbers:
#'
#' pbdbivgeo(1, 2, c(0.2, 0.4, 0.7), lower.tail = TRUE)
#' # [1] 0.79728
#'
#' # If x is a matrix:
#'
#' matr 	<- 	 matrix(c(1,2,3,5), ncol = 2)
#' pbdbivgeo(matr, y = NULL, c(0.2,0.4,0.7), lower.tail = TRUE)
#' # [1] 0.8424384 0.9787478
#'
#' # If lower.tail = FALSE:
#'
#' pbdbivgeo(1, 2, c(0.2, 0.4, 0.7), lower.tail = FALSE)
#' # [1] 0.01568
#'
#' matr 	<- 	 matrix(c(1,2,3,5), ncol = 2)
#' pbdbivgeo(matr, y = NULL, c(0.2,0.4,0.7), lower.tail = FALSE)
#' # [1] 4.390400e-03 6.884147e-05
#'
#' @source
#'
#' \code{\link[BivGeo]{pbdbivgeo}} is calculated directly from the definition.
#'
#' @seealso
#'
#' \code{\link[stats]{Geometric}} for the univariate geometric distribution.
#'
#' @rdname BDjointCDF
#' @export

pbdbivgeo <- function(x, y = NULL, theta = c(), lower.tail = TRUE)
{
	if(is.matrix(x))
	{
		x0	<-	x[,1]
		y0	<-	x[,2]
	}
	else if(is.vector(x) & is.vector(y))
	{
		if(length(x)==length(y))
		{
			x0	<-	x
			y0	<-	y
		}
		else
		{
			stop('lengths of x and y are not equal')
		}
	}
	else
	{
		stop('x is not a matrix or x and y are not vectors')
	}

  	# Max values between x and y
  	z 	  <-  pmax(x0, y0)

  	sf    <-  theta[1]^x0 * theta[2]^y0 * theta[3]^z

  	sfx   <-  (theta[1] * theta[3])^x0
  	sfy   <-  (theta[2] * theta[3])^y0

  	cdf   <-  1 - sfx - sfy + sf

  	if(lower.tail) return(cdf) else return(sf)
}
