#' @importFrom stats runif
#'
#' @name pbivgeocure
#' @aliases pbivgeocure
#'
#' @title Joint Cumulative Function for the Basu-Dhar Bivariate Geometric Distribution in Presence of Cure Fraction
#'
#' @description This function computes the joint cumulative function of the Basu-Dhar bivariate geometric distribution assuming arbitrary parameter values in presence of cure fraction.
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
#' de Oliveira, R. P., Achcar, J. A., Peralta, D., & Mazucheli, J. (2018). Discrete and continuous bivariate lifetime models in presence of cure rate: a comparative study under Bayesian approach. \emph{Journal of Applied Statistics}, 1-19.
#'
#' @param x matrix or vector containing the data. If x is a matrix then it is considered as x the first column and y the second column (y argument need be setted to NULL). Additional columns and y are ignored.
#'
#' @param y vector containing the data of y. It is used only if x is also a vector. Vectors x and y should be of equal length.
#'
#' @param theta vector (of length 3) containing values of the parameters \eqn{\theta_1, \theta_2} and \eqn{\theta_{3}} of the Basu-Dhar bivariate Geometric distribution. The parameters are restricted to \eqn{0 < \theta_i < 1, i = 1,2} and \eqn{0 < \theta_{3} \le 1}.

#' @param phi vector (of length 4) containing values of the cure fraction incidence parameters \eqn{\phi_{11}, \phi_{10}, \phi_{01}} and \eqn{\phi_{00}}. The parameters are restricted to \eqn{\phi_{11} + \phi_{10} + \phi_{01} + \phi_{00}= 1}.
#'
#' @param lower.tail logical; If TRUE (default), probabilities are \eqn{P(X \le x, Y \le y)} otherwise \eqn{P(X > x, Y > y)}.
#'
#' @return \code{\link[BivGeo]{pbivgeocure}} gives the values of the cumulative function in presence of cure fraction.
#'
#' @return Invalid arguments will return an error message.
#'
#' @usage
#'
#' pbivgeocure(x, y, theta, phi, lower.tail = TRUE)
#'
#' @details
#'
#' The joint cumulative function for a random vector (\eqn{X}, \eqn{Y}) following a Basu-Dhar bivariate geometric distribution in presence of cure fraction could be written as:
#' \deqn{P(X \le x, Y \le y) = 1 - (\phi_{11} + \phi_{10}) (\theta_1 \theta_3)^x - (\phi_{01} + \phi_{00}) - (\phi_{11} + \phi_{01}) (\theta_2 \theta_3)^y - (\phi_{10} + \phi_{00})}
#' \deqn{+ \phi_{11} (\theta_{1}^{x} \theta_{2}^{y}\theta_{3}^{\max(x,y)}) + \phi_{10} (\theta_1 \theta_{3})^x + \phi_{01} (\theta_2 \theta_{3})^y + \phi_{00}}
#' and the joint survival function is given by:
#' \deqn{P(X > x, Y > y) = \phi_{11} (\theta_{1}^{x} \theta_{2}^{y}\theta_{3}^{\max(x,y)}) + \phi_{10} (\theta_1 \theta_{3})^x + \phi_{01} (\theta_2 \theta_{3})^y + \phi_{00}}
#'
#' @examples
#'
#' # If lower.tail = TRUE:
#'
#' pbivgeocure(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), phi = c(0.2, 0.3, 0.3, 0.2), lower.tail = TRUE)
#' # [1] 0.159456
#'
#' matr 	<- 	 matrix(c(1,2,3,5), ncol = 2)
#' pbivgeocure(x=matr,y=NULL,theta=c(0.2, 0.4, 0.7),phi=c(0.2, 0.3, 0.3, 0.2),lower.tail = TRUE)
#' # [1] 0.1684877 0.1957496
#'
#' # If lower.tail = FALSE:
#'
#' pbivgeocure(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), phi = c(0.2, 0.3, 0.3, 0.2), lower.tail = FALSE)
#' # [1] 0.268656
#'
#' matr 	<- 	 matrix(c(1,2,3,5), ncol = 2)
#' pbivgeocure(x=matr,y=NULL,theta=c(0.2, 0.4, 0.7),phi=c(0.2, 0.3, 0.3, 0.2),lower.tail = FALSE)
#' # [1] 0.2494637 0.2064101
#'
#' @source
#'
#' \code{\link[BivGeo]{pbivgeocure}} is calculated directly from the definition.
#'
#' @seealso
#'
#' \code{\link[stats]{Geometric}} for the univariate geometric distribution.
#'
#' @rdname pbivgeocure
#' @export

pbivgeocure <- function(x, y = NULL, theta = c(), phi = c(), lower.tail = TRUE)
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

	if(theta[1] <= 0 | theta[1] >= 1) return('theta1 out of bounds')
	if(theta[2] <= 0 | theta[2] >= 1) return('theta2 out of bounds')
	if(theta[3] <= 0 | theta[3] > 1)  return('theta3 out of bounds')

	if(phi[1] <= 0 | phi[1] >= 1) return('phi11 out of bounds')
	if(phi[2] <= 0 | phi[2] >= 1) return('phi10 out of bounds')
	if(phi[3] <= 0 | phi[3] >= 1) return('phi01 out of bounds')
	if(phi[4] <= 0 | phi[4] >= 1) return('phi00 out of bounds')

	if(phi[1] + phi[2] + phi[3] + phi[4] != 1) return('wrong values for the parameters phi11, phi10, phi01 and phi00')

	# Max values between x and y
	z 	  <-  pmax(x0, y0)

	sf    <-  phi[1] * (theta[1]^x0 * theta[2]^y0 * theta[3]^z) + phi[2] * (theta[1] * theta[3])^x0 + phi[3] * (theta[2] * theta[3])^y0 + phi[4]

	sfx   <-  (phi[1] + phi[2]) * (theta[1] * theta[3])^x0 + (phi[3] + phi[4])
	sfy   <-  (phi[1] + phi[3]) * (theta[2] * theta[3])^y0 + (phi[2] + phi[4])

	cdf   <-  1 - sfx - sfy + sf

	if(lower.tail) return(cdf) else return(sf)
}
