% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbivgeocure.R
\name{dbivgeocure}
\alias{dbivgeocure}
\title{Joint Probability Mass Function for the Basu-Dhar Bivariate Geometric Distribution in Presence of Cure Fraction}
\source{
\code{\link[BivGeo]{dbivgeocure}} is calculated directly from the definition.
}
\usage{
dbivgeocure(x, y, theta, phi11, log = FALSE)
}
\arguments{
\item{x}{matrix or vector containing the data. If x is a matrix then it is considered as x the first column and y the second column (y argument need be setted to NULL). Additional columns and y are ignored.}

\item{y}{vector containing the data of y. It is used only if x is also a vector. Vectors x and y should be of equal length.}

\item{theta}{vector (of length 3) containing values of the parameters \eqn{\theta_1, \theta_2} and \eqn{\theta_{3}} of the Basu-Dhar bivariate Geometric distribution. The parameters are restricted to \eqn{0 < \theta_i < 1, i = 1,2} and \eqn{0 < \theta_{3} \le 1}.}

\item{phi11}{real number containing the value of the cure fraction incidence parameter \eqn{\phi_{11}} restricted to \eqn{0 < \phi_{11} < 1} and \eqn{\phi_{11} + \phi_{10} + \phi_{01} + \phi_{00}= 1} where \eqn{\phi_{10}, \phi_{01}} and \eqn{\phi_{00}} are the complementary cure fraction incidence parameters for the joint cdf and sf functions.}

\item{log}{logical argument for calculating the log probability or the probability function. The default value is FALSE.}
}
\value{
\code{\link[BivGeo]{dbivgeocure}} gives the values of the probability mass function in presence of cure fraction.

Invalid arguments will return an error message.
}
\description{
This function computes the joint probability mass function of the Basu-Dhar bivariate geometric distribution assuming arbitrary parameter values in presence of cure fraction.
}
\details{
The joint probability mass function for a random vector (\eqn{X}, \eqn{Y}) following a Basu-Dhar bivariate geometric distribution in presence of cure fraction could be written as:
\deqn{P(X = x, Y = y) = \phi_{11}(\theta_{1}^{x - 1} \theta_{2}^{y - 1} \theta_{3}^{z_1} - \theta_{1}^{x} \theta_{2}^{y - 1} \theta_{3}^{z_2} - \theta_{1}^{x - 1} \theta_{2}^{y} \theta_{2}^{z_3} + \theta_{1}^{x} \theta_{2}^{y} \theta_{3}^{z_4})}
where \eqn{x,y > 0} are positive integers and \eqn{z_1 = \max(x - 1, y - 1),z_2 = \max(x, y - 1), z_3 = \max(x - 1, y), z_4 = \max(x, y)}.
}
\examples{

# If log = FALSE:

dbivgeocure(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), phi11 = 0.4, log = FALSE)
# [1] 0.064512

# If log = TRUE:

dbivgeocure(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), phi11 = 0.4, log = TRUE)
# [1] -2.740904

}
\references{
Basu, A. P., & Dhar, S. K. (1995). Bivariate geometric distribution. \emph{Journal of Applied Statistical Science}, \bold{2}, 1, 33-44.

Achcar, J. A., Davarzani, N., & Souza, R. M. (2016). Basu–Dhar bivariate geometric distribution in the presence of covariates and censored data: a Bayesian approach. \emph{Journal of Applied Statistics}, \bold{43}, 9, 1636-1648.

de Oliveira, R. P., & Achcar, J. A. (2018). Basu-Dhar's bivariate geometric distribution in presence of censored data and covariates: some computational aspects. \emph{Electronic Journal of Applied Statistical Analysis}, \bold{11}, 1, 108-136.

de Oliveira, R. P., Achcar, J. A., Peralta, D., & Mazucheli, J. (2018). Discrete and continuous bivariate lifetime models in presence of cure rate: a comparative study under Bayesian approach. \emph{Journal of Applied Statistics}, 1-19.
}
\seealso{
\code{\link[stats]{Geometric}} for the univariate geometric distribution.
}
\author{
Ricardo P. Oliveira \email{rpuziol.oliveira@gmail.com}

Jorge Alberto Achcar \email{achcar@fmrp.usp.br}
}
