\name{brea_mcmc}

\alias{brea_mcmc}

\title{Bayesian Discrete Survival Inference}

\description{
This function performs Metropolis-Hastings exploration of the posterior distribution for Bayesian discrete time survival models. The models may have competing risks, semiparametric covariate effects (like arbitrary baseline hazards), and random effects shared across repeated events and correlated across competing risks.
}

\usage{
brea_mcmc(x, y, S = 1000L, priors = NULL, n = NULL, K = NULL, store_re = FALSE)
}

\arguments{
  \item{x}{an \code{integer} matrix or dataframe all columns of which are factors specifying the values of the (discretized) predictors at each person-time point; specifically, the (i,j) entry is the value of predictor j at discrete person-time observation i}
  \item{y}{an \code{integer} matrix whose (i,j) entry counts the number of events of type j occurring at discrete person-time observation i}
  \item{S}{the number of MCMC iterations to perform}
  \item{priors}{a list with one element for each predictor variable (column of \code{x}) specifying the prior type to use for that predictor; see Details for more information}
  \item{n}{a vector of positive integers with length equal to the number of person-time observations whose entries equal the number of replicated observations that row stands for (defaults to 1 for each observation)}
  \item{K}{a vector of positive integers with length equal to the number of predictors giving the number of distinct values each discretized predictor may assume; this is not used if a dataframe of factors is provided for \code{x}}
  \item{store_re}{if \code{TRUE}, the random effects are stored from each MCMC iteration, and if \code{FALSE} they are not stored}
}

\details{
The data provided to the \code{brea_mcmc} function is specified at the person-time level: there is one row in \code{x} and \code{y} for each discrete time point each person or thing was at risk for event occurrence. All predictors in \code{x} must be encoded as factors (or their corresponding integer codes in the case that \code{x} is an \code{integer} matrix). The underlying type of predictor is specified in the \code{priors} argument, which is a list with one element for each predictor variable which specifies both the type of that predictor and the prior distribution to use. The allowed predictor types are:
\itemize{
  \item \code{"cat"} for categorical variables. The first element of the prior specification list is the string \code{"cat"}, and second element is a positive conditional prior standard deviation parameter.
}
\itemize{
  \item \code{"gmrf"} for underlying continuous predictors; continuous predictors should be \code{cut()} before being included into \code{x}; Gaussian Markov random field (GMRF) priors are then used to smooth the effects of adjacent categories of the discretized continuous predictor. The first element of the prior specification list is the string \code{"gmrf"}, the second element is a prior degrees of freedom for the scaled inverse chi-squared prior for the random walk increment variance, and the third element is a prior scale for the scaled inverse chi-squared.
}
\itemize{
  \item \code{"re"} for variables (like subject id numbers) that represent random effects. The first element of the prior specification list is the string \code{"re"}, the second element is a prior degrees of freedom for an inverse Wishart prior for the random effects covariance, and the third element is a prior scale matrix for the random effects covariance.
}
\itemize{
  \item \code{"zero"} for predictors that are not used in the current MCMC run. This is provided as a convenient way to exclude predictors from certain runs. The first and only element of the prior specification list is the string \code{"zero"}.
}
}

\value{
A list providing the values from each of the \code{S} MCMC iterations for the intercept, the other linear predictor parameters, standard deviations of the GMRF random walk increments, and covariances of the random effects:
\item{b_0_s}{an \code{R} by \code{S} matrix (where \code{R} is the number of competing risks) whose (r,s) entry is the intercept for risk r in iteration s}
\item{b_m_s}{a list whose m^th element is the \code{R} by \code{K[m]} by \code{S} array whose (r,k,s) entry is the value of the linear predictor contribution of value k for predictor m on risk r at MCMC iteration s}
\item{s_m_s}{a list whose m^th element is, in the case of a GMRF prior, the \code{R} by \code{S} matrix whose (r,s) entry is the random walk standard deviation for predictor m and competing risk r at MCMC iteration s, or, in the case of a random effects priors, the \code{R} by \code{R} by \code{S} array whose (,,s) entry is the random effects covariance matrix for predictor m at MCMC iteration s}
\item{b_m_a}{a list whose m^th element is the length \code{K[m]} vector giving the number of accepted Metropolis proposals for each value of predictor m across all \code{S} MCMC iterations}
}

\author{Adam King}

\examples{

# leukemia remission times data from Table 1 of "Regression Models and
# Life-Tables" (Cox, 1972)

# times of event occurrence or right censoring:
study_time <- c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35,
                1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23)

# event/censoring indicator (1 for event occurrence, 0 for right censoring):
cens <- c(0,1,1,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,
          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# treatment assignment (1 for treatment, 2 for control):
treat <- c(rep(1L,21),rep(2L,21))

# total number of person-time observations among all subjects:
N <- sum(study_time)

# create and fill in person-time predictor and event indicator matrices:
x <- matrix(0L,N,2)  # first column is time, second column is treatment
y <- matrix(0L,N,1)
next_row <- 1  # next available row in the person-time matrices
for (i in seq_along(study_time)) {
  x[next_row:(next_row + study_time[i] - 1),] <- cbind(seq_len(study_time[i]),
                                                    rep(treat[i],study_time[i]))
  y[next_row:(next_row + study_time[i] - 1),] <- c(rep(0L,study_time[i] - 1),
                                                   as.integer(cens[i]))
  next_row <- next_row + study_time[i]
}

# group the time variable into 6-week intervals:
x[,1] <- cut(x[,1],seq(0,36,6),labels=FALSE)

# use GMRF prior for time, and categorical prior for treatment group:
priors <- list(list("gmrf",3,.01),list("cat",2))

# run 1,100 MCMC iterations:
fit <- brea_mcmc(x, y, 1100, priors)

# look at the structure of the returned posterior samples and acceptance counts:
str(fit)

# approximate posterior mean hazard ratio:
ss <- 101:1100  # use last 1,000 samples, discarding first 100 as burn-in
exp(2*mean(fit$b_m_s[[2]][1,2,ss]))
# hazard ratio of approximately 5

}

\keyword{survival}
\keyword{regression}
