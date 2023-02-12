# Bayesian_RelativeRisk_FirstAid
Code for an R function in the style of the `BayesianFirstAid' package to analyse contingency tables, and estimate relative risk of a given outcome.

#'Relative Risk test for a contingency table, created for compatibility with BayesianFirstAid
#'
#'\code{bayes.rr.test} estimates relative risk of a given outcome against a reference condition 
#' for a given 2-dimensional contingency table. 
#' This is intended as a replacement to the frequentist \code{\link{mosaic::relrisk()}}.
#'
#'Given a contingency table, \code{bayes.rr.test} 
#'estimates \ifelse{latex}{\eqn{\theta_{1...m}}}{\eqn{\theta}[1...m]}, the relative
#'frequencies of success for each of the \eqn{m} groups. The following model is
#'assumed for each group:
#'
#'\deqn{x \sim \mathrm{Binom}(\theta, n)}{x ~ Binomial(\theta, n)} \deqn{\theta 
#'\sim \mathrm{Beta}(1, 1)}{\theta ~ Beta(1, 1)}
#'
#'
#'Here the prior on the \eqn{\theta}s is a non-informative \eqn{\mathrm{Beta}(1,
#'1)}{Beta(1, 1)} distribution which is identical to a \eqn{\mathrm{Uniform}(0, 
#'1)}{Uniform(0, 1)} distribution. By \code{plot}ing and looking at a 
#'\code{summary} of the object returned by \code{bayes.prop.test} you can get 
#'more information about the shape of the posterior and the posterior predictive
#'distribution. \code{\link{model.code}} prints out the corresponding R code 
#'underlying \code{bayes.prop.test} which can be copy-n-pasted into an R script 
#'and modified, for example, changing the prior on \eqn{\theta}.
#'
#'The \code{print} and \code{plot} function will only work well with a small 
#'number of groups (2 to 6). If you have more groups you might want to run the 
#'model with a small number of groups and then print the model code using 
#'\code{\link{model.code}} and fit the model using that code. The current model does
#'not assume any dependency between the groups, if this is an unreasonable assumption
#'you might want to modify the model code (from \code{\link{model.code}}) to 
#'include a dependency between the groups (see 
#'\href{http://lingpipe-blog.com/2009/09/23/bayesian-estimators-for-the-beta-binomial-model-of-batting-ability/}{here}
#'for an example).
#'
#'@param x a vector of counts of successes, a one-dimensional table with two 
#'  entries, or a two-dimensional table (or matrix) with 2 columns, giving the 
#'  counts of successes and failures, respectively.
#'@param n a vector of counts of trials; ignored if x is a matrix or a table.
#'@param comp.rr a vector of fixed relative frequencies of success to compare
#'  with the estimated relative frequency of success. The length of 
#'  \code{comp.rr} must be the same as the number of groups specified by 
#'  \code{x}, and its elements must be greater than 0 and less than 1. This 
#'  argument fills a similar role as \code{p} in \code{\link{prop.test}}.
#'@param alternative ignored and is only retained in order to mantain 
#'  compatibility with \code{\link{prop.test}}.
#'@param cred.mass the amount of probability mass that will be contained in 
#'  reported credible intervals. This argument fills a similar role as 
#'  \code{cred.mass} in \code{\link{prop.test}}.
#'@param correct ignored and is only retained in order to mantain compatibility 
#'  with \code{\link{prop.test}}
#'@param n.iter The number of iterations to run the MCMC sampling.
#'@param progress.bar The type of progress bar. Possible values are "text", 
#'  "gui", and "none".
#'@param p same as \code{comp.rr} and is only retained in order to mantain 
#'  compatibility with \code{\link{prop.test}}.
#'@param cred.mass same as \code{cred.mass} and is only retained in order to 
#'  mantain compatibility with \code{\link{prop.test}}.
#'  
#'  
#'@return A list of class \code{bayes_prop_test} that contains information about
#'  the analysis. It can be further inspected using the functions 
#'  \code{summary}, \code{plot}, \code{\link{diagnostics}} and 
#'  \code{\link{model.code}}.
#'  
#' @examples
#' 
#' 
#' # Data from Muller, F. H., Tobakmissbrauch und Lungencarcinom,
#' # Zeit. f. Krebsforsch. 49, 57-85, 1939. One of the early papers
#' # investigating the relation between smoking and lung cancer.
#'
#' # Number of heavy smokers in one group of 86 lung cancer patients
#' # and one group of 86 healthy individuals.
#' no_heavy_smokers <- c(56, 31)
#' no_cases <- c(86, 86)
#'
#' bayes.prop.test(no_heavy_smokers, no_cases)
#' 
#' # Save the return value in order to inspect the model result further.
#' fit <- bayes.prop.test(no_heavy_smokers, no_cases)
#' summary(fit)
#' plot(fit)
#' 
#' # MCMC diagnostics (should not be necessary for such a simple model)
#' diagnostics(fit)
#' 
#' # Print out the R code to run the model. This can be copy'n'pasted into
#' # an R-script and further modified.
#' model.code(fit)
#' 
#' 
#'@seealso \code{\link{bayes.binom.test}} for when you want to estimate the 
#'  relative frequency for only one group.
