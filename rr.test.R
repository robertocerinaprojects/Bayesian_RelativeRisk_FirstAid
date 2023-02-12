#' Relative Risk test for a contingency table, created for compatibility with 
#' Bayesian First Aid. Much of this documentation and code borrows from that
#' library, accessible here: https://github.com/rasmusab/bayesian_first_aid
#'
#' The main function associated with this repositorty is \code{bayes.rr.test}.
#' 
#' This estimates relative risk of a given outcome against a reference condition 
#' for a given 2-dimensional contingency table. 
#' This is intended as a replacement to the frequentist \code{mosaic::relrisk()}
#' , though it currently does not support an odds-ratio test.
#'
#' Given a contingency table, an outcome of interest, and a reference condition, 
#' \code{bayes.rr.test} estimates:
#'
#' 1) the relative risk (RR) associated with developing the given outcome, 
#' conditional on each possible condition, relative to the reference condition ; 
#' 
#' 2) the conditional risk (CR) of developing the outcome given outcome given 
#' each possible condition ;
#' 
#' 3) the relative frequency of each outcome-condition pair (\eqn{\theta}) ; 
#' 
#' 4) the posterior predictive distribution of each outcome-condition pair 
#' (x_pred).
#' 
#' The following joint model is assumed for the outcome-combination pairs in the 
#' contingency table:
#'
#'\deqn{x_{1,...,J} \sim Multinomial(\theta_{1,...,J},n)}{x_{1,...,J} ~ Multinom
#'ial(\theta,n)}
#'\deqn{\theta_{1,...,J} \sim Dirichlet(\alpha_{1,...,J})}{\tehta_{1,...,J} ~ Di
#'richlet(\alpha_{1,...,J})}
#'\deqn{\alpha_j = 1 \forall j}{\alpha_j = 1 \forall j}
#'
#'
#'Here the prior on the \eqn{\theta}s is a non-informative \eqn{\mathrm{Dirichle
#'t}(1,...,1)}{Dirichlet(1,...,1)}. By \code{plot}ing and looking at a 
#'\code{summary} of the object returned by \code{bayes.rr.test} you can get 
#'more information about the shape of the posterior and the posterior predictive
#'distribution. The JAGS code for this model is available in the object returned
#'by \code{bayes.rr.test}. This can be copy-n-pasted into an R script 
#'and modified, for example, changing the prior on \eqn{\theta}.
#'
#'
#'@param con.table is a two-dimensional contingency table.
#'@param outcome is a character-string corresponding to one of the row names 
#'  of the contingency table, and indicating the outcome of interest.
#'@param reference.condition is a character-string corresponding to one of the
#'  column names of the contingency table, and indicating the reference
#'  condition against which we want to compare other conditions' risk.
#'@param comp.rr a real number against which we want to compare the relative
#'  risk metric. This defaults to 1.
#'@param cred.mass the amount of probability mass that will be contained in 
#'  reported credible intervals. .
#'@param n.iter The number of iterations to run the MCMC sampling.
#'@param progress.bar The type of progress bar. Possible values are "text", 
#'  "gui", and "none".
#'  
#'  
#'@return A list of class \code{bayes_rr_test} that contains information about
#'  the analysis. It can be further inspected using the functions 
#'  \code{summary}, \code{plot} and \code{\link{diagnostics}}.
#'  
#'  
#'@export
bayes.relative.risk.test <- 
  function (con.table,
            comp.rr = 1, 
            cred.mass = 0.95, 
            n.iter=15000, 
            progress.bar="none",
            outcome,
            reference.condition
            ) {
  
  if (length(dim(con.table))<2L | any(dim(con.table)==0)) {
      stop("contingency table should have 2 dimensions (rows and columns)")
  }

  if ((length(cred.mass) != 1L) || is.na(cred.mass) ||
      (cred.mass <= 0) || (cred.mass >= 1)) 
    stop("'cred.mass' must be a single number between 0 and 1")
  
  # format data
    data_list = 
    list(
      x = as.numeric(con.table), 
      n = sum(con.table), 
      J = dim(con.table)[1], 
      E = dim(con.table)[2],
      J_id = which(rownames(con.table)==outcome),
      E_id = which(colnames(con.table)==reference.condition)
    )
    
  # model string written in the JAGS language
  model_string <-
    "model {

      x[1:(J*E)] ~ dmulti(theta[1:(J*E)], n)
      theta ~ ddirch(alpha[1:(J*E)])
      
      for(j in 1:(J*E)){
        alpha[j] <- 1 
      }
      
      x_pred ~ dmulti(theta, n)
      
      for(e in 1:E){
        tau[e] <- sum(theta[((e-1)*J + 1) : ( e*J )]) # exposure margin 
        CR[e] <- theta[ (e -1)*J + J_id ]/tau[e] # conditional risk 
        RR[e] <- CR[e]/CR[E_id] # risk relative to reference condition 
      }
      
   }"
  
  # compile model
  model <-
    jags.model(
      textConnection(model_string), 
      data = data_list, 
      n.chains = 3,
      n.adapt = 0
    )
  
  mcmc_samples <- coda.samples(model,c("theta","x_pred",'CR',"RR"),n.iter = ceiling(n.iter / 3), progress.bar=progress.bar)
  
  temp_comp_val <- comp.rr
  stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = temp_comp_val)
  
  #diff_stats <- mcmc_stats(create_theta_diff_matrix(as.matrix(mcmc_samples)), cred_mass = cred.mass, comp_val = 0)
  #stats <- rbind(stats, diff_stats)
  bfa_object <- list(data_list = data_list, 
                     comp_rr = comp.rr, 
                     cred_mass = cred.mass,
                     mcmc_samples = mcmc_samples, 
                     stats = stats,
                     reference.condition = reference.condition,
                     outcome = outcome,
                     con.table = con.table,
                     jags_model_string = model_string) 
  
  class(bfa_object) <- c("bayes_rr_test", "bayesian_first_aid")
  bfa_object
  }
#' @method plot bayes_rr_test
#' @export
plot.bayes_rr_test <- function(x, ...) {

  samples <- as.matrix(x$mcmc_samples)
 
  # plot estimated proportion for each cell of the contingency table 
  
  theta.samples <- samples[,grepl('theta',colnames(samples))]
  
  par(mfrow = c(3,prod(dim(x$con.table))))

  label = c()  
    for(i in 1:prod(dim(x$con.table)[2])){
      for(j in 1:dim(x$con.table)[1]){
    label = 
      c( label,
      paste(
        'outcome =',rownames(x$con.table)[j],
        '&',
        'condition = ',colnames(x$con.table)[i]
      ) )
  } }
    
    
  for(i in 1:prod(dim(x$con.table))){
    
    main = ''
    if(i == 1){
      main = 'Relative Frequency of Outcome'
    }
    
    plotPost(
      theta.samples[,i], 
      cex.lab = 1.5, 
      main = main, 
      xlab = bquote(theta[ .(label[i]) ]),  
      cred_mass = x$cred_mass, 
      xlim = quantile(theta.samples,c(0,1)),
      col="skyblue" , 
      show_median=TRUE
    )
    
  }
  
  # plot posterior predictive check
  
  x_pred.samples <- samples[,grepl('x_pred',colnames(samples))]
  
  for(i in 1:prod(dim(x$con.table))){
    
    main = ''
    if(i == 1){
      main = 'Data w. Post. Pred.'
    }
      
 hist.ppc <- hist(x_pred.samples[,i],plot = F)
 hist.ppc$p.counts <- hist.ppc$counts/sum(hist.ppc$counts)
 
 plot(hist.ppc,
      cex.lab = 1.5, 
      main = main, 
      xlab = bquote(n[ .(label[i]) ]),  
      xlim = quantile(x_pred.samples,c(0,1)),
      col="skyblue",
      border = NA,
      yaxt = 'n',
      ylab = 'Probability'
    )
 axis(side = 2,
      at = seq(0,max(hist.ppc$counts),length.out = 10),
      labels = round(seq(0,max(hist.ppc$p.counts),length.out = 10),2)
      )
 abline( v = x$data_list$x[i],lwd = 2,col = 'red')

 text(x =  x$data_list$x[i], 
      y = max(hist.ppc$counts),
      labels = x$data_list$x[i],
      cex.lab = 1.5)
 
  }
  
  # plot relative risk 
  alt.conditions = colnames(x$con.table)[-which(colnames(x$con.table)==x$reference.condition)]
  
  rr.id <- grep("^RR\\[",colnames(samples))[which(colnames(x$con.table) %in% alt.conditions)]
  rr.samples <- as.matrix(samples[,rr.id])
  
  for(i in 1:length(alt.conditions)){
    main = ''
    if(i == 1){main = paste("Relative Risk of Outcome =",x$outcome)}
    
    plotPost(
      rr.samples[,i], 
      cex.lab = 1.5, 
      xlab= paste('RR:', alt.conditions[i],'v.',x$reference.condition), 
      main= main,  
      cred_mass = x$cred_mass, 
      col="skyblue" , 
      show_median=TRUE, 
      comp_val=x$comp_rr
    )
    
  }
  
  
}
#' @method print bayes_rr_test
#' @export
print.bayes_rr_test <- function(x, ...) {
  
  s <- format_stats(x$stats)
  
  cat("\n")
  cat("\tBaysian estimation of Relative Risk\n")
  cat("\n")
  cat(paste('contingency table size:',dim(x$con.table)[1],'x',dim(x$con.table)[1], sep=" "))
  cat(paste('\noutcome:',x$outcome))
  cat(paste('\nreference condition:',x$reference.condition))
  cat(paste('\nother conditions:',paste0(colnames(x$con.table)[!colnames(x$con.table) %in% x$reference.condition],collapse = ',')))
  cat("\n")
  print(addmargins(x$con.table))
  cat("\n")

  for(i in which(!colnames(x$con.table) %in% x$reference.condition)){
    cat('Estimated relative risk conditional on', colnames(x$con.table)[i],', relative to',x$reference.condition,
        ':',s[paste('RR[',i,']',sep=''),'mean']
        )
    cat("\n")
    cat(x$cred.mass*100,'% credible interval: (',
        paste0(c(s[paste('RR[',i,']',sep=''),'HDIlo'],s[paste('RR[',i,']',sep=''),'HDIup']),collapse = ' , '),
        ')')
    cat("\n")
    cat('The relative risk of the outcome is larger than',x$comp.rr,'with a probability of',s[paste('RR[',i,']',sep=''),'%>comp'])
    cat("\n")
    cat('The relative risk of the outcome is smaller than',x$comp.rr,'with a probability of',s[paste('RR[',i,']',sep=''),'%<comp'])
  }
  cat("\n")
  cat("\n")

  cat("Estimated relative frequency of outcome [", s[1, "HDI%"] ,"% credible interval]:\n", sep="")
  grep("theta\\[",rownames(s))
  for(param_i in 1:prod(dim(x$con.table))) {
    param <- paste("theta[", param_i, "]", sep="")
    cat(paste("Group ", param_i," (",
        apply(expand.grid(
          paste('outcome:',rownames(x$con.table)),
          paste('condition:',colnames(x$con.table))
                          ),1,paste0,collapse = ' & ')[param_i],
        "): " ,s[param, "median"], " [", paste(s[param, c("HDIlo", "HDIup")], collapse = ", "),"]", sep = "")
        )
    cat("\n")
    invisible(NULL)
  }
}      
#' @method summary bayes_rr_test
#' @export
summary.bayes_rr_test <- function(object, ...) {
  
  s <- round(object$stats, 3)
  
  cat("  Data\n")
  cat("\n")
  cat(paste('contingency table size:',dim(object$con.table)[1],'x',dim(object$con.table)[1], sep=" "))
  cat(paste('\noutcome:',object$outcome))
  cat(paste('\nreference condition:',object$reference.condition))
  cat(paste('\nother conditions:',paste0(colnames(object$con.table)[!colnames(object$con.table) %in% object$reference.condition],collapse = ',')))
  cat("\n")
  print(addmargins(object$con.table))
  cat("\n")
  
  cat("  Model parameters and generated quantities\n")
  cat(paste("RR[i]: the relative risk of the outcome developing relative to the ",object$reference.condition,' condition\n',sep=''))
  cat(paste("CR[i]: the conditional risk of the outcome given a condition j\n"))
  cat("theta[i]: the relative frequency of a given combination of condition j and outcomes i\n")
  cat("x_pred[i]: predicted number of a given combination of condition j and outcomes i\n")
  
  cat("  Measures\n" )
  print(s[, c("mean", "sd", "HDIlo", "HDIup", "%<comp", "%>comp")])
  cat("\n")
  cat("'HDIlo' and 'HDIup' are the limits of a ", s[1, "HDI%"] ,"% HDI credible interval.\n", sep="")
  cat("'%<comp' and '%>comp' are the probabilities of the respective parameter being\n")
  cat("smaller or larger than ", s[1, "comp"] , sep="")
  
  cat("\n")
  cat("  Quantiles\n" )
  print(s[, c("q2.5%", "q25%", "median","q75%", "q97.5%")] )
  invisible(object$stats)
}
#' @method diagnostics bayes_rr_test
#' @export
diagnostics.bayes_prop_test <- function(fit) {
  print_mcmc_info(fit$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(fit$stats, 3))
  cat("\n")
  
  cat("  Model parameters and generated quantities\n")
  cat("RR: The Relative Risk of developing the Outcome\n")
  cat("CR: The Conditional Risk of developing the outcome\n")
  cat("theta: The relative frequency of success\n")
  cat("x_pred: Predicted number of successes in a replication\n")
  old_par <- par( mar=c(3.5,2.5,2.5,0.6) , mgp=c(2.25,0.7,0) )
  plot(fit$mcmc_samples)
  par(old_par)
  invisible(NULL)
}
