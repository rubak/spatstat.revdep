pval<-function (x, alternative=c("two.sided", "less", "greater")) {
       alternative <- match.arg(alternative)
       if (is.na(x[1])) 
          return(NA)
        nas <- is.na(x[-1])
        if (any(nas)) {
             warning("some values were NA and have been ignored")
             x <- c(x[1], x[-1][!nas])
         }
         if (sum(x) == 0)  return(1)
    mx <- mean(x)
    if (x[1] < mx) 
        psign <- -1
    else psign <- 1
    pplus <- (1 + sum(x[-1] >= x[1]))/(1 + length(x[-1]))
    pmini <- (1 + sum(x[-1] <= x[1]))/(1 + length(x[-1]))
   
    
    switch(alternative,
         less = return(pmini),
         greater = return(pplus),
	 two.sided =  if (pplus>=0.5 & pmini>=0.5) return(1.0*psign) else return(2 * min(pplus, pmini)* psign) 
         )
}
    
