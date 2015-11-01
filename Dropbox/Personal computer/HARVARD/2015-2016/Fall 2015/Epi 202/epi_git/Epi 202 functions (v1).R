 
# barfs out all possible information for a person-time 2X2 table
person_time_inf = function(a, b, N1, N0) {
  
  ##### Raw Number of Cases #####
    M1 = a+b
    T = N1+N0  
    exp = (N1*M1)/T
    SE = sqrt( (N1*N0*M1)/(T^2) )
   
    Z = (a-exp)/SE
    p = 2*pnorm( -abs(Z) )  # assume 2-sided
    
    # confidence interval
    lo = a - qnorm(0.025)*SE
    hi = a + qnorm(0.025)*SE
    
    cat("\nX (observed cases) =", a)
    cat("\nExpected cases =", exp)
    cat("\nSE =", SE)
    cat("\nZ =", Z)
    cat("\np =", p)
    cat("\n95% CI:", lo, ",", hi)
  
  ##### IRR #####
    IRR = (a/N1) / (b/N0)
    SE = sqrt( (1/a) + (1/b) )
  
    Z = ( log(IRR) - 1) / SE
    p = 2*pnorm( -abs(Z) )  # assume 2-sided
    
    # confidence interval
    lo = exp( log(IRR) - qnorm(0.025)*SE )
    hi = exp( log(IRR) + qnorm(0.025)*SE )
  
    cat("\n\nIRR =", IRR)
    cat("\nSE of ln(IRR) =", SE)
    cat("\nZ =", Z)
    cat("\np =", p)
    cat("\n95% CI on unlogged scale:", lo, ",", hi)
  
  ##### IRD #####
    IRD = (a/N1) - (b/N0)
    SE = sqrt( (1/a^2) + (1/b^2) )
    
    Z = IRD / SE  # null is 0
    p = 2*pnorm( -abs(Z) )  # assume 2-sided
    
    # confidence interval
    lo = IRD - qnorm(0.025)*SE
    hi = IRD + qnorm(0.025)*SE
    
    cat("\n\nIRD =", IRD)
    cat("\nSE =", SE)
    cat("\nZ =", Z)
    cat("\np =", p)
    cat("\n95% CI:", lo, ",", hi)
}

person_time_inf(a=532, b=374, N1=4.51*3438, N0=4.29*3438)
