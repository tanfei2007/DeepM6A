makePWM <- function(pwm, alphabet="DNA"){

  if (is.data.frame(pwm)) pwm <- as.matrix(pwm)
  if (!is.matrix(pwm)) stop("pwm must be a matrix or a dataframe")
  
  if (!alphabet %in% c("DNA","AA"))
    stop("alphabet must be either DNA or AA")
  if (alphabet == "DNA" & nrow(pwm) != 4){
    stop("PWM for DNA motifs must have 4 rows")
  }else if (alphabet == "AA" && nrow(pwm) != 21){
    stop("PWM for amino acid motifs must have 21 rows")
  }
  if (any(abs(1 - apply(pwm,2,sum)) > 0.01)){
    print(apply(pwm,2,sum))
    warning("Columns of PWM must add up to 1.0")
  }

  width <- ncol(pwm)
  colnames(pwm) <- 1:width
  rownames(pwm) <- c("A","C","G","T")

  cons <- pwm2cons(pwm)
  ic <- pwm2ic(pwm)
  
  new("pwm", pwm=pwm, consensus=cons, ic=ic, width=width, alphabet=alphabet)
}



## get information content profile from PWM
pwm2ic<-function(pwm) {
    npos<-ncol(pwm)
    ic<-numeric(length=npos)
    for (i in 1:npos) {
        ic[i]<-2 + sum(sapply(pwm[, i], function(x) { 
            if (x > 0) { x*log2(x) } else { 0 }
        }))
    }    
    ic
}

## get consensus sequence from PWM
pwm2cons<-function(pwm) {
    if (class(pwm)!="matrix") {warning("pwm argument must be of class matrix")}
    letters <- c("A", "C", "G", "T")
    paste(apply(pwm, 2, function(x){letters[rev(order(x))[1]]}), collapse="")
}
