
library(survival)
library(survcomp)
library(ROCR)

calculate_C <- function(true_survival_time, status, cov, beta, tied = T){
  #true_survival_time: true survival time
  #stauts: 1--> occurrence 0 --> censored
  #cov: predicting attributes/features/covariates 
  #beta: the estimated coefficients associated with different covariates
  status <- as.numeric(status)
  status[status != 0] = 1
  if (length(beta) > 1){
    score = cov %*% beta
  }else{
    score = cov * beta
  }
  
  if(tied){
    C = survConcordance(Surv(true_survival_time, status) ~ score)
    return(C$concordance)
  }else{
    C = concordance.index(x=score, surv.time=true_survival_time, surv.event=status, method="noether")
    return(C$c.index)
  }
}

eval <- function(post, ground.truth, thres = 0.5, alternative = "greater"){
  
  #browser()
  ind <- !is.na(post)
  post <- post[ind]
  ground.truth <- ground.truth[ind]
  
  if(alternative == "greater"){
    p.label <- as.numeric(post > thres)
  }else if(alternative == "less"){
    p.label <- as.numeric(post < thres)
  }
  
  
  TP <- sum(p.label * ground.truth)
  
  TN <- sum((1 - p.label) * (1 - ground.truth))
  
  FP <- sum(p.label * (1 - ground.truth))
  
  FN <- sum((1 - p.label) * ground.truth)
  
  
  #sensitivity recall  TP/conditional POSITIVE
  sen <- sum(TP)/sum(ground.truth)
  
  #specifity: TN/conditional NEGATIVE
  spe <- sum(TN)/sum(1 - ground.truth)
  
  #false discorvery rate: FP/predictive POSITIVE
  fdr <- sum(FP)/sum(p.label)
  
  #precision: TP/predictive POSITIVE = 1 - FDR
  pre <- sum(TP)/sum(p.label)
  
  #accuraccy TP+TN/(conditional POSITIVE + conditioal Negative)
  acc <- (sum(TP) + sum(TN))/(length(ground.truth))
  
  #balanced accuracy
  bacc <- 1/2 * (sum(TP)/sum(ground.truth) + sum(TN)/sum(ground.truth == 0))
  
  
  
  #F1 score  harmonic mean of sensitivity and precision (2TP/(2TP + FN + FP))
  f1 <- 2*sum(TP)/(2*sum(TP) + sum(FN) + sum(FP))
  
 
  
  ##mcc
  mcc1 = sum(TP) * sum(TN) - sum(FP) * sum(FN)
  mcc2 = sqrt((sum(TP) + sum(FP)) * (sum(TP) + sum(FN)) * (sum(TN) + sum(FP)) * (sum(TN) + sum(FN)))
  mcc = mcc1/mcc2
  
  
  ##
  ground_positive = sum(ground.truth)
  ground_negative = sum(ground.truth==0)
  predict_positive = sum(p.label)
  predict_negative = sum(1 - p.label)
  true_positive = TP
  true_negative = TN
  percentage_predict_positve = predict_positive/length(p.label)
  
  retval <- c(sen, spe, fdr, pre, acc, bacc, f1, mcc, 
	         ground_positive, ground_negative, predict_positive, 
			 predict_negative, true_positive, true_negative, 
			 percentage_predict_positve, thres)
			 
  names(retval) <- c("sensitivity", "specifity", "false discovery rate", 
  					"precision", "accuracy", "balanced accuracy", "F1 score","MCC",
	   	            "ground_positive", "ground_negative", "predict_positive", 
	   			    "predict_negative", "true_positive", "true_negative", 
	   			    "percentage_predict_positve", "cutoff")
 
  retval <- data.frame(t(retval))
  
  return(retval)	
}


plot.ROC <- function(post, ground.truth){
  
  pred = prediction(post, ground.truth)
  roc = performance(pred, "tpr", "fpr")
  
  setEPS()
  postscript("ROC.eps")
  plot(roc, lwd=2)
  lines(x=c(0, 1), y=c(0, 1), col="black", lwd=1)
  dev.off()
}

#precision@N
prec_n <- function(post, ground.truth, topN){
  eval.ind <- order(post, decreasing = T)
  post.order <- post[eval.ind]
  
  ground.truth <- ground.truth[eval.ind]
  
  precision <- sapply(topN, function(x){
    sum((ground.truth[1:x] == 1))/x
  })
  
  names(precision) <- as.character(topN)
  return(precision)
}


#recall@N
rcal_n <- function(post, ground.truth, topN){
  eval.ind <- order(post, decreasing = T)
  post.order <- post[eval.ind]
  
  ground.truth <- ground.truth[eval.ind]
  
  recall <- sapply(topN, function(x){
    sum((ground.truth[1:x] == 1))/sum(ground.truth)
  })
  
  names(recall) <- as.character(topN)
  return(recall)
}

#s.precision@N
s.prec_n <- function(post, ground.truth, thres = 0.05, topN){
  ind <- !is.na(post)
  post <- post[ind]
  ground.truth <- ground.truth[ind]
  
  eval.ind <- order(post, decreasing = F)
  post.order <- post[eval.ind]
  ground.truth <- ground.truth[eval.ind]
  
  precision <- sapply(topN, function(x){
    sum(!xor(post.order[1:x] < thres, ground.truth[1:x]))/x
  })
  
  names(precision) <- as.character(topN)
  return(precision)
}


#recall@N
s.rcal_n <- function(post, ground.truth, thres = 0.05, topN){
  
  N <- sum(ground.truth)
  
  ind <- !is.na(post)
  post <- post[ind]
  ground.truth <- ground.truth[ind]
  
  
  eval.ind <- order(post, decreasing = F)
  post.order <- post[eval.ind]
  ground.truth <- ground.truth[eval.ind]
  
  recall <- sapply(topN, function(x){
    sum(!xor(post.order[1:x] < thres, ground.truth[1:x]))/N
  })
  
  names(recall) <- as.character(topN)
  return(recall)
}


#school conformatiy
post.with.school <- function(post, school.id){
  #browser()
  #agg.post <- aggregate(post, list(school.id), mean)
  
  agg.post <- aggregate(post, list(school.id), function(x){
    if(sum(x > 0.5) < sum(x < 0.5)){
      val <- 0
    }
    else{
      val <- 1
    }
  })
  
  
  colnames(agg.post) <- c("id", "post")
  
  
  return(agg.post$post)
}

truth.with.school <- function(ground.truth, school.id){
  agg.truth <- aggregate(ground.truth, list(school.id), function(x){
    if(sum(x > 0.5) < sum(x < 0.5)){
      val <- 0
    }
    else{
      val <- 1
    }
  })
  
  colnames(agg.truth) <- c("id", "truth")
  
  return(agg.truth$truth)
}


ndcg <- function(post, ground.truth, topN){
  
  eval.ind <- order(post, decreasing = T)
  
  
  post.order <- ground.truth[eval.ind]
  
  ground.truth.order <- ground.truth[order(ground.truth, decreasing = T)]
  
  ind <- 2:(length(post.order) + 1)
  
  pred.dg <- (2^post.order - 1)/log2(ind)
  truth.dg <- (2^ground.truth.order - 1)/log2(ind)
  ndcg <- sapply(topN, function(x){
    val.pred <- sum(pred.dg[1:x])
    val.truth <- sum(truth.dg[1:x])
    val <- val.pred/val.truth
    return(val)
  })
  
  names(ndcg) <- as.character(topN)
  
  return(ndcg)
}

school.from.top.com <- function(post, school, topN){
  
  eval.ind <- order(post, decreasing = T)
  
  school.order <- school[eval.ind, ]
  school.order <- school.order[1:topN, ]
  
  val <- unique(school.order[, -1])
  
  val.lst <- list(school.info = val, n = nrow(val), prec = sum(val$label)/nrow(val))
  
  return(val.lst)	
}

school.t.test <- function(post, school){
  #options(digits=22)
  options(error = recover)
  
  
  rslt <- aggregate(post, list(school), function(x){
    #x <- log(x)
    if(length(x) < 2 || sd(x) < 10^(-14)){
      return(NA)
    }else{
      val <- t.test(x, mu = 0.5, alternative = "greater")
      p.value <- val$p.value
      return(p.value)
    }
    
  })
  #browser()
  colnames(rslt) <- c("SchoolID", "p.value")
  return(rslt)
}

# mean average precision @k
#apk
apk <- function(k, actual, predicted)
{
  score <- 0.0
  cnt <- 0.0
  
  
  ind <- order(predicted, decreasing = T)
  predicted <- predicted[ind]
  actual <- actual[ind]
  
  predicted <- as.numeric(predicted > 0.5)
  
  for (i in 1:min(k,length(predicted)))
  {
    if (actual[i] == 1)
    {
      cnt <- cnt + 1
      score <- score + cnt/i 
    }
  }
  if(sum(actual) == 0)
    score <- 0
  else
    score <- score / min(sum(actual), k)
  
  return(score)
}


#mapk
mapk <- function (K, actual, predicted)
{
  rslt <- sapply(K, function(k){
    
    if( length(actual)==0 || length(predicted)==0 ) 
    {
      return(0.0)
    }
    
    #browser()
    scores <- rep(0, length(actual))
    for (i in 1:length(scores))
    {
      scores[i] <- apk(k, actual[[i]], predicted[[i]])
    }
    score <- mean(scores)
    return(score)
  })
  names(rslt) <- as.character(K)
  return(rslt)
}



# mapkNZ
mapkNZ <- function (K, actual, predicted){
  
  rslt <- sapply(K, function(k){
    if( length(actual)==0 || length(predicted)==0 ) 
    {
      return(0.0)
    }
    
    count = 0
    #browser()
    scores <- rep(0, length(actual))
    for (i in 1:length(scores))
    {
      scores[i] <- apk(k, actual[[i]], predicted[[i]])
      if(sum(actual[[i]]) != 0)
        count = count + 1
    }
    score <- sum(scores)/count
    return(score)
  })	
  names(rslt) <- as.character(K)
  return(rslt)
}


#mean reciprocal rank
mrr <- function(actual, predicted){
  
  if( length(actual)==0 || length(predicted)==0 ) 
  {
    return(0.0)
  }
  
  scores <- rep(0, length(actual))
  
  
  
  score <- mapply(function(x,y){
    
    ind <- order(y, decreasing = T)
    x_ord <- x[ind]
    x_ord_ind <- which(x_ord == 1)
    if(length(x_ord_ind) == 0)
      rslt <- 0
    else
      rslt <- 1/x_ord_ind[1]
    
    return(rslt)	
  }, actual, predicted)
  
  retval <- mean(score)
  
  return(retval)
}


mrrNZ <- function(actual, predicted){
  
  if( length(actual)==0 || length(predicted)==0 ) 
  {
    return(0.0)
  }
  
  scores <- rep(0, length(actual))
  
  score <- mapply(function(x,y){
    
    ind <- order(y, decreasing = T)
    y_ord <- x[ind]
    y_ord_ind <- which(y_ord == 1)
    
    if(length(y_ord_ind) == 0)
      rslt <- 0
    else
      rslt <- 1/y_ord_ind[1]
    
    return(rslt)	
  }, actual, predicted)
  
  count_vec <- sapply(actual, function(x){
    rslt <- any(x == 1)
    return(rslt)
  })
  
  retval <- sum(score)/sum(count_vec)
  
  return(retval)
}
library(caTools)

#' Area under the ROC curve
#'
#' improt ROCR
auc_roc <- function(obs, pred) {
  pred <- prediction(pred, obs)
  auc  <- performance(pred, "auc")@y.values[[1]]
  return(auc)
}

#' Area under Precision-recall curve
#'
#' import ROCR
#' import caTools
auc_pr <- function(obs, pred) {
  xx.df <- prediction(pred, obs)
  perf  <- performance(xx.df, "prec", "rec")
  xy    <- data.frame(recall=perf@x.values[[1]], precision=perf@y.values[[1]])
  
  # take out division by 0 for lowest threshold
  xy <- subset(xy, !is.nan(xy$precision))
  
  # Designate recall = 0 as precision = x...arbitrary
  xy <- rbind(c(0, 0), xy)
  #xy <- xy[!(rowSums(xy)==0), ]
  
  res   <- trapz(xy$recall, xy$precision)
  res
}

# Function to create raw data needed to plot Precision against recall
#
# For a vector of observed and predicted, creates x-y coordinates for a ROC
# or PR curve.
rocdf <- function(pred, obs, data=NULL, type=NULL) {
  # plot_type is "roc" or "pr"
  if (!is.null(data)) {
    pred <- eval(substitute(pred), envir=data)
    obs  <- eval(substitute(obs), envir=data)
  }
  
  rocr_xy <- switch(type, roc=c("tpr", "fpr"), pr=c("prec", "rec"))
  rocr_df <- prediction(pred, obs)
  rocr_pr <- performance(rocr_df, rocr_xy[1], rocr_xy[2])
  xy <- data.frame(rocr_pr@x.values[[1]], rocr_pr@y.values[[1]])
  
  # If PR, designate first (undefined) point as recall = 0, precision = x
  if (type=="pr") {
    xy[1, 2] <- 0
    #xy <- xy[!(rowSums(xy)==0), ]
  }
  
  colnames(xy) <- switch(type, roc=c("tpr", "fpr"), pr=c("rec", "prec"))
  return(xy)
}


