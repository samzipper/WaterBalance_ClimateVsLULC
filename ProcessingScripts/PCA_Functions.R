## PCA_Functions.R
# These are functions associated with principle component analysis as described in the
# blog post here: http://www.win-vector.com/blog/2016/05/pcr_part1_xonly/
# 
# They are copied from the markdown document found here:
# https://github.com/WinVector/Examples/blob/master/PCR/XonlyPCA.Rmd

barbell_plot = function(frame, xvar, ymin, ymax, colorvar=NULL) {
  if(is.null(colorvar)) {
    gplot = ggplot(frame, aes_string(x=xvar))
  } else {
    gplot = ggplot(frame, aes_string(x=xvar, color=colorvar))
  }
  
  gplot + geom_point(aes_string(y=ymin)) + 
    geom_point(aes_string(y=ymax)) +
    geom_linerange(aes_string(ymin=ymin, ymax=ymax)) +
    ylab("value")
}

dotplot_identity = function(frame, xvar, yvar, colorvar=NULL) {
  if(is.null(colorvar)) {
    gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar))
  } else {
    gplot = ggplot(frame, 
                   aes_string(x=xvar, y=yvar, ymax=yvar, 
                              color=colorvar))
  }
  gplot + geom_point() + geom_linerange(aes(ymin=0))
}

extractProjection <- function(ndim,princ) {
  # pull off the rotation.  
  proj <- princ$rotation[,1:ndim] 
  # sign was arbitrary, so flip in convenient form
  for(i in seq_len(ndim)) {
    si <- sign(mean(proj[,i]))
    if(si!=0) {
      proj[,i] <- proj[,i]*si
    }
  }
  proj
}

rsq <- function(x,y) {
  1 - sum((y-x)^2)/sum((y-mean(y))^2)
}