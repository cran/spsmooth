########################################################################
#
#   Extensions to mgcv package allowing a projection Slepian basis
#   to be used. 
#
########################################################################


########################################################################
#
#   Smooth Constructor
#
#   Forms model matrix from basis vectors consisting of orthogonal
#   Slepian sequences (DPSSs), using dpss code originally written
#   for package:multitaper.
#
########################################################################
smooth.construct.sp.smooth.spec<-function(object,data,knots) {

  # p.order is inapplicable; we are not using a polynomial
  if (!is.na(object$p.order[1])) warning("Specifying 'p' is meaningless for sp smoothers.")

  # is bandwidth W not specified?
  if (! "W" %in% names(object[['xt']]) ) {
    warning("Bandwidth W not specified as sub-parameter of xt.")
    W <- 7/365  # default for smooth functions of time in GAMs, epidemiology
  } else { # W specified, check to be sure correctly ...
    W <- object[['xt']][['W']]
    if (!is.numeric(W)) stop("Bandwidth W must be numeric.")
    if (W > 0.5 || W < 0) stop("Bandwidth W is strictly bounded on (0,0.5).")
  }

  # use the data
  x <- data[[object$term]]

  # no knots for sp(), so this is the dimension of the projection subspace
  nk <- object$bs.dim 
  if (nk >= 0 & nk<=4) stop("Dimension too small for sp smoother")

  # number of input points; in mgcv, the passed data is post-na.action,
  # so we require a mask to determine the actual structure
  if(!is.null(object[['xt']][['mask']])) {
    mask <- object[['xt']][['mask']]
 
    # sanity check 1: is mask TRUE/FALSE?
    if(length(c(which(mask==TRUE), which(mask==FALSE))) != length(mask) ) {
      stop("'mask' must be populated with TRUE/FALSE elements.")
    }

    # sanity check 2: mask specified, does it match up with object$term?
    if(length(which(mask==TRUE)) != length(data[[object$term]])) {
      stop("Mask must correspond to missing data.")
    }

    nx <- length(mask)
    nxF <- length(x)
  } else {
    # assume the user knows what they are doing ... 
    mask <- NULL
    nx <- length(x)
  }

  # k specified ==> use that number of basis vectors
  if (nk > 4) {
    nw <- nx * W
    dof <- floor(2 * nw) - 1
    if (nk >= dof+2 ) {
      warning("k provided is too high! Manually decreased to 2NW-2.")
      nk <- dof
    }
    if(nk <= (floor(2 * nw) - 4)) {
      warning("k provided is smaller than 2NW-3; this will not provide full
          coverage on the index ends. Are you sure you want to do this?")
    }
  } else { # default to 2NW - 1
    nw <- nx * W
    dof <- floor(2 * nw) - 1
    nk <- dof 
  }

  # augmented with a mean value term ... 
  # (if not augmented, mgcv kills the first column and eliminates the
  # highest-concentration Slepian, i.e. the one with the lowest variance) 
  # See smoothCon() and its calls for more details.
  X <- matrix(data = 0, nrow = nx, ncol = (nk+1))
  X[, 1] <- rep(1, nx)
  X[, 2:(nk+1)] <- .dpss(n = nx, k = nk, nw = nw)$v

  if(!is.null(mask)) {
    object$X <- X[mask==TRUE,] 
    object$size <- c(nxF,nk)
  } else {
    object$X <- X 
    object$size <- c(nx,nk)
  }
  object$rank <- nk  # penalty rank

  # *** not sure this is correct ... is it actually used?
  object$null.space.dim <- nx - nk  # dim. of unpenalized space

  # store "sp" specific stuff ...
  object$mask <- mask
  object$dof <- dof     # maximum DoF (if unconstrained)
  object$df <- nk
  object$bs.dim <- object$size[2]
  object$k <- nk
  object$N <- nx
  object$W <- W
  object$v <- X[,2:(nk+1)]  # save only the tapers

  # centering constraints; set to a positive number to remove the
  # 'zero-meaning' of the basis vectors, which we do not want!
  object$C <- 1

  class(object)<-"sp.smooth"  # Give object a class
  object
}


########################################################################
#
#   Predictor
#
#   Forms prediction (for summary, predict, etc.) from model matrix
#   using provided (sub-set) X. 
#
########################################################################
Predict.matrix.sp.smooth<-function(object,data)
{ 
  # die if freq=TRUE is not set
  x <- data[[object$term]]
  if(length(x) != object$size[1]) {
    stop("sp requires that summary.gam be called with freq=TRUE and p.type=5.")
  }

  # grab objects
  nx <- object$size[1]
  nk <- object$size[2]
  mask <- object[['xt']][['mask']]
  if(is.null(mask)) {
    X <- cbind(rep(1,nx), object$v)
  } else {
    nxF <- length(which(mask==TRUE))
    X <- cbind(rep(1,nxF), object$v[mask==TRUE, ])
  }

  X # return the prediction matrix
}

