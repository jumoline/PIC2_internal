#' Function to calculate the area under the curve
#'
#' \code{oAUCurve} is supposed to calculate the area under the curve.
#'
#' @param x numeric vectors giving the x-coordinates
#' @param y numeric vectors giving the y-coordinates
#' @param method the medthod used to calculate the area under curve. It can be 'trapezoid', 'step', 'linear' (linear interpolation) and 'spline' (spline interpolation)
#' @param subdivisions the maximum number of subintervals
#' @return 
#' a numeric value
#' @note none
#' @export
#' @seealso \code{\link{oAUCurve}}
#' @include oAUCurve.r
#' @examples
#' set.seed(825)
#' y <- abs(rnorm(1000))
#' x <- 1:1000/1000
#' oAUCurve(x,y)

oAUCurve <- function(x, y, method=c("trapezoid","step","linear","spline"), subdivisions=100L)
{
    
    method <- match.arg(method)
    
  	if(length(x) != length(y)){
    	stop("length x must equal length y")
    }
    
  	idx <- order(x)
  	x <- x[idx]
  	y <- y[idx]
    
	if(method=='trapezoid'){
		a <- cbind(y[-length(y)], y[-1])
		b <- x[-1] - x[-length(x)]
		res <- sum(apply(a, 1, mean) * b)
	}else if(method=='step'){
		b <- x[-1] - x[-length(x)]
		res <- sum(y[-length(y)] * b)
	}else if(method=='linear'){
		from <- min(x, na.rm=TRUE)
		to <- max(x, na.rm=TRUE)
		xout <- c(from, to, x[x>from & x<to]) %>% unique() %>% sort()
		values <- stats::approx(x, y, xout=xout)
      	res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
	}else if(method=='spline'){
		lower <- min(x, na.rm=TRUE)
		upper <- max(x, na.rm=TRUE)
		myfunction <- stats::splinefun(x, y, method="natural")
		res <- stats::integrate(myfunction, lower=lower, upper=upper, subdivisions=subdivisions)$value
	}
    
    return(res)
}
