print.varComp <-
function(x, ...)
{
	cat("Variance component model fit", '\n')
	cat("\nFixed:", '\n')
	print(coef(x, 'fixed'))
	cat("\nVariance components:", '\n')
	print(coef(x, 'varComp'))
	cat(sprintf("\nNumber of observations: %d\n", nobs(x)))
	invisible(x)
}

summary.varComp = function(object, ...)
{## FIXME: add more summary information
	object$fixef = fixef(object)
	class(object) = c('summary.varComp', class(object))
	object
}
