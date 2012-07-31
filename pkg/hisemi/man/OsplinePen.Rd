\name{OsplinePen      }
\alias{OsplinePen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
O-spline penalty matrix

}
\description{
This function returns the penalty matrix for smoothing spline of any order. 
}
\usage{
OsplinePen(Boundary.knots, knots, ord=1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Boundary.knots}{ A length 2 numeric vector, giving the boundary knot values. }
  \item{knots}{ A nuemric vector of internal knots. }
  \item{ord}{ A numeric integer, which is the order of the derivatives on which squared integral will become the smoothness penalty.}
}
\value{
        A symmetrix penalty matrix. 
}
\references{
Wand, M. P. and Ormerod, J. T. (2008) On semiparametric regression with O'Sullivan penalized splines. Aust. N. Z. J. Stat. 50(2), 179--198.
}
\author{ Long Qu [\email{longorognol@hotmail.com}] }
\note{
When knots are selected by all distinct x values, this returns the penalty matrix of smoothing splines. 

}
\seealso{ 
\code{\link[fda]{bsplinepen}}
}
\examples{
    b.k=c(0,1)
    br=seq(.1,.9,by=.1)
    O1=OsplinePen(b.k, br, 1)
    O2=OsplinePen(b.k, br, 2)
    O3=OsplinePen(b.k, br, 3)
    O4=OsplinePen(b.k, br, 4)
    O5=OsplinePen(b.k, br, 5)
    O6=OsplinePen(b.k, br, 6)

    library(fda)
\dontrun{
    des1=create.bspline.basis(c(0,1),norder=2, breaks=br)
    P1=bsplinepen(des1, 1) # ERROR
    max(abs(P1-O1))
}
    des2=create.bspline.basis(c(0,1),norder=4, breaks=c(b.k[1], br, b.k[2]))
    P2=bsplinepen(des2, 2)
    max(abs(P2-O2))

    des3=create.bspline.basis(c(0,1),norder=6, breaks=c(b.k[1], br, b.k[2]))
    P3=bsplinepen(des3, 3)
    max(abs(P3-O3))

    des4=create.bspline.basis(c(0,1),norder=8, breaks=c(b.k[1], br, b.k[2]))
    P4=bsplinepen(des4, 4, c(0,1))
    max(abs((P4-O4)/(P4+O4)*2),na.rm=TRUE)

    des5=create.bspline.basis(c(0,1),norder=10, breaks=c(b.k[1], br, b.k[2]))
    P5=bsplinepen(des5, 5, c(0,1))
    max(abs((P5-O5)/(P5+O5)*2),na.rm=TRUE)

    des6=create.bspline.basis(c(0,1),norder=12, breaks=c(b.k[1], br, b.k[2]))
    P6=bsplinepen(des6, 6, c(0,1))
    max(abs((P6-O6)/(P6+O6)*2),na.rm=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
