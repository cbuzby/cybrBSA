\name{slope_change}
\alias{slope_change}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Slope change in vector by index
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Find the change in slope of a line with y value as vector values and x as the
index within the vector. Uses lm(y~x) to calculate.
}
\usage{
slope_change(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Vector of numeric values in meaningful order
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Useful for finding QTL peaks based on rolling slopes, when paired with frollapply.

}
\value{
Returns the slope of the line from lm(y~x)
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cassandra Buzby
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
airquality \%>\%
mutate(Day = as.numeric(paste(Month, Day, sep = "."))) \%>\%
  reframe(Day = Day,
          Smooth = frollapply(Temp, n = 5, FUN = slope_change, align = "center"))

## The function is currently defined as
function (x)
{
    df = data.frame(x = x, index = 1:length(x))
    slope = lm(x ~ index, df)$coefficients[2]
    return(slope)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
