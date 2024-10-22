\name{findpeak}
\alias{findpeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Difference between first and last values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds the difference between the first and last values of a numeric vector. This
differs from subtract(), subtract2(), and findchange() in that it calculates
based on index rather than min/max (as in subtract) or first two values (as in
subtract2).
}
\usage{
findpeak(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Numeric vector
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function finds the difference of last - first values; if order is to be
reversed, multiply by -1.
}
\value{
Numeric vector of length 1
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x)
{
    t <- length(x)
    diff <- x[t] - x[1]
    return(diff)
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
