\name{findchange}
\alias{findchange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find if there is a difference between first and last values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds IF the difference between first and last values of a vector are identical;
returns 0 if true and 1 if false. Can be used for determining when the slope of
a vector shifts from positive to negative.
}
\usage{
findchange(x)

myvector \%>\% summarize(frollapply(FUN = findchange, x = mycolumn))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Numeric vector of length 2 or greater
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For difference between two values, use subtract() or subtract2()
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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


## The function is currently defined as
function (x)
{
    t <- length(x)
    diff <- x[t] - x[1]
    if (is.na(diff)) {
        return(NA)
    }
    else if (diff == 0) {
        return(0)
    }
    else {
        return(diff > 0)
    }
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
