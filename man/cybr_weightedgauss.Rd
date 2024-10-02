\name{cybr_weightedgauss}
\alias{cybr_weightedgauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smoothing by weighted gaussian
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
For use with dplyr's summarize() and frollapply(), to roll allele counts from
short-read sequencing as a weighted gaussian rather than mean.
}
\usage{
cybr_weightedgauss(myx)

tibble %>%
  summarize(frollapply(FUN = cybr_weightedgauss, n = 100, myx))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myx}{
%%     ~~Describe \code{myx} here~~
Numeric vector in order of position, as column of a tibble or as a vector.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
require(dplyr)
require(data.table)

iris %>%
  summarize(frollapply(Sepal.Length, FUN = cybr_weightedgauss, n = 10))

## The function is currently defined as
function (myx)
{
    myy <- dnorm(1:length(myx), mean = length(myx)/2, sd = 10)
    return(weighted.mean(x = myx, y = myy, na.rm = TRUE))
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
