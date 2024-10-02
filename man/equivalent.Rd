\name{equivalent}
\alias{equivalent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determines if two values are equivalent.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to find equivalent values using dplyr's piping, and summarizing using
frollapply(). Can be used to find where derivative shifts from negative to
positive.
}
\usage{
equivalent(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A vector of at least 2 values. The first two will be used by the function.
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
Boolean, TRUE if x[1] and x[2] are equal, FALSE if not.
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
  summarize(frollapply(Sepal.Length, FUN = equivalent, n = 2))

## The function is currently defined as
function (x)
{
    x[1] == x[2]
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
