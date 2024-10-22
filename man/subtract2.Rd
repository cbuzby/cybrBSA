\name{subtract2}
\alias{subtract2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subtract two values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Subtracts the second value in a vector from the first value in a vector.
}
\usage{
subtract2(POS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{POS}{
Numeric vector containing at least two values
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Only the first two values will be used; subtract2() will ignore all other elements.
}
\value{
Vector of length 1 containing numeric value of the difference
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
If calculating the range of values, use cybrBSA::subtract()
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
subtract2(mtcars$mpg) #no difference in first 2 values
subtract2(wt) #slight difference in first 2 values

## The function is currently defined as
function (POS)
{
    return(POS[1] - POS[2])
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
