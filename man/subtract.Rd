\name{subtract}
\alias{subtract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find absolute value of difference between two values, or the max and min values
in a vector.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds span of data by subtracting the maximum value and minimum value of a
vector. If only two values are supplied, the linear distance will be calculated.
}
\usage{
subtract(POS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{POS}{
%%     ~~Describe \code{POS} here~~
Vector of numeric values
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Distinction: if the sign of the difference is important, use subtract2(), which
does not calculate the absolute value and instead finds the difference between
the first two values.
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
Cassandra Buzby
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
subtract2
}
\examples{
subtract(mtcars$mpg)

## The function is currently defined as
function (POS)
{
    return(max(POS) - min(POS))
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
