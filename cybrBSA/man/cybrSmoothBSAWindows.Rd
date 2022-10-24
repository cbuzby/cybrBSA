\name{cybrSmoothBSAWindows}
\alias{cybrSmoothBSAWindows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smooth Z scores from logistic regression
%%  ~~function to do ... ~~
}
\description{Smooth the outputs of cybrBSA_GLM()
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybrSmoothBSAWindows(Results, windowsize = 100, chr = unique(Results$CHROM)[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Results}{Dataframe output of cybrBSA_GLM, which includes Z-scores for each factor of a logistic regression.
%%     ~~Describe \code{x} here~~
}
  \item{windowsize}{Numeric value for number of positions which should be included in each sliding window average of Z scores
%%     ~~Describe \code{x} here~~
}
  \item{chr}{Must be a single chromosome, which should correspond to the only chromosome included in output of cybrBSA_GLM. If smoothing multiple chromosomes, run for each separately.
%%     ~~Describe \code{x} here~~
}
}
\details{
Function will take in all Z scores and output a smoothed score for each.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Dataframe }{Input dataframe with additional columns for smoothed Z-scores by window size}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cassandra Y Buzby
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Run each function individually
GLMoutput <- cybrBSA_GLM(inputdata, "II")
cybrSmoothBSAWindows(GLMoutput, windowsize = 200)


#Running the entire genome in a loop
wholegenomeBSA <- foreach(i=unique(testmerge$CHROM), .combine=rbind) %dopar%{
  cybrSmoothBSAWindows(cybrBSA_GLM(testmerge, chr = i))
}


}
