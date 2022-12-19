\name{cybrBSA_GLM_window}
\alias{cybrBSA_GLM_window}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate z-scores of positional window for BSA}
\description{
Within a window of positions, add all read counts together to calculate effects}
\usage{
cybrBSA_GLM_window(lrP, chr = "II", windowsize = 5000,
                    formula = "PAllele~Bulk*Parent")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lrP}{
Input data frame, stands for logistic regression Parent. Must include all columns within formula}
  \item{chr}{
REQUIRED: Chromosome to separate by}
  \item{windowsize}{
Number of base pairs to calculate within}
  \item{formula}{
Input for logistic regression, WITHOUT spaces}

}
\details{
}
\value{
Data frame containing columns CHROM, POS, and coefficient names
\item{CHROM }{chromosome as character of roman numeral}
\item{POS }{position of SNP at central locus within window}
\item{coefficients}{Z score of the effect listed in logistic regression. There will be one coefficient per factor added and interaction terms for each listed as * in formula}

%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cassandra Buzby}

\note{
This is an alternative to the cybrBSA_GLM() script, which calculates without smoothing for every position. This script should not be used wtih cybrBSA_GLM_Smooth() because it is already an averaging of positions within a window.  Note that this does NOT average by number of SNPs, but by total position, so if there is a SNP desert or the window is smaller than those around, it will calculate only a single position.}

\section{Warning }{If spaces are included in the formula input, the code to extract those variables will break.}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cybrBSA_GLM_window()

## The function is currently defined as
function (lrP, chr = "II", windowsize = 5000, formula = "PAllele~Bulk*Parent",
    resultscol = c("Intercept", "Bulk", "Parent", "Interaction"))
{
    require(stringr)
    if (identical(grep(" ", formula), integer(0)) == FALSE) {
        return(print("Remove spaces from formula or use cybrBSA_GLM()"))
    }
    lrP <- subset(lrP, CHROM == chr)
    AllResults <- list()
    for (c in unique(lrP$CHROM)) {
        lrP <- subset(lrP, CHROM == c)
        Results <- foreach(i = unique(lrP$POS), .combine = rbind) \%dopar\%
            {
                windowdata <- lrP[lrP$POS < (i + windowsize) &
                  lrP$POS > (i - windowsize), ]
                mycols <- unlist(str_split(unlist(str_split(formula,
                  pattern = c("~"))), pattern = "\\*"))
                windowdata <- windowdata \%>\% group_by(.[mycols]) \%>\%
                  summarize(ReadCount = sum(ReadCount))
                res <- suppressWarnings(glm(as.formula(formula),
                  weights = ReadCount, family = binomial, data = windowdata))
                c(c, i, summary(res)$coefficients[((length(summary(res)$coefficients)/2) +
                  1):(length(summary(res)$coefficients) - length(summary(res)$coefficients)/4)])
            }
        resnames <- suppressWarnings(glm(as.formula(formula),
            weights = ReadCount, family = binomial, data = lrP[lrP$POS ==
                unique(lrP$POS)[1], ]))
        Results <- as.data.frame(Results)
        colnames(Results) <- c("CHROM", "POS", names(resnames$coefficients))
        for (i in 2:length(colnames(Results))) {
            Results[, i] <- as.numeric(Results[, i])
        }
        Results <- Results \%>\% arrange(POS)
    }
    return(Results)
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
