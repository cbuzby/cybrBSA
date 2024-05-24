\name{cybr_callpeaks}
\alias{cybr_callpeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to identify the peaks in a QTL map, either by the max in a chromosome or all that pass a threshold.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes in the data output of a glm and peaks to produce a plot of the traces as well as peaks from one locus to the chromosome it interacts wtih. Peaks can be called using cybr_callpeaks().
}
\usage{
cybr_callpeaks(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{Dataset of effects, in the format of CHROM | POS | Summary | Label}
  \item{param}{Default NULL; must be a label in the datset to return. Excludes intercept.}
  \item{threshold}{Numeric, default NULL; if present, will be the threshold at which above a peak is called. If NULL, will use the maximum value of Chr 3 as the threshold.}
  \item{include_all}{Logical, default NULL; if TRUE, will return all loci which pass the threshold, rather than just the highest point on a chromosome}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Data frame of peaks which include the chromosome and position of the QTL peak.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Cassandra Buzby
%%  ~~who you are~~
}
\note{This is very specific to my analyses since Chr III can only be used with the MATalpha locus adjustment that we have.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
cybr_callpeaks <- function(dataset, param = NULL, threshold = NULL, include_all = NULL){
  dataset %>% mutate(summary = abs(summary)) %>% filter(label != "intercept",
                                                        label != "Intercept",
                                                        label != "(Intercept)",
                                                        label != "NA") %>%
    na.omit() -> dataset
  #Set cutoff

  Cutoff <- data.frame(X = NA)
  if(is.null(threshold) == FALSE){
    Cutoff$X <- threshold
  }else{
    dataset %>% filter(CHROM == "III") %>%
      ungroup() %>%
      reframe(X = max(summary)) -> Cutoff
  }

  #Find peaks
  if(is.null(include_all)){
    dataset %>%
      filter(summary > Cutoff$X) %>%
      group_by(CHROM, label) %>%
      summarize(summary = max(summary)) %>%
      ungroup() %>%
      merge(dataset) -> output

  }else{
    dataset %>%
      filter(summary > Cutoff$X) -> output
  }

  if(is.null(param) == FALSE){
    output %>% filter(label == param) -> output
  }
  return(output)
}
###############################
