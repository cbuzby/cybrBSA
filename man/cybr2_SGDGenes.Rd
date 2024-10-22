\name{cybr2_SGDGenes}
\alias{cybr2_SGDGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find genes nearby
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybr2_SGDGenes(peaks, GeneList = SGD_Genes, mywindow = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{
%%     ~~Describe \code{peaks} here~~
Data frame of peaks with columns POS for position and CHROM for chromosome.
}
  \item{GeneList}{
%%     ~~Describe \code{GeneList} here~~
Data frame of genes from SGD containing POS_Start, POS_End, and chr
}
  \item{mywindow}{
%%     ~~Describe \code{mywindow} here~~
Numeric vector of length 1 describing the distance from which a peak can be called
as part of a gene.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Function to find the genes near QTL peaks, based on chromosome and position data.
Can be adapted to any species as long as the chromosomes and positions are in
the same format.
}
\value{
Data frame containing genes of each gene that is near an input peak.
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
To download data, see \link{https://www.alliancegenome.org/bluegenes/alliancemine/lists}

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
function (peaks, GeneList = SGD_Genes, mywindow = 1000)
{
    peaks <- peaks \%>\% select(CHROM, POS)
    filtered_df <- GeneList \%>\% mutate(CHROM = gsub(pattern = "chr",
        replacement = "", x = Chromosome)) \%>\% select(CHROM.x = CHROM,
        POS_Start, POS_End) \%>\% crossing(peaks) \%>\% filter(CHROM.x ==
        CHROM) \%>\% rowwise() \%>\% filter(length(intersect(seq(POS -
        mywindow, POS + mywindow), seq(POS_Start, POS_End))) >
        0) \%>\% select(CHROM, POS, POS_Start, POS_End) \%>\% distinct() \%>\%
        merge(., GeneList) \%>\% arrange(CHROM, POS) \%>\% select(-POS) \%>\%
        distinct()
    rm(peaks)
    filtered_df \%>\% return()
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
