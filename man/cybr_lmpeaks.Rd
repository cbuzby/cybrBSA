\name{cybr_lmpeaks}
\alias{cybr_lmpeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate peaks based on troughs in slope
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cybr_lmpeaks(Data, cutoff = 2, width = 700)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
Dataframe containing columns for: zscore, CSS, CHROM, and POS
}
  \item{cutoff}{
Absolute value of cutoff for which |z-scores| should be excluded; 2 represents
likely false discovery rate of 5\%
}
  \item{width}{
Smoothing width for data before finding peaks}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Data frame of peak locations and |z-scores| above the set threshold.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cassandra Buzby}
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
function (Data, cutoff = 2, width = 700)
{
    require(dplyr)
    CrossoverPoints <- Data \%>\% mutate(abs_zscore = abs(zscore)) \%>\%
        mutate_at(vars(abs_zscore), funs(replace(., . < cutoff,
            -2))) \%>\% arrange(POS) \%>\% group_by(CSS, CHROM) \%>\%
        summarize(POS = POS, abs_zscore = abs_zscore, smooth_abs_z = frollapply(abs_zscore,
            mean, n = 1, align = "center")) \%>\% na.omit() \%>\%
        group_by(CHROM, CSS) \%>\% arrange(POS) \%>\% summarize(POS = POS,
        abs_zscore = abs_zscore, smooth_abs_z = smooth_abs_z,
        slope = frollapply(smooth_abs_z, FUN = slope_change,
            n = width, align = "center")) \%>\% mutate(negative = slope <
        0) \%>\% na.omit() \%>\% summarize(POSc = POS, crossover = frollapply(negative,
        FUN = subtract2, n = 2)) \%>\% filter(crossover == 1)
    tempStartEnd <- data.frame(CHROM = as.factor(as.character(as.roman(1:16)))) \%>\%
        merge(data.frame(CSS = c("VIII", "I"))) \%>\% mutate(POSc = 1,
        crossover = 0) \%>\% rbind(CrossoverPoints) \%>\% group_by(CHROM,
        CSS) \%>\% arrange(POSc) \%>\% mutate(order = paste("A",
        row_number(POSc), sep = "_")) \%>\% select(-crossover) \%>\%
        pivot_wider(names_from = order, values_from = POSc) \%>\%
        pivot_longer(cols = starts_with("A"), names_to = "segment",
            values_to = "value") \%>\% filter(!is.na(value)) \%>\%
        arrange(CHROM, CSS, value) \%>\% group_by(CHROM, CSS) \%>\%
        mutate(End = lead(value, default = Inf)) \%>\% ungroup() \%>\%
        rename(Start = value) \%>\% select(CHROM, CSS, Start, End)
    peakdata <- data.frame(CHROM = NA, CSS = NA, zscore = NA)
    for (i in unique(tempStartEnd$CHROM)) {
        for (c in unique(tempStartEnd$CSS)) {
            newtemp <- tempStartEnd \%>\% filter(CHROM == i, CSS ==
                c)
            for (k in 1:length(newtemp$Start)) {
                something <- Data \%>\% filter(CHROM == i, CSS ==
                  c, label == "Bulk") \%>\% filter(POS > newtemp$Start[k],
                  POS < newtemp$End[k]) \%>\% ungroup() \%>\% group_by(CHROM,
                  CSS) \%>\% summarize(zscore = max(abs(zscore)))
                peakdata <- rbind(peakdata, something)
            }
            rm(newtemp)
        }
    }
    peaks <- Data \%>\% mutate(zscore = abs(zscore)) \%>\% merge(peakdata) \%>\%
        filter(zscore > cutoff)
    return(peaks)
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
