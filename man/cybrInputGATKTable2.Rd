\name{cybrInputGATKTable2}
\alias{cybrInputGATKTable2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert GATK VCF Table to R Dataframe
}
\description{
Converts the gatk output into a data frame, renaming the AD, PL, GQ, and DP fields
for each sample in the file. Will create a column for dataset based on each unique
filename, and will split the AD and PL by alternate and reference alleles.
}
\usage{
cybrInputGATKTable2(rawData, yeast = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawData}{
Output of gatk VariantsToTable, required. See Details for how to create table.
Columns required: CHROM (chromosome), POS (Position), REF (reference allele
frequency), ALT (alternate allele frequency), and PER SAMPLE AD (overall allele
depth), DP (?), GQ (genome quality), PL (?).
}
  \item{yeast}{
If TRUE, converts chromosome names to Roman Numerals; if FALSE, does no
conversion

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Update from cybrInputGATKTable : also includes the number of alleles called for
each position; if more than two, can filter as proxy for mapping quality.

}
\value{
Data farme containing columns CHROM, POS, Dataset (unique names), Allele
(ALT/REF), AD (overall allele depth), DP (?), GQ (genome quality), and PL
(Phred-scaled likelihoods of the genotypes).
}
\references{
%% ~put references to the literature/web site here ~
McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res, 20:1297-303. DOI: 10.1101/gr.107524.110.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (rawData, yeast = TRUE)
{
    require(dplyr)
    require(doParallel)
    require(foreach)
    HNGLCDRXY <- read.table(rawData, header = TRUE)
    Samples <- gsub(".AD", "", gsub(".GQ", "", gsub(".DP", "",
        gsub(".PL", "", colnames(select(HNGLCDRXY, -CHROM, -POS,
            -REF, -ALT)))))) \%>\% unique()
    resultscdf <- foreach(i = Samples, .combine = rbind) \%dopar\%
        {
            mydf <- HNGLCDRXY \%>\% select(CHROM, POS, REF, ALT) \%>\%
                mutate(Dataset = i)
            AD <- select(HNGLCDRXY, matches(c(i), ignore.case = FALSE)) \%>\%
                select(., contains("AD"))
            GQ <- select(HNGLCDRXY, matches(c(i), ignore.case = FALSE)) \%>\%
                select(., contains("GQ"))
            DP <- select(HNGLCDRXY, matches(c(i), ignore.case = FALSE)) \%>\%
                select(., contains("DP"))
            PL <- select(HNGLCDRXY, matches(c(i), ignore.case = FALSE)) \%>\%
                select(., contains("PL"))
            mydftotal <- cbind(mydf, AD, GQ, DP, PL)
            colnames(mydftotal) <- c(colnames(mydf), "AD", "GQ",
                "DP", "PL")
            mycdf <- mydftotal \%>\% separate(AD, c("AD.REF", "AD.ALT"),
                extra = "drop") \%>\% separate(PL, c("PL.REF",
                "PL.ALT"), extra = "drop") \%>\% select(CHROM,
                POS, REF, ALT, Dataset, AD.REF, AD.ALT, GQ, DP,
                PL.REF, PL.ALT)
            doublecdf <- mycdf \%>\% filter(grepl(",", ALT)) \%>\%
                separate(ALT, c("A1", "A2"), extra = "merge") \%>\%
                separate(AD.ALT, c("AD1", "AD2"), extra = "merge") \%>\%
                separate(PL.ALT, c("P1", "P2"), extra = "merge") \%>\%
                pivot_longer(c(A1, A2), names_to = "NumAlt",
                  values_to = "ALT") \%>\% pivot_longer(c(AD1,
                AD2), names_to = "NumADAlt", values_to = "AD.ALT") \%>\%
                pivot_longer(c(P1, P2), names_to = "NumPL", values_to = "PL.ALT") \%>\%
                mutate(NumAlt = gsub("A", "", NumAlt), NumADAlt = gsub("AD",
                  "", NumADAlt), NumPL = gsub("P", "", NumPL)) \%>\%
                filter(NumAlt == NumPL, NumPL == NumADAlt) \%>\%
                select(CHROM, POS, REF, ALT, Dataset, AD.REF,
                  AD.ALT, GQ, DP, PL.REF, PL.ALT)
            triplecdf <- doublecdf \%>\% filter(grepl(",", ALT)) \%>\%
                separate(ALT, c("A1", "A2"), extra = "merge") \%>\%
                separate(AD.ALT, c("AD1", "AD2"), extra = "merge") \%>\%
                separate(PL.ALT, c("P1", "P2"), extra = "merge") \%>\%
                pivot_longer(c(A1, A2), names_to = "NumAlt",
                  values_to = "ALT") \%>\% pivot_longer(c(AD1,
                AD2), names_to = "NumADAlt", values_to = "AD.ALT") \%>\%
                pivot_longer(c(P1, P2), names_to = "NumPL", values_to = "PL.ALT") \%>\%
                mutate(NumAlt = gsub("A", "", NumAlt), NumADAlt = gsub("AD",
                  "", NumADAlt), NumPL = gsub("P", "", NumPL)) \%>\%
                filter(NumAlt == NumPL, NumPL == NumADAlt) \%>\%
                select(CHROM, POS, REF, ALT, Dataset, AD.REF,
                  AD.ALT, GQ, DP, PL.REF, PL.ALT)
            newcdf <- rbind(mycdf, doublecdf, triplecdf)
            newcdf
        }
    if (yeast == TRUE) {
        ChromKey <- data.frame(chromosomes = c("I", "II", "III",
            "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI",
            "XII", "XIII", "XIV", "XV", "XVI", "M"), CHROM = c("NC_001133.9",
            "NC_001134.8", "NC_001135.5", "NC_001136.10", "NC_001137.3",
            "NC_001138.5", "NC_001139.9", "NC_001140.6", "NC_001141.2",
            "NC_001142.9", "NC_001143.9", "NC_001144.5", "NC_001145.3",
            "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))
        results <- resultscdf \%>\% left_join(., ChromKey) \%>\%
            select(-CHROM) \%>\% mutate(CHROM = chromosomes) \%>\%
            select(-chromosomes)
    }
    else {
        results <- resultscdf
    }
    return(results)
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
