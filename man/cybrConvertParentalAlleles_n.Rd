\name{cybrConvertParentalAlleles_n}
\alias{cybrConvertParentalAlleles_n}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify parental alleles from Oak and Wine VCF Tables
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Converts output table files of Oak and Wine snps into data frame assigning each to the reference or alternate
}
\usage{
cybrConvertParentalAlleles(P1 = Wine, P2 = Oak,
          P1l = "Wine", P2l = "Oak",
          Truncate = TRUE, yeast = TRUE
    )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P1, P2}{Character vector of text files to call; must be a loaded data
  frame. Oak and Wine (defaults) are loaded in with package automatically for
  R versions 3.5 or higher.
%%     ~~Describe \code{x} here~~
}
  \item{P1l}{Character vector of shorthand parental names
%%     ~~Describe \code{x} here~~
}
  \item{P2l}{Character vector of shorthand parental names
%%     ~~Describe \code{x} here~~
}
  \item{Truncate}{Boolean, if TRUE removes positions which do not align to one or the other parental genome; if FALSE keeps all positions. Truncate recommended for downstream functions.
%%     ~~Describe \code{x} here~~
}
  \item{yeast}{Boolean, if TRUE uses the chromosome names for yeast, if FALSE does nothing
%%     ~~Describe \code{x} here~~
}
}
\details{If you are using different parents or have text files for the output tables,
use cybrConvertParentalAlleles_text() instead.
%%  ~~ If necessary, more details than the description above ~~

}
\value{
%%  ~Describe the value returned
\item{Dataframe }{Locations of all SNPs corresponding to included parents}
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
mydatatotest = "HGV.SortedCat.vcf.output.table"

mydf <- cybrInputGATKTable(mydatatotest)

qualitydf <- cybrQualityFilter(mydf)

parentSNPids <- cybrConvertParentalAlleles_n(Truncate = TRUE)

testmerge <- cybrIDAlleles(BSAdfstart = qualitydf, Parentdf = parentSNPids, yeast = TRUE) %>% na.omit()
}

