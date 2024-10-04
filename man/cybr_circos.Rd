\name{cybr_circos}
\alias{cybr_circos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot interactions between chromosomes I and VIII from BSA experiments.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes in the data output of a glm and peaks to produce a plot of the traces as well as peaks from one locus to the chromosome it interacts wtih. Peaks can be called using cybr_callpeaks().
}
\usage{
cybr_circos(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d1}{Data of interaction peaks with chromosome VIII}
  \item{d8}{Data of interaction peaks with chromosome VIII}
  \item{peaklist1}{Default NULL; if present, will plot the lines from peaks to chromosome I}
  \item{peaklist8}{Default NULL; if present, will plot the lines from peaks to chromosome I}
  \item{maxy}{Default NULL, sets the maximum y axis value for traces of interaction scores}
  \item{color1}{Optional, hex code without opacity. Opacity will be adjusted within for the colors of trace 1 and interactions with Chr I.}
  \item{color8}{Optional, hex code without opacity. Opacity will be adjusted within for the colors of trace 1 and interactions with Chr VIII.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Circular plot faceted by chromosome, showing a separate trace for interaction scores of chromosomes 1 and 8.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
circlize package
%% ~put references to the literature/web site here ~
}
\author{Cassandra Buzby
%%  ~~who you are~~
}
\note{This is very specific to my analyses
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
cybr_circos <- function(d1, d8, peaklist1 = NULL, peaklist8 = NULL, maxy = NULL, color1 = "#7DB0B0", color8 = "#ED7B01"){

  color1fade <- paste(color1, "50", sep = "")
  color8fade <- paste(color8, "50", sep = "")

  #SET UP DATA
  df8 <- data.frame(sectors = as.character(d8$CHROM),
                    x = d8$POS,
                    y = abs(d8$summary),
                    label = d8$label)

  df1 <- data.frame(sectors = as.character(d1$CHROM),
                    x = d1$POS,
                    y = abs(d1$summary),
                    label = d1$label)

  #Take out only interactions - could do this at the prior step anyways?
  df8int <- subset(df8, label == "Interaction")
  df1int <- subset(df1, label == "Interaction")

  #Merge the other chromosomes in, so bind chr1 8 and chr8 1
  # df8int <- rbind(df8int, df1int[df1int$sectors == "VIII",])
  # df1int <- rbind(df1int, df8int[df8int$sectors == "I",])

  #Remove the data from those
  # df8int$y[df8int$sectors == "VIII"] <- 0
  # df1int$y[df1int$sectors == "I"] <- 0

  #REORDER THE CHROMOSOMES
  df8int$sectors <- factor(df8int$sectors, levels = as.character(as.roman(1:16)))
  df1int$sectors <- factor(df1int$sectors, levels = as.character(as.roman(1:16)))

  ##############################################################################
  dfall <- rbind(df1int, df8int)
  circos.par("track.height" = 0.3, start.degree = 90, cell.padding = c(0,0))
  circos.initialize(dfall$sectors, x = dfall$x)

  if(is.null(maxy)){
    #I think this makes the sizes?
    circos.track(ylim = c(max(c(dfall$y)), 0), dfall$sectors, y = dfall$y,
                 panel.fun = function(x, y) {
                   circos.text(CELL_META$xcenter,
                               0 - mm_y(5),
                               CELL_META$sector.index,
                               niceFacing = FALSE)
                   circos.axis(labels.cex = 0.1)
                 })
  }else{
    circos.track(ylim = c(maxy, 0), dfall$sectors, y = dfall$y,
                 panel.fun = function(x, y) {
                   circos.text(CELL_META$xcenter,
                               0 - mm_y(5),
                               CELL_META$sector.index,
                               niceFacing = FALSE)
                   circos.axis(labels.cex = 0.1)
                 })
  }
  #Makes the chromosome overlap parts
  #CHROMOSOME I
  draw.sector(83.5, #RIGHT
              90, #LEFT
              rou1 = 0.99, rou2 = 0.69, clock.wise = FALSE, col = color1fade, border = FALSE)
  #CHROMOSOME VIII
  draw.sector(289.5, #LEFT
              305.4, #RIGHT
              rou1 = 0.99, rou2 = 0.69, clock.wise = FALSE, col = color8fade, border = FALSE)

  #Makes the lines
  circos.trackPoints(df8int$sectors, df8int$x, abs(df8int$y), col = color8, pch = 16, cex = 0.1)
  circos.trackPoints(df1int$sectors, df1int$x, abs(df1int$y), col = color1, pch = 16, cex = 0.1)

  if(is.null(peaklist8) == FALSE){
    for(i in 1:length(peaklist8$POS)){
      circos.link(peaklist8$CHROM[i],
                  peaklist8$POS[i],
                  #Add 8 after
                  "VIII", c(0, max(df8int$x[df8int$sectors == "VIII"])),

                  col = color8fade,
                  h.ratio = 0.3,
                  border = color8fade,
                  lwd = 2)
    }
  }

  if(is.null(peaklist1) == FALSE){
    for(i in 1:length(peaklist1$POS)){
      circos.link("I",c(0, max(df1int$x[df1int$sectors == "I"])),
                  #add 1 first
                  peaklist1$CHROM[i],
                  peaklist1$POS[i],
                  col = color1fade,
                  h.ratio = 0.3,
                  border = color1fade,
                  lwd = 2)
    }
  }

  circos.clear()
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
