\name{cybrPurple}
\alias{cybrPurple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to plot BSA glm result and QTL peaks
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes in the data output of a glm and produces a plot of the traces as well as peaks and specific loci that are known to be adjusted. Peaks can be called using cybr_callpeaks().
}
\usage{
cybrPurple(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{Dataframe in the format of CHROM | POS | Summary | Label, where summary is the statistic and label is the effect (Bulk, Parent, Interaction, etc). Because intercept tends to be on a different scale than the rest of the data, it is removed in newer versions and will not be plotted.
}
  \item{cutoff}{Numeric, default is character "III", indicating that the maximum value on Chr III is the cutoff. Will use a numeric argument instead to draw the cutoff line (cannot be replaced by other characters.)}
  \item{peakslist}{Optional. Dataframe containing columns CHROM and POS, which can be used to plot peaks. If NULL, will not plot peaks.}
  \item{includeFixedGenes}{Logical, default is TRUE. Will plot vertical lines for fixed loci at URA3 and MATalpha}
  \item{peakcolor}{Character hex code. Determines the color of vertical lines representing peaks.}
  \item{mylim}{Optional, Numeric. Y-axis limit; single number, to be used if not all peaks need to be in the frame}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Horizontal ggplot faceted by chromosome, showing a separate trace for up to 4 labels (Bulk, Parent, Interaction, Rep) and vertical lines for peaks and loci of note. Best if plotted in width 18 by height 6 rather than the default.
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
\note{This is very specific to my analyses
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
cybrPurple <- function(dataset,
                       cutoff = "III",
                       peakslist = NULL,
                       includeFixedGenes = TRUE,
                       peakcolor = "#24588A90",
                       mylim = NULL){
  #Make start and end points
  ChromosomeScale <- data.frame(CHROM = factor(as.character(as.roman(1:16)),
                                               levels = as.character(as.roman(1:16))),
                                start = rep(1, 16),
                                end = c(.23,.81,.32, 1.53, .58, .27, 1.09, .56, .44, .75, .67, 1.08, .92, .78, 1.09, .95)*1000000) %>%
    pivot_longer(c(start, end), names_to = "delete", values_to = "POS") %>%
    mutate(Summary = NA, Label = NA) %>% select(-delete)

  #Make Fixed Genes
  FG <- data.frame(Gene = c("Ura3", "MATalphastart", "MATalpha_end"),
                   CHROM = factor(c("V", "III", "III"), levels = as.character(as.roman(1:16))),
                   POS = c(116167, 198671,201177))

  #Make factors for coloring
  dataset$label <- factor(dataset$label, levels = c("Bulk", "Parent", "Interaction", "Rep", "Intercept"))

  dataset %>% select(CHROM, POS, Summary = summary, Label = label) %>%
    rbind(ChromosomeScale) %>%
    ggplot(aes(x = POS, y = abs(Summary), color = Label, linetype = Label %in% c("Parent", "Rep"))) +
    facet_grid(cols = vars(CHROM), scales = "free", space = "free") -> BasePlot

  #Add fixed genes to plot
  if(includeFixedGenes == TRUE){
    BasePlot <- BasePlot + geom_vline(data = FG, aes(xintercept = POS), color = "gray", size = 2, linetype = "dashed")
  }

  #Add peaks to plot
  if(is.null(peakslist) == FALSE){

    peakslist$CHROM <- factor(peakslist$CHROM, levels = as.character(as.roman(1:16)))
    BasePlot <- BasePlot + geom_vline(data = peakslist, aes(xintercept = POS), color = peakcolor, size = 2)
  }


  if(cutoff == "III"){
    FinalPlot <- BasePlot +
      geom_line(size = 1.2) +
      #Cutoff
      geom_hline(aes(yintercept = max(abs(dataset$summary[dataset$CHROM == "III"]))), linetype = "dashed") +
      ylab("") + xlab("")+ ggtitle("")+
      scale_color_manual(values = c("black",  "gray40","#7030A0", "lightpink", "red")) +
      theme(legend.position = "none", axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

  }else if(is.null(cutoff)){
    FinalPlot <- BasePlot +
      geom_line(size = 1.2) +
      ylab("") + xlab("")+ ggtitle("")+
      # facet_grid(~CHROM, scales = "free", space = "free")  +
      scale_color_manual(values = c("black",  "gray40","#7030A0", "lightpink", "red")) +
      theme(legend.position = "none", axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  }else{
    cutoff <- as.numeric(cutoff)
    FinalPlot <- BasePlot +
      geom_line(size = 1.2) +
      geom_hline(aes(yintercept = cutoff)) +
      ylab("") + xlab("")+ ggtitle("")+
      # facet_grid(~CHROM, scales = "free", space = "free")  +
      scale_color_manual(values = c("black",  "gray40","#7030A0", "lightpink", "red")) +
      theme(legend.position = "none", axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  }

  if(is.null(mylim) == FALSE){
    FinalPlot <- FinalPlot + ylim(0, mylim)
  }

  return(FinalPlot)
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
