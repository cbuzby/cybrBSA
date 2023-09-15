
#Load packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(cowplot)
library(dplyr)
require(circlize)

library(foreach) #no longer used
library(doParallel) #no longer used


################################################################################
# Themes and saved data
################################################################################
ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                       "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                       CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                 "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                 "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))


cybr_BSAcolors = c("#345F6F","#D7335C","#FFB05C",
                   "black",
                   "gray",
                   "#f7dcfa",
                   "#d0eaf5",
                   "#cce3d6",
                   "#fcd7df",
                   "#a32a02")

cybr_CHROMcolors =c("#F26430", "#0A369D", "#7EA3CC","#FF9F1C", "#C14953","#92DCE5", "#8FC93A","#4C4C47",
                      "#F26430", "#0A369D", "#7EA3CC","#FF9F1C", "#C14953","#92DCE5", "#8FC93A","#4C4C47",
                      "#F26430", "#0A369D", "#7EA3CC","#FF9F1C", "#C14953","#92DCE5", "#8FC93A","#4C4C47")


theme_cybr <- function(base_size = 11,
                       base_family = "",
                       base_line_size = base_size/22,
                       base_rect_size = base_size/22){

  #theme_minimal
  theme_bw(base_size = base_size, base_family = base_family,
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(axis.ticks = element_blank(), legend.background = element_blank(),
          legend.key = element_blank(), panel.background = element_blank(),
          panel.border = element_blank(), strip.background = element_blank(),
          plot.background = element_blank(),
          #FROM MY PLOTS
          legend.position = "bottom", #axis.text.x=element_blank(),
          axis.ticks.x=element_blank())}

################################################################################
## FUNCTIONS - New
################################################################################

#Make sure that the parameters to look at (Zscore, etc) are under "coefficient"
cybr_callpeaks_chr3 <- function(dataset,
                                statistic = "Effect",
                                statistic_col = "label",
                                value_col = "GLMResult",
                                exclude = "(Intercept)"){

  #Make dataset
  dataset %>% pivot_wider(names_from = statistic_col, values_from = value_col) %>%
    select(CHROM, POS, coefficient, statistic = statistic) -> piv_glm

  #Set cutoff
  piv_glm %>% filter(CHROM == "III", coefficient != "(Intercept)") %>%
    summarize(X = max(abs(statistic))) -> Cutoff

  #Find peaks
  piv_glm%>%
    filter(coefficient %in% exclude == FALSE) %>%
    filter(abs(statistic) > Cutoff$X) %>%
    group_by(CHROM, coefficient) %>%
    summarize(statistic = max(abs(statistic))) %>%
    right_join(piv_glm) -> output

  colnames(output)[colnames(output) == "statistic"] = statistic

  return(output)
}
################################################################################
cybr_callpeaks <- function(dataset, param = NULL, threshold = NULL){
  dataset %>% mutate(summary = abs(summary)) -> dataset
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
  dataset %>%
    filter(summary > Cutoff$X) %>%
    group_by(CHROM, label) %>%
    summarize(summary = max(summary)) %>%
    ungroup() %>%
    merge(dataset) -> output

  if(is.null(param) == FALSE){
    output %>% filter(label == param) -> output
  }
  return(output)
}
################################################################################
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
  df8int <- rbind(df8int, df1int[df1int$sectors == "VIII",])
  df1int <- rbind(df1int, df8int[df8int$sectors == "I",])

  #Remove the data from those
  df8int$y[df8int$sectors == "VIII"] <- 0
  df1int$y[df1int$sectors == "I"] <- 0

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

################################################################################
## Functions - STILL USED
################################################################################

cybrInputGATKTable <- function(rawData, yeast = TRUE){

  require(dplyr)
  require(doParallel)
  require(foreach)

  HNGLCDRXY <- read.table(rawData, header = TRUE)

  #Identify the unique values besides AD/DP/GQ/PL
  gsub(".AD", "",
       gsub(".GQ", "",
            gsub(".DP","",
                 gsub(".PL","",
                      colnames(select(HNGLCDRXY, -CHROM, -POS, -REF, -ALT)))))) %>% unique() -> Samples
  #i <- Samples[1]

  resultscdf <- foreach(i=Samples,.combine=rbind) %dopar% {
    mydf <- HNGLCDRXY %>% select(CHROM, POS, REF, ALT) %>% mutate(Dataset = i)
    AD <- select(HNGLCDRXY,matches(c(i), ignore.case = FALSE)) %>% select(., contains("AD"))
    GQ <- select(HNGLCDRXY,matches(c(i), ignore.case = FALSE)) %>% select(., contains("GQ"))
    DP <- select(HNGLCDRXY,matches(c(i), ignore.case = FALSE)) %>% select(., contains("DP"))
    PL <- select(HNGLCDRXY,matches(c(i), ignore.case = FALSE)) %>% select(., contains("PL"))
    cbind(mydf, AD , GQ , DP, PL) -> mydftotal
    colnames(mydftotal) <- c(colnames(mydf), "AD", "GQ", "DP", "PL")

    mydftotal %>% separate(AD, c('AD.REF','AD.ALT'), extra='drop') %>%
      separate(PL, c('PL.REF','PL.ALT'), extra='drop') -> mycdf

    mycdf
  }

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                         "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                         CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                   "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                   "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                   "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

    resultscdf %>% left_join(.,ChromKey) %>% select(-CHROM) %>% mutate(CHROM = chromosomes) %>% select(-chromosomes) -> results
  }else{
    results <- resultscdf
  }
  return(results)

}


################################################################################
### Filter by quality

cybrQualityFilter <- function(gatkdf, GQcutoff = 98, cleandata = TRUE){

  #Filter by quality
  gatkdf %>% filter(GQ > GQcutoff) %>%
    select(-DP, -GQ, -PL.ALT, -PL.REF) %>%
    pivot_longer(c(AD.REF, AD.ALT), names_to = "AltRef_Allele", values_to = "ReadCount") %>%
    na.omit() -> filteredgatkdf

  if(cleandata == TRUE){
    #REMOVE POSITIONS WHERE THERE ISN'T ALL BULKS
filteredgatkdf %>% group_by(CHROM, POS) %>%
  summarize(uniqueDatasets = length(unique(Dataset))) -> testingforglm

filteredgatkdf %>% left_join(testingforglm, by = c("CHROM", "POS")) %>%
  filter(uniqueDatasets == max(uniqueDatasets)) %>% select(-uniqueDatasets) -> cleanedgatkdf
}else{
  cleanedgatkdf <- filteredgatkdf
}

return(cleanedgatkdf)
}

################################################################################
### Convert Parental VCFs to Data Frame

cybrConvertParentalAlleles <- function(ParentFiles = c("Wine_VCF.txt", "Oak_VCF.txt"),
                                       Parents = gsub("_VCF.txt","", ParentFiles), Truncate = TRUE, yeast = TRUE){
  temparent <- list()
  mergeparents <- foreach(i=1:length(ParentFiles), .combine=rbind) %dopar% {
    read.table(ParentFiles[i], header = TRUE) %>% mutate(parent = Parents[i])
  }

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                           "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                           CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                     "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                     "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

    rbind(mergeparents) %>% arrange(CHROM, POS) %>%
      select(CHROM, POS, REF, ALT, parent) %>%
      merge(ChromKey) %>% select(-CHROM) %>%
      mutate(CHROM = chromosomes) %>% select(-chromosomes) -> ParentalVCF

  }else{
    rbind(mergeparents) %>% arrange(CHROM, POS) %>%
      select(CHROM, POS, REF, ALT, parent) -> ParentalVCF
  }

  ParentalVCF %>% pivot_wider(names_from = parent, values_from = ALT) -> SNPids

  SNPids$Type <- 0
  for(i in Parents){

    #filter rows in which all values of columns of the parent NOT selected are NA
    select(SNPids,-i, -CHROM, -POS, -REF) -> tempdf
    tempdf$Any_NA <- apply(tempdf, 1, function(x) anyNA(x))
    SNPids$Type[which(tempdf$Any_NA)] <- i
    rm(tempdf)
  }


  #Collect it to output
  if(Truncate == TRUE){
    SNPids %>% select(CHROM, POS,  Type) %>% filter(Type != 0) -> SNPids
  }

  return(SNPids)

}


################################################################################
### Combine Parental and Experimental Variants

cybrIDAlleles <- function(BSAdfstart = finaldf, Parentdf = test, yeast = TRUE){

  Parentdf %>% na.omit()
  BSAdf <- left_join(BSAdfstart, Parentdf)
  BSAdf$PAllele <- NA

  Parents <- unique(Parentdf$Type)
  if(length(Parents) == 2){
    BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == Parents[2]] <- Parents[1]
    BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == Parents[1]] <- Parents[2]

    #Run if only two-parent cross
    BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.ALT" & BSAdf$Type == Parents[2]] <- Parents[2]
    BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.ALT" & BSAdf$Type == Parents[1]] <- Parents[1]

  }else{
    for(i in Parents){
      BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == i] <- i
    }
  }

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                           "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                           CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                     "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                     "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))
    BSAdf$CHROM <- factor(BSAdf$CHROM, levels = ChromKey$chromosomes)
  }


  #Convert all to factors for glm
  BSAdf$Dataset <- factor(BSAdf$Dataset)
  BSAdf$AltRef_Allele <- factor(BSAdf$AltRef_Allele)
  BSAdf$PAllele <- factor(BSAdf$PAllele)

  return(BSAdf)
}

#### Reformat data so that it has bulk etc included




################################################################################
## FUNCTIONS NO LONGER INCLUDED
################################################################################

################################################################################
### Run GLM Function

cybrBSA_GLM <-  function(lrP, chr = "II", windowsize = 100, formula = "PAllele ~ Bulk*Parent",
                         resultscol = c("Intercept", "Bulk", "Parent", "Interaction")){

  lrP <- subset(lrP, CHROM == chr)

  AllResults <- list()
  for(c in unique(lrP$CHROM)){
    lrP <- subset(lrP, CHROM == c)
    #Run the glm on that chromosome
    Results <- foreach (i=unique(lrP$POS), .combine=rbind) %dopar% {
      res <- suppressWarnings(glm(as.formula(formula), weights = ReadCount, family = binomial, data = lrP[lrP$POS == i,]))

      #Output of foreach automatically binds rows of what is printed
      c(c, i, summary(res)$coefficients[((length(summary(res)$coefficients)/2)+1):(length(summary(res)$coefficients) - length(summary(res)$coefficients)/4)])
    }

    resnames <- suppressWarnings(glm(as.formula(formula), weights = ReadCount, family = binomial, data = lrP[lrP$POS == unique(lrP$POS)[1],]))

    #Format the results for the next function
    Results <- as.data.frame(Results)
    colnames(Results) <- c("CHROM", "POS", names(resnames$coefficients))
    #colnames(Results) <- c("CHROM", "POS", resultscol)

    for(i in 2:length(colnames(Results))){
      Results[,i] <- as.numeric(Results[,i])
    }
    Results %>% arrange(POS) -> Results
  }
  return(Results)
}

################################################################################
### New Window Script Function

cybrSmoothBSAWindows <- function(Results, windowsize = 100, chr = unique(Results$CHROM)[1]){

  #extract parameters
  Results %>% select(-CHROM, -POS) %>% names() -> params

  #arrange by position
  Results %>% filter(CHROM == chr) %>% arrange(POS) -> Results

  #run window
  WResult <- foreach(i=windowsize+1:(length(Results$POS) - (2*windowsize)), .combine=rbind) %dopar% {

    smoothedparams <- vector()
    for(p in 1:length(params)){
      smoothedparams[p] <- mean(Results[[params[p]]][(i-windowsize):(i+windowsize)])
    }

    #print CHROM, index, POS, and mean
    unlist(c(chr, i, Results[i,2:(2+length(params))],smoothedparams))
  }

  #rename columns
  WResult <- as.data.frame(WResult)

  colnames(WResult) <- c("CHROM", "Index", "POS",
                         paste(params, "Z", sep = "_"),
                         paste(params, "Zprime", sep = "_"))

  #convert to numeric
  for(i in 2:length(colnames(WResult))){
    WResult[,i] <- as.numeric(WResult[,i])
  }

  return(WResult)
}


cybrSmoothBSAWindows_b <- function(Results, windowsize = 100, chr = unique(Results$CHROM)[1]){

  #extract parameters
  Results %>% select(-CHROM, -POS) %>% names() -> params

  #arrange by position
  Results %>% filter(CHROM == chr) %>% arrange(POS) -> Results

  Results %>% summarise(SNPs = length(unique(POS))) %>%
    distinct() %>%
    mutate(maxW = floor(SNPs/2)) -> tablecounts

  if(windowsize < tablecounts$maxW){
    #run window
    WResult <- foreach(i=windowsize+1:(length(Results$POS) - (2*windowsize)), .combine=rbind) %dopar% {

      smoothedparams <- vector()
      for(p in 1:length(params)){
        smoothedparams[p] <- mean(Results[[params[p]]][(i-windowsize):(i+windowsize)])
      }

      #print CHROM, index, POS, and mean
      unlist(c(chr, i, Results[i,2:(2+length(params))],smoothedparams))
    }

    #rename columns
    WResult <- as.data.frame(WResult)
    colnames(WResult) <- c("CHROM", "Index", "POS",
                           paste(params, "Z", sep = "_"),
                           paste(params, "Zprime", sep = "_"))

    #convert to numeric
    for(i in 2:length(colnames(WResult))){
      WResult[,i] <- as.numeric(WResult[,i])
    }
    return(WResult)
  }else{
    warning(print(paste("Window size of chromosome ", chr, "are larger than number of data points")),
            call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
            domain = NULL)
  }
}

cybrSmoothBSAWindows_Med <- function(Results, windowsize = 100, chr = unique(Results$CHROM)[1]){

  #extract parameters
  Results %>% select(-CHROM, -POS) %>% names() -> params

  #arrange by position
  Results %>% filter(CHROM == chr) %>% arrange(POS) -> Results

  Results %>% summarise(SNPs = length(unique(POS))) %>%
    distinct() %>%
    mutate(maxW = floor(SNPs/2)) -> tablecounts

  if(windowsize < tablecounts$maxW){
    #run window
    WResult <- foreach(i=windowsize+1:(length(Results$POS) - (2*windowsize)), .combine=rbind) %dopar% {

      smoothedparams <- vector()
      for(p in 1:length(params)){
        smoothedparams[p] <- median(Results[[params[p]]][(i-windowsize):(i+windowsize)])
      }

      #print CHROM, index, POS, and mean
      unlist(c(chr, i, Results[i,2:(2+length(params))],smoothedparams))
    }

    #rename columns
    WResult <- as.data.frame(WResult)
    colnames(WResult) <- c("CHROM", "Index", "POS",
                           paste(params, "Z", sep = "_"),
                           paste(params, "Zprime", sep = "_"))

    #convert to numeric
    for(i in 2:length(colnames(WResult))){
      WResult[,i] <- as.numeric(WResult[,i])
    }
    return(WResult)
  }else{
    return(NA)
    warning(paste("Window size of chromosome ", chr, "are larger than number of data points"),
            call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
            domain = NULL)
  }
}

################################################################################
# Check number of loci per chromosome for samples to get window sizes

checkWindowLim <- function(Dataset, includechr = TRUE, exceptchr = NULL){
  if(is.null(exceptchr) == FALSE){
    includechr = FALSE
  }
  if(includechr != TRUE){
    if(is.null(exceptchr)){
      Dataset %>% filter(CHROM %in% includechr) -> Dataset
    }else{
      Dataset %>% filter(CHROM %in% exceptchr == FALSE | CHROM %in% includechr) -> Dataset
    }
  }

  Dataset %>% group_by(CHROM) %>%
    summarise(CHROM = CHROM, SNPs = length(unique(POS))) %>%
    distinct() %>%
    mutate(maxW = floor(SNPs/2)) -> tablecounts

  return(tablecounts)
}

################################################################################
## Make function for plotting

cybrPlotZPrime <- function(zprimedf,
                           columns = c("Bulk_Zprime", "Parent_Zprime", "Interaction_Zprime"),
                           chromosomes = "All",
                           title = "Smoothed Z Scores",
                           yeast = TRUE,
                           colvalues = c("#345F6F", "#D7335C", "#FFB05C")){

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                           "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                           CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                     "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                     "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

    zprimedf$CHROM <- factor(zprimedf$CHROM, levels = ChromKey$chromosomes)
  }


  if(chromosomes == "All"){
    zprimedf %>% filter(CHROM != "M") %>%
      pivot_longer(columns, names_to = "Factor", values_to = "Zprime") %>%
      ggplot(aes(x = POS, y = Zprime, color = Factor)) +

      #Add Bulk, Parent, Day, and Interaction Traces
      geom_line(size = 0.75) +
      scale_color_manual(values = colvalues) +

      #Add the breaks for the facets
      scale_x_continuous(breaks = seq(from = 0, to = max(zprimedf$POS), by = 1e5), name = "Genomic Position") +
      facet_grid(~CHROM, scales = "free_x", space = "free_x") + theme_minimal() +
      theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Genomic Position") + ylab("Zprime")+
      ggtitle(paste("Whole Genome ", title, sep = "")) -> finalplot

  }else{

    zprimedf %>% filter(CHROM %in% chromosomes) %>%
      pivot_longer(columns, names_to = "Factor", values_to = "Zprime") %>%
      ggplot(aes(x = POS, y = Zprime, color = Factor)) +

      #Add Bulk, Parent, Day, and Interaction Traces
      geom_line(size = 0.75) +
      scale_color_manual(values = colvalues) +

      #Add the breaks for the facets
      scale_x_continuous(breaks = seq(from = 0, to = max(zprimedf$POS), by = 1e5), name = "Genomic Position") +
      facet_grid(~CHROM, scales = "free_x", space = "free_x") + theme_minimal() +
      theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Genomic Position") + ylab("Zprime")+
      ggtitle(title) -> finalplot
  }
  return(finalplot)
}

cybrPlotZScore <- function(zprimedf = CuSO4_wholegenomeBSA,
                           column = "Bulk_Z",
                           chromosomes = "All",
                           title = "Z Scores by Position",
                           yeast = TRUE,
                           dotcolor = "black"){

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                           "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                           CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                     "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                     "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

    zprimedf$CHROM <- factor(zprimedf$CHROM, levels = ChromKey$chromosomes)
  }


  if(chromosomes == "All"){

    zprimedf %>% pivot_longer(column, names_to = "Factor", values_to = "Zscore") %>%
      filter(CHROM != "M") %>% filter(Factor == column) %>%
      ggplot(aes(x = POS, y = Zscore)) +
      geom_point(size = 0.5, alpha = 0.08, color = dotcolor) +
      #Add the breaks for the facets
      scale_x_continuous(breaks = seq(from = 0, to = max(zprimedf$POS), by = 1e5), name = "Genomic Position") +
      facet_grid(~CHROM, scales = "free_x", space = "free_x") + theme_minimal() +
      theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Genomic Position") + ylab("Zscore") +
      ggtitle(title) -> finalplot

  }else{

    zprimedf %>% pivot_longer(column, names_to = "Factor", values_to = "Zscore") %>%
      filter(CHROM %in% chromosomes) %>% filter(Factor == column) %>%
      ggplot(aes(x = POS, y = Zscore)) +
      geom_point(size = 1, alpha = 0.08, color = dotcolor) +
      #Add the breaks for the facets
      scale_x_continuous(breaks = seq(from = 0, to = max(zprimedf$POS), by = 1e5), name = "Genomic Position") +
      facet_grid(~CHROM, scales = "free_x", space = "free_x") + theme_minimal() +
      theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Genomic Position") + ylab("Zscore") +
      ggtitle(title) -> finalplot
  }
  return(finalplot)
}

################################################################################
# Adding the loci within a window instead of smoothing

cybrBSA_GLM_window <-  function(lrP, chr = "II", windowsize = 5000, formula = "PAllele~Bulk*Parent"){
  require(stringr)
  if(identical(grep(" ", formula), integer(0)) == FALSE){return(print("Remove spaces from formula or use cybrBSA_GLM()"))}


  lrP <- subset(lrP, CHROM == chr)

  AllResults <- list()
  for(c in unique(lrP$CHROM)){
    lrP <- subset(lrP, CHROM == c)
    #Run the glm on that chromosome
    Results <- foreach (i=unique(lrP$POS), .combine=rbind) %dopar% {
      windowdata <- lrP[lrP$POS < (i+windowsize) & lrP$POS > (i-windowsize),]

      #Convert to sum of counts
      mycols <- unlist(str_split(unlist(str_split(formula, pattern = c("~"))), pattern = "\\*"))
      windowdata %>% group_by(.[mycols]) %>% summarize(ReadCount = sum(ReadCount)) -> windowdata

      res <- suppressWarnings(glm(as.formula(formula),
                                  weights = ReadCount,
                                  family = binomial,
                                  data = windowdata))

      #Output of foreach automatically binds rows of what is printed
      c(c, i, summary(res)$coefficients[((length(summary(res)$coefficients)/2)+1):(length(summary(res)$coefficients) - length(summary(res)$coefficients)/4)])
    }

    resnames <- suppressWarnings(glm(as.formula(formula), weights = ReadCount, family = binomial, data = lrP[lrP$POS == unique(lrP$POS)[1],]))

    #Format the results for the next function
    Results <- as.data.frame(Results)
    colnames(Results) <- c("CHROM", "POS", names(resnames$coefficients))

    for(i in 2:length(colnames(Results))){
      Results[,i] <- as.numeric(Results[,i])
    }
    Results %>% arrange(POS) -> Results
  }
  return(Results)

}

