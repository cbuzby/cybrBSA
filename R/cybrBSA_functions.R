
#Load packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(cowplot)
library(dplyr)
library(circlize)
library(foreach)
library(doParallel)
library(data.table)


################################################################################
# Themes and saved data
################################################################################
ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                       "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                       CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                 "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                 "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

ChromosomeScale <- data.frame(CHROM = factor(as.character(as.roman(1:16)),
                                             levels = as.character(as.roman(1:16))),
                              start = rep(1, 16),
                              end = c(.23,.81,.32, 1.53, .58, .27, 1.09, .56, .44, .75, .67, 1.08, .92, .78, 1.09, .95)*1000000) %>%
  pivot_longer(c(start, end), names_to = "delete", values_to = "POS") %>%
  mutate(Summary = NA, Label = NA) %>% select(-delete)

ChromosomeScale2 <- data.frame(CHROM = factor(as.character(as.roman(1:16)),
                                              levels = as.character(as.roman(1:16))),
                               start = rep(1, 16),
                               end = c(.23,.81,.32, 1.53, .58, .27, 1.09, .56, .44, .75, .67, 1.08, .92, .78, 1.09, .95)*1000000) %>%
  pivot_longer(c(start, end), names_to = "delete", values_to = "POS") %>%
  mutate(summary = 0, label = "Bulk") %>% select(-delete)

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
## FUNCTIONS - Misc
################################################################################
subtract <- function(POS){
  return(max(POS) - min(POS))
}

subtract2 <- function(POS){
  return(POS[1] - POS[2])
}

findchange <- function(x){
  t <- length(x)
  diff <- x[t] - x[1]
  if(is.na(diff)){
    return(NA)
  }else if(diff == 0){
    return(0)
  }else{
    return(diff > 0)
  }
}

findpeak <- function(x){
  t <- length(x)
  diff <- x[t] - x[1]
  return(diff)
}

cybr_weightedgauss <- function(myx){
  myy <- dnorm(1:length(myx), mean = length(myx)/2, sd = 10)
  return(weighted.mean(x = myx, y = myy, na.rm = TRUE))
}

equivalent <- function(x){
  x[1] == x[2]
}

slope_change <- function(x){
  df = data.frame(x = x, index = 1:length(x))
  slope = lm(x ~ index, df)$coefficients[2]
  return(slope)
}

cybr2_SGDGenes <- function(peaks, GeneList = SGD_Genes, mywindow = 1000){

  peaks %>% select(CHROM, POS) -> peaks

  #ChatGPT code adapted for filtering (what is this crossing function?)
  filtered_df <- GeneList %>%
    mutate(CHROM = gsub(pattern = "chr", replacement = "", x = Chromosome)) %>%
    select(CHROM.x = CHROM, POS_Start, POS_End) %>%
    crossing(peaks) %>%                               # Create all combinations of df and numbers
    filter(CHROM.x == CHROM) %>%                      # Match category
    rowwise() %>%
    #filter(POS %in% seq((POS_Start - window), (POS_End + window))) %>%      # Check if number is within the range
    #OR
    filter(length(intersect(seq(POS - mywindow, POS + mywindow), seq(POS_Start, POS_End))) > 0) %>%

    select(CHROM, POS, POS_Start, POS_End) %>%        # Remove redundant columns
    distinct() %>%                                    # In case there are duplicate rows
    merge(., GeneList)%>%
    arrange(CHROM, POS) %>%
    select(-POS) %>% distinct()

  rm(peaks)

  filtered_df %>%
    return()
}

################################################################################
## Data loading and processing
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
      separate(PL, c('PL.REF','PL.ALT'), extra='drop') %>%
      #Added 10/18/23:
      select(CHROM, POS,REF,ALT,Dataset,AD.REF,AD.ALT,GQ,DP,PL.REF,PL.ALT) -> mycdf

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

cybrInputGATKTable2 <- function(rawData, yeast = TRUE){

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
      separate(PL, c('PL.REF','PL.ALT'), extra='drop') %>%
      #Added 10/18/23:
      select(CHROM, POS,REF,ALT,Dataset,AD.REF,AD.ALT,GQ,DP,PL.REF,PL.ALT) -> mycdf

    mycdf %>% filter(grepl(",", ALT)) %>%
      separate(ALT, c("A1", "A2"), extra = 'merge') %>%
      separate(AD.ALT, c("AD1", "AD2"), extra = 'merge') %>%
      separate(PL.ALT, c("P1", "P2"), extra = 'merge') %>%

      pivot_longer(c(A1, A2), names_to = "NumAlt", values_to = "ALT") %>%
      pivot_longer(c(AD1, AD2), names_to = "NumADAlt", values_to = "AD.ALT") %>%
      pivot_longer(c(P1, P2), names_to = "NumPL", values_to = "PL.ALT") %>%
      mutate(NumAlt = gsub("A", "", NumAlt),
             NumADAlt = gsub("AD", "", NumADAlt),
             NumPL = gsub("P", "", NumPL)) %>%
      filter(NumAlt == NumPL,
             NumPL == NumADAlt) %>%
      select(CHROM, POS,REF,ALT,Dataset,AD.REF,AD.ALT,GQ,DP,PL.REF,PL.ALT) -> doublecdf

    doublecdf %>% filter(grepl(",", ALT)) %>%
      separate(ALT, c("A1", "A2"), extra = 'merge') %>%
      separate(AD.ALT, c("AD1", "AD2"), extra = 'merge') %>%
      separate(PL.ALT, c("P1", "P2"), extra = 'merge') %>%

      pivot_longer(c(A1, A2), names_to = "NumAlt", values_to = "ALT") %>%
      pivot_longer(c(AD1, AD2), names_to = "NumADAlt", values_to = "AD.ALT") %>%
      pivot_longer(c(P1, P2), names_to = "NumPL", values_to = "PL.ALT") %>%
      mutate(NumAlt = gsub("A", "", NumAlt),
             NumADAlt = gsub("AD", "", NumADAlt),
             NumPL = gsub("P", "", NumPL)) %>%
      filter(NumAlt == NumPL,
             NumPL == NumADAlt) %>%
      select(CHROM, POS,REF,ALT,Dataset,AD.REF,AD.ALT,GQ,DP,PL.REF,PL.ALT) -> triplecdf

    rbind(mycdf, doublecdf, triplecdf) -> newcdf

    newcdf
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
## GLM Analysis Functions
################################################################################

glmer_cb2_short <- function (..., W, formula, numgroups = FALSE, outputlength = 4,
                             return = c("Z"))
{
  data <- list(...)

  require(lme4)
  if (is.null(W) || is.null(formula)) {
    stop("Weights (W) and formula must be provided")
  }
  glm_formula <- as.formula(formula)
  if (!all(names(data) %in% all.vars(glm_formula))) {
    stop("One or more variables in the formula are not provided as arguments")
  }
  for (i in all.vars(glm_formula)) {
    if (length(unique(as.data.frame(data)[, i])) < 2) {
      output <- rep(NA, outputlength)
      return(output)
    }
  }
  glm_fit <- glmer(glm_formula, data = as.data.frame(data), weights = W,
                   family = binomial)
  if (return %in% "Z") {
    output <- summary(glm_fit)$coefficients[((length(summary(glm_fit)$coefficients) *
                                                0.5) + 1):((length(summary(glm_fit)$coefficients) *
                                                              0.75))]
  }
  if (length(output) == outputlength) {
    return(output)
  }
  else {
    return(rep(NA, outputlength))
  }
}


glm_cb2_all <- function(..., W, formula, numgroups = FALSE, outputlength = 8) {
  data <- list(...)

  #Ensure that there is a formula and W parameter
  if (is.null(W) || is.null(formula)) {
    stop("Weights (W) and formula must be provided")
  }
  #Set formula
  glm_formula <- as.formula(formula)
  #Ensure that formula works for the data provided
  if (!all(names(data) %in% all.vars(glm_formula))) {
    stop("One or more variables in the formula are not provided as arguments")
  }

  #########################
  for(i in all.vars(glm_formula)){
    if(length(unique(as.data.frame(data)[,i])) < 2){
      output <- rep(NA, outputlength)
      #print("Not enough levels within groups")

      return(output)
    }
  }

  glm_fit <- glm(glm_formula, data = as.data.frame(data), weights = W, family = binomial)

  output <- summary(glm_fit)$coefficients[c(((length(summary(glm_fit)$coefficients)*0)+1):
                                              ((length(summary(glm_fit)$coefficients)*0.25)),
                                            ((length(summary(glm_fit)$coefficients)*0.5)+1):
                                              ((length(summary(glm_fit)$coefficients)*0.75)))]


  if(length(output) == outputlength){
    return(output)
  }else{
    return(rep(NA, outputlength))
  }

}

#running glm for Z scores within summarize or reframe()
glm_cb2_short <- function(..., W, formula, numgroups = FALSE, outputlength = 4, return = c("Z")) {
  data <- list(...)

  #Ensure that there is a formula and W parameter
  if (is.null(W) || is.null(formula)) {
    stop("Weights (W) and formula must be provided")
  }
  #Set formula
  glm_formula <- as.formula(formula)
  #Ensure that formula works for the data provided
  if (!all(names(data) %in% all.vars(glm_formula))) {
    stop("One or more variables in the formula are not provided as arguments")
  }

  #########################
  #MAYBEWORKS
  for(i in all.vars(glm_formula)){
    if(length(unique(as.data.frame(data)[,i])) < 2){
      output <- rep(NA, outputlength)
      #print("Not enough levels within groups")

      return(output)
    }
  }

  glm_fit <- glm(glm_formula, data = as.data.frame(data), weights = W, family = binomial)

  #Return specific ones like z-score only
  if(return %in% "Z"){
    output <- summary(glm_fit)$coefficients[((length(summary(glm_fit)$coefficients)*0.5)+1):((length(summary(glm_fit)$coefficients)*0.75))]
  }

  if(length(output) == outputlength){
    return(output)
  }else{
    return(rep(NA, outputlength))
  }

}

################################################################################
## Peak-Calling Analysis Functions
################################################################################

cybr_lmpeaks <- function(Data, cutoff = 2, width = 700){
  require(dplyr)

  Data %>% mutate(abs_zscore = abs(zscore)) %>%
    mutate_at(vars(abs_zscore), funs(replace(., .< cutoff, -2))) %>%
    arrange(POS) %>%
    group_by(CSS, CHROM) %>%
    summarize(POS = POS,
              abs_zscore = abs_zscore,
              smooth_abs_z = frollapply(abs_zscore, mean, n = 1, align = "center")) %>%
    na.omit() %>%
    group_by(CHROM, CSS) %>% arrange(POS) %>%
    summarize(POS = POS, abs_zscore = abs_zscore, smooth_abs_z = smooth_abs_z,
              slope = frollapply(smooth_abs_z, FUN = slope_change, n = width, align = "center")) %>%
    mutate(negative = slope < 0) %>% na.omit() %>%
    summarize(POSc = POS, crossover = frollapply(negative, FUN = subtract2, n = 2)) %>%
    filter(crossover == 1) -> CrossoverPoints

  data.frame(CHROM = as.factor(as.character(as.roman(1:16)))) %>%
    merge(data.frame(CSS = c("VIII", "I"))) %>%
    mutate(POSc = 1, crossover = 0) %>%
    rbind(CrossoverPoints) %>%
    group_by(CHROM, CSS) %>%
    arrange(POSc) %>%
    mutate(order = paste("A", row_number(POSc), sep = "_")) %>%
    select(-crossover) %>%
    pivot_wider(names_from = order, values_from = POSc) %>%
    pivot_longer(cols = starts_with("A"), names_to = "segment", values_to = "value") %>%
    filter(!is.na(value)) %>%
    arrange(CHROM, CSS, value) %>%
    group_by(CHROM, CSS) %>%
    mutate(End = lead(value, default = Inf)) %>%
    ungroup() %>%
    rename(Start = value) %>%
    select(CHROM, CSS, Start, End) -> tempStartEnd

  peakdata <- data.frame(CHROM = NA, CSS = NA, zscore = NA)

  for(i in unique(tempStartEnd$CHROM)){
    for(c in unique(tempStartEnd$CSS)){
      tempStartEnd %>% filter(CHROM == i, CSS == c) -> newtemp

      # newtemp$End[is.na(newtemp$End)] <- Inf
      # newtemp$End[is.na(newtemp$Start)] <- Inf

      for(k in 1:length(newtemp$Start)){
        Data %>% filter(CHROM == i, CSS == c, label == "Bulk") %>%
          filter(POS > newtemp$Start[k], POS < newtemp$End[k]) %>%
          ungroup() %>%
          group_by(CHROM, CSS) %>%
          summarize(zscore = max(abs(zscore))) -> something

        peakdata <- rbind(peakdata, something)
      }

      rm(newtemp)

    }
  }
  Data %>%
    mutate(zscore = abs(zscore)) %>% merge(peakdata) %>%
    filter(zscore > cutoff) -> peaks

  return(peaks)

}

################################################################################
## Visualization Functions
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

################################################################################
## OBSOLETE
################################################################################

# #MULTIPLE PEAKS
# cybr_multpeaks <- function(dataset,
#                            param = NULL,
#                            threshold = NULL,
#                            width = 1000,
#                            buffer = 200){
#   dataset %>% mutate(summary = abs(summary)) %>% filter(label != "intercept",
#                                                         label != "Intercept",
#                                                         label != "(Intercept)",
#                                                         label != "NA") %>%
#     mutate(Window = POS %/% width) %>%
#     mutate(Check = 1) %>%
#     na.omit() -> dataset
#
#
#   #Set cutoff
#
#   Cutoff <- data.frame(X = NA)
#   if(is.null(threshold) == FALSE){
#     Cutoff$X <- threshold
#   }else{
#     dataset %>% filter(CHROM == "III") %>%
#       ungroup() %>%
#       reframe(X = max(summary)) -> Cutoff
#   }
#
#   #Select only the parameter in question
#   if(is.null(param) == FALSE){
#     dataset %>% filter(label == param) -> dataset
#   }
#
#   for(i in buffer:(nrow(dataset)-buffer)){
#     if(dataset$Window[i] == dataset$Window[(i-(buffer - 1))]){
#       dataset$Check[(i-(buffer - 1)):(i+(buffer - 1))] <- FALSE
#     }else{
#       dataset$Check[(i-(buffer - 1)):((buffer - 1))] <- TRUE
#     }
#   }
#
#   #Find peaks
#   dataset %>%
#     filter(summary > Cutoff$X) %>%
#     group_by(CHROM, label, Window) %>%
#     summarize(summary = max(summary),
#               Check = Check) %>%
#     ungroup() %>%
#     merge(dataset) %>%
#     distinct() -> output
#
#
#   return(output)
# }
#
# #Make sure that the parameters to look at (Zscore, etc) are under "coefficient"
# cybr_callpeaks_chr3 <- function(dataset,
#                                 statistic = "Effect",
#                                 statistic_col = "label",
#                                 value_col = "GLMResult",
#                                 exclude = "(Intercept)"){
#
#   #Make dataset
#   dataset %>% pivot_wider(names_from = statistic_col, values_from = value_col) %>%
#     select(CHROM, POS, coefficient, statistic = statistic) -> piv_glm
#
#   #Set cutoff
#   piv_glm %>% filter(CHROM == "III", coefficient != "(Intercept)") %>%
#     summarize(X = max(abs(statistic))) -> Cutoff
#
#   #Find peaks
#   piv_glm%>%
#     filter(coefficient %in% exclude == FALSE) %>%
#     filter(abs(statistic) > Cutoff$X) %>%
#     group_by(CHROM, coefficient) %>%
#     summarize(statistic = max(abs(statistic))) %>%
#     right_join(piv_glm) -> output
#
#   colnames(output)[colnames(output) == "statistic"] = statistic
#
#   return(output)
# }
#
# #New version:
# cybr_callpeaks <- function(dataset, param = NULL, threshold = NULL, include_all = NULL){
#   dataset %>% mutate(summary = abs(summary)) %>% filter(label != "intercept",
#                                                         label != "Intercept",
#                                                         label != "(Intercept)",
#                                                         label != "NA") %>%
#     na.omit() -> dataset
#   #Set cutoff
#
#   Cutoff <- data.frame(X = NA)
#   if(is.null(threshold) == FALSE){
#     Cutoff$X <- threshold
#   }else{
#     dataset %>% filter(CHROM == "III") %>%
#       ungroup() %>%
#       reframe(X = max(summary)) -> Cutoff
#   }
#
#   #Find peaks
#   if(is.null(include_all)){
#     dataset %>%
#       filter(summary > Cutoff$X) %>%
#       group_by(CHROM, label) %>%
#       summarize(summary = max(summary)) %>%
#       ungroup() %>%
#       merge(dataset) -> output
#
#   }else{
#     dataset %>%
#       filter(summary > Cutoff$X) -> output
#   }
#
#   if(is.null(param) == FALSE){
#     output %>% filter(label == param) -> output
#   }
#   return(output)
# }
#
# #Plotting
# cybrPurple <- function(dataset,
#                        cutoff = "III",
#                        peakslist = NULL,
#                        includeFixedGenes = TRUE,
#                        peakcolor = "#24588A90",
#                        mylim = NULL){
#   #Make start and end points
#   ChromosomeScale <- data.frame(CHROM = factor(as.character(as.roman(1:16)),
#                                                levels = as.character(as.roman(1:16))),
#                                 start = rep(1, 16),
#                                 end = c(.23,.81,.32, 1.53, .58, .27, 1.09, .56, .44, .75, .67, 1.08, .92, .78, 1.09, .95)*1000000) %>%
#     pivot_longer(c(start, end), names_to = "delete", values_to = "POS") %>%
#     mutate(Summary = NA, Label = NA) %>% select(-delete)
#
#   #Make Fixed Genes
#   FG <- data.frame(Gene = c("Ura3", "MATalphastart", "MATalpha_end"),
#                    CHROM = factor(c("V", "III", "III"), levels = as.character(as.roman(1:16))),
#                    POS = c(116167, 198671,201177))
#
#   #Make factors for coloring
#   dataset$label <- factor(dataset$label, levels = c("Bulk", "Parent", "Interaction", "Rep", "Intercept"))
#
#   dataset %>% select(CHROM, POS, Summary = summary, Label = label) %>%
#     rbind(ChromosomeScale) %>%
#     ggplot(aes(x = POS, y = abs(Summary), color = Label, linetype = Label %in% c("Parent", "Rep"))) +
#     facet_grid(cols = vars(CHROM), scales = "free", space = "free") -> BasePlot
#
#   #Add fixed genes to plot
#   if(includeFixedGenes == TRUE){
#     BasePlot <- BasePlot + geom_vline(data = FG, aes(xintercept = POS), color = "gray", size = 2, linetype = "dashed")
#   }
#
#   #Add peaks to plot
#   if(is.null(peakslist) == FALSE){
#
#     peakslist$CHROM <- factor(peakslist$CHROM, levels = as.character(as.roman(1:16)))
#     BasePlot <- BasePlot + geom_vline(data = peakslist, aes(xintercept = POS), color = peakcolor, size = 2)
#   }
#
#
#   if(cutoff == "III"){
#     FinalPlot <- BasePlot +
#       geom_line(size = 1.2) +
#       #Cutoff
#       geom_hline(aes(yintercept = max(abs(dataset$summary[dataset$CHROM == "III"]))), linetype = "dashed") +
#       ylab("") + xlab("")+ ggtitle("")+
#       scale_color_manual(values = c("black",  "gray40","#7030A0", "lightpink", "red")) +
#       theme(legend.position = "none", axis.text.x=element_blank(),
#             axis.ticks.x=element_blank())
#
#   }else if(is.null(cutoff)){
#     FinalPlot <- BasePlot +
#       geom_line(size = 1.2) +
#       ylab("") + xlab("")+ ggtitle("")+
#       # facet_grid(~CHROM, scales = "free", space = "free")  +
#       scale_color_manual(values = c("black",  "gray40","#7030A0", "lightpink", "red")) +
#       theme(legend.position = "none", axis.text.x=element_blank(),
#             axis.ticks.x=element_blank())
#   }else{
#     cutoff <- as.numeric(cutoff)
#     FinalPlot <- BasePlot +
#       geom_line(size = 1.2) +
#       geom_hline(aes(yintercept = cutoff)) +
#       ylab("") + xlab("")+ ggtitle("")+
#       # facet_grid(~CHROM, scales = "free", space = "free")  +
#       scale_color_manual(values = c("black",  "gray40","#7030A0", "lightpink", "red")) +
#       theme(legend.position = "none", axis.text.x=element_blank(),
#             axis.ticks.x=element_blank())
#   }
#
#   if(is.null(mylim) == FALSE){
#     FinalPlot <- FinalPlot + ylim(0, mylim)
#   }
#
#   return(FinalPlot)
# }
#
#
# ### Combine Parental and Experimental Variants
# cybrIDAlleles <- function(BSAdfstart = finaldf, Parentdf = test, yeast = TRUE){
#
#   Parentdf %>% na.omit()
#   BSAdf <- left_join(BSAdfstart, Parentdf)
#   BSAdf$PAllele <- NA
#
#   Parents <- unique(Parentdf$Type)
#   if(length(Parents) == 2){
#     BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == Parents[2]] <- Parents[1]
#     BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == Parents[1]] <- Parents[2]
#
#     #Run if only two-parent cross
#     BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.ALT" & BSAdf$Type == Parents[2]] <- Parents[2]
#     BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.ALT" & BSAdf$Type == Parents[1]] <- Parents[1]
#
#   }else{
#     for(i in Parents){
#       BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == i] <- i
#     }
#   }
#
#   if(yeast == TRUE){
#     ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
#                                            "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
#                            CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
#                                      "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
#                                      "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
#                                      "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))
#     BSAdf$CHROM <- factor(BSAdf$CHROM, levels = ChromKey$chromosomes)
#   }
#
#
#   #Convert all to factors for glm
#   BSAdf$Dataset <- factor(BSAdf$Dataset)
#   BSAdf$AltRef_Allele <- factor(BSAdf$AltRef_Allele)
#   BSAdf$PAllele <- factor(BSAdf$PAllele)
#
#   return(BSAdf)
# }
#
# cybr2_rollmean <- function(dataframe){
#   dataframe %>% pivot_longer(c(-CHROM, -POS), names_to = "label") %>% group_by(CHROM, label) %>% arrange(POS) %>%
#     summarize(POS = POS, CHROM = CHROM, SmoothCount = ceiling(frollmean(value, n = 100))) %>% na.omit() %>%
#     pivot_wider(names_from = label,values_from = SmoothCount)
# }
#
# # No replicates, fix Parent
# cybrPermute_byParent <- function(dataset){
#   start.time <- Sys.time()
#
#   print("Make sure that dilute bulk is labeled D")
#
#   dataset %>%
#     distinct() %>% ungroup() %>%
#     group_by(CHROM, POS, Allele,
#              Bulk,
#              #Rep, #might not have these
#              Parent) %>%
#     summarize(culprits = length((SmoothCount))) %>%
#     merge(dataset) %>%
#     filter(culprits == 1) %>%
#     ungroup() %>%
#     distinct() %>% #THIS IS IMPORTANT
#     pivot_wider(names_from = Allele, values_from = SmoothCount) -> newnewtest
#
#   #these are now all of the ones that can be used to permute
#   newnewtest %>% filter(Bulk == "D",
#                         CHROM %in% c("I", "III", "V", "VIII", "M") == FALSE) %>%
#     select(CHROM, POS, Parent, Oak, Wine) %>%
#     group_by(Parent) %>% mutate(Loc = paste(CHROM, POS)) %>%
#     summarize(Parent = Parent, Oak = Oak, Wine = Wine,
#               Loc = sample(Loc)) %>%
#     mutate(Bulk = "A") -> shuffled_DiluteA2
#
#   newnewtest %>% filter(Bulk == "D",
#                         CHROM %in% c("I", "III", "V", "VIII", "M") == FALSE) %>%
#     select(CHROM, POS, Parent, Oak, Wine) %>%
#     group_by(Parent) %>% mutate(Loc = paste(CHROM, POS)) %>%
#     summarize(Parent = Parent, Oak = Oak, Wine = Wine,
#               Loc = sample(Loc)) %>%
#     mutate(Bulk = "B") -> shuffled_DiluteB2
#
#   rbind(shuffled_DiluteA2, shuffled_DiluteB2) %>% pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "SmoothCount") -> shuffletoglm2
#
#   #Trying this again
#
#   shuffletoglm2 %>% na.omit() %>%
#     #Original Script
#     group_by(Loc) %>%
#     mutate_if(is.character, as.factor) %>%
#     summarize(Summary = glm_cb2_short(Allele = Allele,
#                                       Bulk = Bulk,
#                                       Parent = Parent,
#                                       #Rep = Rep,
#                                       W = SmoothCount,
#                                       formula = "Allele ~ Bulk * Parent",
#                                       outputlength = 4),
#               Factor = (c("Intercept", "Bulk", "Parent", "Interaction"))) -> glmresult
#   end.time = Sys.time()
#   print(end.time - start.time)
#   return(glmresult)
# }
#
# # No replicates, fix Parent
# cybrPermute_cb2_all <- function(dataset,
#                                 R = 0, perp = 1, outlength = 4,
#                                 glmform = "Allele ~ Bulk * Parent",
#                                 inputlabels = c("Intercept", "Bulk", "Parent", "Interaction")){
#   start.time <- Sys.time()
#
#   print("Make sure that dilute bulk is labeled aD")
#   print(paste("Your labels are:", inputlabels))
#
#   if(R > 0){ #INCLUDE REPLICATES
#     dataset %>%
#       distinct() %>% ungroup() %>%
#       group_by(CHROM, POS, Allele, Bulk, Rep, Parent) %>%
#       summarize(culprits = length((SmoothCount))) %>%
#       merge(dataset) %>%
#       filter(culprits == 1) %>%
#       ungroup() %>%
#       distinct() %>% #THIS IS IMPORTANT
#       pivot_wider(names_from = Allele, values_from = SmoothCount) -> newnewtest
#   }else{ #DO NOT INCLUDE REPLICATES
#     dataset %>%
#       distinct() %>% ungroup() %>%
#       group_by(CHROM, POS, Allele, Bulk, Parent) %>%
#       summarize(culprits = length((SmoothCount))) %>%
#       merge(dataset) %>%
#       filter(culprits == perp) %>%
#       ungroup() %>%
#       distinct() %>% #THIS IS IMPORTANT
#       pivot_wider(names_from = Allele, values_from = SmoothCount) -> newnewtest
#   }
#
#   #PERMUTE TWICE
#   if(R > 0){
#     newnewtest %>% filter(Bulk == "aD",
#                           CHROM %in% c("I", "III", "V", "VIII", "M") == FALSE) %>%
#       select(CHROM, POS, Rep,
#              Parent, Oak, Wine) %>%
#       group_by(Parent) %>% mutate(Loc = paste(CHROM, POS)) %>%
#       summarize(Parent = Parent, Rep = Rep,
#                 Oak = Oak, Wine = Wine, Loc = sample(Loc)) %>%
#       mutate(Bulk = "A") -> shuffled_DiluteA2
#
#     newnewtest %>% filter(Bulk == "aD",
#                           CHROM %in% c("I", "III", "V", "VIII", "M") == FALSE) %>%
#       select(CHROM, POS,
#              Rep, Parent, Oak, Wine) %>%
#       group_by(Parent) %>% mutate(Loc = paste(CHROM, POS)) %>%
#       summarize(Parent = Parent, Rep = Rep,
#                 Oak = Oak, Wine = Wine, Loc = sample(Loc)) %>%
#       mutate(Bulk = "B") -> shuffled_DiluteB2
#
#   }else{
#     newnewtest %>% filter(Bulk == "aD",
#                           CHROM %in% c("I", "III", "V", "VIII", "M") == FALSE) %>%
#       select(CHROM, POS, Parent, Oak, Wine) %>%
#       group_by(Parent) %>% mutate(Loc = paste(CHROM, POS)) %>%
#       summarize(Parent = Parent, Oak = Oak, Wine = Wine, Loc = sample(Loc)) %>%
#       mutate(Bulk = "A") -> shuffled_DiluteA2
#
#     newnewtest %>% filter(Bulk == "aD",
#                           CHROM %in% c("I", "III", "V", "VIII", "M") == FALSE) %>%
#       select(CHROM, POS, Parent, Oak, Wine) %>%
#       group_by(Parent) %>% mutate(Loc = paste(CHROM, POS)) %>%
#       summarize(Parent = Parent, Oak = Oak, Wine = Wine, Loc = sample(Loc)) %>%
#       mutate(Bulk = "B") -> shuffled_DiluteB2
#
#   }
#
#   rbind(shuffled_DiluteA2, shuffled_DiluteB2) %>% pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "SmoothCount") -> shuffletoglm2
#
#   #RUN THE GLM
#   if(R > 0) {
#     shuffletoglm2 %>% na.omit() %>%
#       group_by(Loc) %>%
#       mutate_if(is.character, as.factor) %>%
#       summarize(Summary = glm_cb2_all(Allele = Allele,
#                                       Bulk = Bulk,
#                                       Parent = Parent,
#                                       Rep = Rep,
#                                       W = SmoothCount,
#                                       formula = glmform,
#                                       outputlength = length(inputlabels)*2),
#                 Factor = rep(inputlabels, 2),
#                 d = c(rep("Effect", length(inputlabels)),
#                       rep("Z", length(inputlabels)))) -> glmresult
#   }else{
#     shuffletoglm2 %>% na.omit() %>%
#       #Original Script
#       group_by(Loc) %>%
#       mutate_if(is.character, as.factor) %>%
#       summarize(Summary = glm_cb2_all(Allele = Allele,
#                                       Bulk = Bulk,
#                                       Parent = Parent,
#                                       W = SmoothCount,
#                                       formula = glmform,
#                                       outputlength = length(inputlabels)*2),
#                 Factor = rep(inputlabels, 2),
#                 d = c(rep("Effect", length(inputlabels)),
#                       rep("Z", length(inputlabels)))) -> glmresult
#   }
#
#
#   end.time = Sys.time()
#   print(end.time - start.time)
#   return(glmresult)
# }
