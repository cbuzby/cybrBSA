# cybrBSA
Package for analyzing sequencing data for multilevel Bulk Segregant Analysis experiment

### Convert VCF to GATK table

Script is in bash

```{bash, eval = FALSE}
module load gatk/4.2.0.0

gatk VariantsToTable \
     -V ${myfile} \
     -F CHROM -F POS -F REF -F ALT \
     -GF AD -GF DP -GF GQ -GF PL \
     -O ${myfile}.output.table

```
### Running Pipeline of Functions

```{r, eval = FALSE, warning=FALSE, message=FALSE}
mydatatotest = "HGV.SortedCat.vcf.output.table"

mydf <- cybrInputGATKTable(mydatatotest)
head(mydf)

qualitydf <- cybrQualityFilter(mydf)
head(qualitydf)

parentSNPids <- cybrConvertParentalAlleles(Truncate = TRUE)
head(parentSNPids)

testmerge <- cybrIDAlleles(BSAdfstart = qualitydf, Parentdf = parentSNPids, yeast = TRUE) %>% na.omit()
head(testmerge)

```

### Piping

```{r, warning=FALSE, message=FALSE}
mydatatotest = "HGV.SortedCat.vcf.output.table"

parentSNPids <- cybrConvertParentalAlleles(Truncate = TRUE)

testmerge <- cybrInputGATKTable(mydatatotest) %>% 
  cybrQualityFilter() %>% 
  cybrIDAlleles(BSAdfstart = ., Parentdf = parentSNPids, yeast = TRUE) %>% 
  na.omit()

```

### Reformat data so that it has bulk etc included

This part needs to be done manually to ensure that the formula matches the data. The two factors added here (Bulk and Parent) are then the inputs in the formula for the GLM.

```{r}
unique(testmerge$Dataset)

testmerge$Bulk <- NA
testmerge$Bulk[testmerge$Dataset == "HGVMVDRX2_n01_OakI_Dilute_B.fastq" |
              testmerge$Dataset == "HGVMVDRX2_n01_WineI_Dilute_B.fastq" |
              testmerge$Dataset == "HGVMVDRX2_n01_WineI_Dilute_C.fastq" | 
              testmerge$Dataset == "HGVMVDRX2_n01_WineI_Dilute_D.fastq"] <- "Dilute"

testmerge$Bulk[testmerge$Dataset == "HGVMVDRX2_n01_OakI_Fluc_B.fastq" |
              testmerge$Dataset == "HGVMVDRX2_n01_OakI_Fluc_C.fastq" |
              testmerge$Dataset == "HGVMVDRX2_n01_OakI_Fluc_D.fastq" | 
              testmerge$Dataset == "HGVMVDRX2_n01_WineI_Fluc_C.fastq"] <- "Fluconazole"

testmerge$Parent <- NA
testmerge$Parent[testmerge$Dataset == "HGVMVDRX2_n01_OakI_Dilute_B.fastq" |
              testmerge$Dataset == "HGVMVDRX2_n01_OakI_Fluc_B.fastq" |
              testmerge$Dataset == "HGVMVDRX2_n01_OakI_Fluc_C.fastq" | 
              testmerge$Dataset == "HGVMVDRX2_n01_OakI_Fluc_D.fastq"] <- "OakI"

testmerge$Parent[testmerge$Dataset == "HGVMVDRX2_n01_WineI_Dilute_B.fastq" |
              testmerge$Dataset == "HGVMVDRX2_n01_WineI_Dilute_C.fastq" |
              testmerge$Dataset == "HGVMVDRX2_n01_WineI_Dilute_D.fastq" | 
              testmerge$Dataset == "HGVMVDRX2_n01_WineI_Fluc_C.fastq"] <- "WineI"

testmerge$Parent <- factor(testmerge$Parent)
testmerge$Bulk <- factor(testmerge$Bulk)
testmerge$ReadCount <- as.numeric(testmerge$ReadCount)

head(testmerge)

```
```{r, warning=FALSE, message=FALSE}

#Running the entire genome in a loop instead of just chr II
wholegenomeBSA <- foreach(i=unique(testmerge$CHROM), .combine=rbind) %dopar%{
  cybrBSA_GLM(testmerge, chr = i) %>% cybrSmoothBSAWindows()
}

```

### Plot

```{r}
cybrPlotZPrime(CuSO4_wholegenomeBSA, columns = c("Bulk_Zprime"))

cybrPlotZScore(CuSO4_wholegenomeBSA)

cybrPlotZPrime(wholegenomeBSA, chromosomes = ChromKey$chromosomes[2:16])

cybrPlotZScore(wholegenomeBSA)

```

# Updated R Pipeline

### Analysis
***
#### Combine oak and wine parental alleles, and define the bulks, parents, and replicates
```
setwd("../../../GitHub/Sequencing/Analysis/")

parentSNPids <- cybrConvertParentalAlleles(ParentFiles = c("Wine_VCF.txt", "Oak_VCF.txt"), Truncate = TRUE)

parentSNPids %>% group_by(CHROM, POS) %>% summarize(Type = Type, Unique = length(unique(Type))) %>% filter(Unique == 1) -> pSNPs

#CHANGE THIS
mydatatotest = "Data/HKTFTDRX2.SortedCat.vcf.output.table"

FilteredData <- cybrInputGATKTable(mydatatotest) %>% 
  cybrQualityFilter() %>% 
  cybrIDAlleles(BSAdfstart = ., Parentdf = pSNPs, yeast = TRUE) %>% 
  na.omit()

#Using Gsub for this
gsub(FilteredData$Dataset, "HKTFTDRX2_n01_", "") #CHANGE THIS

FilteredData %>% mutate(DShort = gsub("HKTFTDRX2_n01_", "", Dataset),
                       DS = gsub(".fastq", "", DShort)) %>% select(-Dataset, -DShort) -> tempFilteredData

tempFilteredData$Bulk <- NA
tempFilteredData$Parent <- NA
tempFilteredData$Rep <- NA

tempFilteredData$Bulk[grep("C", tempFilteredData$DS)] <- "CuSO4" #CHANGE THIS
tempFilteredData$Bulk[grep("D", tempFilteredData$DS)] <- "Dilute"

tempFilteredData$Rep[grep("a", tempFilteredData$DS)] <- "A"
tempFilteredData$Rep[grep("b", tempFilteredData$DS)] <- "B"

tempFilteredData$Parent[grep("O", tempFilteredData$DS)] <- "Oak"
tempFilteredData$Parent[grep("W", tempFilteredData$DS)] <- "Wine"

tempFilteredData$ReadCount <- as.numeric(tempFilteredData$ReadCount)
  
# #THIS IGNORES REPLICATES
tempFilteredData %>% select(CHROM, POS, PAllele, ReadCount, Bulk, Parent, Rep) %>% distinct %>% 
 pivot_wider(names_from = c(Bulk, Parent, Rep, PAllele), values_from = ReadCount) -> cybr2Data
 ```
 
#### Check log alleles per flask and coverage across sequencing run
```
cybr2Data %>% pivot_longer(c(-CHROM, -POS), names_to = c("Bulk", "Parent", "Rep", "Allele"), names_sep = "_") %>%
  pivot_wider(names_from = Allele, values_from = value) %>% mutate(Coverage = Wine + Oak, logWineOak = log(Wine/Oak)) -> RawCountSummary
```
 
#### Smooth data by rolling mean or median

```
#Use rolling average of 100 SNPs, finding the mean
cybr2Data %>% cybr2_rollmean() -> rollmeanData

#Find the rolling median or change n instead
cybr2Data %>% pivot_longer(c(-CHROM, -POS), names_to = "label") %>% 
  group_by(CHROM, label) %>% arrange(POS) %>% 
  summarize(POS = POS, 
            CHROM = CHROM, 
            SmoothCount = ceiling(frollapply(value, n = 100, FUN = median))) %>% #CHANGE THIS LINE
  na.omit() %>% 
  pivot_wider(names_from = label,values_from = SmoothCount) -> rollData
```

#### Caluclate GLM of rolling data

```
#Change for different datasets
mydata <- rollData

#Run GLM - options are glmfixed_rep(), glmfixed_rep3(), and glmfixed()
#IF NOT USING REPS, change 1:5 to 1:4, and remove "Rep" from labels
mydata %>% group_by(CHROM, POS) %>% 
  summarise(summary = glmfixed_rep(HOOa = CuSO4_Oak_A_Oak, 
                           HOWa = CuSO4_Oak_A_Wine, 
                           HWOa = CuSO4_Wine_A_Oak,
                           HWWa = CuSO4_Wine_A_Wine,
                           LOOa = Dilute_Oak_A_Oak,
                           LOWa = Dilute_Oak_A_Wine,
                           LWOa = Dilute_Wine_A_Oak, 
                           LWWa = Dilute_Wine_A_Wine,
                           
                           HOOb = CuSO4_Oak_B_Oak, 
                           HOWb = CuSO4_Oak_B_Wine, 
                           HWOb = CuSO4_Wine_B_Oak,
                           HWWb = CuSO4_Wine_B_Wine,
                           LOOb = Dilute_Oak_B_Oak,
                           LOWb = Dilute_Oak_B_Wine,
                           LWOb = Dilute_Wine_B_Oak, 
                           LWWb = Dilute_Wine_B_Wine)[1:5],
                                                   label = c("intercept", "Bulk", "Parent", "Rep", "Interaction")) -> GLMdata
```

### Visualizing
***
#### Single GLM Plot for this Data for reference

```
GLMdata %>% 
  filter(label != "intercept", CHROM != "I", CHROM != "M") %>% 
  
  ggplot(aes(x = POS, y = summary, color = label)) + geom_line() + 
  
  facet_grid(~CHROM, scales = "free", space = "free") +
  theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) + ggtitle("GLM of Data")
```

#### Log Odds of Alleles Plot for reference

```
RawCountSummary %>% 
  filter(CHROM != "I", CHROM != "M") %>% 
  
  ggplot(aes(x = POS, y = logWineOak, shape = paste(Bulk, Parent, Rep, sep = "_"), color = Bulk)) + 
  geom_point(alpha = 0.3) + 
  
  facet_grid(~CHROM, scales = "free", space = "free") +
  theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
  scale_color_manual(values = c("Violet", "Black"))
```
