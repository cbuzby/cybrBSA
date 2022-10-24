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
