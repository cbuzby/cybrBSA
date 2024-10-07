# cybrBSA
R Package for analyzing sequencing data for multilevel Bulk Segregant Analysis experiment. For more information, please see pre-print at doi: or our [epicQTL](https://github.com/Siegallab/epicQTL) github.

*** 
## Pre-processing
### Convert VCF to GATK table (pre-processing for cybrBSA)

_cybrBSA_ takes in a data frame of variants (REF and ALT alleles) and their counts per sequencing bulk. 

Bash script to convert the VCF, keeping specific fields listed after -GF. AD (allele depth) is necessary for all analysis, as this is the read count to analyze in the glm. GQ (genome quality) can be used as a filter in the `cybrQualityFilter()` function. DP and PL are optional. myfile is the name of the VCF to convert.

```{bash, eval = FALSE}
module load gatk/4.2.0.0

gatk VariantsToTable \
     -V ${myfile} \
     -F CHROM -F POS -F REF -F ALT \
     -GF AD -GF DP -GF GQ -GF PL \
     -O ${myfile}.output.table

```

***  
## Analysis

### Variant Table Processing

The following file is a truncated version of the raw data which includes only Chr XI. For the full dataset, please visit [epicQTL](https://github.com/Siegallab/epicQTL) and replace `AllCuSO4.SortedCat.vcf.output.table` with unzipped [rawdata](https://github.com/Siegallab/epicQTL/blob/main/BSA_Analysis/Input/AllCuSO4.REF_.SortedCat.vcf.output.zip)
```{r}
mydatatotest = "AllCuSO4.SortedCat.vcf.output.table"

cybrInputGATKTable(mydatatotest) %>% mutate(Coverage = as.numeric(AD.REF) + as.numeric(AD.ALT)) %>%
  select(POS, CHROM, Dataset, GQ, AD.REF, AD.ALT, Coverage) -> rawdata
```
Ensure that parental alleles Oak.txt and Wine.txt are in the same folder as this script, and then run the following to call the parent of origin for each ALT allele.
```
parentSNPids <- cybrConvertParentalAlleles(Truncate = TRUE)

rawdata %>% 
  merge(parentSNPids) %>% 
  filter(grepl("HNGLVDRXY", Dataset)) %>% #Filter for smaller dataset
  mutate(REFW = as.numeric(Type == "Wine"), REFO = as.numeric(Type == "Oak")) %>%
  group_by(Dataset, CHROM, POS) %>%
  mutate(Wine = max(REFW * as.numeric(AD.REF), REFO * as.numeric(AD.ALT)),
         Oak = max(REFW * as.numeric(AD.ALT), REFO * as.numeric(AD.REF))) %>%
  select(Dataset, POS, CHROM, Coverage, Wine, Oak) %>%
  pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "Reads") -> rawdata_called
```

### Smooth Data
Next, smooth data using `cybr_weightedgass()`, which will apply a weighted gaussian that uses n/10 as the standard deviation.

```{r}
rawdata_called %>% group_by(Dataset, CHROM, Allele) %>% arrange(POS) %>%
  reframe(POS = POS, 
          SmoothCount = ceiling(frollapply(Reads, n = 200, FUN = cybr_weightedgauss, align = "center"))
          ) -> rawdata_smoothed 

```

### Convert factors and contrasts
Formatting output.table data frame for consistent analysis; convert coefficients to factors. Here we select only the **HNGLVDRXY** sequencing bulk for demonstration.
```{r}
# Pick out the specific groups of the dataset
rawdata_smoothed %>% 
    mutate(Dataset = gsub("HNGLVDRXY_n01_CuSO4_CSSI_", "", Dataset)) %>%
    mutate(Dataset = gsub(".fastq", "", Dataset)) %>%
    mutate(Dataset = gsub("Unselected", "Unselected_", Dataset)) %>%
    mutate(Dataset = gsub("Selected", "Selected_", Dataset)) %>%
  separate(Dataset, into = c("Bulk", "Parent"), sep = "_") %>% 
  mutate_if(is.character, as.factor) -> rawdata_glm_prep
```
Often, the contrasts of the coefficient factors might need to be adjusted, depending on the question interrogated. We ensure that the unselected bulk is labeled 0 and the selected 1, for ease of interpretation of our results.
```{r}
contrasts(rawdata_glm_prep$Parent) <- matrix(c(-0.5, 0.5))
contrasts(rawdata_glm_prep$Bulk) <- matrix(c(1,0))
```

### Calculate z-scores for glm

First, isolate a single position and run the glm or glmer to verify the number of output coefficients and labels to include:
```
test <- rawdata_glm_prep[rawdata_glm_prep$CHROM == "I" & rawdata_glm_prep$POS == min(rawdata_glm_prep$POS)),]
W = SmoothCount
F = "Allele ~ Bulk * Parent"

testg <- glm(data = test, formula = F, weights = W, family = "binomial")
summary(testg)
```

To calculate logistic regression, use either `glm_cb2_short()` or `glmer_cb2_short()` (requires [lme4](https://cran.r-project.org/web/packages/lme4/index.html) package), which will conduct the logistic regression for EACH individual position. Note that we use dplyr's `summarize` (or `reframe`) to parallelize for all positions. Using a mixed model will take significantly more time. The outputlength parameter of `glm_cb2_short()` should match the labels (here shown as "Factor") within `summarize()`.
```
# Run full dataset
rawdata_glm_prep %>% na.omit() %>% 
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  mutate_if(is.character, as.factor) %>%
  summarize(summary = glm_cb2_short(Allele = Allele,
                             Bulk = Bulk,
                             Parent = Parent,
                             #Rep = Rep, #any parameters not used should NOT be included
                             W = SmoothCount,
                             formula = "Allele ~ Bulk * Parent",
                             outputlength = 4),
            #MAKE SURE THIS IS THE SAME LENGTH AS OUTPUT LENGTH
            Factor = (c("Intercept", "Bulk", "Parent", "Interaction"))) -> processed_glm_all
```

### Permutations for False Discovery Rate
Complete the following TWICE, one for each new pseudoreplicate "Bulk". 
```{r}
rawdata_glm_prep %>% na.omit() %>% 
  filter(Bulk == "Unselected") %>%
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  pivot_wider(names_from = Allele, values_from = SmoothCount) %>% unnest() %>%
  mutate(label = paste(CHROM, POS, sep = "_")) %>%
  group_by(Bulk, Parent) %>%
  summarize(NewPOS = sample(label),
            Oak = Oak,
            Wine = Wine,
            Parent = Parent) %>%
  separate(NewPOS, into = c("CHROM", "POS"), sep = "_") %>%
  mutate(Bulk = "Selected") -> rd_shuffled_selected

rawdata_glm_prep %>% na.omit() %>% 
  filter(Bulk == "Unselected") %>%
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  pivot_wider(names_from = Allele, values_from = SmoothCount) %>% unnest() %>%
  mutate(label = paste(CHROM, POS, sep = "_")) %>%
  group_by(Bulk, Parent) %>%
  summarize(NewPOS = sample(label),
            Oak = Oak,
            Wine = Wine,
            Parent = Parent) %>%
  separate(NewPOS, into = c("CHROM", "POS"), sep = "_") %>%
  mutate(Bulk = "Unselected") -> rd_shuffled_unselected
```

Combine the two permuted datasets and set contrasts to match original data analysis.
```
rbind(rd_shuffled_selected,rd_shuffled_unselected) %>%
  ungroup() %>%
  pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "SmoothCount") %>% unnest() %>%
  mutate(POS = as.numeric(POS)) %>%
  mutate_if(is.character, as.factor) -> Perm1contrasts

#Set contrasts for experiment
contrasts(Perm1contrasts$Parent) <- matrix(c(-0.5, 0.5))
contrasts(Perm1contrasts$Bulk) <- matrix(c(1,0))
```
Conduct the analysis for permuted data as done above.
```
Perm1contrasts %>% na.omit() %>% 
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  mutate_if(is.character, as.factor) %>%
  summarize(summary = glm_cb2_short(Allele = Allele,
                             Bulk = Bulk,
                             Parent = Parent,
                             W = SmoothCount,
                             formula = "Allele ~ Bulk * Parent",
                             outputlength = 4),
            #MAKE SURE THIS IS THE SAME LENGTH AS OUTPUT LENGTH
            Factor = (c("Intercept", "Bulk", "Parent", "Interaction"))) -> permuted_glm_all
```

Find 5% quantile for each coefficient (labeled Factor):
```
permuted_glm_all %>% ungroup() %>% group_by(Factor) %>% quantile(summary, 0.95)
```

### Call Peaks
Renaming columns for `cybr_lmpeaks()` as multiple chromosomes and fixed chromosomes (CSS) might be used:
```{r}
processed_glm_all$zscore <- processed_glm_all$summary
processed_glm_all$CSS <- "I"
processed_glm_all$label <- processed_glm_all$Factor
```

Script for finding Bulk peaks within 300bp windows. This can be parallelized using dplyr's `rollapply()`.
```
testpeaks <- cybr_lmpeaks(processed_glm_all[processed_glm_all$Factor == "Bulk",], width = 300)
```

## Visualizing Data
### Manhattan Plot
Data from this example will include one chromosome, but it is important to facet by chromosome for datasets containing multiple chromosomes.
```{r}
processed_glm_all %>% 
  ggplot(aes(x = POS, y = abs(summary), color = Factor)) + geom_line() +
  facet_grid(~CHROM, space = "free", scales = "free") + ggtitle("BSA Bulk Z Scores")
```

### Circos Plot
Interactions with a fixed chromosome can be visualized using `cybr_circos()`, which takes in data frames for 
```{r}
#Ensure that data frame contains the correct column names for this function:
processed_glm_all$label <- processed_glm_all$Factor
processed_glm_all <- processed_glm_all %>% select(-Factor)

Plot:
cybr_circos2(d1 = processed_glm_all, 
            d8 = processed_glm_all,
            peaklist8 = testpeaks)
```
