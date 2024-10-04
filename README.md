# cybrBSA
R Package for analyzing sequencing data for multilevel Bulk Segregant Analysis experiment. For more information, please see pre-print at doi: 

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
### Variant Table Processing

#### Processing by Individual Functions

```{r, eval = FALSE, warning=FALSE, message=FALSE}
mydatatotest = "HGV.SortedCat.vcf.output.table"

mydf <- cybrInputGATKTable(mydatatotest)

#quality can be filtered but we've found that this is not useful for our dataset and have not used this step
qualitydf <- cybrQualityFilter(mydf)

parentSNPids <- cybrConvertParentalAlleles(Truncate = TRUE)

rawdata_called <- cybrIDAlleles(BSAdfstart = qualitydf, Parentdf = parentSNPids, yeast = TRUE) %>% na.omit()

```

#### Processing by Piping

```{r, warning=FALSE, message=FALSE}
mydatatotest = "HGV.SortedCat.vcf.output.table"

cybrInputGATKTable(mydatatotest) %>% mutate(Coverage = as.numeric(AD.REF) + as.numeric(AD.ALT)) %>%
  select(POS, CHROM, Dataset, GQ, AD.REF, AD.ALT, Coverage) -> rawdata

parentSNPids <- cybrConvertParentalAlleles(Truncate = TRUE)

rawdata %>% 
  merge(parentSNPids) %>% 
  mutate(REFW = as.numeric(Type == "Wine"), REFO = as.numeric(Type == "Oak")) %>%
  group_by(Dataset, CHROM, POS) %>%
  mutate(Wine = max(REFW * as.numeric(AD.REF), REFO * as.numeric(AD.ALT)),
         Oak = max(REFW * as.numeric(AD.ALT), REFO * as.numeric(AD.REF))) %>%
  select(Pool, Dataset, POS, CHROM, Coverage, Wine, Oak) %>%
  pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "Reads") -> rawdata_called

```

### Smooth data by median

For dense, non-intercross data, use the data.table function frollapply to summarize within each experimental flask (ie unique bulk, parent, etc) within a sliding window. We use 200 SNPs of data as our window width, and align to the center so that the peak is likely at the specific gene that we are looking for. Smoothing takes advantage of linkage within the genome to even out sequencing noise.
```
rawdata_called %>% group_by(Pool, Dataset, CHROM, Allele) %>% arrange(POS) %>%
  reframe(POS = POS, SmoothCount = ceiling(frollapply(Reads, n = 200, FUN = median, align = "center"))) -> rawdata_smoothed 
```

To smooth by weighted Gaussian, instead use the cybr_weightedgauss() function:

```
rawdata_called %>% group_by(Pool, Dataset, CHROM, Allele) %>% arrange(POS) %>%
  reframe(POS = POS, SmoothCount = ceiling(frollapply(Reads, n = 200, FUN = cybr_weightedgauss, align = "center"))) -> rawdata_smoothed 
```

### Separate out columns to run the GLM

This part needs to be done for each experiment to ensure that the formula matches the data. The two factors added here (Bulk and Parent) are then the inputs in the formula for the GLM.

This sample has a Dataset column which contains the annotations from the sequencing run. By using `mutate()`, we can create underscores between areas to separate, and then use dplyr's `separate()` to make separate columns for each. Once these columns are separate, the data is ready for the `glm()` function.

```{r}
# Pick out the specific groups of the dataset
rawdata_smoothed %>% 
    mutate(Dataset = gsub("CuSO4_CSSI_", "", Dataset)) %>%
    mutate(Dataset = gsub("Unselected", "Unselected_", Dataset)) %>%
    mutate(Dataset = gsub("Selected", "Selected_", Dataset)) %>%
  separate(Dataset, into = c("Bulk", "Parent"), sep = "_") %>% 
  mutate_if(is.character, as.factor) -> rawdata_glm_prep
```

### Calculate GLM of rolling data

#### Test with a single locus to identify output length

```
#Test once using a position that has enough levels to run the glm
testdata <- rawdata_glm_prep %>% na.omit() %>% filter(CHROM == "II", POS == 65975) %>% arrange(Bulk, Parent)

#Does it have enough levels, ie at least two different values for each element in the formula?
table(testdata$Bulk, testdata$Parent)

#Define what will be used in the function
formula = "Allele ~ Bulk * Parent"

#Run the test glm
testglm <- glm(formula = formula, family = "binomial", data = testdata, weights = SmoothCount)

#print to view
testglm

#actual outputs
summary(testglm)$coefficients
```

#### Parallelize using dplyr and glm_cb2_short()

This uses the glm_cb2_short() or glmer_cb2_short() function to take in a list of columns which should be used in the formula and returns the Z score column of the glm (or glmer if mixed effect model) summary coefficients
```
# Run full dataset
rawdata_glm_prep %>% na.omit() %>% 
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  mutate_if(is.character, as.factor) %>%
  summarize(Summary = glm_cb2_short(Allele = Allele,
                             Bulk = Bulk,
                             Parent = Parent,
                             #Rep = Rep,
                             W = SmoothCount,
                             formula = "Allele ~ Bulk * Parent",
                             outputlength = 4),
            #MAKE SURE THIS IS THE SAME LENGTH AS OUTPUT LENGTH
            Factor = (c("Intercept", "Bulk", "Parent", "Interaction"))) -> processed_glm_all
```

### Peak Calling

False discovery rate can be determined by permuting non-selected bulks; remove selected bulks and all altered chromosomes (fixed or containing selection markers), then on each bulk individually, shuffle positions while keeping allele ratios consistent. Conduct twice so that there are two (selected, unselected) pseudo-bulks to compare, and run analysis as before. Function for permutations (cybrPermute_byParent_Rep) can be used if data is specifically in the correct format.

```
rawdata_glm_prep %>% na.omit() %>% 
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  pivot_wider(names_from = Allele, values_from = SmoothCount) %>% unnest() %>%
  mutate(label = paste(CHROM, POS, sep = "_")) %>%
  group_by(Pool, Dataset) %>%
  summarize(NewPOS = sample(label),
            Oak = Oak,
            Wine = Wine,
            Background = Background) %>%
  separate(NewPOS, into = c("CHROM", "POS"), sep = "_") %>%
  mutate(Bulk = "Selected") -> rd_shuffled_selected

rawdata_glm_prep %>% na.omit() %>% 
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  pivot_wider(names_from = Allele, values_from = SmoothCount) %>% unnest() %>%
  mutate(label = paste(CHROM, POS, sep = "_")) %>%
  group_by(Pool, Dataset) %>%
  summarize(NewPOS = sample(label),
            Oak = Oak,
            Wine = Wine,
            Background = Background) %>%
  separate(NewPOS, into = c("CHROM", "POS"), sep = "_") %>%
  mutate(Bulk = "Selected") -> rd_shuffled_unselected

rbind(rd_shuffled_selected,rd_shuffled_unselected) %>%
  ungroup() %>%
  pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "SmoothCount") %>% unnest() %>%
  mutate(POS = as.numeric(POS)) %>%
  mutate_if(is.character, as.factor) -> Perm1contrasts

contrasts(Perm1contrasts$Background) <- matrix(c(-0.5, 0.5))
contrasts(Perm1contrasts$Bulk)

#Then run as if perm1contrasts is rawdata_glm_prep in analysis, and calculate quantile for each coefficient individually

```

Peaks from z-scores can be called by piping through two functions based on the derivative: 
```
processed_glm_all %>% 
   na.omit() %>% group_by(CHROM, CSS) %>% arrange(POS) %>%
  summarize(POS = POS, abs_zscore = abs_zscore, smooth_abs_z = smooth_abs_z,
            slope = frollapply(smooth_abs_z, FUN = slope_change, n = 200, align = "center")) %>% 
  mutate(negative = slope < 0) %>% 
  na.omit() %>%
  summarize(POS = POS, crossover = frollapply(negative, FUN = subtract2, n = 2)) %>%
  filter(crossover == 1) -> temp4

```

### Visualize Data

Data can be plotted in ggplot() by faceting by chromosome (CHROM) so that positions do not overlap from different chromosomes. 
```
processed_glm_all %>% 
  ggplot(aes(x = POS, y = abs(Summary), color = Factor)) + geom_line() +
  facet_grid(~CHROM, space = "free", scales = "free") + ggtitle("BSA Bulk Z Scores")
```

Peaks for interactions between chromosomes can be visualized using cybr_circos(), which will use the z-scores to map traces and then connect chromosomes with the peaks from a list of chromosome and position.

```
cybr_circos(d1 = d1AllC, 
            d8 = d8AllC,
            peaklist8 = peaks8)
```