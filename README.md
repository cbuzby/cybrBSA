# cybrBSA
R Package for analyzing sequencing data for multilevel Bulk Segregant Analysis experiment. 

### Convert VCF to GATK table (pre-processing for cybrBSA)

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
  select(POS, CHROM, Dataset, GQ, AD.REF, AD.ALT, Coverage) %>% mutate(Pool = RawFiles$Pool[i])-> rawdata

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

This uses the glm_cb2_short() function to take in a list of columns which should be used in the formula and returns the Z score column of the glm summary coefficients
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

### Visualize Data

Data can be plotted in ggplot() by faceting by chromosome (CHROM) so that positions do not overlap from different chromosomes. 
```
processed_glm_all %>% 
  ggplot(aes(x = POS, y = abs(Summary), color = Factor)) + geom_line() +
  facet_grid(~CHROM, space = "free", scales = "free") + ggtitle("BSA Bulk Z Scores")
```
