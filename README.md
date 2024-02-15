# CNFPD 
Cumulative Normalized Fractional Population Density: a biodiversity index for large camera trap networks

CNFPD is a biodiversity index designed to provide intuitive and ecologically informative measures of biodiversity when working with large camera trap networks.
Camera trap networks have a number of inherent biases which make typical diversity indexes like Shannon diversity less informative. E.g. not all species are equally detectable, detections follow a poisson distribution rather than a normal distribution, and species of different body sizes or trophic levels have different maximum population densities. 
CNFPD accounts for these idiosycrasies of camera trap data.



## The CNFPD diversity index behaves differently than other diversity indexes (e.g. Shannon diversity) in several key ways:
1. **All species are valued equally**. With Shannon Diversity, low-density or rare species have less impact on the overall diversity score. With CNFPD, all species contribute equally to the diversity score.
2. **Diversity increases linearly with species richness**. Shannon Diversity has a plateauing curvilinear distribution, so that the resolution is reduced at higher species richness. CNFPD has a linear relationship with richness, and performs equally well at high and low species richness. When all species are at their average population densities, CNFPD is equal to ½ species richness. When all species are at their maximum population densities, CNFPD is equal to species richness. 
3. **An increase in richness or abundance always results in an increase in diversity**. With Shannon diversity, if the abundances of species are very uneven, then a decrease in population or a complete loss of the most common species can counterintuitively result in a higher measured diversity. This is never the case with CNFPD.
4. **Population decline results in a loss of diversity**. If all species become scarce but evenness remains the same, CNFPD decreases while Shannon diversity remains constant. 
5. **Species detectability does not bias the index**. Practically speaking, the measured abundance of a species, or its detection rate, is never exactly the same as the true abundance. Some species are more easily detectable than others. Because CNFPD values all species equally, and dsite,i and dmax,i are affected by the same set of biases, species detectability bias is effectively canceled out and does not affect the index, so long as a species that is present is detected. 
6. **Multiple sites are required**. The main drawback of CNFPD when compared to Shannon diversity is that it requires more information. Shannon diversity can be calculated for a single site, while CNFPD requires either multiple sites or multiple time points in order to calculate  dmax and mean fractional population density.



## Data requirements
To calculate CNFPD, detection rates for multiple species at multiple sites are required. 
Each row should contain data for one site and each column should contain detection rates for one species.

For example:

ID    | Species 1  | Species 2  | Species 3
----- | ---------- | ---------- | ---------
Site1 |  0.2       |  8         |  0.001
Site2 |  0.3       |  7         |  0.007
Site3 |  0.1       |  9         |  0.004




## Equation

Creating the following 4 functions in R allows you to calculate CNFPD:

```{r}
FPD <- function(species_col) {species_col / max(species_col)}
meanFPD <- function(species_col) {mean(FPD(species_col)[FPD(species_col)>0])}
NFPD <- function(species_col) {(FPD(species_col))^(log(0.5)/log(meanFPD(species_col)))}

CNFPD <- function(data, species_cols, ID_col){
  data.frame(ID = data[[ID_col]],
             CNFPD = rowSums(sapply(data[species_cols], NFPD), na.rm=TRUE))
}
```

**FPD = Fractional Population Density**

The detection rate of a species at a site as a fraction of the maximum detection rate of that species across all sites. 



**meanFPD = mean Fractional Population Density**

A single value for each species. A measure of the skewness of the data distribution. 



**NFPD = Normalized Fractional Population Density**








