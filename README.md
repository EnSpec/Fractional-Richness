# CNFPD 
Cumulative Normalized Fractional Population Density: a biodiversity index for large camera trap networks

CNFPD is a biodiversity index designed to provide intuitive and ecologically informative measures of biodiversity when working with large camera trap networks.
Camera trap networks have a number of inherent biases which make typical diversity indexes like Shannon diversity less informative. E.g. not all species are equally detectable, detections follow a poisson distribution rather than a normal distribution, and species of different body sizes or trophic levels have different maximum population densities. 
CNFPD accounts for these idiosycrasies of camera trap data.

**Language**: R



## The CNFPD diversity index behaves differently than other diversity indexes (e.g. Shannon diversity) in several key ways:
1. **All species are valued equally**. With Shannon Diversity, low-density or rare species have less impact on the overall diversity score. With CNFPD, all species contribute equally to the diversity score.
2. **Diversity increases linearly with species richness**. Shannon Diversity has a plateauing curvilinear distribution, so that the resolution is reduced at higher species richness. CNFPD has a linear relationship with richness, and performs equally well at high and low species richness. When all species are at their average population densities, CNFPD is equal to ½ species richness. When all species are at their maximum population densities, CNFPD is equal to species richness. 
3. **An increase in richness or abundance always results in an increase in diversity**. With Shannon diversity, if the abundances of species are very uneven, then a decrease in population or a complete loss of the most common species can counterintuitively result in a higher measured diversity. This is never the case with CNFPD.
4. **Population decline results in a loss of diversity**. If all species become scarce but evenness remains the same, CNFPD decreases while Shannon diversity remains constant. 
5. **Species detectability does not bias the index**. Practically speaking, the measured abundance of a species, or its detection rate, is never exactly the same as the true abundance. Some species are more easily detectable than others. Because CNFPD values all species equally, and dsite,i and dmax,i are affected by the same set of biases, species detectability bias is effectively canceled out and does not affect the index, so long as a species that is present is detected. 
6. **Multiple sites are required**. The main drawback of CNFPD when compared to Shannon diversity is that it requires more information. Shannon diversity can be calculated for a single site, while CNFPD requires either multiple sites or multiple time points in order to calculate  dmax and mean fractional population density.



# Data requirements
To calculate CNFPD, detection rates for multiple species at multiple sites are required. 
Each row should contain data for one site and each column should contain detection rates for one species.

For example:

ID    | Species 1  | Species 2  | Species 3
----- | ---------- | ---------- | ---------
Site1 |  0.2       |  8         |  0.001
Site2 |  0.3       |  7         |  0.007
Site3 |  0.1       |  9         |  0.004




# Equation

$CNFPD_{site} = \sum_{i=1}^{n}\left ( \frac{d_{site,i}}{d_{max,i}} \right )^{\frac{log0.5}{log\overline{\left ( \frac{d_{site,i}}{d_{max,i}} \right )}}}$

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
The mean FPD of all sites where the species is present, excluding sites where the species has not been detected.



**NFPD = Normalized Fractional Population Density**

Fractional population density normalized so that the distribution of detection rates centers around 0.5. 
Normalization is species specific. For each species, sites with an average detection rate will have a NFPD of 0.5. 



**CNFPD = Cumulative Normalized Population Density**

The sum of NFPD for each species. At sites where all species are at their typical population densities, CNFPD = 1/2 species richness.



**Arguments**
 + data:             a dataframe containing species detection rates
 + species_cols:     columns which contain species detection rates, one column per species
 + ID_col:           column containing site ID



# Usage

Calculate CNFPD for a selected community of species and visualize using mapview:

```{r}
library(sf)
library(mapview)

#load dataset
DBC <- read.csv("DBC_clean.csv")

#Calculate CNFPD for only sensitive species
Sensitive.CNFPD <- CNFPD(data=DBC, 
                         species_cols = c(5:7, 15, 16, 17, 26, 33), 
                         ID_col = "CAMERA_ID")

#bind CNFPD values to original dataframe
Sensitive.CNFPD <- cbind(DBC, Sensitive.CNFPD)

#visualize
Sensitive.CNFPD <- st_as_sf(Sensitive.CNFPD, coords = c("LON", "LAT"), crs=4326)
mapview(Sensitive.CNFPD, zcol="CNFPD")
```

Intermediate steps can also be calculated independently:

```{r}
x <- FPD(DBC$BearsPerDay)

x <- meanFPD(DBC$BearsPerDay)

x <- NFPD(DBC$BearsPerDay)
```


# Provided dataset

An example dataset is available in this github which can be used to calculate CNFPD: Snapshot_DBC.csv

Snapshot_DBC.csv is data from the Snapshot Wisconsin camera trap network (https://dnr.wisconsin.gov/topic/research/projects/snapshot). 
This dataset includes detection rates of 33 species across 2218 camera locations in Wisconsin.
Detection rates are measured as the number of detection events / the number of days the camera location was active.
LAT/LON coordinates are coarsened for the privacy of community scientists in this publicly available dataset.

**Variables**
  + CAMERA_ID: Unique ID number of each camera location
  +	LAT: Latitude in decimal degrees, coarsened to 2 decimal places to maintain privacy of community scientists
  +	LON: Longitude in decimal degrees, coarsened to 1 decimal place to maintain privacy of community scientists
  + SpeciesPerDay: One column per species. The number of times that species was detected / the number of days that camera site was active.


# Citation

For more details, or if you use CNFPD in your work, please cite:

Berman, L.; Schneider, F.; Stenglein, J.; Bemowski, R.; Dean, M.; Townsend, P. (in preparation). Cumulative Normalized Fractional Population Density: a diversity index for camera trap networks. 




