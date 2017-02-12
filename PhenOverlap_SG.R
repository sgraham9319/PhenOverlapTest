##################################
# Quantifying phenological overlap
##################################

# Load data
chaut.dat <- read.csv("../Chautdata.csv")

# Check for mispelled species names
table(chaut.dat$Species)

# Restrict ordinal date to 6/1 to 9/16
chaut <- chaut.dat[chaut.dat$OrdinalDate >= 152 &
                     chaut.dat$OrdinalDate <= 259, ]

# Subset to adults
chaut <- droplevels(chaut[which(!is.na(chaut$Stage6)),])

# Find out which species had more than 50 adults sampled
x <- tapply(chaut$Stage6, chaut$Species, FUN = sum)
x1 <- x[!is.na(x) & x > 50]
names(x1)

# Subset data to these common species
chaut <- droplevels(chaut[chaut$Species %in% names(x1),])
table(chaut$Species)

#############################
# First try doing just 1 year
#############################

# Subset to 1 year
chaut <- chaut[chaut$year == 2008, ]

########################
# Combinations both ways
########################

unique(chaut$Species)
sps <- sort(unique(chaut$Species))
sps[3]

results <- matrix(nrow = length(sps), ncol = length(sps))
colnames(results) <- sps
rownames(results) <- sps

for(sp1 in 1:length(sps)) {
  # Find ordinal dates when focal species present
  dates <- chaut[chaut$Species == sps[sp1] & chaut$Stage6 > 0, "OrdinalDate"]
  # Subset to these dates
  dat <- chaut[chaut$OrdinalDate %in% dates,]
  # Calculate total number adults of focal species in year
  tot.sp1 <- tapply(dat$Stage6, dat$Species, FUN = sum)[sps[sp1]]
  # Create empty vector for overlap by ordinal date
  measures <- rep(NA, times = length(dates))
  
  # Loop through each competitor species
  for(comp in 1:length(sps)) {
    if(comp != sp1) {
      # Subset to the 2 species of interest
      dat1 <- dat[dat$Species == sps[sp1] | dat$Species == sps[comp], ]
      
      # Loop through dates, calculating overlap for each one
      for (ord in 1:length(dates)) {
        # Extract number of adults for focal species
        focal <- dat1[dat1$OrdinalDate == dates[ord] & dat1$Species == sps[sp1],
                   "Stage6"]
        # Extract number of adults for competitor species
        competitor <- dat1[dat1$OrdinalDate == dates[ord] & dat1$Species == sps[comp],
                   "Stage6"]
        # Check for a row representing competitor species on this date
        if(length(competitor) == 0){
          # If no row, make number of competitors zero
          competitor <- 0
        }
        # Calculate ratio
        ratio <- competitor/focal
        # Weight the ratio by the proportion of focal adults for the year that
        # were sampled on this date
        measures[ord] <- ratio*(focal/tot.sp1)
      }
      # Take weighted mean of overlap measures as measure of phenological overlap
      results[sp1, comp] <- sum(measures)
    }
  }
}
