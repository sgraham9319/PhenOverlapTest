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

# Create vector of species names in alphabetical order
sps <- sort(unique(chaut$Species))

# Create a matrix with a cell for each species combination
results <- matrix(nrow = length(sps), ncol = length(sps))

# Attach species names to the columns and rows
colnames(results) <- sps
rownames(results) <- sps

#################################
# Evaluating phenological overlap
#################################

# Loop through the focal species (sp1)
for(sp1 in 1:length(sps)) {
  # Find ordinal dates when focal species present
  dates <- chaut[chaut$Species == sps[sp1] & chaut$Stage6 > 0, "OrdinalDate"]
  # Subset to these dates
  sub.dat <- chaut[chaut$OrdinalDate %in% dates,]
  # Calculate total number adults of focal species in year
  tot.sp1 <- tapply(sub.dat$Stage6, sub.dat$Species, FUN = sum)[sps[sp1]]
  # Create empty vector for overlap by ordinal date
  measures <- rep(NA, times = length(dates))
  
  # Loop through each competitor species (sp2)
  for(sp2 in 1:length(sps)) {
    # Avoid calculating overlap between a focal species and itself
    if(sp2 != sp1) {
      # Subset to the 2 species of interest
      sub.sp <- sub.dat[sub.dat$Species == sps[sp1] | 
                          sub.dat$Species == sps[sp2], ]
      
      # Loop through dates, calculating overlap for each date
      for (ord in 1:length(dates)) {
        # Extract number of adults for focal species
        focals <- sub.sp[sub.sp$OrdinalDate == dates[ord] & 
                        sub.sp$Species == sps[sp1], "Stage6"]
        # Extract number of adults for competitor species
        competitors <- sub.sp[sub.sp$OrdinalDate == dates[ord] & 
                             sub.sp$Species == sps[sp2], "Stage6"]
        # Check that a row representing competitor species on this date exists
        if(length(competitors) == 0){
          # If no row exists make number of competitors zero
          competitors <- 0
        }
        # Calculate ratio of competitors to focals
        ratio <- competitors/focals
        # Weight the ratio by the proportion of focal adults for the year that
        # were present on this date
        measures[ord] <- ratio*(focals/tot.sp1)
      }
      # Sum measures of overlap to obtain weighted mean and add to results matrix
      results[sp1, sp2] <- sum(measures)
    }
  }
}
