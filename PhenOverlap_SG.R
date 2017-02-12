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
chaut <- chaut[which(!is.na(chaut$Stage6)),]

# Find out which species had more than 50 adults sampled
x <- tapply(chaut$Stage6, chaut$Species, FUN = sum)
x1 <- x[!is.na(x) & x > 50]
names(x1)

# Subset data to these common species
chaut <- chaut[chaut$Species %in% names(x1),]
chaut <- droplevels(chaut)
table(chaut$Species)

#############################
# First try doing just 1 year
#############################

# Subset to 1 year
chaut <- chaut[chaut$year == 2008, ]

# Find all species combinations
sp.combs <- combn(unique(chaut$Species), m = 2)
ncol(sp.combs)

results <- rep(NA, times = ncol(sp.combs))

for(com in 1:ncol(sp.combs)) {
  # Subset to the 2 species of interest
  dat <- chaut[chaut$Species == sp.combs[1, com] | chaut$Species == sp.combs[2, com], ]
  # Find ordinal dates when focal species present
  dates <- dat[dat$Species == sp.combs[1, com] & dat$Stage6 > 0, "OrdinalDate"]
  # Subset to these dates
  dat <- dat[dat$OrdinalDate %in% dates,]
  
  # Create empty vector for overlap by ordinal date
  measures <- rep(NA, times = length(dates))
  # Calculate total number adults of focal species in year
  tot.sp1 <- tapply(dat$Stage6, dat$Species, FUN = sum)[sp.combs[1, com]]
  
  # Loop through dates, calculating overlap for each one
  for (i in 1:length(dates)) {
    # Extract number of adults for focal species
    sp1 <- dat[dat$OrdinalDate == dates[i] & dat$Species == sp.combs[1, com],
                 "Stage6"]
    # Extract number of adults for competitor species
    sp2 <- dat[dat$OrdinalDate == dates[i] & dat$Species == sp.combs[2, com],
                 "Stage6"]
    # Check for a row representing competitor species on this date
    if(length(sp2) == 0){
      # If no row, make number of competitors zero
      sp2 <- 0
    }
    # Calculate ratio
    ratio <- sp2/sp1
    # Weight the ratio by the proportion of focal adults for the year that
    # were sampled on this date
    measures[i] <- ratio*(sp1/tot.sp1)
  }
  # Take weighted mean of overlap measures as measure of phenological overlap
  results[com] <- sum(measures)
}
results
