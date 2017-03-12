##################################
# Quantifying phenological overlap
##################################

# Load grasshopper sampling datasets
chaut <- read.csv("../Chautdata.csv")
a1 <- read.csv("../A1data.csv")
b1 <- read.csv("../B1data.csv")
c1 <- read.csv("../C1data.csv")

###############
# Cleaning data
###############

# Remove unwanted column from B1.dat
unwanted.col <- which(names(b1) == "Genu_spec")
b1 <- b1[,-unwanted.col]

# Check for mispelled species names
all.dat <- rbind(chaut, a1, b1, c1)
sort(unique(as.character(all.dat$Species)))

# Correct mistakes in species names
all.dat$Species <- gsub(" $", "", all.dat$Species) # Remove any spaces accidentally
                                                   # added at end of species names
all.dat$Species <- sub("Cratypedes neglectus", "Cratypledes neglectus", all.dat$Species)
all.dat$Species <- sub("Melanoplus bivitattus", "Melanoplus bivatattus", all.dat$Species)
all.dat$Species <- sub("Melanoplus bivattatus", "Melanoplus bivatattus", all.dat$Species)
all.dat$Species <- sub("Melanoplus bivittatus", "Melanoplus bivatattus", all.dat$Species)

# Check that all mistakes are corrected
sort(unique(as.character(all.dat$Species)))

# Restrict ordinal date to June 1st to September 16th
all.dat <- all.dat[all.dat$OrdinalDate >= 152 &
                     all.dat$OrdinalDate <= 259, ]

# Subset to adults
all.dat <- all.dat[which(!is.na(all.dat$Stage6)),]

# Subset to common species
com.sp <- c("Camnula pellucida", "Melanoplus dawsoni", "Melanoplus sanguinipes", 
            "Melanoplus dodgei", "Melanoplus bivatattus", "Melanoplus packardii",
            "Melanoplus fasciatus")
all.dat <- all.dat[all.dat$Species %in% com.sp, ]

# Split data by site for analysis
a1 <- all.dat[all.dat$Site == "A1",]
b1 <- all.dat[all.dat$Site == "B1",]
c1 <- all.dat[all.dat$Site == "C1",]
chaut <- all.dat[all.dat$Site == "Chaut",]

##############################################################
# Calculating phenological overlap for all years within a site
##############################################################

# Select dataset
dat <- chaut

# Create vector of years included in sampling
years <- unique(dat$year)

# Create vector of species names in alphabetical order
sps <- sort(unique(dat$Species))

# Create an array with a cell for each species combination in each year
results <- array(dim = c(length(sps), length(sps), length(years)), 
                 dimnames = list(sps, sps, years))

# Loop through the years
for(year in 1:length(years)) {
  # Subset to 1 year
  sub.year <- dat[dat$year == years[year], ]
  
  # Loop through the focal species (sp1)
  for(sp1 in 1:length(sps)) {
    # Find ordinal dates when focal species present
    dates <- sub.year[sub.year$Species == sps[sp1] & sub.year$Stage6 > 0, 
                      "OrdinalDate"]
    # Subset to these dates
    sub.dat <- sub.year[sub.year$OrdinalDate %in% dates,]
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
        for(ord in 1:length(dates)) {
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
          # Calculate product of competitors and focals
          product <- competitors*focals
          # Divide this product by the total number of focal adults sampled in the 
          # season to make measures comparable between species of differing abundance
          measures[ord] <- product/tot.sp1
        }
        # Sum measures of overlap to obtain weighted mean and add to results matrix
        results[sp1, sp2, year] <- sum(measures)
      } # Close if statement loop
    } # Close competitor species loop
  } # Close focal species loop
} # Close years loop


############
# Make plots
############

# Load GDD data
chautGDD <- read.csv("../ChautauquaGDD.csv")
a1GDD <- read.csv("../A1GDD.csv")
b1GDD <- read.csv("../B1GDD.csv")
c1GDD <- read.csv("../C1GDD.csv")

# Define standard error function
se <- function(x) {
  sd(x, na.rm = TRUE)/sqrt(length(!is.na(x)))
}

# Calculate GDDs accumulated by September 16th (ordinal date = 259) for each year
GDDs <- rep(NA, times = 4)
#years1 <- c("X1959", "X1960", "X2006", "X2007", "X2008")
years1 <- c("X1959", "X1960", "X2007", "X2008")
for(a in 1:length(years1)) {
  GDDs[a] <- chautGDD[chautGDD$OrdinalDate == 259, years1[a]]
}


# All of species 1 as focal within a year
Species <- row.names(results[,,1])
Means <- apply(results[,,1], MARGIN = 1, FUN = se)
SEs <- apply(results[,,1], MARGIN = 1, FUN = mean, na.rm = TRUE)
upp <- Means + SEs
low <- Means - SEs
GDD <- rep(GDDs[1], times = length(low))
plotData <- data.frame(Species, Means, SEs, upp, low, GDD)

for(b in 2:4) {
  Species <- row.names(results[,,b])
  Means <- apply(results[,,b], MARGIN = 1, FUN = se)
  SEs <- apply(results[,,b], MARGIN = 1, FUN = mean, na.rm = TRUE)
  upp <- Means + SEs
  low <- Means - SEs
  GDD <- rep(GDDs[b], times = length(low))
  newData <- data.frame(Species, Means, SEs, upp, low, GDD)
  plotData <- rbind(plotData, newData)
}
View(plotData)

plot(y = plotData$Means, x = plotData$GDD, xlab = "GDDs accumulated by 9/16",
     main = "Chautauqua site",
     ylab = "PhenOver", col = 1:3, pch = 1:3)
legend(x = "top", legend = row.names(results[,,1]), col = 1:3, 
       cex = 0.8, pch = 1:3, bty = "n")


######################################################
# Determining whether species are early or late season
######################################################

# Select dataset
dat <- b1

# Create vector of years included in sampling
years <- unique(dat$year)

# Create vector of species names in alphabetical order
sps <- sort(unique(dat$Species))

for(i in 1:length(sps)) {
  dat1 <- dat[dat$Species == sps[i], ]
  plot(x = dat1$OrdinalDate, y = dat1$Stage6)
}

sps

