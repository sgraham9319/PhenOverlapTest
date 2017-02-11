##################################
# Quantifying phenological overlap
##################################

# Load data
chaut.dat <- read.csv("Chautdata.csv")

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

# First try calculation with just one species pair (column 10 in sp.combs)

# Subset to the 2 species of interest
chaut <- chaut[chaut$Species == "Melanoplus dawsoni" | 
                 chaut$Species == "Melanoplus sanguinipes", ]
# Find ordinal dates when focal species present
dates <- chaut[chaut$Species == "Melanoplus dawsoni" & chaut$Stage6 > 0, 
               "OrdinalDate"]

# Subset to these dates
chaut <- chaut[chaut$OrdinalDate %in% dates,]

measures <- rep(NA, times = length(dates))
tot.sp1 <- tapply(chaut$Stage6, chaut$Species, FUN = sum)["Melanoplus dawsoni"]
for (i in 1:length(dates)) {
  sp1 <- chaut[chaut$OrdinalDate == dates[3] & chaut$Species == "Melanoplus dawsoni",
               "Stage6"]
  sp2 <- chaut[chaut$OrdinalDate == dates[3] & chaut$Species == "Melanoplus sanguinipes",
               "Stage6"]
  if(length(sp2) == 0){
    sp2 <- 0
  }
  ratio <- sp2/sp1
  measures[3] <- ratio*(sp1/tot.sp1)
}
phenOver <- sum(measures)
phenOver
