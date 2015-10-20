## All screened compound data can be downloaded from https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-60+Growth+Inhibition+Data
## the GI50 Data (Sept 2014) release contains a file CANCER60GI50.csv
## We will use this file to extract a few compounds of interest for prediction models

Screen <- read.csv(file = "./Data/CANCER60GI50.csv", header = TRUE, colClasses = c("numeric", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "character"), strip.white = TRUE)

# restrict to approved and investigational agents subset, and only log high concentration of -4
AOD_IOA <- read.csv("AOD_IOA_inSep2014.csv")

## use only the subset of NSC values and only the -4 high concentration screens
Screen_sub <- Screen[which(Screen$NSC %in% AOD_IOA$NSC & Screen$LCONC == -4), ]

OUT <- merge(Screen_sub[, c("NSC", "PANEL", "CELL", "NLOGGI50", "INDN", "TOTN")], AOD_IOA[, c("NSC", "NAME")])
colnames(OUT) <- sub("NAME", "DRUGNAME", colnames(OUT))

write.csv(OUT, "AOD_IOA_GI50.csv")
