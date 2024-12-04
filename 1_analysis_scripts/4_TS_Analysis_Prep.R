## DIVIDING SAMPLES ------------------------------------------------------------

rm(list = ls(all.names = T))

# Setting Format for Targets for TS analysis

TS_target <- read.delim("2_rawdata/2_targets/TEMP_targets.txt", sep = ",")
names(TS_target)[names(TS_target) == 'X'] <- 'samples'

TS_target$replicate <- c("infected_1","infected_2","infected_3",
                        "infected_1","infected_2","infected_3",
                        "infected_1","infected_2","infected_3",
                        "mock_1","mock_2","mock_3",
                        "mock_1","mock_2","mock_3",
                        "mock_1","mock_2","mock_3")
TS_target$replicate <- as.factor(TS_target$replicate)

write.csv(TS_target, "2_rawdata/4_TS_Data/sample_file.csv", row.names = F)

# Setting Format for the Counts for TS Analysis

TS_counts <- read.delim("2_rawdata/1_counts/TEMP_counts.txt", sep = ",", row.names = 1)
rows <- rownames(TS_counts)
ids <- colnames(TS_counts)

for (i in ids){
  dfname <- paste0("2_rawdata/4_TS_Data/raw_counts_TS/",i,".counts")
  dfsubset <- data.frame(TS_counts[,i])
  dfsubset <- cbind(rows, dfsubset)
  names(dfsubset) <- NULL
  write.table(dfsubset, dfname, row.names = F, quote = F, sep = "\t")
}

# Now we should have inside of 4_TS_Data/...
#       * A sample_file.csv - this contains the metadata and is COMMA separated
#       * A file called raw_counts_TS with all our SAMPLEID.counts files inside and are TAB separated
#
# These files can now be copied and used within the TiSA pipeline (https://github.com/shelleydan/TimeSeriesAnalysis)