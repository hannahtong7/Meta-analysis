
library(dplyr)

#Filter cocolithophores for accepted files only
setwd("~/Desktop/Meta-analysis/Cocos")
cocos <- read.csv("Cocos.csv_Screened.csv", header=TRUE)
cocos_accepted <- cocos %>% 
  filter(Screen == "Accept")
head(cocos_accepted)
tail(cocos_accepted)
nrow(cocos_accepted)
#8 accepted files

#Filter diatoms for accepted files only
setwd("~/Desktop/Meta-analysis/Diatoms")
Diatoms <- read.csv("Diatoms.csv_Screened.csv", header=TRUE)
Diatoms_accepted <- Diatoms %>% 
  filter(Screen == "Accept")
head(Diatoms_accepted)
tail(Diatoms_accepted)
nrow(Diatoms_accepted)
#12 accepted files

#Filter cyanobacteria for accepted files only
setwd("~/Desktop/Meta-analysis/Cyanobacteria")
Cyano <- read.csv("Cyano.csv_Screened.csv", header=TRUE)
Cyanos_accepted <- Cyano %>% 
  filter(Screen == "Accept")
head(Cyanos_accepted)
tail(Cyanos_accepted)
nrow(Cyanos_accepted)
#11 accepted files

#Filter Dinoflagellates for accepted files only
setwd("~/Desktop/Meta-analysis/Dinos")
Dinos <- read.csv("Dinos.csv_Screened.csv", header=TRUE)
Dinos_accepted <- Dinos %>% 
  filter(Screen == "Accept")
head(Dinos_accepted)
tail(Dinos_accepted)
nrow(Dinos_accepted)
#14 accepted files


library(dplyr)

#Add a new column identifying taxa
cocos_accepted <- cocos_accepted %>% mutate(Taxa = "Coccolithophores")
Diatoms_accepted <- Diatoms_accepted %>% mutate(Taxa = "Diatoms")
Cyanos_accepted <- Cyanos_accepted %>% mutate(Taxa = "Cyanobacteria")
Dinos_accepted <- Dinos_accepted %>% mutate(Taxa = "Dinoflagellates")

#Combine all accepted datasets into one data frame
all_accepted <- bind_rows(cocos_accepted, Diatoms_accepted, Cyanos_accepted, Dinos_accepted)

#Save the combined dataset
write.csv(all_accepted, "All_Accepted_Studies.csv", row.names = FALSE)

#Check the result
head(all_accepted)
tail(all_accepted)

#Save combined data frame as CSV
write.csv(all_accepted, "~/Desktop/Meta-analysis/accepted_papers.csv", row.names = FALSE)

