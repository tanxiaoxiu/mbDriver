rm(list = ls(all = TRUE)) 
library(dplyr)
library(rhdf5)
### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

setwd("~/mbDriver/real_data/other/p15/Fiber")
data_fiber_diet <- read.delim("data_fiber.txt",  sep = '\t', check.names = FALSE)

threshold <- 0.8 * nrow(data_fiber_diet)
filtered_data <- data_fiber_diet %>%
  dplyr::select(Row.names, 6:ncol(data_fiber_diet)) %>%
  select_if(~sum(. != 0) > threshold)
filtered_data <- bind_cols(data_fiber_diet %>% dplyr::select(1:5), filtered_data)

filtered_data <- filtered_data[,-6]

####TOP15
data_fiber_diet_top <- filtered_data[,1:20]

####################
####################Control
Con <- data_fiber_diet_top %>% filter(Group == "Control")
Con_sorted <- Con %>%
  arrange(Subject, Time) %>%
  mutate(subjectID = as.integer(factor(Subject, levels = unique(Subject)))) %>%  
  relocate(subjectID, .after = names(Con_sorted)[3]) %>% 
  mutate(sampleID = row_number()) %>% 
  relocate(sampleID, .after = names(Con_sorted)[1])
write.table(Con_sorted,file = "Control/Con_meta.txt",row.names = F,col.names = T, sep = "\t",quote = F)

Con_counts <- Con_sorted %>%
  dplyr::select(-c(1, 3:7))
Con_counts <-t(Con_counts)
write.table(Con_counts,file = "Control/counts.txt",row.names = T,col.names = F, sep = "\t",quote = F)

meta <- Con_sorted %>%
  dplyr::select(c(2, 5:6))
meta$isIncluded <- "1"
meta$perturbID <- "0"
meta$exptblock <- "1"
meta$intv <- "0"

metadata <- meta %>%
  rename(measurementID = Time) 
metadata <- metadata %>%
  relocate(isIncluded, .after = names(metadata)[1])
write.table(metadata,file = "Control/metadata.txt",row.names = F,col.names = T, sep = "\t",quote = F)


####################RS
Rs <- data_fiber_diet_top %>% filter(Group == "Resistant starch")
Rs_sorted <- Rs %>%
  arrange(Subject, Time) %>%
  mutate(subjectID = as.integer(factor(Subject, levels = unique(Subject))))

Rs_sorted <- Rs_sorted %>% 
relocate(subjectID, .after = names(Rs_sorted)[3]) %>% 
mutate(sampleID = row_number()) %>% 
relocate(sampleID, .after = names(Rs_sorted)[1])
write.table(Rs_sorted,file = "Rs/Rs_meta.txt",row.names = F,col.names = T, sep = "\t",quote = F)

Rs_counts <- Rs_sorted %>%
  dplyr::select(-c(1, 3:7))
Rs_counts <-t(Rs_counts)
write.table(Rs_counts,file = "Rs/counts.txt",row.names = T,col.names = F, sep = "\t",quote = F)

meta <- Rs_sorted %>%
  dplyr::select(c(2, 5:6))
meta$isIncluded <- "1"
meta$perturbID <- "0"
meta$exptblock <- "1"
meta$intv <- "0"

metadata <- meta %>%
  rename(measurementID = Time) 
metadata <- metadata %>%
  relocate(isIncluded, .after = names(metadata)[1])
write.table(metadata,file = "Rs/metadata.txt",row.names = F,col.names = T, sep = "\t",quote = F)


####################IN
In <- data_fiber_diet_top %>% filter(Group == "Inulin")
In_sorted <- In %>%
  arrange(Subject, Time) %>%
  mutate(subjectID = as.integer(factor(Subject, levels = unique(Subject))))
In_sorted <- In_sorted %>% 
  relocate(subjectID, .after = names(In_sorted)[3]) %>% 
  mutate(sampleID = row_number()) %>% 
  relocate(sampleID, .after = names(In_sorted)[1])
write.table(In_sorted,file = "In/In_meta.txt",row.names = F,col.names = T, sep = "\t",quote = F)

In_counts <- In_sorted %>%
  dplyr::select(-c(1, 3:7))
In_counts <-t(In_counts)
write.table(In_counts,file = "In/counts.txt",row.names = T,col.names = F, sep = "\t",quote = F)

meta <- In_sorted %>%
  dplyr::select(c(2, 5:6))
meta$isIncluded <- "1"
meta$perturbID <- "0"
meta$exptblock <- "1"
meta$intv <- "0"
metadata <- meta %>%
  rename(measurementID = Time) 
metadata <- metadata %>%
  relocate(isIncluded, .after = names(metadata)[1])
write.table(metadata,file = "In/metadata.txt",row.names = F,col.names = T, sep = "\t",quote = F)


