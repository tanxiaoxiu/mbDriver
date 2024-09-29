rm(list = ls(all = TRUE)) 
library(dplyr)

### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

setwd("~/mbDriver/real_data/other/p15/UC")
data_UC <- read.delim("data_UC.txt",  sep = '\t', check.names = FALSE)

threshold <- 0.5 * nrow(data_UC)
filtered_data <- data_UC %>%
  dplyr::select(Row.names, 5:ncol(data_UC)) %>%
  select_if(~sum(. != 0) > threshold)
filtered_data <- bind_cols(data_UC %>% dplyr::select(1:4), filtered_data)
filtered_data <- filtered_data[,-5]

####TOP15
data_UC_top <- filtered_data[,1:19]

####################
####################H1
H1 <- data_UC_top %>% filter(Group %in% c("H1", "H2"))
H1_sorted <- H1 %>%
  arrange(Subject, Time)


time_sets <- H1_sorted %>%
  group_by(Subject) %>%
  summarise(Times = list(unique(Time)), .groups = 'drop')
common_times <- Reduce(intersect, time_sets$Times)
H1_common <- H1_sorted %>%
  filter(Time %in% common_times)

H1_common_sorted <- H1_common %>%
  mutate(subjectID = as.integer(factor(Subject, levels = unique(Subject)))) %>%
  relocate(subjectID, .after = names(H1_common)[2]) %>% 
  mutate(sampleID = row_number()) %>% 
  relocate(sampleID, .after = names(H1_common)[1])
write.table(H1_common_sorted,file = "H/H_meta.txt",row.names = F,col.names = T, sep = "\t",quote = F)

H1_counts <- H1_common_sorted %>%
  dplyr::select(-c(1, 3:6))
H1_counts <-t(H1_counts)
write.table(H1_counts,file = "H/counts.txt",row.names = T,col.names = F, sep = "\t",quote = F)

meta <- H1_common_sorted %>%
  dplyr::select(c(2, 4:6))
meta$isIncluded <- "1"

meta <- meta %>%
  mutate(perturbID = ifelse(Group == "H1", 0, ifelse(Group == "H2", 1, NA)))

meta$exptblock <- "1"
meta$intv <- "0"

metadata <- meta %>%
  rename(measurementID = Time) 
metadata <- metadata %>%
  relocate(isIncluded, .after = names(metadata)[1])

metadata <- metadata %>%
  select(-Group)

write.table(metadata,file = "H/metadata.txt",row.names = F,col.names = T, sep = "\t",quote = F)


####################UC1
UC1 <- data_UC_top %>% filter(Group %in% c("UC1", "UC2"))
UC1_sorted <- UC1 %>%
  arrange(Subject, Time)


time_sets <- UC1_sorted %>%
  group_by(Subject) %>%
  summarise(Times = list(unique(Time)), .groups = 'drop')
common_times <- Reduce(intersect, time_sets$Times)
UC1_common <- UC1_sorted %>%
  filter(Time %in% common_times)

UC1_common_sorted <- UC1_common %>%
  mutate(subjectID = as.integer(factor(Subject, levels = unique(Subject)))) %>%
  relocate(subjectID, .after = names(UC1_common)[2]) %>% 
  mutate(sampleID = row_number()) %>% 
  relocate(sampleID, .after = names(UC1_common)[1])
write.table(UC1_common_sorted,file = "UC/UC_meta.txt",row.names = F,col.names = T, sep = "\t",quote = F)

UC1_counts <- UC1_common_sorted %>%
  dplyr::select(-c(1, 3:6))
UC1_counts <-t(UC1_counts)
write.table(UC1_counts,file = "UC/counts.txt",row.names = T,col.names = F, sep = "\t",quote = F)

meta <- UC1_common_sorted %>%
  dplyr::select(c(2, 4:6))
meta$isIncluded <- "1"

meta <- meta %>%
  mutate(perturbID = ifelse(Group == "UC1", 0, ifelse(Group == "UC2", 1, NA)))

meta$exptblock <- "1"
meta$intv <- "0"

metadata <- meta %>%
  rename(measurementID = Time) 
metadata <- metadata %>%
  relocate(isIncluded, .after = names(metadata)[1])

metadata <- metadata %>%
  select(-Group)
write.table(metadata,file = "UC/metadata.txt",row.names = F,col.names = T, sep = "\t",quote = F)


