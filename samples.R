library(magrittr)
library(dplyr)
library(purrr)
library(stringr)
library(openxlsx)
library(readr)
library(vadr)

batches.v4 <- unique(targets.v4.reduce$Batch) %>% paste(collapse = "|")
targets.v4a.only <- filter(targets.v4a.reduced, !grepl(batches.v4, Batch))

targets.combined <- rbind(targets.v1.reduce, targets.v2.reduce, targets.v3.reduce, targets.v4.reduce, targets.v4a.only)
targets.ucsf <- filter(targets.combined, !str_detect(Code, "UCLA"))

targets.ucsf$Code %<>% str_replace("\\_2", "")
targets.ucsf$Code <- str_replace(targets.ucsf$Code, "\\_2", "")

#PIDN.temp <- str_replace(targets.ucsf$Code, "\\_2", "") 
#PIDN.duped <- PIDN.temp[duplicated(PIDN.temp)] %>% paste(collapse = "|")
#duped <- filter(targets.ucsf, grepl(PIDN.duped, Code)) %>% arrange(Code)
write_csv(targets.ucsf, "targets.ucsf.csv")

needed.PIDNs <- read.xlsx("~/Downloads/60arrayed.xlsx", colNames = FALSE)$X1
missing.PIDNs <- needed.PIDNs[is.na(match(needed.PIDNs, targets.ucsf$Sample))]
found.PIDNs <- needed.PIDNs[!is.na(match(needed.PIDNs, targets.ucsf$Sample))] %>% paste(collapse = "|")
needed.targets <- filter(targets.ucsf, grepl(found.PIDNs, Sample))

rna <- read.xlsx("~/Documents/_Research/Dementia Project/phenotypedata/dan_dementia_rna_20151201T161608.xlsx", detectDates = TRUE)
test.group <- by(rna, factor(rna$PIDN), select, Sample.Num, Date.Received, Date.Processed)

test.merge <- merge()

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

load("./workspace.RData")

rna[is.na(rna$Notes),]$Notes <- "None"
rna.nobiomarker <- filter(rna, !(grepl("Biomarker", Notes)))

rna.count <- by(rna.nobiomarker, factor(rna.nobiomarker$PIDN), nrow) %>% reduce(c)
rna.df <- data.frame(PIDN = levels(factor(rna.nobiomarker$PIDN)), RNA = rna.count) %>% filter(!grepl("Easton", PIDN)) %>% arrange(desc(RNA)) 
write.xlsx(rna.df, "rna.counts.xlsx")

rna.twoplus <- filter(rna.df, RNA > 2)$PIDN %>% match.exact

targets.twoplus <- filter(targets.ucsf, grepl(rna.twoplus, Code))
rna.table.twoplus <- filter(rna, grepl(rna.twoplus, PIDN))
write.xlsx(targets.twoplus, "targets.twoplus.xlsx")
write.xlsx(rna.table.twoplus, "rna.twoplus.xlsx") 
