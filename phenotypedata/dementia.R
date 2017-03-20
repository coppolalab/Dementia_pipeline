library(openxlsx)
library(readr)
library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(purrr)
library(functional)
library(lubridate)
library(vadr)
library(broom)
library(ggplot2)
library(Cairo)
library(RColorBrewer)

load("./workspace.RData")
rna.match <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\rna.match.match.xlsx")


Sys.setenv(R_ZIPCMD= "C:\\Rtools\\bin\\zip")

ringman <- read.xlsx("C:\\Users\\Youngjun Park\\Downloads\\parkyj_ringman_20160112T162623.xlsx")
ucsf <- read.xlsx("C:\\Users\\Youngjun Park\\Downloads\\parkyj_ucsf_20160112T171146.xlsx")

targets.v4 <- read_csv("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\targets_v4_updated_052012.csv")
targets.v1 <- read_tsv("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\targets_platform_v1.tsv") %>% select(-pID)
targets.v2 <- read_tsv("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\targets_platform_v2.tsv")
targets.v3 <- read_tsv("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\targets_platform_v3.tsv")
targets.v4a <- read_csv("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\targets_v4_additional.csv")

colnames(targets.v4)

targets.v1.reduce <- select(targets.v1, -pID)
select(targets.v1, -pID)
#select sorts the columns for you too

temp.match <- match(colnames(targets.v3), colnames(targets.v4))
missing <- colnames(targets.v3)[is.na(temp.match)]
present <- colnames(targets.v3)[!is.na(temp.match)] #%>% paste(collapse = ",")
present.n <- present[!str_detect(present, "^Status$")]

colnames(targets.v1)[12] <- "Platform"

targets.v4a.reduce <- select_(targets.v4a, .dots = present.n)
targets.v1.reduce <- select_(targets.v1, .dots = present.n)
targets.v2.reduce <- select_(targets.v2, .dots = present.n)
targets.v3.reduce <- select_(targets.v3, .dots = present.n)
targets.v4.reduce <- select_(targets.v4, .dots = present.n)
targets.merged <- rbind(targets.v1.reduce, targets.v2.reduce, targets.v3.reduce, targets.v4.reduce, targets.v4a.reduce)
colnames(targets.merged)[2] <- "PIDN"

#rbind to combine rows to put in one table join command, merge command for
#rowwise merging code(targets) and PIDN(patient)
#change code in targets to PIDN
ucsf.merged <- rbind(targets.merged, ucsf)

etargets.1 <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\2011-211_Sample Key.xlsx")
etargets.2 <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\2012-067_Sample Key.xlsx")
etargets.3 <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\2014-050 sample key CORRECTED.xlsx")
etargets.merged <- rbind(etargets.1, etargets.2, etargets.3)


duplicated(targets.merged$Slide)
#table of logic:duplicated or not
dup <- duplicated(targets.v4a.reduce[5])
#make a subset with no duplicates
targets.v4a.r.u <- subset(targets.v4a.reduce, !dup)
#check for duplicates
duplicated(rbind(targets.v4a.r.u[5], targets.v1.reduce[5]))

dup2 <- duplicated(rbind(targets.v4a.r.u[5], targets.v4.reduce[5]))
targets.v4.merged <- rbind(targets.v4a.r.u, targets.v4.reduce)
targets.v4.merged.unique <- subset(targets.v4.merged, !dup2)

#join the targets and the ucsf data together
targets.merged.unique <- rbind(targets.v1.reduce, targets.v2.reduce, targets.v3.reduce, targets.v4.merged.unique)
colnames(targets.merged.unique)[2] <- "PIDN"
ucsf.merged.unique <- inner_join(targets.merged.unique, ucsf)

#make the "Slide" consist of all unique IDs and merge.
targets.4 <- targets.v4a$Slide[!is.na(targets.v4a$position)]
targets.4p <- targets.v4a$position[!is.na(targets.v4a$position)]
targets.j <- paste(targets.4, targets.4p, sep = "_")
targets.v4a$Slide[!is.na(targets.v4a$position)] <- targets.j

targets.v4a.reduced <- select_(targets.v4a, .dots = present.n)
dup <- duplicated(c(targets.v4a.reduced, targets.v4.reduce))
targets.v4.merged <- rbind(targets.v4a.reduced, targets.v4.reduce)
targets.4.merged.unique <- subset(targets.v4.merged, !dup)

#grel, paste, filter out samples
batches.v4 <- unique(targets.v4.reduce$Batch) %>% paste(collapse = "|")
batches.key <- paste("^", batches.v4, sep = "") %>% paste("$", sep = "") %>% paste(collapse = "|")
targets.v4a.only <- filter(targets.v4a.reduced, !grepl(batches.key, Batch))
targets.merged <- rbind(targets.v1.reduce, targets.v2.reduce, targets.v3.reduce, targets.v4.reduce, targets.v4a.only)

#remove the sameples with Code that contains UCLA and replace "_2"
duplicated(targets.merged$Slide) %>% which() %>% length()
targets.merged.ucsf <- filter(targets.merged, !str_detect(targets.merged$Code, "UCLA"))
targets.merged.ucsf$Code %<>% str_replace("\\_2", "")
#factor
#summary

PIDNs <- targets.merged.ucsf$Code
pidn.key <- paste("^", PIDNs, sep = "") %>% paste("$", sep ="") %>% paste(collapse = "|")
grepl(pidn.key, ucsf$PIDN)
ucsf.match <- filter(ucsf, grepl(pidn.key, PIDN))
ucsf.nomatch <- filter(ucsf, !grepl(pidn.key, PIDN))

rna <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\dan_dementia_rna_20151201T161608.xlsx")

rna.filtered <- filter(rna, !str_detect(rna$PIDN, "Easton"))

pidn <- ucsf.match$PIDN
pidn.key <- paste("^", pidn, sep = "") %>% paste("$", sep = "") %>% paste(collapse = "|")
rna.match.match <- filter(rna.filtered, grepl(pidn.key, PIDN))
rna.match.nomatch <- filter(rna.filtered, !grepl(pidn.key, PIDN))

#pidn <- ucsf.nomatch$PIDN
#pidn.key <- paste("^", pidn, sep = "") %>% paste("$", sep = "") %>% paste(collapse = "|")
#rna.nomatch.match <- filter(rna.filtered, grepl(pidn.key, PIDN))
#rna.nomatch.nomatch <- filter(rna.filtered, !grepl(pidn.key, PIDN))

write.xlsx(rna.match.match, "rna.match.match.xlsx")
write.xlsx(rna.match.nomatch, "rna.match.nomatch.xlsx")
write.xlsx(rna.nomatch.nomatch, "rna.nomatch.nomatch.xlsx")
write.xlsx(rna.nomatch.match, "rna.nomatch.match.xlsx")

typeof()

split(rna.match.match, rna.match.match$PIDN) %>% map_dbl(nrow)
data.frame(PIDN=name(number), number=number)


filter(duplicates, duplicates$number > 2)
morethantwo <- filter(duplicates, duplicates$number > 2)
morethantwo.l <- morethantwo$names.number.
key <- paste("^", morethantwo.l, sep = "") %>% paste("$", sep = "") %>% paste(collapse = "|")
rna.morethantwo <- filter(rna.match.match, grepl(key, PIDN))


split <- split(ucsf.match, ucsf.match$Categorized.Diagnosis)
#compose: kind of like pipeline for functions, unlist: turns it into array
gender <- map(split, select, Gender) %>% map(unlist)
unique.values <- reduce(gender, c) %>% unique
gender.summary <- map(gender, Compose(factor, summary))

add.missing <- function(missing.array, all.values)
{
  if (length(missing.array) < length(all.values))
  {
     paste.exact <- paste %<<<% "^" %<<% c("$", sep = "")   
     found.key <- map_chr(names(missing.array), paste.exact) %>% paste(collapse = "|")
     missing.key <- !grepl(found.key, all.values)
     add.values <- all.values[missing.key]
     add.vector <- rep(0, length(add.values))
     names(add.vector) <- add.values
     new.vector <- c(missing.array, add.vector) %>% sort
     return(new.vector)
  }
  else
  {
    return(missing.array)
  }
}

new.list <- map(gender.summary, add.missing, unique.values) 

#get every row in gender.summary to have Male, Female, and Unknown
#then use reduce(gender.summary, rbind) to make it into a table

#get the list of names of diagnosis and for each diagnosis filter for MALE, FEMALE, UNKNOWN. if there is none, add one

diagnosis <- names(new.list)
new.list.reduce <- reduce(new.list, rbind)
Total <- rowSums(new.list.reduce)
new.df <- data.frame(diagnosis, new.list.reduce, Total)
write.xlsx(new.df, "patient-totals.xlsx")

age <- rna.match$AgeAtDraw
pidn <- rna.match$PIDN
age.df <- data.frame(pidn, age)
split <- split(age.df, age.df$pidn) %>% map_dbl(nrow)
split <- map(split, unlist)
split[split > 1]

age.n <-by(age.df, age.df$pidn, select, age)
#apply unique to each value in age.n
laply(age.n, unique)
#sapply is exactly like map()....
age.nn <- sapply(age.n, unique)
age.nn <- map(age.n, unlist)

#map_int makes it into a vector
age.count <- map_int(age.nn, length)
#separate the ones with only one value for each pidn in age.nn
age.count.df <- data.frame(names(age.count), age.count)
age.onedrawdate <- filter(age.count.df, age.count.df$age.count == 1)
#list of patients with more than one draw
age.morethanonedraw <- filter(age.count.df, age.count.df$age.count != 1)


age.nn.v <- unlist(age.nn)
age.names <- age.onedrawdate$names.age.count.
age.names.key <- paste(age.names, collapse = "|")

age.nn.df <- data.frame(names(age.nn.v), age.nn.v)
age.new.df <- filter(age.nn.df, grepl(age.names.key, age.nn.df$PIDN))

age.test <- filter(age.nn.df, !grepl(age.names.key, age.nn.df$PIDN))

#inner_join(age.onedrawdate, age.new.df, by = )


#on the microarray table, merge sex and diagnosis and age at drawn
#targets.merged.ucsf
colnames(age.new.df)[1] <- "Code"
age.names <- age.new.df$Code
targets.filtered <- filter(targets.merged.ucsf, is.element(targets.merged.ucsf$Code, age.names))
targets.age <- inner_join(targets.filtered, age.new.df, by = "Code")

####################################################################################
####### patients with more than one age ############################################
####################################################################################
colnames(age.df)[1] <- "PIDN"
colnames(age.nn.df)[1] <- "PIDN"
ucsf.morethanoneage <- filter(ucsf.match, is.element(ucsf.match$PIDN, age.morethanonedraw$PIDN))
age.morethanoneage.df <- filter(age.df, is.element(age.df$PIDN, ucsf.morethanoneage$PIDN))
patients.morethanoneage <- inner_join(ucsf.morethanoneage, age.df)
colnames(patients.morethanoneage)[17] <- "Age"
patients.morethanoneage.data <- select(patients.morethanoneage, PIDN, Gender, Age)
patients.key <- patients.morethanoneage.data$PIDN %>% unique()
temp <- select(patients.morethanoneage.data, PIDN, Age)
temp0 <- filter(temp, temp$Age != -1)
temp1 <- split(temp0, temp0$PIDN)
temp2 <- map(temp1, select, Age)
temp3 <- map(temp2, min)
temp4 <- unlist(temp3)
patients.plusage <- data.frame(patients.key, temp4)
colnames(patients.plusage)[1] <- "PIDN"
colnames(patients.plusage)[2] <- "Age"
ucsf.plusage <- filter(ucsf.morethanoneage, !duplicated(ucsf.morethanoneage$PIDN))
patients.ageage <- inner_join(ucsf.plusage, patients.plusage)
patients.age <- select(patients.ageage, PIDN, Gender, Age)
diagnosis.dupage <- filter(diagnosis.df, is.element(diagnosis.df$PIDN, patients.age$PIDN))
diagnosis.dupage2 <- filter(diagnosis.dupage, diagnosis.dupage$DX.Type != "MOLECULAR")
patients.age.diagnosis <- inner_join(patients.age, diagnosis.dupage2)
patients.2.data <- select(patients.age.diagnosis, PIDN, Gender, Age, Diagnosis)

##############################################################################################################
ucsf.match.filtered <- filter(ucsf.match, is.element(ucsf.match$PIDN, age.names))
colnames(targets.age)[2] <- "PIDN"
microarray.df <- inner_join(targets.age, ucsf.match.filtered, by = "PIDN")

microarray.dups <- filter(microarray.df, duplicated(microarray.df$PIDN))

duplicated.pidn <- microarray.dups$PIDN
microarray.duplicates <- filter(microarray.df, is.element(microarray.df$PIDN, duplicated.pidn))
rna.duplicates <- filter(rna.match, is.element(rna.match$PIDN, duplicated.pidn))

write.xlsx(microarray.duplicates, "microarray.duplicates.xlsx")
write.xlsx(rna.duplicates, "rna.duplicates.xlsx")

microarray.morethantworna <- filter(microarray.df, RNA > 2)
microarray.plustwo <- filter(microarray.morethantworna, !duplicated(PIDN))

write.xlsx(microarray.plustwo, "microarray.plustwo.xlsx")
morethantwo.pidn <- microarray.plustwo$PIDN
rna.plustwo <- filter(rna.match, is.element(rna.match$PIDN, morethantwo.pidn))
write.xlsx(rna.plustwo, "rna.plustwo.xlsx")

age.tobefixed <- filter(microarray.df, microarray.df$Age == -1 | microarray.df$Age > 150)
write.xlsx(age.tobefixed.rna, "nodrawdate.xlsx")
write.xlsx(age.tobefixed, "age.tobefixed.xlsx")

#to change date
install.packages("lubridate")
library(lubridate)

#get average age for each diagnosis
age.exclude <- age.tobefixed$PIDN
microarray.filtered <- filter(microarray.df, !is.element(microarray.df$PIDN, age.exclude))
map2_lgl(toupper(microarray.filtered$Dx_Status), toupper(microarray.filtered$Categorized.Diagnosis), identical)
unmatched.diagnosis <- filter(microarray.filtered, !map2_lgl(toupper(microarray.filtered$Dx_Status), toupper(microarray.filtered$Categorized.Diagnosis), identical))

unknown1 <- toupper(unmatched.diagnosis$Categorized.Diagnosis) != "UNCATEGORIZED" & toupper(unmatched.diagnosis$Categorized.Diagnosis) != "UNSPECIFIED" & toupper(unmatched.diagnosis$Dx_Status) == "UNKNOWN"
unmatched.diagnosis[unknown1,]$Dx_Status <- unmatched.diagnosis[unknown1,]$Categorized.Diagnosis

unknown2 <- (toupper(unmatched.diagnosis$Categorized.Diagnosis) == "UNCATEGORIZED" | toupper(unmatched.diagnosis$Categorized.Diagnosis) == "UNSPECIFIED") & toupper(unmatched.diagnosis$Dx_Status) != "UNKNOWN"
unmatched.diagnosis[unknown2,]$Categorized.Diagnosis <- unmatched.diagnosis[unknown2,]$Dx_Status

unmatched.diagnosis.unmatchcheck <- filter(unmatched.diagnosis, !map2_lgl(toupper(unmatched.diagnosis$Dx_Status), toupper(unmatched.diagnosis$Categorized.Diagnosis), identical))
filter(unmatched.diagnosis.unmatchcheck, grepl("UNCATEGORIZED|UNSPECIFIED", Categorized.Diagnosis))$Dx_Status
diagnosis.conflict <- filter(unmatched.diagnosis.unmatchcheck, Dx_Status != "Unknown")
diagnosis.conflict$Categorized.Diagnosis %<>% str_replace(" ", "")

diagnosis.conflict1 <- filter(diagnosis.conflict, !map2_lgl(toupper(diagnosis.conflict$Dx_Status), toupper(diagnosis.conflict$Categorized.Diagnosis), identical))
write.xlsx(diagnosis.conflict1, "diagnosis.conflict.xlsx")

#diagnosis table
diagnosis.df <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\parkyj_diagnosis_20160222T122536.xlsx")
diagnosis.filtered <- filter(diagnosis.df, is.element(diagnosis.df$PIDN, microarray.filtered$PIDN))

#table with molecular diagnosis
diagnosis.molecular <- filter(diagnosis.filtered, diagnosis.filtered$DX.Type == "MOLECULAR")

#table with clinical diagnosis
diagnosis.clinical <- filter(diagnosis.filtered, diagnosis.filtered$DX.Type != "MOLECULAR")

#filter clinical diagnosis for duplicate PIDNs
diagnosis.dups <- filter(diagnosis.clinical, duplicated(diagnosis.clinical$PIDN))
diagnosis.duplicated <- filter(diagnosis.clinical, is.element(diagnosis.clinical$PIDN, diagnosis.dups$PIDN))
write.xlsx(diagnosis.duplicated, "diagnosis.duplicated.xlsx")

#diagnosis table without any duplicate PIDNs //made mistake so same as diagnosis.clinical
diagnosis.nodups <- filter(diagnosis.clinical, !is.element(diagnosis.clinical$PIDN, diagnosis.dups))
#microarray.notindiagnosis <- filter(microarray.filtered, !is.element(microarray.filtered$PIDN, diagnosis.filtered$PIDN))
diagnosis.withdrawn <- filter(diagnosis.filtered, diagnosis.filtered$Diagnosis == "Withdrawn" | diagnosis.filtered$Diagnosis == "Withdrew/Exclude")

#got rid of diagnosis that had Withdrawn
diagnosis.use <- filter(diagnosis.nodups, !is.element(diagnosis.nodups$PIDN, diagnosis.withdrawn$PIDN))
withdrawn.patients <- filter(diagnosis.nodups, is.element(diagnosis.nodups$PIDN, diagnosis.withdrawn$PIDN))
diagnosis.duplicated <- filter(diagnosis.use, duplicated(diagnosis.use$PIDN))

write.xlsx(diagnosis.use, "diagnosis.xlsx")
diagnosis.nodups <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\diagnosis.xlsx")



#filter microarray for the ones that is in diagnosis table and join with diagnosis table
microarray.n <- filter(microarray.filtered, is.element(microarray.filtered$PIDN, diagnosis.nodups$PIDN))
microarray.nn <- filter(microarray.n, !duplicated(microarray.n$PIDN))

microarray.nomatch <- filter(microarray.filtered, !is.element(microarray.filtered$PIDN, diagnosis.nodups$PIDN))
microarray.nomatch.merge <- inner_join(microarray.nomatch, diagnosis.nodups, "PIDN")

microarray.diagnosis <- inner_join(diagnosis.nodups, microarray.nn, "PIDN")
patients.1.data <- select(microarray.diagnosis, PIDN, Gender, Age, Diagnosis)
patients.data <- rbind(patients.1.data, patients.2.data)
by(patients.data, "Diagnosis")

######################################################################################################
### split by Diagnosis ###############################################################################
######################################################################################################

split <- split(patients.data, patients.data$Diagnosis)
#compose: kind of like pipeline for functions, unlist: turns it into array
gender <- map(split, select, Gender) %>% map(unlist)
unique.values <- reduce(gender, c) %>% unique
gender.summary <- map(gender, Compose(factor, summary))

##################################################
### use if there are missing values ##############
##################################################

add.missing <- function(missing.array, all.values)
{
  if (length(missing.array) < length(all.values))
  {
    paste.exact <- paste %<<<% "^" %<<% c("$", sep = "")   
    found.key <- map_chr(names(missing.array), paste.exact) %>% paste(collapse = "|")
    missing.key <- !grepl(found.key, all.values)
    add.values <- all.values[missing.key]
    add.vector <- rep(0, length(add.values))
    names(add.vector) <- add.values
    new.vector <- c(missing.array, add.vector) %>% sort
    return(new.vector)
  }
  else
  {
    return(missing.array)
  }
}

gender.list <- map(gender.summary, add.missing, unique.values)
###################################################

diagnosis.gender <- reduce(gender.list, rbind)
Total <- rowSums(diagnosis.gender)
gender.diagnosis <- data.frame(names(gender.list), diagnosis.gender, Total)
colnames(gender.diagnosis)[1] <- "Diagnosis"

###############################################################################################
######### get average age #####################################################################
###############################################################################################

split <- split(patients.data, patients.data$Diagnosis)
split.age <- map(split, select, Age)
age.mean <- map(split.age, colMeans)
diagnosis.age <- reduce(age.mean, rbind)
age.diagnosis <- data.frame(names(age.mean), diagnosis.age)
colnames(age.diagnosis)[1] <- "Diagnosis"
colnames(age.diagnosis)[2] <- "Mean Age"
patients.summary <- inner_join(gender.diagnosis, age.diagnosis)

write.xlsx(patients.summary, "PatientSummary.xlsx")




#ucsf.age <- inner_join(ucsf.match.filtered, age.new.df)
############################################################################################################
#####use microarray with one age. age sex and diagnosis, batch##############################################
#targets.age include only microarrays with one age##########################################################
############################################################################################################

targets.1 <- targets.v4a.batch$Slide[!is.na(targets.v4a.batch$position)]
targets.2 <- targets.v4a.batch$position[!is.na(targets.v4a.batch$position)]
targets.join <- paste(targets.1, targets.2, sep = "_")
targets.v4a.batch$Slide[!is.na(targets.v4a.batch$position)] <- targets.join

targets.v4a.batch <- filter(targets.v4a.batch, targets.v4a.batch$Batch > 33)
filter(targets.age, is.na(targets.age$Batch))
#fill in missing data by comparing two dataframes: FillIn
#devtools::source_gist("4959237")

targets <- select(targets.age, PIDN, Slide, Batch, Age)
temp <- filter(targets, is.na(targets.age$Batch))
temp0 <- select(temp, PIDN, Slide, Age)   #changed PIDN to Slide

temp.other <- filter(targets, !is.na(targets.age$Batch))
temp1 <- filter(targets.v4a.batch, is.element(targets.v4a.batch$Slide, temp$Slide))
temp2 <- select(temp1, Slide, Batch)
#colnames(temp2)[1] <- "PIDN"

temp3 <- map(temp2, as.character, temp2$Slide)
temp4 <- data.frame(temp3)


#temp$Batch <- temp2$Batch
temp.merged <- inner_join(temp0, temp4)

#temp0 <- filter(targets, !is.element(targets$PIDN, temp$PIDN))
targets.batch <- rbind(temp.other, temp.merged)

data.1 <- filter(patients.data, is.element(patients.data$PIDN, targets.batch$PIDN))
data.2 <- select(data.1, PIDN, Gender, Diagnosis)
microarray.data <- inner_join(targets.batch, data.2)
dups <- filter(microarray.data, duplicated(microarray.data$PIDN))
microarray.dup.pidn <- filter(microarray.data, is.element(microarray.data$PIDN, dups$PIDN))

write.xlsx(microarray.data, "microarray_data.xlsx")
write.xlsx(microarray.dup.pidn, "microarray_duplicate_pidn.xlsx")

temp <- filter(targets.batch, duplicated(targets.batch$PIDN))
tmp <- filter(targets.batch, is.element(targets.batch$PIDN, temp$PIDN))


###################################################################################################
######### add new microarray data #################################################################
###################################################################################################
targets.new <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\2015-9253A_Sample Key.xlsx")

######################
## add batch number ##
######################
targets.new$batch <- "37"

######################
### join slide #######
######################
targets.4 <- targets.new$general.array[!is.na(targets.new$genexstripe.controling.stripe)]
targets.4p <- targets.new$genexstripe.controling.stripe[!is.na(targets.new$genexstripe.controling.stripe)]
targets.j <- paste(targets.4, targets.4p, sep = "_")

targets.new$Slide[!is.na(targets.new$genexstripe.controling.stripe)] <- targets.j

colnames(targets.new)[8] <- "Batch"
targets.new$External.ID %<>% str_replace(" ", "") %>% str_replace("_.*$", "")

#filter out reKO stuff
targets.reko <- filter(targets.new, grepl("reKO", targets.new$External.ID))
targets.new <- filter(targets.new, !grepl("reKO", targets.new$External.ID))

#fix notations..
targets.new$External.ID %<>% str_replace("-.*$", "")f
colnames(targets.new)[5] <- "PIDN"
targets.new$Platform <- "v4"


###################################################################################################
######### add new rna data ########################################################################
###################################################################################################
rna.ucsf.new <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\new_ucsf.xlsx")          #3117
rna.noage <- filter(rna.ucsf.new, rna.ucsf.new$AgeAtDraw < 1 | rna.ucsf.new$AgeAtDraw > 150)          #43
rna.ucsf.validage <- filter(rna.ucsf.new, !is.element(rna.ucsf.new$PIDN, rna.noage$PIDN))             #3070
rna.ucsf.noage <- filter(rna.ucsf.new, is.element(rna.ucsf.new$PIDN, rna.noage$PIDN))                 #47
rna.ucsf.noage$AgeAtDraw <- NA
rna.use <- rbind(rna.ucsf.validage, rna.ucsf.noage)

###################################################################################################
######### add new patients data ###################################################################
###################################################################################################
patients.ucsf <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\patients_ucsf.xlsx")

###################################################################################################
######### add new diagnosis data ##################################################################
###################################################################################################
diagnosis.filtered <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\new_diagnosis.xlsx")
#fix the date:
#make a column with the formula =TEXT(K2,"yyyy-mm-dd")
#make third column and only paste the text
#paste third column into original

#table with molecular diagnosis
diagnosis.molecular <- filter(diagnosis.filtered, diagnosis.filtered$DX.Type == "MOLECULAR")

#table with clinical diagnosis
diagnosis.clinical <- filter(diagnosis.filtered, diagnosis.filtered$DX.Type != "MOLECULAR")

#filter clinical diagnosis for duplicate PIDNs
diagnosis.dups <- filter(diagnosis.clinical, duplicated(diagnosis.clinical$PIDN))
diagnosis.duplicated <- filter(diagnosis.clinical, is.element(diagnosis.clinical$PIDN, diagnosis.dups$PIDN))
write.xlsx(diagnosis.duplicated, "diagnosis.duplicated.xlsx")

#diagnosis table without any duplicate PIDNs //made mistake so same as diagnosis.clinical
diagnosis.nodups <- filter(diagnosis.clinical, !is.element(diagnosis.clinical$PIDN, diagnosis.dups))

#manually take out duplicate diagnosis (use the most recent diagnosis)
diagnosis.filtered.duplicated <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\diagnosis.duplicated.xlsx")
diagnosis.data <- rbind(diagnosis.nodups, diagnosis.filtered.duplicated)


####microarray.notindiagnosis <- filter(microarray.filtered, !is.element(microarray.filtered$PIDN, diagnosis.filtered$PIDN))
#diagnosis.withdrawn <- filter(diagnosis.filtered, diagnosis.filtered$Diagnosis == "Withdrawn" | diagnosis.filtered$Diagnosis == "Withdrew/Exclude")

###got rid of diagnosis that had Withdrawn
#diagnosis.use <- filter(diagnosis.nodups, !is.element(diagnosis.nodups$PIDN, diagnosis.withdrawn$PIDN))
#withdrawn.patients <- filter(diagnosis.nodups, is.element(diagnosis.nodups$PIDN, diagnosis.withdrawn$PIDN))
#diagnosis.duplicated <- filter(diagnosis.use, duplicated(diagnosis.use$PIDN))

#write.xlsx(diagnosis.use, "diagnosis.xlsx")
#diagnosis.nodups <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\diagnosis.xlsx")

###################################################################################################
######### merge microarray info ###################################################################
###################################################################################################
#first get diagnosis for each microarray
#diagnosis.targets.new <- filter(diagnosis.data, is.element(diagnosis.data$PIDN, targets.new$PIDN))
#targets.nodiagnosis <- filter(targets.new, !is.element(targets.new$PIDN, diagnosis.data$PIDN))
write.xlsx(targets.nodiagnosis, "targets_nodiagnosis.xlsx")


#diagnosis.temp <- select(diagnosis.targets.new, PIDN, Diagnosis)
#targets.new.diagnosis <- left_join(targets.new, diagnosis.temp)

#patients.temp <- select(patients.ucsf, PIDN, Gender)
#temp <- filter(patients.temp, is.element(patients.temp$PIDN, targets.new.diagnosis$PIDN))
#targets.new.diagnosis.gender <- inner_join(targets.new.diagnosis, temp)

#age.temp <- select(rna.ucsf.validage, PIDN, AgeAtDraw)
#temp <- filter(age.temp, duplicated(age.temp$PIDN))
#temp.dup <- filter(age.temp, is.element(age.temp$PIDN, temp$PIDN))

#temp <- filter(age.temp, is.element(age.temp$PIDN, targets.new.diagnosis.gender$PIDN))


###################################################################################################
##### microarray table unfiltered #################################################################
###################################################################################################

temp <- filter(targets.merged, is.na(targets.merged$Batch))
temp0 <- select(temp, Array, Code, Sample, Slide, Platform)

temp.other <- filter(targets.merged, !is.na(targets.merged$Batch))
temp.other2 <- select(temp.other, Array, Code, Sample, Slide, Platform, Batch)
temp1 <- filter(targets.v4a.batch, is.element(targets.v4a.batch$Slide, temp$Slide))
temp2 <- select(temp1, Slide, Batch)

#temp3 <- map(temp2, as.character, temp2$Slide)
#temp4 <- data.frame(temp3)
temp.merged <- inner_join(temp0, temp2)

targets.first <- rbind(temp.other2, temp.merged)

colnames(targets.first)[2] <- "PIDN"

####merge new targets and the original###########################
## only select Slide, PIDN, Batch, Platform
targets.a <- select(targets.first, Slide, PIDN, Batch, Platform)
targets.b <- select(targets.new, Slide, PIDN, Batch, Platform)
## merge
targets.complete.list <- rbind(targets.a, targets.b)


####targets with invalid age and valid age ######################
targets.validage <- filter(targets.complete.list, is.element(targets.complete.list$PIDN, rna.ucsf.validage$PIDN))    #1763
targets.invalidage <- filter(targets.complete.list, is.element(targets.complete.list$PIDN, rna.ucsf.noage$PIDN))     #29

temp <- rbind(targets.validage, targets.invalidage)

targets.no.rna <- filter(targets.complete.list, !is.element(targets.complete.list$PIDN, temp$PIDN))
targets.no.rna$AgeAtDraw <- NA                                                                                       #31
write.xlsx(targets.no.rna, "targets_no_rna.xlsx")


####rna filtered for targets ####################################
rna.validage.filtered <- filter(rna.ucsf.validage, is.element(rna.ucsf.validage$PIDN, targets.validage$PIDN))
rna.invalidage.filtered <- filter(rna.ucsf.noage, is.element(rna.ucsf.noage$PIDN, targets.invalidage$PIDN))

#############################################################################################
###### filter patients and rna ##############################################################
#############################################################################################
## use the rna table as basis
rna.reduced <- select(rna.use, RNA.ID, PIDN, Sample.Num, AgeAtDraw)
patients.reduced <- select(patients.ucsf, PIDN, DNA, RNA, Gender)
rna.patients <- inner_join(rna.reduced, patients.reduced)                #why did i do this?

diagnosis.reduced <- select(diagnosis.data, PIDN, Diagnosis)
#temp <-filter(rna.patients, is.element(rna.patients$PIDN, diagnosis.reduced$PIDN))
#temp <- left_join(rna.patients, diagnosis.reduced)

#temp <- filter(rna.patients, is.element(rna.patients$PIDN, targets.complete.list$PIDN))

#use the targets.first table because there is no duplicates in the other one
targets.duplicated <- filter(targets.first, duplicated(targets.first$PIDN))
targets.plus2 <- filter(targets.first, is.element(targets.first$PIDN, targets.duplicated$PIDN))
#rna.plus2 <- filter(rna.ucsf.validage, is.element(rna.ucsf.validage$PIDN, targets.duplicated$PIDN))
write.xlsx(targets.plus2, "targets_duplicated.xlsx")
write.xlsx(rna.plus2, "rna_duplicated.xlsx")

#targets.not.duplicated <- filter(targets.first, !is.element(targets.first$PIDN, targets.duplicated$PIDN))
#temp <- filter( rna.ucsf.validage, is.element(rna.ucsf.validage$PIDN, targets.not.duplicated$PIDN))

#Assign NA to microarray with invalid age ####################
targets.invalidage$AgeAtDraw <- NA



### make data frame with number of unique age for each PIDN ######
### using the rna table with valid age ###########################
age.list <- by(rna.validage.filtered, rna.validage.filtered$PIDN, select, AgeAtDraw) %>% map(compose(unique, unlist))
pidn <- names(age.list)
number <- map_int(age.list, length)
age.df <- data.frame(pidn, number)

oneage.df <- filter(list.age, number == 1)
oneage.list <- age.list[names(age.list) %in% oneage.df$pidn] 
oneage.vector <- reduce(oneage.list, c)

targets.validage.oneage <- filter(targets.validage, is.element(targets.validage$PIDN, names(oneage.list)))     #1650
targets.validage.plus2 <- filter(targets.validage, !is.element(targets.validage$PIDN, names(oneage.list)))     #113
write.xlsx(targets.validage.plus2, "targets_more_than_1_unique_age.xlsx")

rna.validage.plus2 <- filter(rna.validage.filtered, !is.element(rna.validage.filtered$PIDN, names(oneage.list)))
write.xlsx(rna.validage.plus2, "rna_more_than_1_unique_age.xlsx")


oneage.df <- data.frame(names(oneage.list), oneage.vector)
colnames(oneage.df)[1] <- "PIDN"
colnames(oneage.df)[2] <- "AgeAtDraw"

targets.validage.oneage.merged <- inner_join(targets.validage.oneage, oneage.df)                        #1650

######### targets with only one age + targets with invalid age + targets with no rna data
targets.1age <- rbind(targets.validage.oneage.merged, targets.invalidage, targets.no.rna)               #1650 + 29 + 31 = 1710
#########

#check the number of microarrays for each pidn and if it matches the number of ages separate them and match
plus2age.df <- filter(list.age, number > 1)
colnames(plus2age.df)[1] <- "PIDN"
colnames(plus2age.df)[2] <- "Ages"

mnumber.list <- by(targets.validage.plus2, targets.validage.plus2$PIDN, select, PIDN) %>% map(length) %>% unlist

mnumber.df <- data.frame(names(mnumber.list), mnumber.list)
colnames(mnumber.df)[1] <- "PIDN"
colnames(mnumber.df)[2] <- "Arrays"

arrays.ages.number <- inner_join(mnumber.df, plus2age.df)
#map(arrays.ages.number, arrays.ages.number$Arrays, arrays.ages.number$Ages)
rownames(arrays.ages.number) <- arrays.ages.number$PIDN
arrays.ages.lgl<- select(arrays.ages.number, -PIDN) %>% apply(1, reduce, identical)

#filter targets table based on the logical table
arrays.ages.identical <- filter(arrays.ages.number, arrays.ages.lgl)
PIDNs <- arrays.ages.identical$PIDN
targets.plus2.identical <- filter(targets.validage.plus2, is.element(targets.validage.plus2$PIDN, PIDNs))
#number of microarrays and ages not identical
arrays.ages.not.identical <- filter(arrays.ages.number, !arrays.ages.lgl)
PIDNs.n <- arrays.ages.not.identical$PIDN
targets.plus2.not.identical <- filter(targets.validage.plus2, is.element(targets.validage.plus2$PIDN, PIDNs.n))
age.not.identical <- age.list[PIDNs.n]
age.not.identical.rna <- filter(rna.validage.plus2, is.element(rna.validage.plus2$PIDN, PIDNs.n))

write.xlsx(targets.plus2.not.identical, "targets_fix.xlsx")
write.xlsx(age.not.identical.rna, "rna_fix.xlsx")


#check if its in order
temp <- arrange(targets.plus2.identical, PIDN, as.numeric(Batch))
temp$Platform %<>% str_replace("v", "") %>% as.numeric
by(temp, temp$PIDN, select, Platform) %>% map_lgl(reduce, `<=`)

#fixxxxxxx
age.identical <- age.list[PIDNs]
age.identical.vector <- map(age.identical, sort) %>% reduce(c)


targets.plus2.identical$Batch %<>% as.integer
targets.plus2.identical %<>% arrange(PIDN, Batch)
targets.plus2.identical$AgeAtDraw <- age.identical.vector

age.difference <- map(age.identical, reduce, `-`) %>% map_dbl(abs)
age.difference.df <- data.frame(names(age.difference), age.difference)
colnames(age.difference.df)[1] <- "PIDN"
colnames(age.difference.df)[2] <- "Age Difference"
write.xlsx(age.difference.df, "age_difference.xlsx")
#17570
#4577


########## targets with one age + targets with invalid age + targets with no rna data + targets with identical number of microarray&rna
targets.1age.plus.identical <- rbind(targets.1age, targets.plus2.identical)           #1710 + 76 = 1786
##########

targets.2age <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\targets_fix.xlsx")




######## complete list of the targets ######################################
targets.data <- rbind(targets.1age.plus.identical, targets.2age)
############################################################################




##################################################################################
############ fixing diagnosis table ##############################################
##################################################################################
#2447
diagnosis.n <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\new_diagnosis.xlsx")

#drop molecular
#2421
diagnosis.new <- filter(diagnosis.n, diagnosis.n$DX.Type != "MOLECULAR")

#filtered for PIDN
#2421
diagnosis.new.filtered <- filter(diagnosis.new, is.element(diagnosis.new$PIDN, targets.data$PIDN))
diagnosis.new.nfiltered <- filter(diagnosis.new, !is.element(diagnosis.new$PIDN, targets.data$PIDN))

diagnosis.sep <- by(diagnosis.new.filtered, diagnosis.new.filtered$PIDN, select, Diagnosis)
diagnosis.count <- map_int(diagnosis.sep, length)
diagnosis.count.df <- data.frame(names(diagnosis.count), diagnosis.count)
colnames(diagnosis.count.df)[1] <- "PIDN"
colnames(diagnosis.count.df)[2] <- "Count"

diagnosis.1 <- filter(diagnosis.count.df, diagnosis.count.df$Count == 1)
diagnosis.1.pidn <- diagnosis.1$PIDN

#contains 1 diagnosis per PIDN
#1663
diagnosis.one <- filter(diagnosis.new.filtered, is.element(diagnosis.new.filtered$PIDN, diagnosis.1.pidn))
diagnosis.one.filtered <- select(diagnosis.one, PIDN, Diagnosis)
#join with targets with one diagnosis
#1717
targets.one.diagnosis <- filter(targets.data, is.element(targets.data$PIDN, diagnosis.one$PIDN))
#1717
targets.one.diagnosis.join <- inner_join(targets.one.diagnosis, diagnosis.one.filtered)

#contains more than 1 diagnosis per PIDN
#25
diagnosis.more <- filter(diagnosis.new.filtered, !is.element(diagnosis.new.filtered$PIDN, diagnosis.1.pidn))
#targets with more than one diagnosis
#8
targets.more.diagnosis <- filter(targets.data, is.element(targets.data$PIDN, diagnosis.more$PIDN))
diagnosis.list.unique <- by(diagnosis.more, diagnosis.more$PIDN, select, Diagnosis) %>% map(unique)
diagnosis.names <- c("Control", "Other", "Withdrawn", "Withdrawn", "Withdrawn", "Withdrawn", "Withdrawn", "Other")
diagnosis.list.unique.df <- data.frame(names(diagnosis.list.unique), diagnosis.names )
colnames(diagnosis.list.unique.df)[1] <- "PIDN"
colnames(diagnosis.list.unique.df)[2] <- "Diagnosis"
targets.more.diagnosis.join <- inner_join(targets.more.diagnosis, diagnosis.list.unique.df)

#Microarray has no diagnosis
#98
#1717 + 8 + 98 = 1823 = targets.data
targets.no.diagnosis <- filter(targets.data, !is.element(targets.data$PIDN, diagnosis.one$PIDN) & !is.element(targets.data$PIDN, diagnosis.more$PIDN))
targets.no.diagnosis$Diagnosis <- "Unknown"


######## complete list of all diagnosis with per microarray data ###########
diagnosis.data <- rbind(targets.one.diagnosis.join, targets.more.diagnosis.join, targets.no.diagnosis)
###########################################################################

## get sex from patients table join by PIDN
patients.ucsf <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\parkyj_new_patients_20160414T182054.xlsx")
patients.filtered <- filter(patients.ucsf, is.element(patients.ucsf$PIDN, targets.data$PIDN))
patients.sex <- select(patients.filtered, PIDN, Gender)

#1798
targets.sex <- filter(diagnosis.data, is.element(diagnosis.data$PIDN, patients.sex$PIDN))
#25
targets.no.sex <- filter(diagnosis.data, !is.element(diagnosis.data$PIDN, patients.sex$PIDN))
targets.no.sex$Gender <- "UNKNOWN"

patients.merged <- inner_join(targets.sex, patients.sex)
##### final microarray data with gender, diagnosis, age ############
patients.data <- rbind(patients.merged, targets.no.sex)
####################################################################
write.xlsx(patients.data, "patients.data.xlsx")



#male, female, unknown, average age per diagnosis#####################################
#for the purpose of this analysis use the max age on arrays with more than one age
patients.data.max.age <- by(patients.data, patients.data$PIDN, select, AgeAtDraw) %>% map(max)

patients.dd <- reduce(patients.data.max.age, rbind)
patients.dd.df <- data.frame(names(patients.data.max.age), patients.dd)
colnames(patients.dd.df)[1] <- "PIDN"
colnames(patients.dd.df)[2] <- "AgeAtDraw"
patients.s <- filter(patients.data, is.element(patients.data$PIDN, patients.dd.df$PIDN))
patients.dd.diagnosis <- select(patients.s, PIDN, Diagnosis, Gender)
patients.dd.join <- inner_join(patients.dd.df, patients.dd.diagnosis)
patients.dd.filtered <- filter(patients.dd.join, !duplicated(patients.dd.join$PIDN))



#sex data
gender <- by(patients.dd.filtered, patients.dd.filtered$Diagnosis, select, Gender) %>% map(unlist)
unique.values <- reduce(gender, c) %>% unique
gender.summary <- map(gender, Compose(factor, summary))

add.missing <- function(missing.array, all.values)
{
  if (length(missing.array) < length(all.values))
  {
    paste.exact <- paste %<<<% "^" %<<% c("$", sep = "")   
    found.key <- map_chr(names(missing.array), paste.exact) %>% paste(collapse = "|")
    missing.key <- !grepl(found.key, all.values)
    add.values <- all.values[missing.key]
    add.vector <- rep(0, length(add.values))
    names(add.vector) <- add.values
    new.vector <- c(missing.array, add.vector) %>% sort
    return(new.vector)
  }
  else
  {
    return(missing.array)
  }
}

gender.summary.fixed <- map(gender.summary, add.missing, unique.values)
diagnosis.names <- names(gender.summary.fixed)
gender.summary.fixed.reduce <- reduce(gender.summary.fixed, rbind)
Total <- rowSums(gender.summary.fixed.reduce)
patients.data.sex <- data.frame(diagnosis.names, gender.summary.fixed.reduce, Total)
colnames(patients.data.sex)[1] <- "Diagnosis"

#average age
patients.not.na <- filter(patients.dd.filtered, !is.na(patients.dd.filtered$AgeAtDraw))

patients.mean.age <- by(patients.not.na, patients.not.na$Diagnosis, select, AgeAtDraw) %>% map(mean)
patients.mean.diagnosis <- names(patients.mean.age)
patients.data.mean.age <- data_frame(patients.mean.diagnosis, patients.mean.age)
colnames(patients.data.mean.age)[1] <- "Diagnosis"
colnames(patients.data.mean.age)[2] <- "Mean Age"


diagnosis.sex.age <- inner_join(patients.data.sex, patients.data.mean.age)

#### get number of ages that are missing per diagnosis
patients.noage.length <- by(patients.dd.filtered, patients.dd.filtered$Diagnosis, select, AgeAtDraw) %>% map(map_lgl, is.na) %>% map_int(Compose(which, length))
patients.noage.df <- data.frame(names(patients.noage.length), patients.noage.length)
rownames(patients.noage.df) <- 1:28
colnames(patients.noage.df)[1] <- "Diagnosis"
colnames(patients.noage.df)[2] <- "Unknown Age"
diagnosis.summary <- inner_join(diagnosis.sex.age, patients.noage.df)
diagnosis.summary$`Mean Age` %<>% unlist

write.xlsx(diagnosis.summary, "diagnosis_summary.xlsx")

############################################################

sequencing.summary <- read.xlsx("C:\\Users\\Youngjun Park\\Desktop\\leb\\dementia\\parkyj_summary.ucsf_20160415T185900.xlsx")
s.filtered <- filter(sequencing.summary, is.element(sequencing.summary$PIDN, patients.data$PIDN))
s.positive <- filter(s.filtered, s.filtered$Mutation == "Positive")
write.xlsx(s.positive, "mutations_summary.xlsx")


##### make summary table with the PIDNs with the "Positive" mutations ##########
################################################################################

patients.positive <- filter(patients.dd.filtered, is.element(patients.dd.filtered$PIDN, s.positive$PIDN))


#sex data
gender <- by(patients.positive, patients.positive$Diagnosis, select, Gender) %>% map(unlist)
unique.values <- reduce(gender, c) %>% unique
gender.summary <- map(gender, Compose(factor, summary))


gender.summary.fixed <- map(gender.summary, add.missing, unique.values)
diagnosis.names <- names(gender.summary.fixed)
gender.summary.fixed.reduce <- reduce(gender.summary.fixed, rbind)
Total <- rowSums(gender.summary.fixed.reduce)
patients.data.sex <- data.frame(diagnosis.names, gender.summary.fixed.reduce, Total)
colnames(patients.data.sex)[1] <- "Diagnosis"

#average age
patients.not.na <- filter(patients.positive, !is.na(patients.positive$AgeAtDraw))

patients.mean.age <- by(patients.not.na, patients.not.na$Diagnosis, select, AgeAtDraw) %>% map(mean)
patients.mean.diagnosis <- names(patients.mean.age)
patients.data.mean.age <- data_frame(patients.mean.diagnosis, patients.mean.age)
colnames(patients.data.mean.age)[1] <- "Diagnosis"
colnames(patients.data.mean.age)[2] <- "Mean Age"


diagnosis.sex.age <- inner_join(patients.data.sex, patients.data.mean.age)

#### get number of ages that are missing per diagnosis
patients.noage.length <- by(patients.positive, patients.positive$Diagnosis, select, AgeAtDraw) %>% map(map_lgl, is.na) %>% map_int(Compose(which, length))
patients.noage.df <- data.frame(names(patients.noage.length), patients.noage.length)
rownames(patients.noage.df) <- 1:13
colnames(patients.noage.df)[1] <- "Diagnosis"
colnames(patients.noage.df)[2] <- "Unknown Age"
diagnosis.positive.summary <- inner_join(diagnosis.sex.age, patients.noage.df)
diagnosis.positive.summary$`Mean Age` %<>% unlist
diagnosis.positive.summary %<>% select(-UNKNOWN)

write.xlsx(diagnosis.positive.summary, "diagnosis_POSITIVE.xlsx")

#####################################################################################
#drop diagnosis with small samples and prepare for linear modeling###################
#####################################################################################

diagnosis.include <- c("AD", "ALS", "bvFTD", "CBS", "clinically normal", "Control", "CONTROL", "FTD/ALS", "MCI", "nfPPA", "PSP", "svPPA") %>% paste(collapse = "|")
patients.data.filtered <- filter(patients.dd.filtered, grepl(diagnosis.include, patients.dd.filtered$Diagnosis))
#patients.data.filtered$Diagnosis %>% unique
#table(patients.data.filtered$Diagnosis)
patients.data.filtered %<>% filter(!is.na(AgeAtDraw))
patients.data.filtered$Diagnosis %<>% str_replace_all("\\/ ", "") %>% toupper

linear.model <- lm(AgeAtDraw ~ Diagnosis, patients.data.filtered) %>% anova %>% tidy
linear.model.aov <- aov(AgeAtDraw ~ Diagnosis, patients.data.filtered) %>% TukeyHSD %>% tidy

linear.model.sex <- lm(as.integer(factor(Gender)) ~ Diagnosis, patients.data.filtered) %>% anova %>% tidy
linear.model.sex.aov <- aov(as.integer(factor(Gender)) ~ Diagnosis, patients.data.filtered) %>% TukeyHSD %>% tidy

#create box plot for AgeAtDraw vs Diagnosis
p <- ggplot(patients.data.filtered, aes(x = Diagnosis, y = AgeAtDraw)) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Age at Draw")
CairoPDF("age_boxplot.pdf", height = 6, width = 11)
plot(p)
dev.off()

p <- ggplot(patients.data.filtered, aes(x = Diagnosis, fill = factor(Gender))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title.x = element_blank()) + ylab("Patients") + scale_fill_discrete(name = "Sex")
CairoPDF("sex_barplot.pdf", height = 6, width = 15)
plot(p)
dev.off()


####################################################################################
######## bar plot for the gene mutations ###########################################
####################################################################################
Genes <- select(s.positive, PIDN, Gene)
patients.mutations <- inner_join(Genes, patients.dd.filtered)
patients.mutations.filtered <- filter(patients.mutations, !duplicated(patients.mutations$PIDN))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#sample goes through and picks random subset, by default sample w/o replacement
colors <- sample(col_vector, 13)


p <- ggplot(patients.mutations.filtered, aes(x = Gene, fill = Diagnosis)) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + scale_fill_manual(values = colors)
CairoPDF("mutations_barplot.pdf", height = 15, width = 10)
plot(p)
dev.off()

