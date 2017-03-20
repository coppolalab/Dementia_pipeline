#For DE
library(lumi)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(limma)
library(annotate)
library(biomaRt)
library(siggenes)
library(sva)
library(broom)
library(WGCNA)
library(tools)

#String operations
library(stringr)

#Reading and writing tables
library(readr)
library(openxlsx)

#Plotting
library(ggplot2)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)

#Data arrangement
library(dplyr)
library(tidyr)

#Functional programming
library(magrittr)
library(purrr)
library(rlist)

table.2007.50 <- read_csv("../raw_data/platform_v2/2007-050 sample probe profile.csv")
table.2008.72 <- read_csv("../raw_data/platform_v2/2008-072 sample probe profile.csv")
write.table(table.2007.50, "../raw_data/platform_v2/2007-050.tsv", sep = "\t")
write.table(table.2008.72, "../raw_data/platform_v2/2008-072.tsv", sep = "\t")

platform.1 <- list.files("../raw_data/platform_v1", full.names = TRUE) 
platform.2 <- list.files("../raw_data/platform_v2", full.names = TRUE) 
platform.3 <- list.files("../raw_data/platform_v3", full.names = TRUE)
platform.4.tsv <- list.files("../raw_data/platform_v4", full.names = TRUE) %>% str_subset("tsv|txt")
platform.4.csv <- list.files("../raw_data/platform_v4", full.names = TRUE) %>% str_subset("csv")

lumi.v1 <- lumiR.batch(platform.1, QC = FALSE, transform = 'vst', convertNuID = TRUE, checkDupId = TRUE, lib.mapping = "lumiHumanIDMapping")
#lumi.v1 <- lumiR.batch(platform.1, QC = FALSE, transform = 'vst')

lumi.v2.new <- lumiR.batch(platform.2[c(1,8)], QC = FALSE, transform = 'vst', convertNuID = TRUE, checkDupId = TRUE, lib.mapping = "lumiHumanIDMapping")
lumi.v2.old <- lumiR.batch(platform.2[c(2:7,9)], QC = FALSE, transform = 'vst', convertNuID = TRUE, checkDupId = TRUE, lib.mapping = "lumiHumanIDMapping")
lumi.v2.new <- lumi.v2.new[featureNames(lumi.v2.old),]
lumi.v2.annot <- nuID2IlluminaID(featureNames(lumi.v2.new), "lumiHumanIDMapping", idType = 'Gene')
lumi.v2 <- BiocGenerics::combine(lumi.v2.old, lumi.v2.new)
#featureNames(lumi.v2) %<>% nuID2probeID
#lumi.v2.annot.df <- data.frame(ProbeID = featureNames(lumi.v2), Symbol = lumi.v2.annot)
write.csv(lumi.v2.annot.df, "v2.annotation.csv", row.names = FALSE)

lumi.v3 <- lumiR.batch(platform.3, QC = FALSE, transform = 'vst', convertNuID = TRUE, checkDupId = TRUE, lib.mapping = "lumiHumanIDMapping")
#lumi.v3 <- lumiR.batch(platform.3, QC = FALSE, transform = 'vst')

lumi.4.tsv.list <- map(platform.4.tsv, read_tsv) 
tsv.4.colnames <- map(lumi.4.tsv.list, colnames) %>% map(toupper) %>% map(str_replace, "PROBEID", "PROBE_ID")
lumi.4.tsv.list <- map2(lumi.4.tsv.list, tsv.4.colnames, set_colnames)
lumi.4.csv.list <- map(platform.4.csv, read_csv)
lumi.4.table <- reduce(lumi.4.tsv.list, left_join) %>% left_join(reduce(lumi.4.csv.list, left_join))
write.table(lumi.4.table, "../raw_data/platform_v4/platform_v4.tsv", row.names = FALSE, sep = "\t")
lumi.v4 <- lumiR("../raw_data/platform_v4/platform_v4.tsv", QC = FALSE, convertNuID = TRUE, checkDupId = TRUE, lib.mapping = "lumiHumanIDMapping") %>% lumiT("vst")
#lumi.v4 <- lumiR("../raw_data/platform_v4/platform_v4.tsv", QC = FALSE) %>% lumiT("vst")

pdata <- read.xlsx("../phenotypedata/final_microarray_data.xlsx") 
pdata$Batch %<>% as.integer
pdata %<>% arrange(Batch, Slide)
pdata$Diagnosis %<>% toupper %>% str_replace("SVPPA", "SD") %>% str_replace("NFPPA\\/ UNSPEC\\. PPA", "PNFA")
diagnosis.vector <- str_c("AD", "BVFTD", "CBS", "CONTROL", "MCI", "PNFA", "PSP", "SD", sep = "|")
#pdata.filter <- filter(pdata, grepl(diagnosis.vector, Diagnosis))

pdata.newer <- read.xlsx("../2016-9088_Sample Key.xlsx") 
pdata.newer$Slide <- str_c(pdata.newer$general.array, pdata.newer$genexstripe.controling.stripe, sep = "_")
pdata.newer.reduce <- select(pdata.newer, Slide, External.ID)
colnames(pdata.newer.reduce)[2] <- "PIDN"
pdata.newer.reduce$Batch <- 38
pdata.newer.reduce$Platform <- "v4"
pdata.newer.reduce$AgeAtDraw <- NA
pdata.newer.reduce$Diagnosis <- NA
pdata.newer.reduce$Gender <- NA

pdata.newest <- read.xlsx("../2016-9110A_Sample Key.xlsx") 
pdata.newest$Slide <- str_c(pdata.newest$general.array, pdata.newest$genexstripe.controling.stripe, sep = "_")
pdata.newest.reduce <- select(pdata.newest, Slide, External.ID)
colnames(pdata.newest.reduce)[2] <- "PIDN"
pdata.newest.reduce$Batch <- 39
pdata.newest.reduce$Platform <- "v4"
pdata.newest.reduce$AgeAtDraw <- NA
pdata.newest.reduce$Diagnosis <- NA
pdata.newest.reduce$Gender <- NA

pdata.all <- rbind(pdata, pdata.newer.reduce, pdata.newest.reduce)
pdata.all$PIDN %<>% str_replace(" ", "")
pdata.filter <- filter(pdata.all, !str_detect(PIDN, "[a-z]|[A-Z]"))
pdata.filter$Sample.Num <- str_split_fixed(pdata.filter$PIDN, "_", 2)[,2]
pdata.filter$PIDN %<>% str_replace("_.*$", "")

hillblom.request <- read.xlsx("../Hillblom RNA_microarray data request reformat.xlsx")
hillblom.shared <- intersect(pdata.filter$PIDN, hillblom.request$PIDN)
hillblom.found <- filter(hillblom.request, PIDN %in% hillblom.shared)
pdata.hillblom <- filter(pdata.filter, PIDN %in% hillblom.shared)

pdata.v1.hillblom <- filter(pdata.hillblom, Platform == "v1") %>% filter(Slide %in% sampleNames(lumi.v1))
pdata.v1.hillblom.missing <- filter(pdata.hillblom, Platform == "v1") %>% filter(!(Slide %in% sampleNames(lumi.v1)) )
lumi.v1.reduce <- lumi.v1[,pdata.v1.hillblom$Slide]
#lumi.v1.symbol <- getSYMBOL(featureNames(lumi.v1.reduce), "lumiHumanAll.db")
lumi.v1.norm <- lumiN(lumi.v1.reduce, method = "rsn")
write.csv(exprs(lumi.v1.norm), "./v1_expression.csv")
write.csv(detection(lumi.v1.norm), "./v1_detection.csv")
write.xlsx(pdata.v1.hillblom, "v1_phenodata.xlsx")

pdata.v2.hillblom <- filter(pdata.hillblom, Platform == "v2") %>% filter(Slide %in% sampleNames(lumi.v2))
pdata.v2.hillblom.missing <- filter(pdata.hillblom, Platform == "v2") %>% filter(!(Slide %in% sampleNames(lumi.v2)) )
lumi.v2.reduce <- lumi.v2[,pdata.v2.hillblom$Slide]
#lumi.v2.symbol <- getSYMBOL(featureNames(lumi.v2.reduce), "lumiHumanAll.db")
lumi.v2.norm <- lumiN(lumi.v2.reduce, method = "rsn")
write.csv(exprs(lumi.v2.norm), "./v2_expression.csv")
write.csv(detection(lumi.v2.norm), "./v2_detection.csv")
write.xlsx(pdata.v2.hillblom, "v2_phenodata.xlsx")

pdata.v3.hillblom <- filter(pdata.hillblom, Platform == "v3") %>% filter(Slide %in% sampleNames(lumi.v3))
pdata.v3.hillblom.missing <- filter(pdata.hillblom, Platform == "v3") %>% filter(!(Slide %in% sampleNames(lumi.v3)) )
lumi.v3.reduce <- lumi.v3[,pdata.v3.hillblom$Slide]
#lumi.v3.symbol <- getSYMBOL(featureNames(lumi.v3.reduce), "lumiHumanAll.db")
lumi.v3.norm <- lumiN(lumi.v3.reduce, method = "rsn")
write.csv(exprs(lumi.v3.norm), "./v3_expression.csv")
write.csv(detection(lumi.v3.norm), "./v3_detection.csv")
write.xlsx(pdata.v3.hillblom, "v3_phenodata.xlsx")

pdata.v4.hillblom <- filter(pdata.hillblom, Platform == "v4") %>% filter(Slide %in% sampleNames(lumi.v4)) 
pdata.v4.hillblom.missing <- filter(pdata.hillblom, Platform == "v4") %>% filter(!(Slide %in% sampleNames(lumi.v4)) )
lumi.v4.reduce <- lumi.v4[,pdata.v4.hillblom$Slide]
#lumi.v4.symbol <- getSYMBOL(featureNames(lumi.v4.reduce), "lumiHumanAll.db")
lumi.v4.norm <- lumiN(lumi.v4.reduce, method = "rsn")
write.csv(exprs(lumi.v4.norm), "./v4_expression.csv")
write.csv(detection(lumi.v4.norm), "./v4_detection.csv")
write.xlsx(pdata.v4.hillblom, "v4_phenodata.xlsx")

pdata.v1.missing <- filter(pdata.filter, Platform == "v1") %>% filter(!(Slide %in% sampleNames(lumi.v1)))
write.xlsx(pdata.v1.missing, "pdata.v1.missing.xlsx")

pdata.v1 <- filter(pdata.filter, Platform == "v1")# %>% filter(Slide %in% sampleNames(lumi.v1))
rownames(pdata.v1) <- pdata.v1$Slide
lumi.v1.all <- lumi.v1[,pdata.v1$Slide]
pData(lumi.v1.all) <- pdata.v1 
lumi.v1.complete <- lumi.v1.all[,!is.na(lumi.v1.all$Gender) & !is.na(lumi.v1.all$AgeAtDraw)]
lumi.v1.norm.all <- lumiN(lumi.v1.complete, method = "rsn")
lumi.v1.symbol <- getSYMBOL(featureNames(lumi.v1.norm.all), "lumiHumanAll.db")

v1.matrix <- model.matrix( ~ AgeAtDraw + Gender, pData(lumi.v1.complete))
v1.combat <- ComBat(dat = exprs(lumi.v1.norm.all), batch = factor(lumi.v1.norm.all$Batch), mod = v1.matrix)
v1.rmcov <- empiricalBayesLM(t(v1.combat), removedCovariates = select(pData(lumi.v1.complete), AgeAtDraw, Gender), robustPriors = FALSE) 
v1.rmcov.expr <- t(v1.rmcov$adjustedData)
lumi.v1.grn <- v1.rmcov.expr[which(str_detect("GRN", lumi.v1.symbol)),]
v1.grn.sd <- (lumi.v1.grn - mean(lumi.v1.grn)) / sd(lumi.v1.grn) 
v1.grn.df <- data.frame(Slide = sampleNames(lumi.v1.complete), GRN = v1.grn.sd, GRN.expr = lumi.v1.grn) %>% left_join(pdata.v1) %>% arrange(GRN)

v1.grn.plot <- arrange(v1.grn.df, PIDN)
v1.grn.plot$Sample <- 1:nrow(v1.grn.plot)
v1.grn.plot$Outlier <- factor(abs(v1.grn.plot$GRN) > 2, levels = c("TRUE", "FALSE"))
p <- ggplot(v1.grn.plot, aes(x = Sample, y = GRN.expr, col = Outlier, label = PIDN)) + geom_point() + theme_bw()
p <- p + geom_text(data = filter(v1.grn.plot, PIDN == 3812), nudge_y = -0.15, col = "black") 
p <- p + geom_point(data = filter(v1.grn.plot, PIDN == 3812), col = "black")
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), legend.position = "none")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept = mean(v1.grn.plot$GRN.expr), col = "black", size = 1)
p <- p + ylab("VST GRN Expression")
CairoPDF("v1.grn.pdf", width = 6, height = 4, bg = "transparent")
print(p)
dev.off()

pdata.v2.missing <- filter(pdata.filter, Platform == "v2") %>% filter(!(Slide %in% sampleNames(lumi.v2)))
write.xlsx(pdata.v2.missing, "pdata.v2.missing.xlsx")

pdata.v2 <- filter(pdata.filter, Platform == "v2") %>% filter(Slide %in% sampleNames(lumi.v2))
rownames(pdata.v2) <- pdata.v2$Slide
lumi.v2.all <- lumi.v2[,pdata.v2$Slide]
pData(lumi.v2.all) <- pdata.v2 
lumi.v2.complete <- lumi.v2.all[,!is.na(lumi.v2.all$Gender) & !is.na(lumi.v2.all$AgeAtDraw)]
lumi.v2.norm.all <- lumiN(lumi.v2.complete, method = "rsn")

v2.matrix <- model.matrix( ~ AgeAtDraw + Gender, pData(lumi.v2.complete))
v2.combat <- ComBat(dat = exprs(lumi.v2.norm.all), batch = factor(lumi.v2.norm.all$Batch), mod = v2.matrix)
v2.rmcov <- empiricalBayesLM(t(v2.combat), removedCovariates = select(pData(lumi.v2.complete), AgeAtDraw, Gender), robustPriors = FALSE) 
v2.rmcov.expr <- t(v2.rmcov$adjustedData)
lumi.v2.grn <- v2.rmcov.expr[which(str_detect("GRN", lumi.v2.symbol)),]
v2.grn.sd <- (lumi.v2.grn - mean(lumi.v2.grn)) / sd(lumi.v2.grn) 
v2.grn.df <- data.frame(Slide = sampleNames(lumi.v2.complete), GRN1 = v2.grn.sd[1,], GRN2 = v2.grn.sd[2,], GRN1.expr = lumi.v2.grn[1,], GRN2.expr = lumi.v2.grn[2,]) %>% left_join(pdata.v2) %>% arrange(GRN1)

v2.grn.plot <- arrange(v2.grn.df, PIDN)
v2.grn.plot$Sample <- 1:nrow(v2.grn.plot)
v2.grn.plot$Outlier <- factor(abs(v2.grn.plot$GRN1) > 2, levels = c("TRUE", "FALSE")) 
p <- ggplot(v2.grn.plot, aes(x = Sample, y = GRN1.expr, col = Outlier, label = PIDN)) + geom_point() + theme_bw()
#p <- p + geom_text(data = filter(v2.grn.plot, Outlier == "TRUE"), nudge_y = 0.05, col = "black") + theme_bw()
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), legend.position = "none")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept = mean(v2.grn.plot$GRN1.expr), col = "black", size = 1)
p <- p + ylab("VST GRN Expression")
CairoPDF("v2.grn1.pdf", width = 10, height = 8, bg = "transparent")
print(p)
dev.off()

v2.grn.plot <- arrange(v2.grn.df, PIDN)
v2.grn.plot$Sample <- 1:nrow(v2.grn.plot)
v2.grn.plot$Outlier <- factor(abs(v2.grn.plot$GRN2) > 2, levels = c("TRUE", "FALSE"))
p <- ggplot(v2.grn.plot, aes(x = Sample, y = GRN2.expr, col = Outlier, label = PIDN)) + geom_point() + theme_bw() 
#p <- p + geom_text(data = filter(v2.grn.plot, Outlier == "TRUE"), nudge_y = 0.05, col = "black") + theme_bw()
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), legend.position = "none")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept = mean(v2.grn.plot$GRN2.expr), col = "black", size = 1)
p <- p + ylab("VST GRN Expression")
CairoPDF("v2.grn2.pdf", width = 10, height = 8, bg = "transparent")
print(p)
dev.off()

pdata.v3 <- filter(pdata.filter, Platform == "v3") %>% filter(Slide %in% sampleNames(lumi.v3))
rownames(pdata.v3) <- pdata.v3$Slide
lumi.v3.all <- lumi.v3[,pdata.v3$Slide]
pData(lumi.v3.all) <- pdata.v3 
lumi.v3.complete <- lumi.v3.all[,!is.na(lumi.v3.all$Gender) & !is.na(lumi.v3.all$AgeAtDraw)]
lumi.v3.norm.all <- lumiN(lumi.v3.complete, method = "rsn")

v3.matrix <- model.matrix( ~ AgeAtDraw + Gender, pData(lumi.v3.complete))
v3.combat <- ComBat(dat = exprs(lumi.v3.norm.all), batch = factor(lumi.v3.norm.all$Batch), mod = v3.matrix)
v3.rmcov <- empiricalBayesLM(t(v3.combat), removedCovariates = select(pData(lumi.v3.complete), AgeAtDraw, Gender), robustPriors = FALSE) 
v3.rmcov.expr <- t(v3.rmcov$adjustedData)
lumi.v3.grn <- v3.rmcov.expr[which(str_detect("GRN", lumi.v3.symbol)),]
v3.grn.sd <- (lumi.v3.grn - mean(lumi.v3.grn)) / sd(lumi.v3.grn) 
v3.grn.df <- data.frame(Slide = sampleNames(lumi.v3.complete), GRN1 = v3.grn.sd[1,], GRN2 = v3.grn.sd[2,], GRN1.expr = lumi.v3.grn[1,], GRN2.expr = lumi.v3.grn[2,]) %>% left_join(pdata.v3) %>% arrange(GRN1)

v3.grn.plot <- arrange(v3.grn.df, PIDN)
v3.grn.plot$Sample <- 1:nrow(v3.grn.plot)
v3.grn.plot$Outlier <- factor(abs(v3.grn.plot$GRN1) > 2, levels = c("TRUE", "FALSE"))
p <- ggplot(v3.grn.plot, aes(x = Sample, y = GRN1.expr, col = Outlier, label = PIDN)) + geom_point() + theme_bw()
#p <- p + geom_text(data = filter(v3.grn.plot, Outlier == "TRUE"), nudge_y = 0.05, col = "black") + theme_bw()
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), legend.position = "none")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept = mean(v3.grn.plot$GRN1.expr), col = "black", size = 1)
p <- p + ylab("VST GRN Expression")
CairoPDF("v3.grn1.pdf", width = 10, height = 8, bg = "transparent")
print(p)
dev.off()

v3.grn.plot <- arrange(v3.grn.df, PIDN)
v3.grn.plot$Sample <- 1:nrow(v3.grn.plot)
v3.grn.plot$Outlier <- factor(abs(v3.grn.plot$GRN2) > 2, levels = c("TRUE", "FALSE"))
p <- ggplot(v3.grn.plot, aes(x = Sample, y = GRN2.expr, col = Outlier, label = PIDN)) + geom_point() + theme_bw()
#p <- p + geom_text(data = filter(v3.grn.plot, Outlier == "TRUE"), nudge_y = 0.05, col = "black") + theme_bw()
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), legend.position = "none")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept = mean(v3.grn.plot$GRN2.expr), col = "black", size = 1)
p <- p + ylab("VST GRN Expression")
CairoPDF("v3.grn2.pdf", width = 10, height = 8, bg = "transparent")
print(p)
dev.off()

pdata.v4.missing <- filter(pdata.filter, Platform == "v4") %>% filter(!(Slide %in% sampleNames(lumi.v4)))
write.xlsx(pdata.v4.missing, "pdata.v4.missing.xlsx")

pdata.v4 <- filter(pdata.filter, Platform == "v4") %>% filter(Slide %in% sampleNames(lumi.v4))
rownames(pdata.v4) <- pdata.v4$Slide
lumi.v4.all <- lumi.v4[,pdata.v4$Slide]
pData(lumi.v4.all) <- pdata.v4 
lumi.v4.norm.all <- lumiN(lumi.v4.all, method = "rsn")
lumi.v4.symbol <- getSYMBOL(featureNames(lumi.v4.norm.all), "lumiHumanAll.db")

v4.combat <- ComBat(dat = exprs(lumi.v4.norm.all), batch = factor(lumi.v4.norm.all$Batch))
lumi.v4.grn <- v4.combat[which(str_detect("GRN", lumi.v4.symbol)),] #%>% exprs
v4.grn.sd1 <- (lumi.v4.grn[1,] - mean(lumi.v4.grn[1,])) / sd(lumi.v4.grn[1,]) 
v4.grn.sd2 <- (lumi.v4.grn[2,] - mean(lumi.v4.grn[2,])) / sd(lumi.v4.grn[2,]) 
v4.grn.df <- data.frame(Slide = sampleNames(lumi.v4.norm.all), GRN1 = v4.grn.sd1, GRN2 = v4.grn.sd2, GRN1.expr = lumi.v4.grn[1,], GRN2.expr = lumi.v4.grn[2,]) %>% left_join(pdata.v4) %>% arrange(GRN1)

v4.grn.plot <- arrange(v4.grn.df, PIDN)
v4.grn.plot$Sample <- 1:nrow(v4.grn.plot)
v4.grn.plot$Outlier <- factor(abs(v4.grn.plot$GRN1) > 2, levels = c("TRUE", "FALSE"))
p <- ggplot(v4.grn.plot, aes(x = Sample, y = GRN1.expr, col = Outlier, label = PIDN)) + geom_point() 
p <- p + geom_text(data = filter(v4.grn.plot, PIDN == 13552 | PIDN == 13655), nudge_y = 0.05, col = "black") + theme_bw()
p <- p + geom_point(data = filter(v4.grn.plot, PIDN == 13552 | PIDN == 13655), col = "black") 
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), legend.position = "none")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept = mean(v4.grn.plot$GRN1.expr), col = "black", size = 1)
p <- p + ylab("VST GRN Expression")
CairoPDF("v4.grn1.pdf", width = 12, height = 10, bg = "transparent")
print(p)
dev.off()

v4.grn.plot <- arrange(v4.grn.df, PIDN)
v4.grn.plot$Sample <- 1:nrow(v4.grn.plot)
v4.grn.plot$Outlier <- factor(abs(v4.grn.plot$GRN2) > 2, levels = c("TRUE", "FALSE"))
p <- ggplot(v4.grn.plot, aes(x = Sample, y = GRN2.expr, col = Outlier, label = PIDN)) + geom_point() 
p <- p + geom_text(data = filter(v4.grn.plot, PIDN == 13552 | PIDN == 13655), nudge_y = 0.05, col = "black") + theme_bw()
p <- p + geom_point(data = filter(v4.grn.plot, PIDN == 13552 | PIDN == 13655), col = "black") 
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), legend.position = "none")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept = mean(v4.grn.plot$GRN2.expr), col = "black", size = 1)
p <- p + ylab("VST GRN Expression")
CairoPDF("v4.grn2.pdf", width = 12, height = 10, bg = "transparent")
print(p)
dev.off()

#Read in sequencing data
all.sequencing <- read.xlsx("../dan_new_summary_20170315T155130.xlsx") %>% filter(Center == 'UCSF' & Gene == 'GRN')
all.patient <- read.xlsx("../dan_ucsf_20170315T165103.xlsx") %>% select(PIDN, Categorized.Diagnosis)
#all.grn.df <- rbind(v1.grn.df, v2.grn.df, v3.grn.df, v4.grn.df)
grn.sequencing.v1 <- left_join(v1.grn.df, all.sequencing) %>% left_join(all.patient) %>% select(PIDN, GRN, Categorized.Diagnosis, Diagnosis, AgeAtDraw, Gender, Summary:Uncategorized, Sequencing.Site, Center, Gene:Sequencing.ID, Aliquot.ID, Sample.Num, Batch, Platform)
grn.sequencing.v2 <- left_join(v2.grn.df, all.sequencing) %>% left_join(all.patient) %>% select(PIDN, GRN1, GRN2, Categorized.Diagnosis, Diagnosis, AgeAtDraw, Gender, Summary:Uncategorized, Sequencing.Site, Center, Gene:Sequencing.ID, Aliquot.ID, Sample.Num, Batch, Platform)
grn.sequencing.v3 <- left_join(v3.grn.df, all.sequencing) %>% left_join(all.patient) %>% select(PIDN, GRN1, GRN2, Categorized.Diagnosis, Diagnosis, AgeAtDraw, Gender, Summary:Uncategorized, Sequencing.Site, Center, Gene:Sequencing.ID, Aliquot.ID, Sample.Num, Batch, Platform)
grn.sequencing.v4 <- left_join(v4.grn.df, all.sequencing) %>% left_join(all.patient) %>% select(PIDN, GRN1, GRN2, Categorized.Diagnosis, Diagnosis, AgeAtDraw, Gender, Summary:Uncategorized, Sequencing.Site, Center, Gene:Sequencing.ID, Aliquot.ID, Sample.Num, Batch, Platform)

wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "v1", gridLines = TRUE)
addWorksheet(wb = wb, sheetName = "v2", gridLines = TRUE)
addWorksheet(wb = wb, sheetName = "v3", gridLines = TRUE)
addWorksheet(wb = wb, sheetName = "v4", gridLines = TRUE)
writeDataTable(wb = wb, sheet = 1, x = grn.sequencing.v1)
writeDataTable(wb = wb, sheet = 2, x = grn.sequencing.v2)
writeDataTable(wb = wb, sheet = 3, x = grn.sequencing.v3)
writeDataTable(wb = wb, sheet = 4, x = grn.sequencing.v4)
sig.pvalues <- createStyle(fontColour = "red")
conditionalFormatting(wb, 1, cols = 2, rows = 1:nrow(grn.sequencing.v1), rule = "<-2.0", style = sig.pvalues)
conditionalFormatting(wb, 2, cols = c(2,3), rows = 1:nrow(grn.sequencing.v1), rule = "<-2.0", style = sig.pvalues)
conditionalFormatting(wb, 3, cols = c(2,3), rows = 1:nrow(grn.sequencing.v1), rule = "<-2.0", style = sig.pvalues)
conditionalFormatting(wb, 4, cols = c(2,3), rows = 1:nrow(grn.sequencing.v1), rule = "<-2.0", style = sig.pvalues)
setColWidths(wb, 1, cols = 1:ncol(grn.sequencing.v1), widths = "auto")
setColWidths(wb, 2, cols = 1:ncol(grn.sequencing.v1), widths = "auto")
setColWidths(wb, 3, cols = 1:ncol(grn.sequencing.v1), widths = "auto")
setColWidths(wb, 4, cols = 1:ncol(grn.sequencing.v1), widths = "auto")
freezePane(wb, 1, firstRow = TRUE)
freezePane(wb, 2, firstRow = TRUE)
freezePane(wb, 3, firstRow = TRUE)
freezePane(wb, 4, firstRow = TRUE)
saveWorkbook(wb, "grn.sequencing.xlsx", overwrite = TRUE) 

#Missing samples
hillblom.missing <- filter(hillblom.request, !(PIDN %in% hillblom.shared)) %>% arrange(RIN)
hillblom.badrin <- filter(hillblom.missing, RIN < 7)
hillblom.goodrin <- filter(hillblom.missing, RIN >= 7)

write.xlsx(hillblom.found, "hillblom.found.xlsx")
write.xlsx(hillblom.badrin, "hillblom.badrin.xlsx")
write.xlsx(hillblom.goodrin, "hillblom.goodrin.xlsx")

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

