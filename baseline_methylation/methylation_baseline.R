library(RnBeads)
#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(RColorBrewer)

#Reading and writing tables
library(readr)
library(openxlsx)

#For DE analysis
library(matrixStats)
library(R.utils)
library(unixtools)

#For batch correction and PEER
library(sva)

#Data arrangement
library(dplyr)
library(parallel)

#Functional programming
library(magrittr)
library(purrr)
library(tidyr)
library(broom)

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

wource("../../FRDA project/common_functions.R")

sample.data <- read.xlsx("./methy_pheno.xlsx")
sample.data$barcode <- paste(sample.data$Slide, sample.data$Array, sep = "_")

write_csv(sample.data, "./methy_pheno.csv")

idat.dir <- "./IDAT_AD"
sample.annotation <- "./methy_pheno.csv"
report.dir <- "./reports"

rnb.initialize.reports(report.dir)
logger.start(fname = NA)
parallel.setup(7)
options(fftempdir="~/tmp/Rtmp")
set.tempdir("~/tmp/Rtmp")
rnb.options(disk.dump.big.matrices=TRUE)
rnb.options(logging.disk=TRUE)
rnb.options(enforce.memory.management=TRUE)

patient.annotation <- read.xlsx("../phenotypedata/dan_all_patients_20160726T190552.xlsx")
patient.PIDNs <- select(patient.annotation, PIDN, Gender)
colnames(patient.PIDNs)[2] <- "Sex"
patient.merge <- join(pheno(rnb.set), patient.PIDNs)

data.source <- c(idat.dir, sample.annotation)
result <- rnb.execute.import(data.source = data.source, data.type = "infinium.idat.dir")
saveRDS.gz(result, "result.rda")
rnb.set <- result
#save.rnb.set(rnb.set, "./save/rnb.set.rda", FALSE)
#rnb.set <- load.rnb.set("./save/rnb.set.rda")

rnb.set <- addPheno(rnb.set, patient.merge$Sex, "Sex") 
remove.diagnosis <- grepl("AD|CONTROL", pheno(rnb.set)$Diagnosis) 
remove.diagnosis.key <- which(!remove.diagnosis)
remove.unknown.sex <- is.na(pheno(rnb.set)$Sex) %>% which
remove.unknown.age <- is.na(pheno(rnb.set)$AGE) %>% which
remove.all <- c(remove.diagnosis.key, remove.unknown.sex, remove.unknown.age) %>% unique
rnb.known <- remove.samples(rnb.set, as.character(remove.all))
pheno.export <- pheno(rnb.known)
write.xlsx(pheno.export, "./pheno.export.xlsx")

pheno.stratify <- read.xlsx("./phenodata.methy.adjusted.xlsx")
pheno.drop <- pheno(rnb.known)$PIDN %in% pheno.stratify$PIDN
rnb.stratify <- remove.samples(rnb.known, !pheno.drop)
saveRDS.gz(rnb.stratify, "./save/rnb_stratify.rda")

rnb.run.qc(rnb.stratify, report.dir)

rnb.filter <- rnb.execute.context.removal(rnb.stratify)$dataset
rnb.filter <- rnb.execute.snp.removal(rnb.filter, snp = "any")$dataset
rnb.filter <- rnb.execute.sex.removal(rnb.filter)$dataset

rnb.greedy <- rnb.execute.greedycut(rnb.filter)
saveRDS.gz(rnb.greedy, "./rnb.greedy.rda")
filter.sites <- rnb.greedy$sites

rnb.filter <- remove.sites(rnb.filter, filter.sites)
#rnb.filter.save <- tempfile(pattern = "filter.save")
#save(rnb.filter, file = rnb.filter.save, compress = FALSE)
#rnb.test <- reload(rnb.filter, save.file = "/home/daniel/tmp/Rtmp")

rnb.filter <- rnb.execute.na.removal(rnb.filter)$dataset
gc()
rnb.filter <- rnb.execute.variability.removal(rnb.filter, 0.005)$dataset
gc()
saveRDS.gz(rnb.filter, "./save/rnb.filter.rda")

rnb.norm <- rnb.execute.normalization(rnb.filter, method = "bmiq", bgcorr.method = "methylumi.noob", verbose = TRUE)
saveRDS.gz(rnb.norm, "./save/rnb.norm.rda")
sites.norm <- meth(rnb.norm)
saveRDS.gz(sites.norm, "./save/sites.norm.rda")
promoters.norm <- meth(rnb.norm, type = "promoters")
saveRDS.gz(promoters.norm, "./save/promoters.norm.rda")
genes.norm <- meth(rnb.norm, type = "genes")
saveRDS.gz(genes.norm, "./save/genes.norm.rda")

pheno.final <- pheno(rnb.norm)
pheno.final$Diagnosis %<>% factor(levels = c("CONTROL", "AD")) %>% droplevels
pheno.final$Sex %<>% factor %>% droplevels
model.design <- model.matrix( ~ Diagnosis + Sex + AGE, pheno.final )
promoters.combat <- ComBat(dat = promoters.norm, batch = pheno.final$Batch, mod = model.design[,-1])
genes.combat <- ComBat(dat = genes.norm, batch = pheno.final$Batch, mod = model.design[,-1])

rownames(promoters.combat) <- rownames(annotation(rnb.norm, type = "promoters"))
rownames(genes.combat) <- rownames(annotation(rnb.norm, type = "genes"))

promoters.pgrn <- promoters.combat["ENSG00000030582",]
genes.pgrn <- genes.combat["ENSG00000030582",]

promoters.pgrn.df <- data.frame(PGRN = promoters.pgrn, select(pheno.final, Diagnosis, Sex, AGE))
write.xlsx(promoters.pgrn.df, "pgrn_promoter.xlsx")
genes.pgrn.df <- data.frame(PGRN = genes.pgrn, select(pheno.final, Diagnosis, Sex, AGE))
write.xlsx(genes.pgrn.df, "pgrn_gene.xlsx")

diagnosis.sex <- lm(as.integer(Diagnosis) ~ Sex, pheno.final) %>% anova %>% tidy
diagnosis.age <- lm(as.integer(Diagnosis) ~ AGE, pheno.final) %>% anova %>% tidy

pgrn.promoter.anova <- lm(PGRN ~ Diagnosis + Sex + AGE, promoters.pgrn.df) %>% anova %>% tidy
pgrn.promoter.aov <- aov(PGRN ~ Sex, promoters.pgrn.df) %>% TukeyHSD %>% tidy
pgrn.gene.anova <- lm(PGRN ~ Diagnosis + Sex + AGE, genes.pgrn.df) %>% anova %>% tidy

STOP
#dred.sites <- rnb.execute.dreduction(rnb.norm)
#dred.promoters <- rnb.execute.dreduction(rnb.norm, target = "promoters")
#dred <- list(sites = dred.sites, promoters = dred.promoters)
#pca.colors <- ifelse(pheno(rnb.norm)$Recovered == TRUE, "red", "blue")

#dred.plot <- data.frame(dred.promoters$mds$euclidean[,1:2])
#colnames(dred.plot) <- c("PCA1", "PCA2")
#dred.plot %<>% mutate(Recovered = pheno(rnb.norm)$Recovered)
##Convert this plot ggplot
#CairoPDF("pca_allsites", width = 6, height = 6)
#p <- ggplot(data = dred.plot, aes(x = PCA1, y = PCA2, col = Recovered)) + geom_point()
#print(p)
#dev.off()

#HBP.fix <- pheno(rnb.norm)$HBP %>% droplevels
#rnb.norm <- addPheno(rnb.norm, HBP.fix, "HBP.fix")
#rnb.options(exploratory.columns = c("Age", "Sex", "Ethnicity", "HBP.fix", "BMI", "Recovered"))
#assoc <- rnb.execute.batcheffects(rnb.norm, pcoordinates = dred)
#assoc.qc <- rnb.execute.batch.qc(rnb.norm, pcoordinates = dred)

#clustering.sites <- rnb.execute.clustering(rnb.norm, region.type = "sites")
#clustering.promoters <- rnb.execute.clustering(rnb.norm, region.type = "promoters")

#promoters.beta <- meth(rnb.norm, type = "promoters")
#promoters.m <- lumi::beta2m(meth(rnb.norm, type = "promoters"))
#sites.beta <- meth(rnb.norm)
#sites.m <- lumi::beta2m(meth(rnb.norm))

#gen.heatmap(promoters.beta, 1000, "promoters_beta_heatmap", clustering.promoters)
#gen.heatmap(promoters.m, 1000, "promoters_mvalue_heatmap", clustering.promoters)
#gen.heatmap(sites.beta, 1000, "sites_beta_heatmap", clustering.sites)
#gen.heatmap(sites.m, 1000, "sites_mvalue_heatmap", clustering.sites)

#gen.heatmap <- function(meths, ngenes, file.name, cluster.object)
#{
    #sites.ordered <-  apply(meths, 1, mad) %>% order(decreasing = TRUE)
    #sites.plot <- meths[sites.ordered[1:ngenes],]
    #cluster.tree <- cluster.object[[7]]@result
    #attr(cluster.tree, "class") <- "hclust"
    #cluster.dendro <- as.dendrogram(cluster.tree)

    #CairoPDF(file.name, width = 10, height = 10)
    #heatmap.2(sites.plot, Rowv = TRUE, Colv = cluster.dendro, dendrogram = "both", scale = "none", trace = "none", labRow = FALSE)
    #dev.off()
#}

#rnb.options("covariate.adjustment.columns" = c("Age", "Sex", "Ethnicity", "HBP.fix", "BMI"))
#comp.cols <- "Recovered"
#reg.types <- c("genes", "promoters")
#diffmeth.adj <- rnb.execute.computeDiffMeth(rnb.norm, comp.cols, region.types = reg.types)

#comparison <- get.comparisons(diffmeth.adj)
#tab.sites <- get.table(diffmeth.adj, comparison, "sites", return.data.frame = TRUE) 
#tab.promoters <- get.table(diffmeth.adj, comparison, "promoters", return.data.frame = TRUE)

#promoters.recovered <- promoters.beta[,pheno(rnb.norm)$Recovered == TRUE] 
#promoters.coef <- promoters.recovered - tab.promoters$mean.mean.g1
#sites.cutoff <- sort(tab.sites$combinedRank)[1000]
#promoters.cutoff <- sort(tab.promoters$combinedRank)[1000]
#promoters.500 <- sort(tab.promoters$combinedRank)[500]

#exportDMRs2regionFile(rnb.norm, diffmeth.adj, "./sites_top1000.bed", comparison, "sites", rank.cut = sites.cutoff)
#exportDMRs2regionFile(rnb.norm, diffmeth.adj, "./promoters_top1000.bed", comparison, "promoters", rank.cut = promoters.cutoff)

#promoters.key <- which(tab.promoters$combinedRank <= promoters.cutoff)
#promoters.top1000 <- promoters.coef[promoters.key,]

#promoters.500.key <- which(tab.promoters$combinedRank <= promoters.500)
#promoters.top500 <- promoters.coef[promoters.500.key,]

##Generate anova heatmaps
#gen.anova.heatmap <- function(file.name, dataset)
#{ 
    #CairoPDF(file.name, width = 10, height = 10)
    #heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 1.0, labRow = "", keysize = 0.9)
    #dev.off()
#}

#gen.anova.heatmap("top1000_DMR", promoters.top1000)

#thresholds <- c(0.01, 0.005, 0.001)
#promoters.threshold <- map(thresholds, get.sizes, tab.promoters)
#names(promoters.threshold) <- thresholds
#threshold.df <- melt(promoters.threshold)
#colnames(threshold.df) <- c("mysum", "Direction", "Test")
#threshold.df$Test %<>% factor
#levels(threshold.df$Test) <- c("0.01", "0.005", "0.001")

#get.sizes <- function(p.val, dataset)
#{
    #dataset.sig <- dataset[dataset$comb.p.val < p.val,]
    #dataset.up <- dataset.sig[dataset.sig$mean.mean.quot.log2 > 0,]
    #dataset.down <- dataset.sig[dataset.sig$mean.mean.quot.log2 < 0,]
    #return(list(positive = dim(dataset.up)[1], negative = -(dim(dataset.down)[1])))
#}

#gen.decideplot("threshold_selection", threshold.df)
#gen.decideplot <- function(file.name, decide.plot)
#{
    #p <- ggplot()
    #p <- p + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = Test, y = mysum), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    #p <- p + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Test, y = mysum, ymax = max(mysum) + 60, hjust = -1.1, label = mysum), position = position_dodge(width = 1))
    #p <- p + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = Test, y = mysum), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    #p <- p + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Test, y = mysum, ymax = min(mysum) - 60, hjust = 1.5, label = abs(mysum)), position = position_dodge(width = 1))
    ##p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    #p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    #p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0)) + ylab("Differentially Methylated Promoters")
    #CairoPDF(file.name, width = 6, height = 7)
    #print(p)
    #dev.off()
#}

#rownames(promoters.beta) <- rownames(annotation(rnb.norm, type = "promoters"))
#saveRDS.gz(promoters.beta, "./save/promoters_beta.rda")
#promoters.annotation <- annotation(rnb.norm, type = "promoters")
#saveRDS.gz(promoters.annotation, "./save/promoters_annotation.rda")

#joined.table <- data.frame(tab.promoters, promoters.annotation)

#top.1001 <- which(tab.promoters$combinedRank <= promoters.cutoff)
#annot.1000 <- annotation(rnb.norm, type = "promoters")[top.1000,] %>% filter(!is.na(symbol)) %>% filter(!str_detect(symbol, ";"))
#annot.1000$symbol %<>% str_replace("-.*$", "")
##trouble.symbols <- filter(annot.1000, str_detect(symbol, "-"))$symbol
#write.xlsx(annot.1000, "annot_1000.xlsx")
##rnb.run.tnt(rnb.norm)
##annotation(rnb.norm) %>% str
