library(affy)
library(oligo)

require(hugene10sttranscriptcluster.db)
library("pd.hugene.1.0.st.v1")

celFiles <- list.celfiles(listGzipped = TRUE)
affyRaw <- read.celfiles(celFiles)
eset <- rma(affyRaw) # preprocessing
write.exprs(eset,file="data.txt")
my_frame <- data.frame(exprs(eset)) # getting the data


# con <- db(pd.hugene.1.0.st.v1)
# types <- dbGetQuery(con, "select * from featureSet;")

# ******************************************************annotation with lookup
Annot <- data.frame(ACCNUM=sapply(contents(hugene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene10sttranscriptclusterGENENAME), paste, collapse=", "))
all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)


Annot2 <- data.frame(ENSEMBL=sapply(contents(hugene10sttranscriptclusterENSEMBL), paste, collapse=", "), 
                    SYMBOL=sapply(contents(hugene10sttranscriptclusterSYMBOL), paste, collapse=", "), 
                    DESC=sapply(contents(hugene10sttranscriptclusterGENENAME), paste, collapse=", "),
                    UNIPROT=sapply(contents(hugene10sttranscriptclusterUNIPROT), paste, collapse=", ")
                    )
write.table(Annot2,file="anotacije.txt",sep="\t") # for ground truth table


write.table(all,file="data_all_annotated.txt",sep="\t")
# ******************************************************annotation with annotateEset
library(affycoretools)
esetAn <- annotateEset(eset, hugene10sttranscriptcluster.db) # 22122 jih je anotiranih
esetAn <- annotateEset(eset, pd.hugene.1.0.st.v1) # 25293 jih je anotiranih
q <- fData(esetAn)$GENENAME
q <- as.character(q)
sum(!is.na(q))


Annot2 <- as.data.frame(fData(esetAn))
my_frame <- cbind(row.names(my_frame), my_frame)
colnames(my_frame)[1] <- "PROBEID"
allwithCore <- merge.data.frame(x = my_frame, y = Annot2, by.x="PROBEID", by.y ="PROBEID")
# *****************************************************annotation with file, gene_assignment=25293
fileIN <- read.csv(file = "GPL6244-17930.txt", sep = "\t", header = TRUE, skip = 12)
q <- as.character(allwithFile$gene_assignment)
withF1 <- which(q == "---")

# drugi file, ta ima manj izpolnjenih mest kot prvi
fileIN <- read.csv(file = "GPL6244.annot", sep = "\t", header = TRUE, skip = 27) # Gene.title = 22195
q <- as.character(allwithFile$Gene.title)
withF2 <- which(q == "")
# ---

my_frame <- data.frame(exprs(eset))
my_frame <- cbind(row.names(my_frame), my_frame)
colnames(my_frame)[1] <- "ID"
allwithFile <- merge.data.frame(x = my_frame, y = fileIN, by.x="ID", by.y ="ID")

# ***************************************************** QC
# output od rma je že v log2 intenziteti, zato uporabimo identity
boxplot(eset,  transfo=identity, xlab="", ylab="", main="") #graf intenzitet èipov in varianc

hist(eset, transfo=identity, xlab="", ylab="", main="") # distribucija po gostoti za izmerjene vrednosti

read.affybatch(celFiles)
CLLbatch <- read.affybatch(celFiles)
dataPLM <- fitPLM(CLLbatch)
boxplot(dataPLM, main = "NUSE", ylim = c(0.95, 1.22), outline = FALSE, col = "lightblue", las = 3, whisklty = 0, staplelty = 0)
Mbox(dataPLM, main = "RLE", ylim = c(-0.4, 0.4), outline = FALSE, col = "mistyrose", las = 3, whisklty = 0, staplelty = 0)

# intensity before/after normalization
oligo::boxplot(affyRaw, target = "core", 
               main = "", col = "#00000000")

oligo::boxplot(eset, target = "core", 
               main = "", col = "#00000000")

library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = affyRaw,
                    outdir = getwd(),
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("Factor.Value.disease.", "Factor.Value.phenotype."))


palmieri_medians <- rowMedians(Biobase::exprs(eset))

hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = 4, col = "coral4", lwd = 2)


heatmap(exprs(eset[1:100,])) # heatmap za prvih 100 sond
newnames <- c()
for (i in 1:length(rownames(my_frame))) {
  newnames[i] <- paste(rownames(my_frame)[i], all$DESC[i])
}
rownames(my_frame) <- newnames
colnames(my_frame) <- 1:48
heatmap(as.matrix(my_frame[c(33287,33288,33290,33294,33296, 33263,33270),]))

plotdata <- data.frame(t(my_frame[21570,]))
colnames(plotdata)[1] <- "data"
ggplot() + geom_point(data = plotdata, mapping = aes(y = data, x = 1:48)) + xlab("Ura") + ylab("Vrednost") + ggtitle("Sonde 21570 v odvisnosti od èasa")
  
# *****************************************************ANALYSIS
vhodAnaliza <- all[,5:52]
require("MetaCycle")
write.csv(cbind(as.character(all[,1]), vhodAnaliza), file = "vhodJTK1.csv", row.names = FALSE)
analizaJTK <- meta2d(infile = "vhodJTK1.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,48)[-1], minper = 20, maxper = 28)
write.csv(analizaJTK$JTK, file = "analizaJTKResults.csv", row.names = FALSE)
analizaJTK24 <- meta2d(infile = "vhodJTK1.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,48)[-1], minper = 24, maxper = 24)

require("rain")
analizaRain <- rain(t(vhodAnaliza), period = 24, period.delta = 4, deltat = 1, verbose = TRUE, adjp.method = "BH")
write.csv(analizaRain, file = "analizaRAINResults.csv", row.names = FALSE)
analizaRain24 <- rain(t(vhodAnaliza), period = 24, period.delta = 0, deltat = 1, verbose = TRUE, adjp.method = "BH")

require("cosinor2")
analizaCosinor <- data.frame()
cosPval <- c()
dolociPrimernoPeriodo <- function(vhodCosinor, cas) {
  analiza <- periodogram(vhodCosinor, time = cas, periods = c(20:30))
  periode <- analiza$data$periodogram
  indeks <- which.max(periode)
  #print(paste(indeks, cas[indeks]))
  return(cas[indeks]+19)
  # vse v data$periodogram
}
for (i in 1:length(vhodAnaliza[,1])) {
  perioda <- dolociPrimernoPeriodo(vhodAnaliza[i,], seq(0,48)[-1])
  mydf <- data.frame(seq(0,48)[-1], t(vhodAnaliza[i,]))
  colnames(mydf) <- c('time', 'y')
  fit.cosinor <- cosinor.lm(y ~ time(time), period = perioda, data = mydf)
  pVal <- cosinor.detect(fit.cosinor)[4]
  cosPval <- rbind(cosPval, pVal)
  analizaCosinor <- rbind(analizaCosinor, c(fit.cosinor$coefficients["amp"], pVal, perioda))
}
colnames(analizaCosinor) <- c("Amplituda", "p-vrednost", "perioda")
pBH <- p.adjust(as.numeric(analizaCosinor$`p-vrednost`), method = "BH")
analizaCosinor <- cbind(analizaCosinor, pBH)
write.csv(analizaCosinor, file = "analizaCosinorResults.csv", row.names = FALSE)
# 24 urna perioda
#
#
analizaCosinor24 <- data.frame()
cosPval24 <- c()
for (i in 1:length(vhodAnaliza[,1])) {
  mydf <- data.frame(seq(0,48)[-1], t(vhodAnaliza[i,]))
  colnames(mydf) <- c('time', 'y')
  fit.cosinor <- cosinor.lm(y ~ time(time), period = 24, data = mydf)
  pVal <- cosinor.detect(fit.cosinor)[4]
  cosPval24 <- rbind(cosPval24, pVal)
  analizaCosinor24 <- rbind(analizaCosinor24, c(fit.cosinor$coefficients["amp"], pVal, perioda))
}
colnames(analizaCosinor24) <- c("Amplituda", "p-vrednost", "perioda")
pBH <- p.adjust(as.numeric(analizaCosinor24$`p-vrednost`), method = "BH")
analizaCosinor24 <- cbind(analizaCosinor24, pBH)
# *****************************************************
library(ggplot2)
ggplot(analizaCosinor[analizaCosinor$pBH <= 0.05,], aes(x=perioda)) + geom_histogram() + xlab("Periode") + ylab("")
ggplot(analizaRain[analizaRain$pVal <= 0.05,], aes(x=period)) + geom_histogram() + xlab("Periode") + ylab("")
ggplot(analizaJTK$JTK[analizaJTK$JTK$BH.Q <= 0.05,], aes(x=PER)) + geom_histogram() + xlab("Periode") + ylab("")

ggplot(analizaCosinor, aes(x=pBH)) + geom_histogram() + xlab("P-vrednosti") + ylab("")
ggplot(analizaCosinor[analizaCosinor$pBH <= 0.05,], aes(x=pBH)) + geom_histogram() + xlab("P-vrednosti") + ylab("")

ggplot(analizaRain, aes(x=pVal)) + geom_histogram() + xlab("P-vrednosti") + ylab("")
ggplot(analizaRain[analizaRain$pVal <= 0.05,], aes(x=pVal)) + geom_histogram() + xlab("P-vrednosti") + ylab("")

ggplot(analizaJTK$JTK, aes(x=BH.Q)) + geom_histogram() + xlab("P-vrednosti") + ylab("")
ggplot(analizaJTK$JTK[analizaJTK$JTK$BH.Q < 1,], aes(x=BH.Q)) + geom_histogram() + xlab("P-vrednosti") + ylab("")

allSelect <- subset(all, DESC != "NA")
write.table(allSelect,file="data.ann.hasDesc.txt",sep="\t")
# -------------------MAKING DIAGRAMS
cosordered <- analizaCosinor
rownames(cosordered) <- all$Row.names # da dobijo ob st serije se st sonde
cosordered <- cosordered[order(cosordered$pBH),]
cosordered <- merge(cosordered, Annot, by.x=0, by.y=0)

rainordered <- analizaRain
rownames(rainordered) <- all$Row.names 
rainordered <- rainordered[order(rainordered$pVal),]
rainordered <- merge(rainordered, Annot, by.x=0, by.y=0)

jtkordered <- analizaJTK$JTK
rownames(jtkordered) <- all$Row.names
jtkordered <- jtkordered[order(jtkordered$BH.Q),]
jtkordered <- merge(jtkordered, Annot, by.x=0, by.y=0)
library(VennDiagram)
sezCos <- (cosordered[cosordered$pBH <= 0.05,])$Row.names
sezRain <- (rainordered[rainordered$pVal <= 0.05,])$Row.names
sezJtk <- (jtkordered[jtkordered$BH.Q <= 0.05, ])$Row.names
write(sezCos, file = "sezCos.txt", sep = "\n")
write(sezJtk, file = "sezjtk.txt", sep = "\n")
write(sezRain, file = "sezrain.txt", sep = "\n")


venn.diagram(
  x = list(sezCos, sezRain, sezJtk),
  category.names = c("Cosinor" , "RAIN" , "JTK_CYCLE"),
  filename = 'VENN2.png',
  imagetype = "png",
  output=TRUE,
  scaled=TRUE
)
overrideTriple=T
ven <- draw.triple.venn(      area1           = 1299,
                              area2           = 4623,
                              area3           = 55,
                              n12             = 1282,
                              n23             = 55,
                              n13             = 55,
                              n123            = 55,
                              category        = c('A', 'B', 'C'),
                              #fill            = c('red', 'blue', 'green'),
                              cat.col         = c('red', 'blue', 'green'),
                              cex = rep(2, 7),
                              #cat.dist = -2,
                              euler.d =         F,
                              scaled          = F
)
jpeg(filename = "Triple_Venn_diagram.jpg", height = 750, width = 800)
grid.draw(ven);
dev.off();

# -------------------
# 24
#
cosordered24 <- analizaCosinor24
rownames(cosordered24) <- all$Row.names # da dobijo ob st serije se st sonde
cosordered24 <- cosordered24[order(cosordered24$pBH),]
cosordered24 <- merge(cosordered24, Annot, by.x=0, by.y=0)

rainordered24 <- analizaRain24
rownames(rainordered24) <- all$Row.names 
rainordered24 <- rainordered24[order(rainordered24$pVal),]
rainordered24 <- merge(rainordered24, Annot, by.x=0, by.y=0)

jtkordered24 <- analizaJTK24$JTK
rownames(jtkordered24) <- all$Row.names
jtkordered24 <- jtkordered24[order(jtkordered24$BH.Q),]
jtkordered24 <- merge(jtkordered24, Annot, by.x=0, by.y=0)
library(VennDiagram)
sezCos24 <- (cosordered24[cosordered24$pBH <= 0.05,])$Row.names
sezRain24 <- (rainordered24[rainordered24$pVal <= 0.05,])$Row.names
sezJtk24 <- (jtkordered24[jtkordered24$BH.Q <= 0.05, ])$Row.names
venn.diagram(
  x = list(sezCos24, sezRain24, sezJtk24),
  category.names = c("Cosinor" , "RAIN" , "JTK_CYCLE"),
  filename = 'VENN24.png',
  imagetype = "png",
  output=TRUE,
  scaled=TRUE
)
# -------------------
# Primerjava z znanimi cirkadialnimi geni
znaniIzWWW <- read.table(file = "ZnaniIzWWW", header = T, sep = ",")
analizaZnani <- data.frame(matrix(0, nrow = length(znaniIzWWW$ID), ncol = 7))
colnames(analizaZnani) <- c("ID", "cosinor", "rain", "jtk", "c24", "r24", "j24")
analizaZnani$ID <- znaniIzWWW$ID
for (i in 1:length(analizaZnani$ID)) {
  if(analizaZnani$ID[i] %in% sezCos)   {analizaZnani$cosinor[i] <- 1}
  if(analizaZnani$ID[i] %in% sezRain)  {analizaZnani$rain[i] <- 1}
  if(analizaZnani$ID[i] %in% sezJtk)   {analizaZnani$jtk[i] <- 1}
  if(analizaZnani$ID[i] %in% sezCos24) {analizaZnani$c24[i] <- 1}
  if(analizaZnani$ID[i] %in% sezRain24){analizaZnani$r24[i] <- 1}
  if(analizaZnani$ID[i] %in% sezJtk24) {analizaZnani$j24[i] <- 1}
  #print(analizaZnani$ID[i] %in% rownames(my_frame))
  #print(which(rownames(my_frame) == analizaZnani$ID[i]))
}

katerigeni <- c()
for (i in analizaZnani$ID) {
  q <- which(all$Row.names == i)
  r <- which(allwithFile$ID == i)
  opisq <- as.character(all$DESC[q])
  opisr <- as.character(allwithFile$gene_assignment[r])
  #print(opisq)
  #print(strsplit(opisr, split = "//")[[1]][3])
  katerigeni <- cbind(katerigeni, paste(opisq, strsplit(opisr, split = "//")[[1]][3], sep = ";"))
}

katerigeni <- t(katerigeni)
analizaZnani <- taRifx::remove.factors(cbind(analizaZnani, katerigeni))
analizaZnani[15,8] <- "control"
analizaZnani[25,8] <- "isocitrate dehydrogenase 1 (NADP+), soluble (IDH1) pseudogene"

write.table(as.character(analizaZnani$ID), file="KRNEKI", quote=T, sep=",", row.names = F, eol = ",")
colnames(my_frame) <- 1:48

heatmap(as.matrix(my_frame[c("8078272","7938563","7897378","7966052","8073422","8014956","8059996","8079693",
                             "7961371","7929438","7993622","8151816","7906433","8019308","7893966","7941743",
                             "8037166","8003204","8118310","8155849","8178086","8179322","8179324","7954789",
                             "8123800","8121043","8102440","7961891","7941457","7939595","8038117","7992049",
                             "8040211","8043909","8162276","8012349","8135069","8124562"),]), 
        Colv = NA, scale = "none", labRow = imena2, margins = c(2,25)) # stevilke vrstic

png(filename = "heatmapZnani.jpg", height = 1500, width = 2500, res = 300)

katerigeni <- analizaZnani$katerigeni
imena2 <- c()
for (i in 1:length(katerigeni)) {
  s <- strsplit(katerigeni[i], split = ";")
  imena2 <- cbind(imena2, trimws(s[[1]][2]))
}
imena2[15] <- "control"
imena2[25] <- "isocitrate dehydrogenase 1 (NADP+), soluble (IDH1) pseudogene"
imena2[37] <- "serpin family E member 1"
imena2[18] <- "GINS complex subunit 2"
imena2[14] <- "MAF bZIP transcription factor G"
imena2[11] <- "ITP receptor interacting protein-like 2"
imena2[37] <- "serpin peptidase inhibitor, clade E"
imena2[16] <- "leucine rich repeat and fibronectin type III domain"
heatmap.2(as.matrix(my_frame[c("8078272","7938563","7897378","7966052","8073422","8014956","8059996","8079693",
                               "7961371","7929438","7993622","8151816","7906433","8019308","7893966","7941743",
                               "8037166","8003204","8118310","8155849","8178086","8179322","8179324","7954789",
                               "8123800","8121043","8102440","7961891","7941457","7939595","8038117","7992049",
                               "8040211","8043909","8162276","8012349","8135069","8124562"),]), trace = "none", 
          Colv = F, dendrogram = "row", key = F, labRow = imena2, margins = c(2,17), cexRow = 1)
dev.off();
library(gplots)

require(xtable)
toPrint <- xtable(analizaZnani[,c(8,2,3,4,5,6,7)], type = "latex", 
                  caption = "Primerjava identificiranih ritmiènih transkriptov metod z zunanjimi rezultati.", 
                  label = "primerjavatrans", digits = 3, 
                  align = "|c|c|c|c|c|c|c|c|", row)
names(toPrint) <- c("Ime transkripta", "Cosinor", "RAIN", "JTK_CYCLE", "Cosinor 24 ur", "RAIN 24 ur", "JTK_CYCLE 24 ur")
print(toPrint, file = "transtab.tex", table.placement = "H", include.rownames=FALSE, tabular.environment="longtable")

# ------------------------------STATISTIKA ZA ANALIZA ZNANI
# TP: oscilirajo v obeh primerih, FP: oscilirajo samo pri meni
# TN: ne oscilirajo v obeh primerih FN: ne oscilirajo samo pri meni

# TP: èe je 1. FP: ne more biti. TN: ne more biti. FN: èe je 0
statZnani <- as.data.frame(matrix(0, 2,6))
colnames(statZnani) <- c("cosinor", "rain", "jtk", "c24", "r24", "j24")
# prva vrstica je TP, druga je FN
for (i in 2:7) {
  print(colnames(analizaZnani)[i])
  statZnani[1,i-1] <- sum(analizaZnani[,i] == 1)
  statZnani[2,i-1] <- sum(analizaZnani[,i] == 0)
}

sum(analizaCosinor[analizaCosinor$pBH <= 0.05,]$perioda)/length(analizaCosinor[analizaCosinor$pBH <= 0.05,]$perioda)
sum(analizaRain[analizaRain$pVal <= 0.05,]$period)/length(analizaRain[analizaRain$pVal <= 0.05,]$period)
sum(analizaJTK$JTK[analizaJTK$JTK$BH.Q <= 0.05,]$PER)/length(analizaJTK$JTK[analizaJTK$JTK$BH.Q <= 0.05,]$PER)