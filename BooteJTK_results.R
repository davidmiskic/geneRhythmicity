# ANALIZA REZULTATOV
# ./BooteJTK-CalcP.py -f example/TestInput4.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x OTHERTEXT -r 2 -z 10
# header tabele mora biti v posebnem formatu, začne se z # ali ID, imena stolpcev pa so številke

require(pROC)
require(xtable)


dolociTNTPboote <- function(analiza, podatki) {
  # določi število pravih pozitivnih in negativnih ter 
  # lažnih pozitivnih in negativnih opravljenih napovedi
  analiza$stat <- 0
  analiza$truePeriod <- podatki$Cosine.Periods
  stetje <- c(0,0,0,0)
  for (i in 1:25) {
    TP <- analiza$GammaBH[i] <= 0.05 & analiza$truePeriod[i] == 24
    FP <- analiza$GammaBH[i] <= 0.05 & analiza$truePeriod[i] != 24
    TN <- analiza$GammaBH[i] > 0.05 & analiza$truePeriod[i] != 24
    FN <- analiza$GammaBH[i] > 0.05 & analiza$truePeriod[i] == 24
    if(TP) {analiza$stat[i] <- "TP"; stetje[1] <- stetje[1] + 1 }
    if(TN) {analiza$stat[i] <- "TN"; stetje[2] <- stetje[2] + 1}
    if(FP) {analiza$stat[i] <- "FP"; stetje[3] <- stetje[3] + 1}
    if(FN) {analiza$stat[i] <- "FN"; stetje[4] <- stetje[4] + 1}
    if(analiza$stat[i] == 0) {print("NAPAKA V TNTPboote: nobena kategorija ustrezna!")}
  }
  analiza <- analiza[,c("GammaBH", "PeriodMean", "truePeriod", "stat")]
  analiza$PeriodMean <- as.character(analiza$PeriodMean)
  rownames(analiza) <- NULL
  return(list(analiza, stetje))
}

izvoziLatexboote <- function(analiza, podatki, imeZaDatoteko) {
  # pretvori in izvozi analizo z znanimi periodami v latex format
  toPrint <- xtable(dolociTNTPboote(analiza, podatki)[[1]], type = "latex", caption = "Rezultati BooteJTK.", label = "bjtkStab", digits = 3, align = "|c|c|c|c|c|")
  names(toPrint) <- c("p-vrednost", "predvidena perioda", "prava perioda", "")
  print(toPrint, file = paste(imeZaDatoteko, "tabela.tex", sep = ""), table.placement = "H")
  
  stats <- t(data.frame(as.character(dolociTNTPboote(analiza, podatki)[[2]])))
  toPrint2 <- xtable(stats, type = "latex", caption = "Povzetek pravilnosti klasifikacij za BooteJTK.", label = "bjtkSk", digits = 3, align = "|c|c|c|c|c|")
  names(toPrint2) <- c("TP", "TN", "FP", "FN")
  print(toPrint2, file = paste(imeZaDatoteko, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
}

drawROC <- function(analiza, podatki) {
  # nariše graf ROC
  x <- as.numeric( podatki$Cosine.Periods == 1)
  analiza$class <- x
  analiza <- analiza[order(analiza$GammaBH),]
  print(analiza)
  roc(data = analiza, response = class, predictor = GammaBH, plot = TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
      print.auc=TRUE, show.thres=TRUE)
}

#
# MAIN
#

# S1
podatki <-read.csv("Scenarij 1.csv")
vhodBJTK <- podatki[,1:25]
stolpciImena <- c("#", as.character(seq(2,48,2)))
vhodBJTK <- rbind(stolpciImena, vhodBJTK)
vhodBJTK[,1] <- as.character(vhodBJTK[,1])
write.table(vhodBJTK, file="s1Tab1.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# preberi izhodno datoteko, in razvrsti v enakem vrstnem redu kot scenarij
s1tab1 <- read.delim("s1Tab1_Vash_SCENARIJ1_boot10-rep1_1_GammaP.txt")
s1ordered <- s1tab1[order(s1tab1$ID),]
drawROC(s1ordered, podatki)
izvoziLatexboote(s1ordered, podatki, "S1-boote-ZP")


# S2
podatki <- read.csv("Scenarij 2.csv")
vhodBJTK <- podatki[,1:13]
stolpciImena <- c("#", as.character(seq(4,48,4)))
vhodBJTK <- rbind(stolpciImena, vhodBJTK)
vhodBJTK[,1] <- as.character(vhodBJTK[,1])
write.table(vhodBJTK, file="s2Tab.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# preberi izhodno datoteko, in razvrsti v enakem vrstnem redu kot scenarij
s1tab1 <- read.delim("s2Tab_Vash_SCENARIJ2_boot10-rep1_GammaP.txt")
s1ordered<- s1tab1[order(s1tab1$ID),]
drawROC(s1ordered, podatki)
izvoziLatexboote(s1ordered, podatki, "S2-boote-ZP")

# S3
podatki <- read.csv("Scenarij 3.csv")
vhodBJTK <- podatki[,1:13]
stolpciImena <- c("#", as.character(seq(2,24,2)))
vhodBJTK <- rbind(stolpciImena, vhodBJTK)
vhodBJTK[,1] <- as.character(vhodBJTK[,1])
write.table(vhodBJTK, file="s3Tab.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# S4
podatki <- read.csv("Scenarij 4.csv")
vhodBJTK <- podatki[,1:7]
stolpciImena <- c("#", as.character(seq(4,24,4)))
vhodBJTK <- rbind(stolpciImena, vhodBJTK)
vhodBJTK[,1] <- as.character(vhodBJTK[,1])
write.table(vhodBJTK, file="s4Tab.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# S5
podatki <- read.csv("Scenarij 5.csv")
vhodBJTK <- podatki[,1:7]
stolpciImena <- c("#", as.character(seq(8,48,8)))
vhodBJTK <- rbind(stolpciImena, vhodBJTK)
vhodBJTK[,1] <- as.character(vhodBJTK[,1])
write.table(vhodBJTK, file="s5Tab.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# preberi izhodno datoteko, in razvrsti v enakem vrstnem redu kot scenarij
s1tab1 <- read.delim("s5Tab_Vash_SCENARIJ5_boot10-rep1_GammaP.txt")
s1ordered<- s1tab1[order(s1tab1$ID),]
drawROC(s1ordered, podatki)
izvoziLatexboote(s1ordered, podatki, "S1-boote-ZP")
