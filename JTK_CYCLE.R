require("MetaCycle")
require("ggplot2")
require(xtable)
require(pROC)

preciznostPriklicNeznanaPeriodaJTK <- function(analiza, podatki) {
  # izracuna preciznost in priklic opravljenih napovedi
  # precision: number of correctly predicted out of all predicted
  # recall: number of correctly predicted out of the number of actual
  analiza <- cbind(analiza, podatki$Cosine.Periods)
  vsiPer1 <- sum(analiza[,7]==1)
  vsiPer24 <- sum(analiza[,7]==24)
  vsiKlasNiPer <- sum(analiza[,2]>0.05)
  vsiKlasJe24Per <- sum(analiza[,2] <= 0.05 & analiza[,4] == 24)
  vsiKlasNi24Per <- sum(analiza[,2] <= 0.05 & analiza[,4] != 24)
  
  prec1 <- sum(analiza[,7]==1 & analiza[,2]>0.05) / vsiKlasNiPer
  rec1  <- sum(analiza[,7]==1 & analiza[,2]>0.05) / vsiPer1
  prec2 <- sum(analiza[,2] <= 0.05 & analiza[,4] == 24 & analiza[,7] == 24) / vsiKlasJe24Per
  rec2  <- sum(analiza[,2] <= 0.05 & analiza[,4] == 24 & analiza[,7] == 24) / vsiPer24
  
  tp1 <- paste(sum(analiza[,7]==1 & analiza[,2]>0.05), "/" ,vsiKlasNiPer, "=", round(prec1,3))
  tr1 <- paste(sum(analiza[,7]==1 & analiza[,2]>0.05), "/", vsiPer1, "=", round(rec1,3))
  tp2 <- paste(sum(analiza[,2] <= 0.05 & analiza[,4] == 24 & analiza[,7] == 24), "/", vsiKlasJe24Per, "=", round(prec2, 3))
  tr2 <- paste(sum(analiza[,2] <= 0.05 & analiza[,4] == 24 & analiza[,7] == 24), "/", vsiPer24, "=", round(rec2,3))
  
  stat <- data.frame()
  stat <- rbind(stat, c(tp1, tr1, tp2, tr2))
  colnames(stat) <- c("preciznostC1", "priklicC1", "preciznostC2", "priklicC2")
  return(stat)
}

dolociTNTPJTK <- function(analiza, podatki) {
  # doloci stevilo pravih pozitivnih in negativnih ter 
  # laznih pozitivnih in negativnih opravljenih napovedi
  analiza$stat <- 0
  analiza$truePeriod <- podatki$Cosine.Periods
  analiza$trueAmplitude <- podatki$Cosine.Amplitude
  stetje <- c(0,0,0,0)
  for (i in 1:25) {
    TP <- analiza$BH.Q[i] <= 0.05 & analiza$truePeriod[i] == 24
    FP <- analiza$BH.Q[i] <= 0.05 & analiza$truePeriod[i] != 24
    TN <- analiza$BH.Q[i] > 0.05 & analiza$truePeriod[i] != 24
    FN <- analiza$BH.Q[i] > 0.05 & analiza$truePeriod[i] == 24
    if(TP) {analiza$stat[i] <- "TP"; stetje[1] <- stetje[1] + 1 }
    if(TN) {analiza$stat[i] <- "TN"; stetje[2] <- stetje[2] + 1}
    if(FP) {analiza$stat[i] <- "FP"; stetje[3] <- stetje[3] + 1}
    if(FN) {analiza$stat[i] <- "FN"; stetje[4] <- stetje[4] + 1}
    if(analiza$stat[i] == 0) {print("NAPAKA V TNTPJTK: nobena kategorija ustrezna!")}
    analiza$ampVstrueAmp[i] <- (analiza$AMP[i] / analiza$trueAmplitude[i])*100
  }
  analiza$CycID <- NULL
  analiza$ADJ.P <- NULL
  analiza$LAG <- NULL
  analiza$PER <- as.character(analiza$PER)
  return(list(analiza, stetje))
}

izvoziLatexPrecPriklic <- function(stat, ime) {
  # ustvari in izvozi tabelo preciznosti in priklica v latex formatu
  precrec <- stat[1:2]
  toPrint <- xtable(precrec, type = "latex", caption = "Povzetek pravilnosti klasifikacij za JTK\\_CYCLE.", label = "RainSk", digits = 3, align = "|c|c|c|")
  names(toPrint) <- c("Preciznost neperiodicnih podatkov", "Priklic neperiodicnih podatkov")
  print(toPrint, file = paste(ime, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
  
  precrec <- stat[3:4]
  toPrint <- xtable(precrec, type = "latex", caption = "Povzetek pravilnosti klasifikacij za JTK\\_CYCLE", label = "RainSk", digits = 3, align = "|c|c|c|")
  names(toPrint) <- c("Preciznost periodicnih podatkov", "Priklic periodicnih podatkov")
  print(toPrint, file = paste(ime, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
  
}

izvoziLatexJTK <- function(analiza, podatki, imeZaDatoteko) {
  # pretvori in izvozi analizo z znanimi periodami v latex format
  procesirano <- dolociTNTPJTK(analiza, podatki)[[1]]
  procesirano <- procesirano[,c(1,2,5,7)]
  procesirano$stat <- NULL
  toPrint <- xtable(procesirano, type = "latex", caption = "Rezultati JTK\\_CYCLE.", label = "JTKStab", digits = 3, align = "|c|c|c|c|C{2cm}|")
  names(toPrint) <- c("p-vrednost", "predvidena perioda", "prava perioda", "razmerje amplitud [%]")
  print(toPrint, file = paste(imeZaDatoteko, "tabela.tex", sep = ""), table.placement = "H")
  
  stats <- t(data.frame(as.character(dolociTNTPJTK(analiza, podatki)[[2]])))
  toPrint2 <- xtable(stats, type = "latex", caption = "Povzetek pravilnosti klasifikacij za JTK\\_CYCLE.", label = "rainSk", digits = 3, align = "|c|c|c|c|c|")
  names(toPrint2) <- c("TP", "TN", "FP", "FN")
  print(toPrint2, file = paste(imeZaDatoteko, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
}

izvoziLatexJTKNP <- function(analiza, podatki, imeZaDatoteko) {
  # pretvori in izvozi analizo z neznanimi periodami v latex format
  procesirano <- dolociTNTPJTK(analiza, podatki)[[1]]
  procesirano <- procesirano[,c(1,2,5,3,6,7)]
  toPrint <- xtable(procesirano, type = "latex", caption = "Rezultati JTK\\_CYCLE.", label = "JTKStab", digits = 3, align = "|C{0.5cm}|C{1.5cm}|C{2cm}|C{1.2cm}|C{2cm}|C{1.5cm}|C{1.8cm}|")
  names(toPrint) <- c("p-vrednost", "predvidena perioda", "prava perioda", "predvidena amplituda", "prava amplituda", "razmerje amplitud [%]")
  print(toPrint, file = paste(imeZaDatoteko, "tabela.tex", sep = ""), table.placement = "H")
  
  stats <- t(data.frame(as.character(dolociTNTPJTK(analiza, podatki)[[2]])))
  toPrint2 <- xtable(stats, type = "latex", caption = "Povzetek pravilnosti klasifikacij za JTK\\_CYCLE.", label = "rainSk", digits = 3, align = "|c|c|c|c|c|")
  names(toPrint2) <- c("TP", "TN", "FP", "FN")
  print(toPrint2, file = paste(imeZaDatoteko, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
}

drawROC <- function(analiza, podatki) {
  # narise graf ROC
  x <- as.numeric( podatki$Cosine.Periods == 1)
  analiza$class <- x
  analiza <- analiza[order(analiza$BH.Q),]
  roc(data = analiza, response = class, predictor = BH.Q, plot = TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
      print.auc=TRUE, show.thres=TRUE)
}

#
# MAIN
#

# S1
podatki <- read.csv("Scenarij 1.csv")
vhodJTK <- podatki[,1:25]
vhodJTK[,1] <- as.character(vhodJTK[,1])
write.csv(vhodJTK, file = "vhodJTK1.csv", row.names = FALSE)
# ZNANA
analizaJTK <- meta2d(infile = "vhodJTK1.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,48, by=2)[-1], minper = 24, maxper = 24)
drawROC(analizaJTK$JTK, podatki)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S1-JTK-ZP")
# NEZNANA
analizaJTK <- meta2d(infile = "vhodJTK1.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,48, by=2)[-1], minper = 20, maxper = 28)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S1-JTK-NP")
statistika <- preciznostPriklicNeznanaPeriodaJTK(analizaJTK$JTK, podatki)
izvoziLatexPrecPriklic(statistika, "S1-JTK-NP")


# S2
podatki <- read.csv("Scenarij 2.csv")
vhodJTK <- podatki[,1:13]
vhodJTK[,1] <- as.character(vhodJTK[,1])
write.csv(vhodJTK, file = "vhodJTK2.csv", row.names = FALSE)
# ZNANA
analizaJTK <- meta2d(infile = "vhodJTK2.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,48, by=4)[-1], minper = 24, maxper = 24)
drawROC(analizaJTK$JTK, podatki)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S2-JTK-ZP")
# NEZNANA
analizaJTK <- meta2d(infile = "vhodJTK2.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,48, by=4)[-1], minper = 20, maxper = 28)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S2-JTK-NP")
statistika <- preciznostPriklicNeznanaPeriodaJTK(analizaJTK$JTK, podatki)
izvoziLatexPrecPriklic(statistika, "S2-JTK-NP")


# S3
podatki <- read.csv("Scenarij 3.csv")
vhodJTK <- podatki[,1:13]
vhodJTK[,1] <- as.character(vhodJTK[,1])
write.csv(vhodJTK, file = "vhodJTK3.csv", row.names = FALSE)
# ZNANA
analizaJTK <- meta2d(infile = "vhodJTK3.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,24, by=2)[-1], minper = 24, maxper = 24)
drawROC(analizaJTK$JTK, podatki)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S3-JTK-ZP")
# NEZNANA
analizaJTK <- meta2d(infile = "vhodJTK3.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,24, by=2)[-1], minper = 20, maxper = 28)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S3-JTK-NP")
statistika <- preciznostPriklicNeznanaPeriodaJTK(analizaJTK$JTK, podatki)
izvoziLatexPrecPriklic(statistika, "S3-JTK-NP")

# S4
podatki <- read.csv("Scenarij 4.csv")
vhodJTK <- podatki[,1:7]
vhodJTK[,1] <- as.character(vhodJTK[,1])
write.csv(vhodJTK, file = "vhodJTK4.csv", row.names = FALSE)
# ZNANA
analizaJTK <- meta2d(infile = "vhodJTK4.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,24, by=4)[-1], minper = 24, maxper = 24)
drawROC(analizaJTK$JTK, podatki)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S4-JTK-ZP")
# NEZNANA
analizaJTK <- meta2d(infile = "vhodJTK4.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,24, by=4)[-1], minper = 20, maxper = 28)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S4-JTK-NP")
statistika <- preciznostPriklicNeznanaPeriodaJTK(analizaJTK$JTK, podatki)
izvoziLatexPrecPriklic(statistika, "S4-JTK-NP")


# S5
podatki <- read.csv("Scenarij 5.csv")
vhodJTK <- podatki[,1:7]
vhodJTK[,1] <- as.character(vhodJTK[,1])
write.csv(vhodJTK, file = "vhodJTK5.csv", row.names = FALSE)
# ZNANA
analizaJTK <- meta2d(infile = "vhodJTK5.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,48, by=8)[-1], minper = 24, maxper = 24)
drawROC(analizaJTK$JTK, podatki)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S5-JTK-ZP")
# NEZNANA
analizaJTK <- meta2d(infile = "vhodJTK5.csv", outdir = "JTK", filestyle = "csv", cycMethod = "JTK", 
                     outputFile = FALSE, timepoints = seq(0,48, by=8)[-1], minper = 20, maxper = 28)
izvoziLatexJTK(analizaJTK$JTK, podatki, "S5-JTK-NP")
statistika <- preciznostPriklicNeznanaPeriodaJTK(analizaJTK$JTK, podatki)
izvoziLatexPrecPriklic(statistika, "S5-JTK-NP")
