require("rain")
require("ggplot2")
require(pROC)
require(xtable)

preciznostPriklicNeznanaPerioda <- function(analiza, podatki) {
  # izracuna preciznost in priklic opravljenih napovedi
  # precision: number of correctly predicted out of all predicted
  # recall: number of correctly predicted out of the number of actual
  analiza <- cbind(analiza, podatki$Cosine.Periods)
  vsiPer1 <- sum(analiza[,5]==1)
  vsiPer24 <- sum(analiza[,5]==24)
  vsiKlasNiPer <- sum(analiza[,1]>0.05)
  vsiKlasJe24Per <- sum(analiza[,1] <= 0.05 & analiza[,4] == 24)
  vsiKlasNi24Per <- sum(analiza[,1] <= 0.05 & analiza[,4] != 24)
  
  prec1 <- sum(analiza[,5]==1 & analiza[,1]>0.05) / vsiKlasNiPer
  rec1  <- sum(analiza[,5]==1 & analiza[,1]>0.05) / vsiPer1
  prec2 <- sum(analiza[,1] <= 0.05 & analiza[,4] == 24 & analiza[,5] == 24) / vsiKlasJe24Per
  rec2  <- sum(analiza[,1] <= 0.05 & analiza[,4] == 24 & analiza[,5] == 24) / vsiPer24
  
  tp1 <- paste(sum(analiza[,5]==1 & analiza[,1]>0.05), "/" ,vsiKlasNiPer, "=", round(prec1,3))
  tr1 <- paste(sum(analiza[,5]==1 & analiza[,1]>0.05), "/", vsiPer1, "=", round(rec1,3))
  tp2 <- paste(sum(analiza[,1] <= 0.05 & analiza[,4] == 24 & analiza[,5] == 24), "/", vsiKlasJe24Per, "=", round(prec2, 3))
  tr2 <- paste(sum(analiza[,1] <= 0.05 & analiza[,4] == 24 & analiza[,5] == 24), "/", vsiPer24, "=", round(rec2,3))
  
  stat <- data.frame()
  stat <- rbind(stat, c(tp1, tr1, tp2, tr2))
  colnames(stat) <- c("preciznostC1", "priklicC1", "preciznostC2", "priklicC2")
  return(stat)
}

dolociTNTPneznana <- function(analiza, podatki) {
  # doloci stevilo pravih pozitivnih in negativnih ter 
  # laznih pozitivnih in negativnih opravljenih napovedi
  analiza$truePeriod <- podatki$Cosine.Periods
  stetje <- c(0,0,0,0)
  analiza$stat <- 0
  for (i in 1:25) {
    if(analiza$pVal[i] <= 0.05) { #pozitivni
      TP <- analiza$truePeriod[i] == 24 & analiza$period[i] == 24
      FP <- analiza$truePeriod[i] != 24 & analiza$period[i] == 24
    }
    else { # negativni
      TN <- analiza$truePeriod[i] != 24
      FN <- analiza$truePeriod[i] == 24 & analiza$period[i] != 24
    }
    if(TP) {analiza$stat[i] <- "TP"; stetje[1] <- stetje[1] + 1 }
    if(TN) {analiza$stat[i] <- "TN"; stetje[2] <- stetje[2] + 1}
    if(FP) {analiza$stat[i] <- "FP"; stetje[3] <- stetje[3] + 1}
    if(FN) {analiza$stat[i] <- "FN"; stetje[4] <- stetje[4] + 1}
  }
  analiza$phase <- NULL
  analiza$peak.shape <- NULL
  analiza$period <- as.character(analiza$period)
  return(list(analiza, stetje))
}

dolociTNTP <- function(analiza, podatki) {
  # doloci stevilo pravih pozitivnih in negativnih ter 
  # laznih pozitivnih in negativnih opravljenih napovedi
  analiza$truePeriod <- podatki$Cosine.Periods
  stetje <- c(0,0,0,0)
  for (i in 1:25) {
    TP <- analiza$pVal[i] <= 0.05 & analiza$truePeriod[i] == 24
    FP <- analiza$pVal[i] <= 0.05 & analiza$truePeriod[i] != 24
    TN <- analiza$pVal[i] > 0.05 & analiza$truePeriod[i] != 24
    FN <- analiza$pVal[i] > 0.05 & analiza$truePeriod[i] == 24
    if(TP) {analiza$stat[i] <- "TP"; stetje[1] <- stetje[1] + 1 }
    if(TN) {analiza$stat[i] <- "TN"; stetje[2] <- stetje[2] + 1}
    if(FP) {analiza$stat[i] <- "FP"; stetje[3] <- stetje[3] + 1}
    if(FN) {analiza$stat[i] <- "FN"; stetje[4] <- stetje[4] + 1}
  }
  analiza$phase <- NULL
  analiza$peak.shape <- NULL
  analiza$period <- as.character(analiza$period)
  return(list(analiza, stetje))
}

izvoziLatexPrecPriklic <- function(stat, ime) {
  # ustvari in izvozi tabelo preciznosti in priklica v latex formatu
  precrec <- stat[1:2]
  toPrint <- xtable(precrec, type = "latex", caption = "Povzetek pravilnosti klasifikacij za Rain.", label = "RainSk", digits = 3, align = "|c|c|c|")
  names(toPrint) <- c("Preciznost neperiodiènih podatkov", "Priklic neperiodiènih podatkov")
  print(toPrint, file = paste(ime, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
  
  precrec <- stat[3:4]
  toPrint <- xtable(precrec, type = "latex", caption = "Povzetek pravilnosti klasifikacij za Rain.", label = "RainSk", digits = 3, align = "|c|c|c|")
  names(toPrint) <- c("Preciznost periodiènih podatkov", "Priklic periodiènih podatkov")
  print(toPrint, file = paste(ime, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)  
}

izvoziLatex <- function(analiza, podatki, imeZaDatoteko) {
  # pretvori in izvozi analizo z znanimi periodami v latex format
  analiza$phase <- NULL
  analiza$peak.shape <- NULL
  analiza$period <- as.character(analiza$period)
  analiza$truePeriod <- podatki$Cosine.Periods
  toPrint <- xtable(analiza, type = "latex", caption = "Povzetek pravilnosti klasifikacij za RAIN.", label = "rainStab", digits = 3, align = "|c|c|c|c|")
  names(toPrint) <- c("p-vrednost", "predvidena perioda", "prava perioda")
  print(toPrint, file = paste(imeZaDatoteko, "tabela.tex", sep = ""), table.placement = "H")
  
  stats <- t(data.frame(as.character(dolociTNTP(analiza, podatki)[[2]])))
  toPrint2 <- xtable(stats, type = "latex", caption = "Rezultati RAIN.", label = "rainSk", digits = 3, align = "|c|c|c|c|c|")
  names(toPrint2) <- c("TP", "TN", "FP", "FN")
  print(toPrint2, file = paste(imeZaDatoteko, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
}

drawROC <- function(analiza, podatki) {
  # narise graf ROC
  x <- as.numeric( podatki$Cosine.Periods == 1)
  analiza$class <- x
  analiza <- analiza[order(analiza$pVal),]
  roc(data = analiza, response = class, predictor = pVal, plot = TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
      print.auc=TRUE, show.thres=TRUE)
}

#
# MAIN
#

# S1
podatki <- read.csv("Scenarij 1.csv")
vhodRain <- podatki[,2:25]
vhodRain <- t(vhodRain)
# ZNANA
analizaRain <- rain(vhodRain, period = 24, deltat = 2, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S1-RAIN-ZP2")
# NEZNANA
analizaRain <- rain(vhodRain, period = 24, period.delta = 4, deltat = 2, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S1-RAIN-NP2")
statistika <- preciznostPriklicNeznanaPerioda(analizaRain, podatki)
izvoziLatexPrecPriklic(statistika, "S1-RAIN-NP2")


# S2
podatki <- read.csv("Scenarij 2.csv")
vhodRain <- podatki[,2:13]
vhodRain <- t(vhodRain)
# ZNANA
analizaRain <- rain(vhodRain, period = 24, deltat = 4, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S2-RAIN-ZP2")
# NEZNANA
analizaRain <- rain(vhodRain, period = 24, period.delta = 4, deltat = 4, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S2-RAIN-NP2")
statistika <- preciznostPriklicNeznanaPerioda(analizaRain, podatki)
izvoziLatexPrecPriklic(statistika, "S2-RAIN-NP2")


# S3
podatki <- read.csv("Scenarij 3.csv")
vhodRain <- podatki[,2:13]
vhodRain <- t(vhodRain)
# ZNANA
analizaRain <- rain(vhodRain, period = 24, deltat = 2, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S3-RAIN-ZP2")
# NEZNANA
analizaRain <- rain(vhodRain, period = 24, period.delta = 4, deltat = 2, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S3-RAIN-NP2")
statistika <- preciznostPriklicNeznanaPerioda(analizaRain, podatki)
izvoziLatexPrecPriklic(statistika, "S3-RAIN-NP2")


# S4
podatki <- read.csv("Scenarij 4.csv")
vhodRain <- podatki[,2:7]
vhodRain <- t(vhodRain)
# ZNANA
analizaRain <- rain(vhodRain, period = 24, deltat = 4, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S4-RAIN-ZP2")
# NEZNANA
analizaRain <- rain(vhodRain, period = 24, period.delta = 4, deltat = 4, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S4-RAIN-NP2")
statistika <- preciznostPriklicNeznanaPerioda(analizaRain, podatki)
izvoziLatexPrecPriklic(statistika, "S4-RAIN-NP2")


# S5
podatki <- read.csv("Scenarij 5.csv")
vhodRain <- podatki[,2:7]
vhodRain <- t(vhodRain)
# ZNANA
analizaRain <- rain(vhodRain, period = 24, deltat = 8, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S5-RAIN-ZP2")
# NEZNANA
analizaRain <- rain(vhodRain, period = 24, period.delta = 4, deltat = 8, verbose = TRUE, adjp.method = "BH")
drawROC(analizaRain, podatki)
izvoziLatex(analizaRain, podatki, "S5-RAIN-NP2")
statistika <- preciznostPriklicNeznanaPerioda(analizaRain, podatki)
izvoziLatexPrecPriklic(statistika, "S5-RAIN-NP2")
