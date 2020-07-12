require("cosinor2")
require(xtable)
require(pROC)

#COSINOR
#
# enacba:
#   MESOR je sredina fitted values
#   AKROFAZA je correct.acrophase(analizaCosinor$single.cos[[1]]), PAZI le z akrofazo je pravilno
#   AMPLITUDA je enaka
#   eq = function(x){0.2255538 + 0.8410458*cos((2*pi*x)/24 - 3.106365)}
# EVALVACIJA
# cosinor.detect(analizaCosinorA$single.cos[[1]])
#  -korekcija akrofaze ne vpliva na to


analizaSingleCosinor <- function(analizaCosinor, perioda = 24, podatki) {
  # za vsak cosinor v zbirki pridobi podatke in vrne v obliki tabele
  analiza <- data.frame(row.names = c("MESOR", "Amplituda", "Akrofaza", "p-vrednost", "prava_perioda", "predvidena_perioda"))
  for (i in 1:length(analizaCosinor$single.cos)) {
    sinCos <- analizaCosinor$single.cos[[i]]
    A <- sinCos$coefficients["amp"]
    M <- mean(sinCos$fit$fitted.values)
    PHI <- correct.acrophase(sinCos)
    pVal <- cosinor.detect(sinCos)[4]
    #print(paste(M, " + ", A,"* cos(2*3.14*x/24 + ", PHI, ")", "; P-value = ", pVal))
    #jpeg(paste(i, "rplot.jpg"), width = 1000, height = 700)
    #plot(y = as.numeric(vhodCosinor[i,]), x = seq(0,48,4)[-1], xlab = paste("vzorec: ", i, ", perioda: ", podatki[i,15]))
    #modeq = function(x){M + A*cos((2*pi*x)/24 + PHI)}
    #lines(modeq(1:48), col = "red")
    #dev.off()
    analiza <- cbind(analiza, c(M, A, PHI, pVal, podatki$Cosine.Periods[i], perioda))
    #readline(prompt="Press [enter] to continue")
  }
  return(analiza)
}

dolociPrimernoPeriodo <- function(vhodCosinor, cas) {
  # za casovno serijo izbere najbolj primerno periodo izmed obsega podanih
  analiza <- periodogram(vhodCosinor, time = cas)
  periode <- analiza$data$periodogram
  indeks <- which.max(periode)
  return(cas[indeks])
}

preciznostPriklicNeznanaPerioda <- function(a, podatki) {
  # izracuna preciznost in priklic opravljenih napovedi
  # precision: number of correctly predicted out of all predicted
  # recall: number of correctly predicted out of the number of actual
  stat <- data.frame()
  vsiPer1 <- sum(a[5,]==1)
  vsiPer24 <- sum(a[5,]==24)
  vsiKlasNiPer <- sum(a[4,]>0.05)
  vsiKlasJe24Per <- sum(a[4,] <= 0.05 & a[6,] == 24)
  vsiKlasNi24Per <- sum(a[4,] <= 0.05 & a[6,] != 24)
  
  prec1 <- sum(a[5,]==1 & a[4,]>0.05) / vsiKlasNiPer
  rec1  <- sum(a[5,]==1 & a[4,]>0.05) / vsiPer1
  prec2 <- sum(a[4,] <= 0.05 & a[6,] == 24 & a[5,] == 24) / vsiKlasJe24Per
  rec2  <- sum(a[4,] <= 0.05 & a[6,] == 24 & a[5,] == 24) / vsiPer24
  tp1 <- paste(sum(a[5,]==1 & a[4,]>0.05), "/" ,vsiKlasNiPer, "=", round(prec1,3))
  tr1 <- paste(sum(a[5,]==1 & a[4,]>0.05), "/", vsiPer1, "=", round(rec1,3))
  tp2 <- paste(sum(a[4,] <= 0.05 & a[6,] == 24 & a[5,] == 24), "/", vsiKlasJe24Per, "=", round(prec2, 3))
  tr2 <- paste(sum(a[4,] <= 0.05 & a[6,] == 24 & a[5,] == 24), "/", vsiPer24, "=", round(rec2,3))
  stat <- rbind(stat, c(tp1, tr1, tp2, tr2))
 
  colnames(stat) <- c("preciznostC1", "priklicC1", "preciznostC2", "priklicC2")
  return(stat)
}

dolociTNTPcos <- function(analizaSingleCos, podatki) {
  # doloci stevilo pravih pozitivnih in negativnih ter 
  # laznih pozitivnih in negativnih opravljenih napovedi
  stetje <- data.frame(t(c(0,0,0,0)))
  colnames(stetje) <- c("TP", "TN", "FP", "FN")
  ampRazmerje <- c()
  for (i in 1:25) {
    TP <- analizaSingleCos[4,i] <= 0.05 & podatki$Cosine.Periods[i] == 24
    FP <- analizaSingleCos[4,i] <= 0.05 & podatki$Cosine.Periods[i] != 24
    TN <- analizaSingleCos[4,i] > 0.05 & podatki$Cosine.Periods[i] != 24
    FN <- analizaSingleCos[4,i] > 0.05 & podatki$Cosine.Periods[i] == 24
    if(TP) {stetje[1] <- stetje[1] + 1}
    if(TN) {stetje[2] <- stetje[2] + 1}
    if(FP) {stetje[3] <- stetje[3] + 1}
    if(FN) {stetje[4] <- stetje[4] + 1}
    if(sum(TP, FP, TN, FN) != 1) {print("NOBENA USTREZNA")}
    ampRazmerje[i] <- (analizaSingleCos[2,i] / podatki$Cosine.Amplitude[i])*100
  }
  return(list(ampRazmerje, stetje))
}

analizaNeznanePeriode <- function(vhodCosinor, podatki, cas) {
  # za vsakega izmed vzorcev casovnih serij neznane periode opravi analizo
  analizaNeznane <- data.frame(row.names = c("MESOR", "Amplituda", "Akrofaza", "p-vrednost", "prava_perioda", "predvidena_perioda"))
  ampRazmerje <- c()
  for (i in 1:25) { # 25 vzorcev
    perioda <- dolociPrimernoPeriodo(vhodCosinor[i,], cas)
    analizaCosinor <- population.cosinor.lm(vhodCosinor, time = cas, period = perioda)
    analizaLoceni <- analizaSingleCosinor(analizaCosinor, perioda, podatki)
    analizaNeznane <- cbind(analizaNeznane, analizaLoceni[,i])
    ampRazmerje[i] <- (analizaLoceni[2,i] / podatki$Cosine.Amplitude[i])*100
  }
  analizaNeznane[4,] <- p.adjust(as.numeric(analizaNeznane[4,]), method = "BH")
  statistika <- preciznostPriklicNeznanaPerioda(analizaNeznane, podatki)
  return(list(analizaNeznane, statistika, ampRazmerje))
}

izvoziLatexNeznane <- function(analizaNeznani, podatki, ime) {
  # pretvori in izvozi analizo z neznanimi periodami v latex format
  tabelaParametrov <- cbind(t(analizaNeznani[[1]]), analizaNeznani[[3]])
  rownames(tabelaParametrov) <- NULL
  toPrint <- xtable(tabelaParametrov, type = "latex", caption = "Rezultati cosinor.", label = "CosStab", digits = 3, align = "|C{0.4cm}|C{1.4cm}|C{1cm}|C{1.2cm}|C{1.5cm}|C{1.3cm}|C{1.9cm}|C{1.7cm}|")
  names(toPrint) <- c("MESOR", "A", "$Phi$","p-vrednost", "prava perioda", "predvidena perioda", "razmerje amplitud [%]")
  print(toPrint, file = paste(ime, "tabela.tex", sep = ""), table.placement = "H")
  
  precrec <- analizaNeznani[[2]]
  toPrint <- xtable(precrec, type = "latex", caption = "Povzetek pravilnosti klasifikacij za cosinor.", label = "CosSk", digits = 3, align = "|c|c|c|c|c|")
  names(toPrint) <- c("Preciznost neperiodiènih podatkov", "Priklic neperiodiènih podatkov", "Preciznost periodiènih podatkov", "Priklic periodiènih podatkov")
  print(toPrint, file = paste(ime, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
}

izvoziLatexZnane <- function(analizaZnani, statistika, podatki, ime) {
  # pretvori in izvozi analizo z znanimi periodami v latex format
  tabelaParametrov <- cbind(t(analizaZnani), statistika[[1]])
  rownames(tabelaParametrov) <- NULL
  toPrint <- xtable(tabelaParametrov, type = "latex", caption = "Rezultati cosinor.", label = "CosStab", digits = 3, align = "|C{0.4cm}|C{1.4cm}|C{1cm}|C{1.2cm}|C{1.5cm}|C{1.3cm}|C{1.9cm}|C{1.7cm}|")
  names(toPrint) <- c("MESOR", "A", "$Phi$","p-vrednost", "prava perioda", "predvidena perioda", "razmerje amplitud [%]")
  print(toPrint, file = paste(ime, "tabela.tex", sep = ""), table.placement = "H")
  
  stats <- statistika[[2]]
  toPrint2 <- xtable(stats, type = "latex", caption = "Povzetek pravilnosti klasifikacij za cosinor.", label = "cosSk", digits = 0, align = "|c|c|c|c|c|")
  names(toPrint2) <- c("TP", "TN", "FP", "FN")
  print(toPrint2, file = paste(ime, "tabela.tex", sep = ""), table.placement = "H", include.rownames=FALSE, append = TRUE)
}

drawROC <- function(analiza, podatki) {
  # narise graf ROC
  analiza <- t(analiza)
  analiza <- cbind(analiza, as.numeric( podatki$Cosine.Periods == 1))
  analiza <- analiza[order(analiza[,4]),]
  analiza <- data.frame(analiza)
  row.names(analiza) <- NULL
  i <- colnames(analiza)
  i[7] = "class"
  colnames(analiza) <- i
  print(analiza)
  roc(data = analiza, response = class, predictor = p.vrednost, plot = TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
      print.auc=TRUE, show.thres=TRUE)
}

#
# MAIN
#

# S1
podatki <- read.csv("Scenarij1/Scenarij 1.csv")
vhodCosinor <- podatki[,2:25]
# ZNANA
analizaCosinor <- population.cosinor.lm(vhodCosinor, time = seq(0,48,2)[-1], period = 24)
analizaLoceni <- analizaSingleCosinor(analizaCosinor, 24, podatki)
analizaLoceni[4,] <- p.adjust(as.numeric(analizaLoceni[4,]), method = "BH")
statistika <- dolociTNTPcos(analizaLoceni, podatki)
izvoziLatexZnane(analizaLoceni, statistika, podatki, "S1-Cos-ZP")
drawROC(analizaLoceni, podatki)
# NEZNANA
analizaNeznani <- analizaNeznanePeriode(vhodCosinor, podatki, seq(0,48,2)[-1])
izvoziLatexNeznane(analizaNeznani, podatki, "S1-Cos-NP")

# S2
podatki <- read.csv("Scenarij2/Scenarij 2.csv")
vhodCosinor <- podatki[,2:13]
# ZNANA
analizaCosinor <- population.cosinor.lm(vhodCosinor, time = seq(0,48,4)[-1], period = 24)
analizaLoceni <- analizaSingleCosinor(analizaCosinor, 24, podatki)
analizaLoceni[4,] <- p.adjust(as.numeric(analizaLoceni[4,]), method = "BH")
statistika <- dolociTNTPcos(analizaLoceni, podatki)
izvoziLatexZnane(analizaLoceni, statistika, podatki, "S2-Cos-ZP")
drawROC(analizaLoceni, podatki)
# NEZNANA
analizaNeznani <- analizaNeznanePeriode(vhodCosinor, podatki, seq(0,48,4)[-1])
izvoziLatexNeznane(analizaNeznani, podatki, "S2-Cos-NP")

# S3
podatki <- read.csv("Scenarij3/Scenarij 3.csv")
vhodCosinor <- podatki[,2:13]
# ZNANA
analizaCosinor <- population.cosinor.lm(vhodCosinor, time = seq(0,24,2)[-1], period = 24)
analizaLoceni <- analizaSingleCosinor(analizaCosinor, 24, podatki)
analizaLoceni[4,] <- p.adjust(as.numeric(analizaLoceni[4,]), method = "BH")
statistika <- dolociTNTPcos(analizaLoceni, podatki)
izvoziLatexZnane(analizaLoceni, statistika, podatki, "S3-Cos-ZP")
drawROC(analizaLoceni, podatki)
# NEZNANA
analizaNeznani <- analizaNeznanePeriode(vhodCosinor, podatki, seq(0,24,2)[-1])
izvoziLatexNeznane(analizaNeznani, podatki, "S3-Cos-NP")

# S4
podatki <- read.csv("Scenarij4/Scenarij 4.csv")
vhodCosinor <- podatki[,2:7]
# ZNANA
analizaCosinor <- population.cosinor.lm(vhodCosinor, time = seq(0,24,4)[-1], period = 24)
analizaLoceni <- analizaSingleCosinor(analizaCosinor, 24, podatki)
analizaLoceni[4,] <- p.adjust(as.numeric(analizaLoceni[4,]), method = "BH")
statistika <- dolociTNTPcos(analizaLoceni, podatki)
izvoziLatexZnane(analizaLoceni, statistika, podatki, "S4-Cos-ZP")
drawROC(analizaLoceni, podatki)
# NEZNANA
analizaNeznani <- analizaNeznanePeriode(vhodCosinor, podatki, seq(0,24,4)[-1])
izvoziLatexNeznane(analizaNeznani, podatki, "S4-Cos-NP")

# S5
podatki <- read.csv("Scenarij5/Scenarij 5.csv")
vhodCosinor <- podatki[,2:7]
# ZNANA
analizaCosinor <- population.cosinor.lm(vhodCosinor, time = seq(0,48,8)[-1], period = 24)
analizaLoceni <- analizaSingleCosinor(analizaCosinor, 24, podatki)
analizaLoceni[4,] <- p.adjust(as.numeric(analizaLoceni[4,]), method = "BH")
statistika <- dolociTNTPcos(analizaLoceni, podatki)
izvoziLatexZnane(analizaLoceni, statistika, podatki, "S5-Cos-ZP")
drawROC(analizaLoceni, podatki)
# NEZNANA
analizaNeznani <- analizaNeznanePeriode(vhodCosinor, podatki, seq(0,48,8)[-1])
izvoziLatexNeznane(analizaNeznani, podatki, "S5-Cos-NP")
