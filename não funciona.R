library(gplots)
library(jpeg)
library(pixmap)
library(wavethresh)
library()
getwd()
setwd(choose.dir())
library(shapeR)
data(FISH)


shape=detect.outline(shape,threshold=0.2,write.outline.w.org=TRUE)

shape=remove.outline(shape,"IC","403_54")

detect.outline(object, threshold=0.2, mouse.click=FALSE,
display.images=FALSE, write.outline.w.org=FALSE)

dir()

shape=shapeR("/home/ivan/Documentos/otólitos paper ivan/","plan_paper.csv")

shape=detect.outline(shape,threshold=0.2,write.outline.w.org=TRUE)

setwd(choose.dir("C:/Desktop/ShapeAnalysis","FISH.csv")) #aqui é só um exemplo, galera
shapeR("FISH.csv")

show.original.with.outline(shape,"IC","403_54")

shape=smoothout(shape,n=109)

shape=generateShapeCoefficients(shape)

shape = enrich.master.list(shape)

tapply(getMeasurements(shape)$otolith.area, getMasterlist(shape)$pop, mean)
plotWaveletShape(shape, "station", show.angle = F, lwd =2,lty = 1)
shape = stdCoefs(shape, classes = "pop", "length_cm", bonferroni = FALSE)
getMeasurements(shape)
dados=getMeasurements(shape)
tapply(getMeasurements(shape)$otolith.area, getMasterlist(shape)$pop, mean)
est.list = estimate.outline.reconstruction(shape)
outline.reconstruction.plot(est.list, max.num.harmonics = 15)
plotWavelet(shape, level = 5, class.name = "pop", useStdcoef = TRUE)

list<-getMasterlist(shape)
write.csv(list, "/home/ivan/Área de Trabalho/shaper tainhas/dados morfomet")

save(shape,file = "test.RData")

##carregar o Vegan => library(vegan)

cap.res = capscale(getStdWavelet(shape) ~ getMasterlist(shape)$pop)
anova(cap.res, by = "terms", step = 1000)


eig = eigenvals(cap.res,constrained = T)
eig.ratio = eig/sum(eig)

cluster.plot(scores(cap.res)$sites[,1:2],getMasterlist(shape)$pop,
xlim = range(scores(cap.res)$sites[,1]),
ylim = range(scores(cap.res)$sites[,2]),
xlab = paste("CAP1 (",round(eig.ratio[1]*100,1),"%)",sep = ""),
ylab = paste("CAP2 (",round(eig.ratio[2]*100,1),"%)",sep = ""), 
plotCI = TRUE,conf.level = 0.95,las = 1)

pop = factor(getMasterlist(shape)$pop)

library(ipred)
mypredict.lda <- function(object, newdata)
predict(object, newdata = newdata)$class
stdw = getStdWavelet(shape)
pop = factor(getMasterlist(shape)$pop)
dd = data.frame(stdw = stdw,pop = pop)
errorest(pop ~., data = dd, model = lda, estimator = "cv", predict = mypredict.lda,est.para = control.errorest(nboot = 1000))