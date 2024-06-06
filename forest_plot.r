library(forestplot)
library(Ipaper)
library(readxl)
windowsFonts(
  # Ӣ??????
  arial = windowsFont(family = "Arial"),             # Arial????
  Helvetica = windowsFont(family = "Helvetica")      # ӡˢ??
)
#### T2D to VTE ####

lmr_forest <- read_excel('T2D_VTE_results.xlsx', sheet = 'T2D_to_VTE',
                         col_names=F)
lmr_forest$...7[4:7] <- format(round(as.numeric(lmr_forest$...7[4:7]), 3), nsmall=3)
lmr_forest$...7[9:12] <- format(round(as.numeric(lmr_forest$...7[9:12]), 3), nsmall=3)


lmr_forest <- as.data.frame(lmr_forest)
View(lmr_forest)
### ????????ʱ??һ??Ҫ??header????ΪFalse????֤??һ?в?????????
p <- forestplot(labeltext=as.matrix(lmr_forest[,c(1:7)]),
           mean=cbind(lmr_forest$...8),
           lower=cbind(lmr_forest$...9),
           upper=cbind(lmr_forest$...10),
           is.summary = c(T,F,T,rep(F,4),T,rep(F,4)),
           fn.ci_norm = fpDrawCircleCI,
           clip = c(0.7, 1.75),
           zero=1,
           boxsize=0.3,
           lineheight = unit(9, 'mm'),
           colgap= unit(6,'mm'),
           lwd.zero = 2,
           lwd.ci=1.5,
           xticks = c(0.7, 1, 1.75),
           col=fpColors(zero='grey70', lines='black', 
                        box=c("#000000")),
           xlab='OR [95% CI] for VTE',
           txt_gp = fpTxtGp(ticks = gpar(cex=0.8, fontfamily="arial"),
                            xlab  = gpar(cex = 1, fontfamily="arial"),
                            summary = gpar(cex= 1.2, fontfamily="arial"),
                            legend = gpar(cex=1.0,fontfamily="arial")),
           ci.vertices = T,
           ci.vertices.height = .2,
           graph.pos=6,
           graphwidth = unit(60, 'mm'),
           hrzl_lines = list("2"=gpar(col="black", lwd=2))
           # legend=c('UVMR', 'Model1 MVMR', 'Model2 MVMR', 'Model3 MVMR'),
           # legend_args = fpLegend(pos=list(x=.5, y=0.85, "align" = "horizontal"),r = unit(.1, "snpc"))
           )
write_fig(
  p,
  file = "T2D_to_VTE.jpg",
  width = 12,
  height = 6,
  devices = 'jpg',
  res = 600,
  show = T
)

lmr_forest <- read_excel('T2D_VTE_results.xlsx', sheet = 'VTE_to_T2D',
                         col_names=F)
lmr_forest$...7[4] <- format(round(as.numeric(lmr_forest$...7[4]), 3), nsmall=3)
lmr_forest$...7[6:9] <- format(round(as.numeric(lmr_forest$...7[6:9]), 3), nsmall=3)


lmr_forest <- as.data.frame(lmr_forest)
View(lmr_forest)

p <- forestplot(labeltext=as.matrix(lmr_forest[,c(1:7)]),
                mean=cbind(lmr_forest$...8),
                lower=cbind(lmr_forest$...9),
                upper=cbind(lmr_forest$...10),
                is.summary = c(T,F,T,rep(F,1),T,rep(F,4)),
                fn.ci_norm = fpDrawCircleCI,
                clip = c(0.95, 1.05),
                zero=1,
                boxsize=0.3,
                lineheight = unit(9, 'mm'),
                colgap= unit(6,'mm'),
                lwd.zero = 2,
                lwd.ci=1.5,
                xticks = c(0.95, 1, 1.05),
                col=fpColors(zero='grey70', lines='black', 
                             box=c("#000000")),
                xlab='OR [95% CI] for T2D',
                txt_gp = fpTxtGp(ticks = gpar(cex=0.8, fontfamily="arial"),
                                 xlab  = gpar(cex = 1, fontfamily="arial"),
                                 summary = gpar(cex= 1.2, fontfamily="arial"),
                                 legend = gpar(cex=1.0,fontfamily="arial")),
                ci.vertices = T,
                ci.vertices.height = .2,
                graph.pos=6,
                graphwidth = unit(60, 'mm'),
                hrzl_lines = list("2"=gpar(col="black", lwd=2))
                # legend=c('UVMR', 'Model1 MVMR', 'Model2 MVMR', 'Model3 MVMR'),
                # legend_args = fpLegend(pos=list(x=.5, y=0.85, "align" = "horizontal"),r = unit(.1, "snpc"))
)
write_fig(
  p,
  file = "VTE_to_T2D.jpg",
  width = 12,
  height = 5,
  devices = 'jpg',
  res = 600,
  show = T
)
