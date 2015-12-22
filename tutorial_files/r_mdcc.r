library(plyr)
library(ggplot2)
library(reshape)
corr.rr <- read.table("res_res.txt", header=T)
corr.rr.m <- melt(corr.rr, id = c("res_num1.int", "res_num2.int"))
 scc.rr.mmdcc <- subset(corr.rr.m, variable=="correlation.float" & res_num1.int!=res_num2.int)
 scc.rr.dcc <- subset(corr.rr.m, variable=="corr_dcc.float" & res_num1.int!=res_num2.int)
 scc.rr.mmdcc <- cbind(scc.rr.mmdcc, rep("max(mDCC)", nrow(scc.rr.mmdcc)))
 scc.rr.dcc <- cbind(scc.rr.dcc, rep("DCC", nrow(scc.rr.dcc)))
 colnames(scc.rr.mmdcc)[5] <- "label"
colnames(scc.rr.dcc)[5] <- "label"
draw.map <- function(scc.rr.1, scc.rr.2){
  tmp <- scc.rr.2[,1]
   scc.rr.2[,1] <- scc.rr.2[,2]
  scc.rr.2[,2] <- tmp
   scc.rr.plot <- rbind(scc.rr.1, scc.rr.2)
   ( p <- ggplot(scc.rr.plot, aes(res_num1.int, res_num2.int)) + geom_tile(aes(fill = value))  +   scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-1.0,1.0)) +theme_minimal()  )
 }
 scc.rr.diff <- data.frame(scc.rr.mmdcc$res_num1.int, scc.rr.mmdcc$res_num2.int, rep("diff", nrow(scc.rr.mmdcc)), scc.rr.mmdcc$value-scc.rr.dcc$value, rep("diff", nrow(scc.rr.mmdcc)))
 colnames(scc.rr.diff) <- c("res_num1.int","res_num2.int","variable","value","label")
 summary(scc.rr.diff$value)
 (p <- draw.map(scc.rr.mmdcc, scc.rr.diff))
 ggsave(file = "mdcc_diff_heatmap.png", plot = p, dpi = 600, width = 4.3, height = 3.3) 