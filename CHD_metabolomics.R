library(ggplot2)
library(ggpubr)
library(dplyr)
data <- read.csv("path/to/matched/data.csv")
#pre-process---
data_clean <- data
data_clean <- data_clean[which(data_clean$children_sex!="0"),]
data_clean <- data_clean[-which(data_clean$children_weight<(mean(data_clean$children_weight)-3*sd(data_clean$children_weight)) |
                                 data_clean$children_weight>(mean(data_clean$children_weight)+3*sd(data_clean$children_weight))),]
#umap---
data_clean$dia_result <- as.factor(data_clean$dia_result)
iris_sumap <- uwot::umap(data_clean[,c(4,7:52)], n_neighbors = 15, min_dist = 0.01,
                         y = data_clean$dia_result, target_weight = 0.5)
head(iris_sumap)
iris_sumap_res <- data.frame(iris_sumap,Group = data_clean$dia_result)
head(iris_sumap_res)
iris_sumap_res$Group <- gsub("0","Normal",iris_sumap_res$Group)
iris_sumap_res$Group <- gsub("1","non-significant CHD",iris_sumap_res$Group)
iris_sumap_res$Group <- gsub("2","CHD",iris_sumap_res$Group)
ggplot(iris_sumap_res,aes(X1,X2,color=Group)) + 
  stat_ellipse(level = 0.95, show.legend = T, linetype = 2)+
  geom_point(size = 0.5) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()) + 
  labs(x="UMAP_1",y="UMAP_2")
#Significantly Altered Metabolites---
df.metab <- data_clean
df.metab$dia_result <- as.character(df.metab$dia_result)
df.metab$dia_result[which(df.metab$dia_result==0)] <- "Normal"
df.metab$dia_result[which(df.metab$dia_result==1)] <- "non-significant CHD"
df.metab$dia_result[which(df.metab$dia_result==2)] <- "CHD"
type1 <- "CHD"
type2 <- "Normal"
df.test <- rbind(df.metab[which(df.metab$dia_result == type1),], 
                 df.metab[which(df.metab$dia_result == type2),])
matab_num <- ncol(df.test)
test_results <- data.frame(Metabolites = colnames(df.test)[7:(matab_num-1)])
test_results$wilcox  <- apply(df.test[,7:(matab_num-1)], 2,
                             function(x) unlist(wilcox.test(as.numeric(x) ~ df.test$dia_result, data = df.test)[3]))
test_results$wilcox_BH <- p.adjust(test_results$wilcox, method = "BH")
cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.test[,7:(matab_num-1)], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.test$dia_result == type1)]))/
                            mean(as.numeric(x[which(df.test$dia_result == type2)])))
test_results$LOG2FC <- log2(test_results[,4])
#Volcano plot---
colnames(test_results)[4] <- "FC"
test_results$compare <- "CHD_vs_Normal"
resOrdered <- test_results
resOrdered$change <- ifelse(resOrdered$wilcox_BH < 0.05 & abs(resOrdered$LOG2FC) >= 0.03,
                            ifelse(resOrdered$LOG2FC > 0.03 , "Up", "Down"),"NS")
resOrdered$logp <- -log10(resOrdered$wilcox_BH)
resOrdered$symbol <- resOrdered$Metabolites
ggscatter(resOrdered,
          x = 'LOG2FC',
          y = 'logp',
          xlab = "logFC",
          ylab = "-log10(P.value)",
          size = 4,
          shape = 21,
          repel = T,
          color = "change",
          fill = "change",
          label = "symbol",
          label.select = resOrdered$symbol[resOrdered$change!="NS"],
) +
  scale_color_manual(values = c("Up" = "#FC4E07",
                                "NS" = "#999999",
                                "Down" = "#00AFBB")) +
  scale_fill_manual(values = c("Up" = "#FC4E0770",
                               "NS" = "#99999970",
                               "Down" = "#00AFBB70")) +
  facet_grid(.~compare) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
#heatmap---
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
data_clean_scale <- data_clean
data_clean_scale[,7:52] <- scale(data_clean_scale[,7:52])
data_clean_scale$dia_result[which(data_clean_scale$dia_result==0)] <- "Normal"
data_clean_scale$dia_result[which(data_clean_scale$dia_result==1)] <- "non-significant CHD"
data_clean_scale$dia_result[which(data_clean_scale$dia_result==2)] <- "CHD"
mat <- t(data_clean_scale[,7:52])
colnames(mat) <- data_clean_scale$id
cluster_info <- data_clean_scale$dia_result
top_anno <- HeatmapAnnotation(
  df = data.frame(Group = cluster_info),
  col = list(Group = c("Normal" = "#66C2A5",
                       "non-significant CHD" = "#9fafd3",
                       "CHD" = "#FC8D62")))
col_fun = colorRamp2(c(-0.5, 0, 1, 2), c("cornflowerblue","white", "yellow", "red"))
Heatmap(as.matrix(mat),
        col = col_fun,
        cluster_rows = T,
        cluster_columns = T,
        show_column_names = FALSE,
        show_row_names = T,
        column_gap = unit(0.5,"mm"),
        heatmap_legend_param = list(title=NULL),
        column_split = cluster_info,
        top_annotation = top_anno,
        column_title = NULL)
#cor_heatmap---
pheatmap::pheatmap(
  cor(data_clean[,c(2,3,6:53)],data_clean[,c(2,3,6:53)]),
  cluster_cols = T,cluster_rows = T,
  border_color = NA,
  cellwidth = 10,cellheight = 10,
  treeheight_row = 8,treeheight_col = 8
)
#Correlation with diameters---
library(tidyverse)
data_VSD <- read.csv("path/to/VSD/data.csv")
data_VSD$VSD_pos <- as.factor(data_VSD$VSD_pos)
iris_sumap <- uwot::umap(data_VSD[,c(6:53,55)], n_neighbors = 15, min_dist = 0.01,
                         y = data_VSD$VSD_pos, target_weight = 0.5)
head(iris_sumap)
iris_sumap_res <- data.frame(iris_sumap,Group = data_VSD$VSD_pos)
head(iris_sumap_res)
iris_sumap_res$d <- data_VSD$VSD_d
iris_sumap_res$Group <- factor(iris_sumap_res$Group,levels = c("Perimembranous VSD","Muscular VSD","Subarterial VSD","Others"))
ggplot(iris_sumap_res,aes(X1,X2,color=Group,fill=Group)) + 
  scale_color_manual(values = brewer.pal(4,"Set2"))+
  stat_ellipse(level = 0.95, show.legend = T, linetype = 2)+
  geom_point(size = 2) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()) + 
  labs(x="UMAP_1",y="UMAP_2")
#Significantly Altered Metabolites between Perimembranous and Muscular VSD---
df.metab <- data_VSD
df.metab <- df.metab %>% filter(VSD_pos %in% c("Perimembranous VSD","Muscular VSD"))
type1 <- "Perimembranous VSD"
type2 <- "Muscular VSD"
df.test <- rbind(df.metab[which(df.metab$VSD_pos == type1),], 
                 df.metab[which(df.metab$VSD_pos == type2),])
matab_num <- ncol(df.test)
test_results <- data.frame(Metabolites = colnames(df.test)[7:(matab_num-4)])
test_results$wilcox  <- apply(df.test[,7:(matab_num-4)], 2,
                              function(x) unlist(wilcox.test(as.numeric(x) ~ df.test$VSD_pos, data = df.test)[3]))
test_results$wilcox_BH <- p.adjust(test_results$wilcox, method = "BH")
cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.test[,7:(matab_num-4)], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.test$VSD_pos == type1)]))/
                            mean(as.numeric(x[which(df.test$VSD_pos == type2)])))
test_results$LOG2FC <- log2(test_results[,4])

#regression---
library(ggpmisc)
for (i in colnames(data_VSD)[c(2,3,6:53)]) {
  p <- ggplot(data_VSD,aes(x = VSD_d,y = eval(as.symbol(i)))) + 
    geom_point(shape=21,color="#6F6F6F",fill="#c3a27680",size=2) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          panel.grid = element_blank()) +
    stat_smooth(aes(x = VSD_d,y = eval(as.symbol(i))),formula = y~x,method = "lm",colour = "#a9376a",se = T,fill = "#ebc2d4") +
    stat_poly_eq(aes(label = paste(..rr.label..,..p.value.label..,sep = "~~~")), 
                 parse = TRUE) +
    xlab("Diameter of VSD defect (mm)")+
    ylab(i)
  p
}