library(ggplot2)
library(gridExtra) 
library(patchwork) 

PathResults <- "D:/simulateReads/GeneralResults/"
Folders <- c("1M", "3M", "5M", "10M", "25M", "60M")
data <- data.frame()
fdr_data <- data.frame()
for (folder in Folders) {
  pathSample <- paste0(PathResults,folder)
  for (file in list.files(path = pathSample, pattern = "*.RData",full.names = T)) {
    load(file)
  }
  TP <- colSums(!is.na(table_pvalues_positives_SE) & table_pvalues_positives_SE<0.05) +
    colSums(!is.na(table_pvalues_positives_ALT) & table_pvalues_positives_ALT<0.05)
  FN <- colSums(!is.na(table_pvalues_positives_SE) & table_pvalues_positives_SE>0.05) +
    colSums(!is.na(table_pvalues_positives_ALT) & table_pvalues_positives_ALT>0.05) +
    colSums(is.na(table_pvalues_positives_ALT)) + colSums(is.na(table_pvalues_positives_SE))
  
  TN <- colSums(!is.na(table_pvalues_negatives_SE) & table_pvalues_negatives_SE>0.05) +
    colSums(!is.na(table_pvalues_negatives_ALT) & table_pvalues_negatives_ALT>0.05) +
    colSums(is.na(table_pvalues_negatives_ALT)) + colSums(is.na(table_pvalues_negatives_SE))
  FP <- colSums(!is.na(table_pvalues_negatives_SE) & table_pvalues_negatives_SE<0.05) +
    colSums(!is.na(table_pvalues_negatives_ALT) & table_pvalues_negatives_ALT<0.05)
  
  resSample <- data.frame(t(rbind(TP,FP)))
  resSample_TP <- cbind(rep(folder,5),
                     c("rMATS", "EPST", "EPBAM", "MAJIQ", "SUPPA2"),
                     c("Non-annotation", "Annotation-based", "Non-annotation", "Non-annotation", "Annotation-based"),
                     rep("True Positives",5),
                     TP)
  resSample_FP <- cbind(rep(folder,5),
                     c("rMATS", "EPST", "EPBAM", "MAJIQ", "SUPPA2"),
                     c("Non-annotation", "Annotation-based", "Non-annotation", "Non-annotation", "Annotation-based"),
                     rep("False Positives",5),
                     FP)
  res <- rbind(resSample_TP,resSample_FP)
  data <- rbind(res,data)
  rownames(resSample) <- NULL
  fdr_data <- rbind(fdr_data,c(folder, as.numeric(round(FP/(FP+TP),4))))
  
}


adj_TP <- table_pvalues_positives_SE
adj_TP[table_pvalues_positives_SE<0.05] <- 1
adj_TP[table_pvalues_positives_SE>0.05] <- 0
adj_TP[is.na(adj_TP)] <- 0

plot(euler(adj_TP, shape = "ellipse",),fill=c("#C9B6E4", "#F4B6C2", "#f7caac", "#AFDDE3", "#6BAED6"))
UpSetR::upset(adj_TP)

colnames(data) <- c("Depth","Tool","Type",
                    "Metric",
                    "Count")
colnames(fdr_data) <- c("Depth","rMATS", "EPST", "EPBAM", "MAJIQ", "SUPPA2")
data$Depth <- factor(data$Depth, levels = c("1M", "3M", "5M", "10M", "25M", "60M"))
data$Count<-as.numeric(data$Count)



p <- ggplot(data, aes(x = Depth, y = Count, group = interaction(Tool, Metric), color = Tool, fill = Tool)) +
  geom_line(aes(linetype = Metric), size = 1) + 
  geom_point(aes(shape = Type, fill=Tool), size = 4, stroke = 0.8, color = "black") + 
  scale_shape_manual(values = c("Annotation-based" = 24,  "Non-annotation" = 21)) +
  scale_fill_manual(
    values = c("EPST" = "#6BAED6", "SUPPA2" = "#C9B6E4", "EPBAM" = "#AFDDE3", 
               "rMATS" = "#F4B6C2", "MAJIQ" = "#f7caac")
  ) +
  scale_color_manual(
    values = c("EPST" = "#6BAED6", "SUPPA2" = "#C9B6E4", "EPBAM" = "#AFDDE3", 
               "rMATS" = "#F4B6C2", "MAJIQ" = "#f7caac")
  ) +
  scale_linetype_manual(values = c("True Positives" = "solid", "False Positives" = "dashed")) +
  labs(
    title = "TP & FP by Tool for Simulated Data",
    x = "Sequencing Depth",
    y = "Number of Events Detected",
    color = "Tool",
    shape = "Type",
    fill = "Tool"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    # legend.position = "none",
    legend.box = "vertical",
    axis.line = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_line(size = 0.5, color = "white"),
    panel.grid.minor = element_blank(),
    legend.key = element_blank() # Opcional, elimina el fondo de las leyendas para una visualización más limpia
  )

p


ggsave(
  filename = "tp_fp_figure.png",
  plot = p, 
  width = 9,   
  height = 5,   
  dpi = 600 
)





p <- ggplot(data, aes(x = Depth, y = Count, group = interaction(Tool, Metric), color = Tool, fill = Tool)) +
  geom_line(aes(linetype = Metric), size = 1) + 
  geom_point(aes(shape = Type, fill = Tool), size = 4, stroke = 0.7) + # Deja que el borde esté controlado por el mapeo aes()
  scale_shape_manual(values = c("Annotation-based" = 24,  "Non-annotation" = 21)) +
  scale_fill_manual(
    values = c("EPST" = "#6BAED6", "SUPPA2" = "#C9B6E4", "EPBAM" = "#AFDDE3", 
               "rMATS" = "#F4B6C2", "MAJIQ" = "#f7caac")
  ) +
  scale_color_manual(
    values = c("EPST" = "#6BAED6", "SUPPA2" = "#C9B6E4", "EPBAM" = "#AFDDE3", 
               "rMATS" = "#F4B6C2", "MAJIQ" = "#f7caac")
  ) +
  scale_linetype_manual(values = c("True Positives" = "solid", "False Positives" = "dashed")) +
  labs(
    title = "TP & FP by Tool for Simulated Data",
    x = "Sequencing Depth",
    y = "Number of Events Detected",
    color = "Tool",
    shape = "Type",
    fill = "Tool"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    # legend.position = "none",
    legend.box = "vertical",
    axis.line = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_line(size = 0.5, color = "white"),
    panel.grid.minor = element_blank(),
    legend.key = element_blank() # Opcional, elimina el fondo de las leyendas para una visualización más limpia
  )

p
