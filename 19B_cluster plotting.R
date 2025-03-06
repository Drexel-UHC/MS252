# Define the labels for the clusters in the desired order
cluster_labels <- c("Hot & dispersed", "Hot & concentrated", "Warm & concentrated", "Warm & dispersed", "Mild", "Cool")
cluster_order <- c(5, 3, 2, 4, 1, 6)  # Map to the correct cluster numbers
cluster_counts <- c(19, 103, 82, 30, 21, 17)  # Corresponding counts for each cluster

# Define colors for the clusters
#col <- c("firebrick2", "orange1", "mediumpurple", "seagreen4", "sienna", "skyblue3")
col <- c(brewer.pal(8, "RdYlBu")[1:4], brewer.pal(8, "RdYlBu")[6:7]) # Yellow to red gradient
# Prepare the data for ggplot2
plot_data <- data %>%
  filter(cluster_ward_std_6 %in% cluster_order) %>%
  select(L1ADtemp_pw, cluster = cluster_ward_std_6) %>%
  mutate(cluster = factor(cluster, levels = cluster_order, labels = cluster_labels)) # Adjust the levels to match the desired order

################################################################################################################
# png("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Submission-LTPH/fig2a.png", width = 5, height = 6.5, units = "in", res = 300)  # 8x8 inches, 300 DPI

# Create ggplot2 visualization
ggplot(plot_data, aes(x = L1ADtemp_pw, color = cluster)) +
  geom_histogram(aes(y = ..density..), bins = 100, alpha = .1, position = "identity") +  # Keep histograms unfilled
  geom_density(aes(fill = cluster), size = 0.6, alpha = .1) +  # Density plot with fill for legend
  annotate("text", x = -Inf, y = Inf, label = "a", hjust = 3.5, vjust = .5, size = 10) +  # Annotation outside plot
  coord_cartesian(clip = "off") +  # Allow annotation outside plot area
  scale_color_manual(values = col, name=NULL) +  # Color for boundaries
  scale_fill_manual(values = col, name=NULL) +  # Color for fill in the legend
  labs(x = "Temperature (°C)", y = "Density") +
    theme_linedraw() +
  theme(
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 16),
    legend.position = c(0.32, 0.85),
    legend.background = element_blank(),
    axis.title.x = element_text(size = 18),  # Adjust x-axis title size
    axis.title.y = element_text(size = 18),   # Adjust y-axis title size
    axis.text.x = element_text(size = 18),   # Adjust x-axis tick font size
    axis.text.y = element_text(size = 18),   # Adjust y-axis tick font size
    panel.grid.major = element_line(color = alpha("black", 0.1), size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = alpha("black", 0.1), size = 0.25), # Minor grid lines
    panel.background = element_blank(),      # Remove background color
    panel.border = element_rect(color = "black", fill = NA)  # Add a border if needed
  ) 

# dev.off()
################################################################################################################
plot_data[, .(
  mean_value = mean(L1ADtemp_pw, na.rm = TRUE),            # Mean
  p1 = quantile(L1ADtemp_pw, probs = 0.01, na.rm = TRUE),  # 1st percentile
  p99 = quantile(L1ADtemp_pw, probs = 0.99, na.rm = TRUE)  # 99th percentile
), by = cluster]



png("/Users/cheng-kaihsu/Library/Mobile Documents/com~apple~CloudDocs/Berkeley/Fall 2023/SALURBAL/Submission-LTPH/fig2b-g.png", width = 6, height = 3.8, units = "in", res = 300)  # 8x8 inches, 300 DPI

# Set up the plotting area for base R (2 rows and 3 columns)
par(mfrow = c(2, 3), mar = c(5, 5, 2, 1), las = 1.5, mgp = c(3.8, 1.2, 0))

# Loop through the clusters based on the defined order for base R plots
for (i in rev(c(seq_along(cluster_order)))) {
  cluster_index <- cluster_order[i]  # Get the correct cluster number
  
  # Plot using the index from pred_reconstruct_pooled
  plot(pred_reconstruct_pooled[[cluster_index]], "overall", lag = 0, cumul = TRUE, ylim = c(0.9, 1.4), 
       ylab = "", col = col[i], lwd = 2, ci = 'area', lty = 1,
       xlab = paste("\n\nTemperature (%tile)\n",cluster_labels[i], " (N=", cluster_counts[i], ")", sep = ""),
       ci.arg = list(col = alpha("gray30", 0.1)),
       cex.lab = 1.,   # Adjust the size of axis labels
       cex.main = 1.2,  # Adjust the size of the main title if you have one
       cex.axis = 1.1,
       las=1)
  
  abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")
  mtext("RR", side = 2, line = 2.8, cex = .9, las=0)  # Adjust `line` for distance
  
  # Add plot labels (a, b, c, d, e, f) on top of each plot
  mtext(letters[8-i], side = 3.2, line = 0.5, adj = -0.45, cex = 1.3)
}

dev.off()
# 
# pred_reconstruct_pooled[[6]]$allRRfit[991]
# pred_reconstruct_pooled[[6]]$allRRlow[991]
# pred_reconstruct_pooled[[6]]$allRRhigh[991]

# 1) mild N=21; 2) warm/con, N=82; 3) hot/con N=103; 4) warm/dis N=30; 5) hot/dis N=19; 6) cool N=17



