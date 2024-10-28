# Define the labels for the clusters in the desired order
cluster_labels <- c("Hot & dispersed", "Hot & concentrated", "Warm & concentrated", "Warm & dispersed", "Mild", "Cool")
cluster_order <- c(5, 3, 2, 4, 1, 6)  # Map to the correct cluster numbers
cluster_counts <- c(21, 82, 103, 30, 19, 17)  # Corresponding counts for each cluster

# Define colors for the clusters
#col <- c("firebrick2", "orange1", "mediumpurple", "seagreen4", "sienna", "skyblue3")
col <- c(brewer.pal(8, "RdYlBu")[1:4], brewer.pal(8, "RdYlBu")[6:7]) # Yellow to red gradient
# Prepare the data for ggplot2
plot_data <- data %>%
  filter(cluster_ward_std_6 %in% cluster_order) %>%
  select(L1ADtemp_pw, cluster = cluster_ward_std_6) %>%
  mutate(cluster = factor(cluster, levels = cluster_order, labels = cluster_labels)) # Adjust the levels to match the desired order

# Create ggplot2 visualization
ggplot(plot_data, aes(x = L1ADtemp_pw, color = cluster)) +
  geom_histogram(aes(y = ..density..), bins = 100, alpha = .1, position = "identity") +  # Keep histograms unfilled
  geom_density(aes(fill = cluster), size = 0.6, alpha = .1) +  # Density plot with fill for legend
  scale_color_manual(values = col, name=NULL) +  # Color for boundaries
  scale_fill_manual(values = col, name=NULL) +  # Color for fill in the legend
  labs(x = "Temperature (°C)", y = "Density") +
    theme_linedraw() +
  theme(
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 13),
    legend.position = c(0.31, 0.84),
    legend.background = element_blank(),
    axis.title.x = element_text(size = 13),  # Adjust x-axis title size
    axis.title.y = element_text(size = 13)   # Adjust y-axis title size
  )
  

# Set up the plotting area for base R (2 rows and 3 columns)
par(mfrow = c(2, 3), mar = c(5, 3, 2, 1.5), las = 0.5, mgp = c(3.6, 1, 0))

# Loop through the clusters based on the defined order for base R plots
for (i in rev(c(seq_along(cluster_order)))) {
  cluster_index <- cluster_order[i]  # Get the correct cluster number
  
  # Plot using the index from pred_reconstruct_pooled
  plot(pred_reconstruct_pooled[[cluster_index]], "overall", lag = 0, cumul = TRUE, ylim = c(0.9, 1.4), 
       ylab = "RR", col = col[i], lwd = 1.5, ci = 'area', lty = 1,
       xlab = paste(cluster_labels[i], "\n(N=", cluster_counts[i], ")", sep = ""),
       ci.arg = list(col = alpha("gray30", 0.1)),
       cex.lab = 1.1,   # Adjust the size of axis labels
       cex.main = 1,  # Adjust the size of the main title if you have one
       cex.axis = 1)
  
  abline(v = c(quan01, quan25, medT, quan75, quan99), lty = 4, col = "gray")
  
  # Add plot labels (a, b, c, d, e, f) on top of each plot
  mtext(letters[8-i], side = 3, line = 0.5, adj = -0.25, cex = 1.2)
}

