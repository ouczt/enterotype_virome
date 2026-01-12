
# Load required libraries
library(corrplot)
library(RColorBrewer)

# Function to create heatmap for a given correlation and p-value matrix
create_heatmap <- function(cor_mat, pval_mat, title) {
  # Create a significance mask (e.g., p-value < 0.05)
  sig_mask <- pval_mat < 0.05
  
  # Apply mask to correlation matrix
  cor_mat[!sig_mask] <- 0
  
  # Define color scheme: red for positive, blue for negative
  col <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Create correlation plot
  corrplot(cor_mat, method = "circle", asp = 1, type = "full", 
           col = col, tl.col = "black", tl.srt = 30,
           diag = FALSE,  
           tl.cex = 0.7, cl.cex = 0.7, 
           addgrid.col = "gray", 
           order = "original")
  
  # Add title
  title(main = title, cex.main = 1)
}

# Set up the plotting area for multiple heatmaps
#par(mfrow = c(1, 3))  # Adjust based on the number of conditions (HC, Remission, Flare-up)
# Reset plotting area
par(mfrow = c(1, 1))

# Add legend (manually, as corrplot doesn't handle this automatically)
legend("bottomright", legend = c("Correlation coefficient", "0.8", "0.4", "0", "-0.4", "-0.8"),
       fill = c("white", "red", "pink", "white", "lightblue", "blue"), border = "black", cex = 0.7)
legend("bottomright", inset = c(0, -0.2), legend = c("p-value (-log10)", "3.0", "1.5", "0"),
       fill = c("white", "darkred", "red", "white"), border = "black", cex = 0.7)

# Create heatmaps for each condition
create_heatmap(E1_bv_cor_subset_numeric, E1_bv_pval_subset_numeric, "E1")
create_heatmap(E2_bv_cor_subset_numeric, E2_bv_pval_subset_numeric, "E2")

