#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

ref     <- args[1]
cohort  <- args[2]
mode    <- args[3]
eigenvec<- args[4]
eigenval<- args[5]
pop     <- args[6]

# Load data
eigen_vec <- readr::read_delim(
    eigenvec,
    col_names = c('FID', 'IID', paste0('PC', 1:20))
)

populations <- read.table(pop, header = FALSE)
names(populations) <- c('FID', 'IID', 'pop', 'superpop')
      
pcs <- dplyr::left_join(eigen_vec, populations)

# Plot
pca_plot <- (
    ggplot2::ggplot(pcs) +
    ggplot2::geom_point(
        data = dplyr::filter(pcs, !is.na(superpop)),
        ggplot2::aes(x = PC1, y = PC2, color = superpop),
        size = .5, alpha = .5
    ) +
    ggplot2::geom_point(
        data = dplyr::filter(pcs, is.na(superpop)),
        ggplot2::aes(x = PC1, y = PC2),
        shape = 3, alpha = .5
    ) +
    ggplot2::labs(color = '') +
    ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.ticks.length = ggplot2::unit(.15, "cm"),
        axis.ticks = ggplot2::element_line(color = 'lightgrey'),
        axis.title.y = ggplot2::element_text(vjust = 2)
    )
)

# Save plot
ggplot2::ggsave(
    plot = pca_plot,
    filename = paste(ref, cohort, mode, 'png', sep = '.'),
    width = 4.5, height = 3.5, units = 'in', dpi = 300
)
