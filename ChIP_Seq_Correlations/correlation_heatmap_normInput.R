#### Plot heatmap ####

library(tidyverse)
library(reshape2)

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/ChIP_Seq_Correlations/'
setwd(wd)

## Norm by input corrs
corr_df <- read_csv('./norm_by_input_corrs.csv') %>%
  rename(Corr_With = `...1`)

x <- corr_df %>%
    select(-Corr_With)

mean(x[lower.tri(x, diag = FALSE)])
max(x[lower.tri(x, diag = FALSE)])
min(x[lower.tri(x, diag = FALSE)])

## Plot correlations
m_corr_df <- pivot_longer(corr_df, cols = !Corr_With)

p <- ggplot(m_corr_df, aes(x=Corr_With, y=name, fill=value))
p <- p + geom_tile()
p


reorder_cormat <- function(cormat){
  ## Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


corr_mtx <- corr_df %>%
  select(-Corr_With) %>%
  as.matrix()

rownames(corr_mtx) <- corr_df$Corr_With

## Reorder the correlation matrix
sort_corr_mtx <- reorder_cormat(corr_mtx)
m_corr_mtx <- melt(sort_corr_mtx)
head(m_corr_mtx)

p <- ggplot(m_corr_mtx, aes(x=Var1, y=Var2, fill=value))
p <- p + geom_tile()
p <- p + theme(
           axis.title.x=element_blank(),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
           axis.title.y=element_blank(),
           legend.position="top"
           )
p <- p + scale_fill_gradient(
           low = "yellow",
           high = "red",
           limits = c(0.7,1),
           name = 'Correlation',
           )
p

ggsave('corr_norm_by_imput.pdf', p, device = 'pdf')
