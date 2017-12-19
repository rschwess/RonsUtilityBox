# Common Plot functions

# Themes ======================================================
require(ggplot2)
require(RColorBrewer)

# Standard Themes
science_theme <- theme(
  panel.grid.major = element_line(size = 0.5, color = "grey"),
  panel.grid.minor = element_blank(),
  plot.title = element_text(hjust = 0.5),
  text = element_text(size = 14),
  axis.line = element_line(color="black", size = 0.7),
  axis.line.x = element_line(color="black", size = 0.7),
  axis.line.y = element_line(color="black", size = 0.7),
  # plot.margin = unit(c(0.7,0.7,0.7,0.7), "lines"),
  panel.border=element_blank(),
  strip.background = element_blank()
)

# Plot Functions ==============================================

# Pie Chart Plot ----------------------------------------------
ggPieChart <- function(t, count_type="count", num_indent=1.25, label_indent=1, text_size=4){
  # Function to plot a pie chart from a data.frame or tible
  # Arguments:
  #   t: data.frame or tibble (2 columns {label / count})
  #   count_type: select if to plot "count" or "percent"
  #   num_indent: indent value 1 - 2 where to place the numerical value (count or percent)
  #   label_indent: same for label
  #   test_size: text size for coutn and label
  #   scale colours and modify background, theme etc as if it where a ggplot (it is :])
  # Returns:
  #   ggplot2 pie chart plot object
  
  # needs scales
  require(scales)
  
  # rename columns
  colnames(t) <- c("label", "count")
  # make factor data frame
  t$label <- factor(t$label, levels=rev(t$label))
  # make cumsum column
  t <- t %>%
    mutate(cumsum = count/2 + c(0, cumsum(count)[-length(count)]))
  
  # init plot
  p <- ggplot(t, aes(x="", y=count, fill=label)) + 
    geom_bar(width = 1, stat = "identity", col="black") + 
    coord_polar("y", start=0)
  # add label and count data
  p <- p + geom_text(data=t, aes(x = label_indent, y = cumsum, label = label), size=text_size)
  if(count_type == "percent"){
    p <- p + geom_text(aes(x = num_indent, y = count/2 + c(0, cumsum(count)[-length(count)]), label = paste0(round((count/sum(t$count))*100, digits=1), "%")), size=text_size)
  }else{
    p <- p + geom_text(aes(x = num_indent, y = count/2 + c(0, cumsum(count)[-length(count)]), label = count), size=text_size)
  }
  # add right colour scheme and theme
  p <- p + scale_fill_brewer(palette = "Spectral") + theme_bw() + science_theme + theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
  
  return(p)
  
}