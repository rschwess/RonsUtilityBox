# R functions for plotting Motifs, ICM matrices with ggplot
# Author: Ron Schwessinger
# Date: 13.12.2016

# Plotting ICM matrices using ggplot2 to combine those plots with other ggplots and cowplot etc.
# The functions implement the plotting of the ICM/Motif by having the letters as polygons
# Pretty much what you see if you download an svg from JASPAR
# Usage:
# plotICM(data.frame) to plot the ICM matrix (4 rows A, C, G, T times the number of positions)
#
# Note:
# all alpha like usable so happy about any comments
# *circle function adapted from Joran from stackoverflow
# * don't tell Hadley ;)

# required
require(ggplot2)
require(RColorBrewer)

# THEME -----------------------------------------------------------------------
motif_theme <- theme(
  panel.grid = element_blank(),
  text = element_text(size = 14),
  axis.title.x = element_blank(),
  axis.line = element_line(color="black", size = 0.7),
  axis.line.x = element_line(color="black", size = 0.7),
  axis.line.y = element_line(color="black", size = 0.7),
  plot.margin = unit(c(0.7,0.7,0.7,0.7), "lines"),
  panel.border=element_blank(),
  strip.background = element_blank()
)

# FUNCTIONS -------------------------------------------------------------------
circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
  # circle plotting helper function modified from "joran" (stackoverflow)
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  xx[xx < .5] <- xx[xx < .5] + (0.5 - xx[xx < .5]) * .2 #squash a bit
  xx[xx > .5] <- xx[xx > .5] - (xx[xx > .5] - .5) * .2
  return(data.frame(x = xx, y = yy))
}

plotBitBase <- function(letter, scale=1, xintercept=0, yintercept=0){
  # Function to plot a Base letter in Bit content format using ggplot2 and polygons
  # Input:
  #   letter: letter [A,C,G,T]
  #   scale: REAL scale factor to multiply with the coordinates (defaults to 1). 
  #           Use to scale the letter size to the bit content (will only scale in y direction to keep column width)
  #   xintercept: REAL adjust the x-postion (center) of the base letter (default=0) use to place the letter in your plot
  #   yintercept: REAL adjust the y-postion (default=0,) use to stack letters
  # Returns: list of geom_polygon(s) that can be added to any ggplot2 bject/plot
  
  # check letter
  if(!letter %in% c("A", "C", "G", "T")){
    warning(paste0(letter, " is not a valid base letter!"))
    return(NA_real_)
  }
  
  letter.colours <- brewer.pal(5, "Set1")[c(3,2,5,1)]
  
  # A
  if(letter == "A"){
    a.pol1 <- data.frame(
      x=c(.1, .25, .35, .65, .75, .9, .6, .4),
      y=c(0, 0, .35, .35, 0, 0, 1, 1)
    )
    a.pol2 <- data.frame(
      x=c(.4, .6, .5),
      y=c(.5, .5, .85)
    )
    # scale
    a.pol1$x <- a.pol1$x + xintercept
    a.pol1$y <- a.pol1$y * scale + yintercept
    a.pol2$x <- a.pol2$x + xintercept
    a.pol2$y <- a.pol2$y * scale + yintercept
    
    l <- list(
      geom_polygon(data=a.pol1, aes(x=x, y=), fill=letter.colours[1]),
      geom_polygon(data=a.pol2, aes(x=x, y=), fill="white")
    )
    
    # C
  }else if(letter == "C"){
    c.pol1 <- circleFun(c(.5,.5), 1, npoints = 100)
    c.pol2 <- circleFun(c(.5,.5), .75, npoints =100)
    c.pol3 <- data.frame(x=c(.75, 1, 1, .75), y=c(.3, .3, .7, .7))
    # scale
    c.pol1$x <- c.pol1$x + xintercept
    c.pol2$x <- c.pol2$x + xintercept
    c.pol3$x <- c.pol3$x + xintercept
    c.pol1$y <- c.pol1$y * scale + yintercept
    c.pol2$y <- c.pol2$y * scale + yintercept
    c.pol3$y <- c.pol3$y * scale + yintercept
    # report list
    l <- list(
      geom_polygon(data=c.pol1, aes(x=x, y=y), fill=letter.colours[2]),
      geom_polygon(data=c.pol2, aes(x=x, y=y), fill="white"),
      geom_polygon(data=c.pol3, aes(x=x, y=y), fill="white")
    )
    # G
  }else if(letter == "G"){
    g.pol1 <- circleFun(c(.5,.5), 1, npoints = 100)
    g.pol2 <- circleFun(c(.5,.5), .75, npoints =100)
    g.pol3 <- data.frame(x=c(.75, 1, 1, .75), y=c(.3, .3, .7, .7))
    g.pol4 <- data.frame(x=c(.75, .9, .9, .5, .5, .75), y=c(0, 0, .45, .45, .35, .35))
    # scale
    g.pol1$x <- g.pol1$x + xintercept
    g.pol2$x <- g.pol2$x + xintercept
    g.pol3$x <- g.pol3$x + xintercept
    g.pol4$x <- g.pol4$x + xintercept
    g.pol1$y <- g.pol1$y * scale + yintercept
    g.pol2$y <- g.pol2$y * scale + yintercept
    g.pol3$y <- g.pol3$y * scale + yintercept
    g.pol4$y <- g.pol4$y * scale + yintercept
    #report list
    l <- list(
      geom_polygon(data=g.pol1, aes(x, y), fill=letter.colours[3]),
      geom_polygon(data=g.pol2, aes(x=x, y=y), fill="white"),
      geom_polygon(data=g.pol3, aes(x=x, y=y), fill="white"),
      geom_polygon(data=g.pol4, aes(x=x, y=y), fill=letter.colours[3])
    )
    # T
  }else if(letter == "T"){
    t.pol <- data.frame(
      x=c(.425, .575, .575, .9, .9, .1, .1, .425),
      y=c(0, 0, .85, .85, 1, 1, .85, .85)
    )
    # scale
    t.pol$x <- t.pol$x + xintercept
    t.pol$y <- t.pol$y * scale + yintercept
    # report list
    l <- list(geom_polygon(data=t.pol, aes(x=x, y=y), fill=letter.colours[4]))
  }
  
  return(l)
  
}

plotICM <- function(icm){
  # Plot function to generate a polygon based base letter representation of the ICM matrix/motif
  # Input:
  #   icm: 4 * X dataframe (rows: A, C, G, T, columns: as many positions as there are in the motif)
  # Returns: ggplot2 object with the motif plot
  
  # space for icm data frame checking
  if(nrow(icm) != 4){
    warnings("Must be a 4 row (A,C,G,T) data frame")
    return(NA_character_)
  }
  
  # convert data frame to data frame listing the bit content of letters stagged on top of each other with the largst on top
  stagged.icm <- apply(icm, 2, function(x){
    temp <- data.frame(from=rep(0,4), to=rep(0,4)) # init from_to table
    row.names(temp) <- c("A", "C", "G", "T")
    ordered <- order(x) # order entries
    cum <- 0 # init cumulative value
    # run over bases and count up cum and set from to values for plot
    for(i in c(1:4)){
      temp[ordered[i], "from"] <- cum
      cum <- cum + x[ordered[i]]
      temp[ordered[i], "to"] <- cum
    }
    temp <- unlist(temp)
    return(temp)
  })
  # Convert those to a from - to plot valued, long dataframe for ggplot2
  stagged.icm <- t(stagged.icm)
  temp.df <- data.frame(
    pos=rep(c(1:ncol(icm)), times=4),
    base=rep(c("A", "C", "G", "T"), each=ncol(icm)),
    from=c(stagged.icm[,"from1"], stagged.icm[,"from2"], stagged.icm[,"from3"], stagged.icm[,"from4"]),
    to=c(stagged.icm[,"to1"], stagged.icm[,"to2"], stagged.icm[,"to3"], stagged.icm[,"to4"])
  ) # convert to dataframe
  
  # lay plot base
  p <- ggplot(data.frame(x=factor(c(0:ncol(icm)), levels=c(0:ncol(icm))), y=c(0, max(temp.df$to))), aes(x=x, y=y)) + 
    labs(x="pos", y="bits") +
    xlim(.5, nrow(temp.df)/4+.5) +  ylim(0,2) + theme_bw() + motif_theme 
  
  # add base letters using the position as xintercept, the bit content as scale 
  # and the summed up bit content of the lower letters as yintercept
  for(i in c(1:nrow(temp.df))){
    p <- p + plotBitBase(temp.df$base[i], xintercept = temp.df$pos[i]-1+.5, yintercept = temp.df$from[i], scale= temp.df$to[i] - temp.df$from[i])  
  }
  return(p)
}
