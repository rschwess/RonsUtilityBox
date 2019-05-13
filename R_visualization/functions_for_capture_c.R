# FUNCTIONS FOR CAPTURE C ANALYSIS
# based on pipeline output  + multi cis normalizer
# utilizing genomic ranges, tidyverse, rtracklayer

# Load Packages ==================================
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(reshape2)
# library(AnnotationHub)
library(rtracklayer)
library(GenomicRanges)

cpal <- rev(brewer.pal(3, "Set1")[1:2])

# Reference Files ===============================


# Helper Functions ============================================================
# Helpe for rounding arbitarily
custom_round <- function(x, base) {
  round(x/base)*base
}
custom_floor <- function(x, base) {
  r <- round(x/base)*base
  if(r > x){
    r <- r - base
  }
  return(r)
}
custom_ceil <- function(x, base) {
  r <- round(x/base)*base
  if(r < x){
    r <- r + base
  }
  return(r)
}

vec_custom_floor <- Vectorize(custom_floor)
vec_custom_ceil <- Vectorize(custom_ceil)
vec_custom_round <- Vectorize(custom_round)


# copied from http://crazyhottommy.blogspot.co.uk/2016/02/compute-averagessums-on-granges-or.html
# need to validedte integrity
binnedSum <- function(bins, numvar, mcolname){
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewSums(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}  # need to check integrity as this is a custom function!

binnedMax <- function(bins, numvar, mcolname){
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewMaxs(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}  # need to check integrity as this is a custom function!




# Chromsome Sizes File to GRange object ------
chromSizesToGRange <- function(sizes.file){
  chrom.sizes <- as.tibble(read.table(sizes.file))
  chrom.sizes$start <- 1
  chrom.sizes$end <- chrom.sizes$V2
  chrom.sizes <- chrom.sizes[,-c(2)]
  names(chrom.sizes)[1] <- 'seqnames'
  chrom.sizes <- chrom.sizes[chrom.sizes$seqnames != 'chrM',]
  gr <- makeGRangesFromDataFrame(chrom.sizes, seqinfo = as.character(chrom.sizes$seqnames))
  return(gr)
}

# GR to tibble DF ----------------------------
GRangeToDataframe <- function(grange){
  d <- as.tibble(as.data.frame(grange))
  d <- d %>%
    mutate(pos=as.integer((end-start)/2+start))
  d <- d[,-c(2:5)]
  d <- d[,c(1,ncol(d), 2:(ncol(d)-1))]
  return(d)
}

# Manipulate Data Frames --------------------
getStateViewpoints <- function(data, vwps, state.col="mes.state", viewp.col="viewp", states.select="off_pc", melt=FALSE){
  # Function to retrieve viewpoints matching to their mES/CD71 state
  # where states.select is a character string or character vector indicating which states to retrieve
  # melt indicates if the data frame is melted --> will treat separately
  vwps <- vwps %>% filter(!!as.name(state.col) %in% states.select)
  if(melt == FALSE){
    data <- data %>% select(1:2, which(names(data) %in% vwps$name))
  }else{
    data <- data %>% filter(!!as.name(viewp.col) %in% vwps$name)
  }
  return (data)
}

# helper function to round to arbitrary numbers
roundTo <- function(x, to=500){
  rounded <- (x %/% to) * to
  r <- x %% to
  if(r >= to/2){
    rounded <- rounded + to
  }
  return(rounded)
}
roundTo <- Vectorize(roundTo)  # vectorize it


# Fit three component distal decay function -----------------------------------
fitTripleDistanceDecay<- function(
  data,  # melted datatframe with cond, distance, value, .. columns,
  bin_size = 500,
  conditions = c("mES", "CD71"),
  close.threshold = 2500,
  near.threshold = 150000,
  intermediate.threshold = 2500000,
  far.threshold = 20000000,
  colours = rev(brewer.pal(3, "Set1")[1:2])
){
  # Fitting a triple component log - log linear decay function
  # Returns three fit formula components, distance thresholds for applying and fit plots
  
  # # 1) Bin Distances (rounded to 500 bp) ---------------------------
  # print("Rounding Bins ...")
  # data <- data %>%
  #   mutate(cond = factor(cond, levels=conditions)) %>%
  #   mutate(distance = roundTo(distance, to=bin_size))
  # 
  
  # 2) ### NEAR FIT ### --------------------------------------------
  print("Fitting Near Cis ...")
  min.dist.for.fit <- close.threshold
  max.dist.for.fit <- near.threshold
  # summarize per bin and plot
  data.mean <- data %>%
    filter(value != 0) %>%
    filter(distance >= min.dist.for.fit & distance <= max.dist.for.fit) %>%
    group_by(cond, distance) %>%
    summarise(value = mean(value), n=n())
  # fit linear model to log values
  # one model per condition?
  dfA.mean <- data.mean %>% 
    filter(cond == conditions[1])
  dfB.mean <- data.mean %>% 
    filter(cond == conditions[2])
  ll.model.A <- lm(log(dfA.mean$value) ~ log(dfA.mean$distance))
  ll.intercept.A <- ll.model.A$coefficients[1]
  ll.slope.A <- ll.model.A$coefficients[2]
  ll.model.B <- lm(log(dfB.mean$value) ~ log(dfB.mean$distance))
  ll.intercept.B <- ll.model.B$coefficients[1]
  ll.slope.B <- ll.model.B$coefficients[2]
  # get fitted values for curve
  data.fit <- data.mean %>% 
    mutate(log.distance = log(distance)) %>%
    mutate(ll.fit = if_else(cond == conditions[1], exp(ll.slope.A * log.distance + ll.intercept.A), exp(ll.slope.B * log.distance + ll.intercept.B)))
  # Make Plot with fit curve
  p.value <- ggplot(data.fit, aes(x = distance, y = value)) + geom_point() +
    geom_line(aes(y = ll.fit, colour=cond), size=1) + 
    scale_color_manual(values = cpal) + 
    ggtitle("Interaction Value ~ Distance (log-log linear fit)") +
    science_theme + facet_wrap(~cond, nrow = 2)
  p.loglog <- ggplot(data.fit, aes(x = log(distance), y = log(value))) + geom_point() + 
    geom_line(aes(y = log(ll.fit), colour=cond), size=1) + 
    scale_color_manual(values = cpal) +
    ggtitle("Linear Log-Log relation and Fit") + 
    science_theme + facet_wrap(~cond, nrow = 2)
  p <- plot_grid(p.value, p.loglog, nrow=1, rel_widths = c(2,1.5))
  # save near cis objects 
  plot.fit.near <- p
  ll.model.slope.A.near <- ll.slope.A
  ll.model.slope.B.near <- ll.slope.B
  ll.model.intercept.A.near <- ll.intercept.A
  ll.model.intercept.B.near <- ll.intercept.B
  
  # 3) ### INTERMEDIATE FIT ### ---------------------------------
  print("Fitting Inter Cis ...")
  min.dist.for.fit <- near.threshold
  max.dist.for.fit <- intermediate.threshold
  # summarize per bin and plot
  data.mean <- data %>%
    filter(value != 0) %>%
    filter(distance >= min.dist.for.fit & distance <= max.dist.for.fit) %>%
    group_by(cond, distance) %>%
    summarise(value = mean(value), n=n())
  # fit linear model to log values
  # one model per condition?
  dfA.mean <- data.mean %>% 
    filter(cond == conditions[1])
  dfB.mean <- data.mean %>% 
    filter(cond == conditions[2])
  ll.model.A <- lm(log(dfA.mean$value) ~ log(dfA.mean$distance))
  ll.intercept.A <- ll.model.A$coefficients[1]
  ll.slope.A <- ll.model.A$coefficients[2]
  ll.model.B <- lm(log(dfB.mean$value) ~ log(dfB.mean$distance))
  ll.intercept.B <- ll.model.B$coefficients[1]
  ll.slope.B <- ll.model.B$coefficients[2]
  # get fitted values for curve
  data.fit <- data.mean %>% 
    mutate(log.distance = log(distance)) %>%
    mutate(ll.fit = if_else(cond == conditions[1], exp(ll.slope.A * log.distance + ll.intercept.A), exp(ll.slope.B * log.distance + ll.intercept.B)))
  # Make Plot with fit curve
  p.value <- ggplot(data.fit, aes(x = distance, y = value)) + geom_point() +
    geom_line(aes(y = ll.fit, colour=cond), size=1) + 
    scale_color_manual(values = cpal) + 
    ggtitle("Interaction Value ~ Distance (log-log linear fit)") +
    science_theme + facet_wrap(~cond, nrow = 2)
  p.loglog <- ggplot(data.fit, aes(x = log(distance), y = log(value))) + geom_point() + 
    geom_line(aes(y = log(ll.fit), colour=cond), size=1) + 
    scale_color_manual(values = cpal) +
    ggtitle("Linear Log-Log relation and Fit") + 
    science_theme + facet_wrap(~cond, nrow = 2)
  p <- plot_grid(p.value, p.loglog, nrow=1, rel_widths = c(2,1.5))
  # save near cis objects 
  plot.fit.intermediate <- p
  ll.model.slope.A.intermediate <- ll.slope.A
  ll.model.slope.B.intermediate <- ll.slope.B
  ll.model.intercept.A.intermediate <- ll.intercept.A
  ll.model.intercept.B.intermediate <- ll.intercept.B
  
  # 4) ### FAR FIT ### --------------------------
  print("Fitting Far Cis ...")
  min.dist.for.fit <- intermediate.threshold
  max.dist.for.fit <- far.threshold
  # summarize per bin and plot
  data.mean <- data %>%
    filter(value != 0) %>%
    filter(distance >= min.dist.for.fit & distance <= max.dist.for.fit) %>%
    group_by(cond, distance) %>%
    summarise(value = mean(value), n=n())
  # fit linear model to log values
  # one model per condition?
  dfA.mean <- data.mean %>% 
    filter(cond == conditions[1])
  dfB.mean <- data.mean %>% 
    filter(cond == conditions[2])
  ll.model.A <- lm(log(dfA.mean$value) ~ log(dfA.mean$distance))
  ll.intercept.A <- ll.model.A$coefficients[1]
  ll.slope.A <- ll.model.A$coefficients[2]
  ll.model.B <- lm(log(dfB.mean$value) ~ log(dfB.mean$distance))
  ll.intercept.B <- ll.model.B$coefficients[1]
  ll.slope.B <- ll.model.B$coefficients[2]
  # get fitted values for curve
  data.fit <- data.mean %>% 
    mutate(log.distance = log(distance)) %>%
    mutate(ll.fit = if_else(cond == conditions[1], exp(ll.slope.A * log.distance + ll.intercept.A), exp(ll.slope.B * log.distance + ll.intercept.B)))
  # Make Plot with fit curve
  p.value <- ggplot(data.fit, aes(x = distance, y = value)) + geom_point() +
    geom_line(aes(y = ll.fit, colour=cond), size=1) + 
    scale_color_manual(values = cpal) + 
    ggtitle("Interaction Value ~ Distance (log-log linear fit)") +
    science_theme + facet_wrap(~cond, nrow = 2)
  p.loglog <- ggplot(data.fit, aes(x = log(distance), y = log(value))) + geom_point() + 
    geom_line(aes(y = log(ll.fit), colour=cond), size=1) + 
    scale_color_manual(values = cpal) +
    ggtitle("Linear Log-Log relation and Fit") + 
    science_theme + facet_wrap(~cond, nrow = 2)
  p <- plot_grid(p.value, p.loglog, nrow=1, rel_widths = c(2,1.5))
  # save near cis objects 
  plot.fit.far <- p
  ll.model.slope.A.far <- ll.slope.A
  ll.model.slope.B.far <- ll.slope.B
  ll.model.intercept.A.far <- ll.intercept.A
  ll.model.intercept.B.far <- ll.intercept.B
  
  # 5) Make overall log log cis plot --------------------------
  print("Wrapping Up ...")
  data.mean <- data %>%
    filter(value != 0) %>%
    filter(distance >= close.threshold) %>%
    group_by(cond, distance) %>%
    summarise(value = mean(value), n=n())
  plot.loglog_all.cis <- ggplot(data.mean, aes(x = log(distance), y = log(value), col=cond)) + geom_point() + 
    geom_vline(xintercept = c(log(near.threshold), log(intermediate.threshold), log(far.threshold)), linetype = "dotted") + 
    facet_wrap(~cond, nrow=2) + scale_color_manual(values = cpal) + science_theme
  
  # 6) Assemble Return --------------------------
  out <- list(
    "near" = list(
      "plot.fit" = plot.fit.near, 
      "slope.A" = ll.model.slope.A.near,
      "intercept.A" = ll.model.intercept.A.near, 
      "slope.B" = ll.model.slope.B.near,
      "intercept.B" = ll.model.intercept.B.near
      ),
    "inter" = list(
      "plot.fit" = plot.fit.intermediate, 
      "slope.A" = ll.model.slope.A.intermediate,
      "intercept.A" = ll.model.intercept.A.intermediate, 
      "slope.B" = ll.model.slope.B.intermediate,
      "intercept.B" = ll.model.intercept.B.intermediate
    ),
    "far" = list(
      "plot.fit" = plot.fit.far, 
      "slope.A" = ll.model.slope.A.far,
      "intercept.A" = ll.model.intercept.A.far, 
      "slope.B" = ll.model.slope.B.far,
      "intercept.B" = ll.model.intercept.B.far
    ),
    "plot.all" = plot.loglog_all.cis
  )
  return(out)
    
    
}


# Predict distance function with three component distal decay function -----------------------------------
predictTripleDistanceDecay<- function(
  data,  # melted datatframe with cond, distance, value, .. columns,
  ll,  # log log triple component fit model
  newcol = "f.d", # name of columns to add 
  conditions = c("mES", "CD71"),
  close.threshold = 2500,
  near.threshold = 150000,
  intermediate.threshold = 2500000,
  far.threshold = 20000000
){
  
  # add new data column
  data <- data %>% mutate(newcol = 0)
  
  # predict separately for conditions
  # predict separately for near, inter and far distances
  
  # A 
  data[(data$cond == conditions[1] & data$distance <= near.threshold), "newcol"] <- 
    exp(ll$near$slope.A * log(data[(data$cond == conditions[1] & data$distance <= near.threshold), "distance"]) + ll$near$intercept.A)
  data[(data$cond == conditions[1] & data$distance > near.threshold & data$distance <= intermediate.threshold), "newcol"] <- 
    exp(ll$inter$slope.A * log(data[(data$cond == conditions[1] & data$distance > near.threshold & data$distance <= intermediate.threshold), "distance"]) + ll$inter$intercept.A)
  data[(data$cond == conditions[1] & data$distance > intermediate.threshold), "newcol"] <- 
    exp(ll$far$slope.A * log(data[(data$cond == conditions[1] & data$distance > intermediate.threshold), "distance"]) + ll$far$intercept.A)
  
  # B
  data[(data$cond == conditions[2] & data$distance <= near.threshold), "newcol"] <- 
    exp(ll$near$slope.B * log(data[(data$cond == conditions[2] & data$distance <= near.threshold), "distance"]) + ll$near$intercept.B)
  data[(data$cond == conditions[2] & data$distance > near.threshold & data$distance <= intermediate.threshold), "newcol"] <- 
    exp(ll$inter$slope.B * log(data[(data$cond == conditions[2] & data$distance > near.threshold & data$distance <= intermediate.threshold), "distance"]) + ll$inter$intercept.B)
  data[(data$cond == conditions[2] & data$distance > intermediate.threshold), "newcol"] <- 
    exp(ll$far$slope.B * log(data[(data$cond == conditions[2] & data$distance > intermediate.threshold), "distance"]) + ll$far$intercept.B)
  
  # rename column
  names(data)[ncol(data)] <- newcol

  # return data
  return(data)

}
  
# roundTo <- Vectorize(roundTo)  # vectorize it



# Fit three component distal decay function SINGLE CONDITION -----------------------------------
fitTripleDistanceDecaySingleCondition <- function(
  data,  # melted datatframe with mean distance pe, value, .. columns,
  bin_size = 1000,
  close.threshold = 2500,
  near.threshold = 150000,
  intermediate.threshold = 2500000,
  far.threshold = 15000000,
  colours = rev(brewer.pal(3, "Set1")[1:2])
){
  # Fitting a triple component log - log linear decay function
  # Returns three fit formula components, distance thresholds for applying and fit plots
  
  # 1) Bin Distances (rounded to 500 bp) ---------------------------
  print("Rounding Bins ...")
  data <- data %>% mutate(bin.distance = vec_custom_floor(distance, bin.size) + bin.size/2)
  
  # 2) ### NEAR FIT ### --------------------------------------------
  print("Fitting Near Cis ...")
  min.dist.for.fit <- close.threshold
  max.dist.for.fit <- near.threshold
  # summarize per bin and plot
  data.mean <- data %>%
    filter(value != 0) %>%
    filter(bin.distance >= min.dist.for.fit & bin.distance <= max.dist.for.fit) %>%
    group_by(bin.distance) %>%
    summarise(value = mean(value))
  
  # fit linear model to log values
  ll.model <- lm(log(data.mean$value) ~ log(data.mean$bin.distance))
  ll.intercept <- ll.model$coefficients[1]
  ll.slope <- ll.model$coefficients[2]
  
  # get fitted values for curve
  data.fit <- data.mean %>% 
    mutate(log.distance = log(bin.distance)) %>%
    mutate(ll.fit = exp(ll.slope * log.distance + ll.intercept))
  
  # Make Plot with fit curve
  p.value <- ggplot(data.fit, aes(x = bin.distance, y = value)) + geom_point() +
    geom_line(aes(y = ll.fit, colour=cond), size=1, col = colours[1]) + 
    scale_color_manual(values = cpal) + 
    ggtitle("Interaction Value ~ Distance (log-log linear fit)") +
    science_theme
  p.loglog <- ggplot(data.fit, aes(x = log(bin.distance), y = log(value))) + geom_point() + 
    geom_line(aes(y = log(ll.fit), colour=cond), size=1, col = colours[1]) + 
    scale_color_manual(values = cpal) +
    ggtitle("Linear Log-Log relation and Fit") + 
    science_theme
  p <- plot_grid(p.value, p.loglog, nrow=1, rel_widths = c(2,1.5))
  # save near cis objects 
  plot.fit.near <- p
  ll.model.slope.near <- ll.slope
  ll.model.intercept.near <- ll.intercept
  
  # 3) ### INTERMEDIATE FIT ### ---------------------------------
  print("Fitting Inter Cis ...")
  min.dist.for.fit <- near.threshold
  max.dist.for.fit <- intermediate.threshold
  # summarize per bin and plot
  data.mean <- data %>%
    filter(value != 0) %>%
    filter(bin.distance >= min.dist.for.fit & bin.distance <= max.dist.for.fit) %>%
    group_by(bin.distance) %>%
    summarise(value = mean(value))
  
  # fit linear model to log values
  ll.model <- lm(log(data.mean$value) ~ log(data.mean$bin.distance))
  ll.intercept <- ll.model$coefficients[1]
  ll.slope <- ll.model$coefficients[2]
  # get fitted values for curve
  data.fit <- data.mean %>% 
    mutate(log.distance = log(bin.distance)) %>%
    mutate(ll.fit = exp(ll.slope * log.distance + ll.intercept))
  # Make Plot with fit curve
  p.value <- ggplot(data.fit, aes(x = bin.distance, y = value)) + geom_point() +
    geom_line(aes(y = ll.fit, colour=cond), size=1, col = colours[1]) + 
    scale_color_manual(values = cpal) + 
    ggtitle("Interaction Value ~ Distance (log-log linear fit)") +
    science_theme
  p.loglog <- ggplot(data.fit, aes(x = log(bin.distance), y = log(value))) + geom_point() + 
    geom_line(aes(y = log(ll.fit), colour=cond), size=1, col = colours[1]) + 
    scale_color_manual(values = cpal) +
    ggtitle("Linear Log-Log relation and Fit") + 
    science_theme
  p <- plot_grid(p.value, p.loglog, nrow=1, rel_widths = c(2,1.5))
  # save intermediate cis objects 
  plot.fit.intermediate <- p
  ll.model.slope.intermediate <- ll.slope
  ll.model.intercept.intermediate <- ll.intercept
  
  # 4) ### FAR FIT ### --------------------------
  print("Fitting Far Cis ...")
  min.dist.for.fit <- intermediate.threshold
  max.dist.for.fit <- far.threshold
  # summarize per bin and plot
  data.mean <- data %>%
    filter(value != 0) %>%
    filter(bin.distance >= min.dist.for.fit & bin.distance <= max.dist.for.fit) %>%
    group_by(bin.distance) %>%
    summarise(value = mean(value))
  # fit linear model to log values
  ll.model <- lm(log(data.mean$value) ~ log(data.mean$bin.distance))
  ll.intercept <- ll.model$coefficients[1]
  ll.slope <- ll.model$coefficients[2]
  
  # get fitted values for curve
  data.fit <- data.mean %>% 
    mutate(log.distance = log(bin.distance)) %>%
    mutate(ll.fit = exp(ll.slope * log.distance + ll.intercept))
  
  # Make Plot with fit curve
  p.value <- ggplot(data.fit, aes(x = bin.distance, y = value)) + geom_point() +
    geom_line(aes(y = ll.fit, colour=cond), size=1, col = colours[1]) + 
    scale_color_manual(values = cpal) + 
    ggtitle("Interaction Value ~ Distance (log-log linear fit)") +
    science_theme
  p.loglog <- ggplot(data.fit, aes(x = log(bin.distance), y = log(value))) + geom_point() + 
    geom_line(aes(y = log(ll.fit), colour=cond), size=1, col = colours[1]) + 
    scale_color_manual(values = cpal) +
    ggtitle("Linear Log-Log relation and Fit") + 
    science_theme
  p <- plot_grid(p.value, p.loglog, nrow=1, rel_widths = c(2,1.5))
  # save far cis objects 
  plot.fit.far <- p
  ll.model.slope.far <- ll.slope
  ll.model.intercept.far <- ll.intercept
  
  # 5) Make overall log log cis plot --------------------------
  print("Wrapping Up ...")
  data.mean <- data %>%
    filter(value != 0) %>%
    filter(bin.distance >= close.threshold) %>%
    group_by(bin.distance) %>%
    summarise(value = mean(value))
  plot.loglog_all.cis <- ggplot(data.mean, aes(x = log(bin.distance), y = log(value))) +
    geom_point() +
    geom_vline(xintercept = c(log(near.threshold), log(intermediate.threshold), log(far.threshold)), linetype = "dotted") +
    science_theme
  
  # 6) Assemble Return --------------------------
  out <- list(
    "near" = list(
      "plot.fit" = plot.fit.near, 
      "slope" = ll.model.slope.near,
      "intercept" = ll.model.intercept.near
    ),
    "inter" = list(
      "plot.fit" = plot.fit.intermediate, 
      "slope" = ll.model.slope.intermediate,
      "intercept" = ll.model.intercept.intermediate
    ),
    "far" = list(
      "plot.fit" = plot.fit.far, 
      "slope" = ll.model.slope.far,
      "intercept" = ll.model.intercept.far
    ),
    "plot.all" = plot.loglog_all.cis
  )
  return(out)
  
  
}




# Predict distance function with three component distal decay function -----------------------------------
predictTripleDistanceDecaySingleCondition<- function(
  data,  # melted datatframe with cond, distance, value, .. columns,
  ll,  # log log triple component fit model
  newcol = "f.d", # name of columns to add 
  close.threshold = 2500,
  near.threshold = 150000,
  intermediate.threshold = 2500000,
  far.threshold = 20000000
){
  
  # add new data column
  data <- data %>% mutate(newcol = 0)
  
  # predict separately for conditions
  # predict separately for near, inter and far distances
  data[(data$distance <= near.threshold), "newcol"] <- 
    exp(ll$near$slope * log(data[(data$distance <= near.threshold), "distance"]) + ll$near$intercept)
  data[(data$distance > near.threshold & data$distance <= intermediate.threshold), "newcol"] <- 
    exp(ll$inter$slope * log(data[(data$distance > near.threshold & data$distance <= intermediate.threshold), "distance"]) + ll$inter$intercept)
  data[(data$distance > intermediate.threshold), "newcol"] <- 
    exp(ll$far$slope * log(data[(data$distance > intermediate.threshold), "distance"]) + ll$far$intercept)
  
  # rename column
  names(data)[ncol(data)] <- newcol
  
  # return data
  return(data)
  
}