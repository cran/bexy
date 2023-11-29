#-------------------------------------------------------------------------------
# Package BeXY
# Visualize and parse the output of 'BeXY'
# Developed by Madleina Caduff
#-------------------------------------------------------------------------------

library(TeachingDemos)
library(Ternary)

#-------------------------------------------------------------------------------
# Internal methods, not exported
#-------------------------------------------------------------------------------

#' Function to open a file and generate error message if it was not found
#' @param path A file path.
#' @param files A vector of character strings corresponding to file names found in the path.
#' @param pattern A pattern to search for within 'files'.
#' @return A file connection.
#' @keywords internal
.openFile.bexy <- function(path, files, pattern){
  filename <- files[grepl(pattern, files)]
  if (length(filename) == 0){
    stop(paste0("No file containing the pattern '", pattern, "' found in directory '", path, "'!"))
  }
  if (length(filename) > 1){
    stop(paste0("Multiple files containing the pattern '", pattern, "' found in directory '", path, "'!"))
  }
  if (file.size(paste0(path, filename)) == 3){ # empty file (happens if z was not inferred)
    return(data.frame())
  }
  f <- read.table(paste0(path, filename), header = TRUE, check.names = FALSE)
  return(f)
}

#' Function to reformat column names
#' @param file A file whose header should be renamed.
#' @param pattern Where to split.
#' @return A new file with re-formatted column names.
#' @keywords internal
.renameHeader.bexy <- function(file, pattern){
  new <- sapply(names(file), function(x) return(strsplit(x, pattern)[[1]][2]))
  names(file) <- new
  return(file)
}

#' Function that returns all possible sex karyotypes used in BeXY
#' @return A vector of characters.
#' @keywords internal
.getAllPossibleSexKaryotypes.bexy <- function(){
  return(c("XY", "XX", "X0", "XXY", "XYY", "XXX", "XXYY"))
}

#' Function that returns all possible scaffold types used in BeXY
#' @return A vector of characters.
#' @keywords internal
.getAllPossibleScaffoldTypes.bexy <- function(){
  return(c("Autosome", "Y-linked", "X-linked", "Different"))
}

#' Function that returns all possible autosomal trisomy states used in BeXY
#' @param scaffoldName A vector of size one with the name of a scaffold.
#' @return A vector of characters.
#' @keywords internal
.getAllPossibleAutosomalTrisomy.bexy <- function(scaffoldName){
  return(c("diploid", paste0("trisomy_", scaffoldName)))
}

#' Function that returns the scaffold names of z
#' @param names_z A character vector with the column names of z
#' @param numSamples The number of samples
#' @return A vector of characters.
#' @keywords internal
.getScaffoldNames.bexy <- function(names_z, numSamples){
  all <- sapply(names_z, function(x) strsplit(strsplit(x, "z_")[[1]][2], "_")[[1]][1])
  uq <- all[seq(1, length(all), by = numSamples)]
  return(as.character(uq))
}

#' Function that returns the index of a given scaffold & checks if it is valid
#' @param object A bexy object.
#' @param scaffoldName A vector of size one with the name of a scaffold.
#' @return A vector of characters.
#' @keywords internal
.getIndexScaffoldAutosomalTrisomy.bexy <- function(object, scaffoldName){
  # check if trisomy is allowed
  if (!object$hasZ){
    stop("BeXY was not run with autosomal trisomy!
         Please re-run BeXY using the argument '--triploidAutosomes A', where A are the autosome names where trisomy is allowed")
  }
  # check if it is only a single scaffold
  if (length(scaffoldName) > 1){
    stop("Please specify only a single scaffold at a time.")
  }

  # get all possible autosome names
  autosomes <- names(object$posterior_t)[getPosteriorModeScaffoldTypes(object) == "Autosome"]
  errMessage <- paste0("Available autosomes are: ", paste0(autosomes, collapse = ","))

  # get index of match in scaffold names
  ix_counts <- which(colnames(object$posterior_t) == as.character(scaffoldName))

  # check if this is an autosome
  if (getPosteriorModeScaffoldTypes(object)[ix_counts] != "Autosome"){
    stop(paste0("Scaffold '", scaffoldName, "' was not classified as an autosome. ", errMessage, "."))
  }

  # get index of match in z
  ix_z      <- which(row.names(object$posterior_mode_z) == as.character(scaffoldName))

  # check if there is a match
  if (length(ix_z) == 0){
    stop(paste0("Could not detect output of autosomal trisomy for scaffold '", scaffoldName, "'.
                Please re-run BeXY with the argument '--triploidAutosomes ", scaffoldName, "'."))
  }


  return(list(ix_z=ix_z, ix_counts=ix_counts))
}

#' Function that transforms state posteriors of sex karyotype to ternary-format
#' @param object A bexy object.
#' @return A matrix with 3 rows: P(aneuploid), P(XX), P(XY).
#' @keywords internal
.getInputTernary.bexy <- function(object){
  return(rbind(colSums(object$posterior_s[3:nrow(object$posterior_s),]), object$posterior_s[2,], object$posterior_s[1,]))
}

#' Function to select samples that should be highlighted
#' @param object A bexy object.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @return A vector of integers corresponding to indices of samples to highlight
#' @keywords internal
.getSamplesToHighlight.bexy <- function(object, sampleNamesToHighlight){
  if (length(sampleNamesToHighlight) == 0){
    # automatically based on probabilities
    return(which(object$posterior_s[1,] < 0.9 & object$posterior_s[2,] < 0.9))
  }
  # get index of match
  return(which(names(object$posterior_s) %in% sampleNamesToHighlight))
}

#' Function to select samples that should be highlighted for autosomal trisomy
#' @param object A bexy object.
#' @param ix A list containing the index of a specific scaffold within z.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @return A vector of integers corresponding to indices of samples to highlight
#' @keywords internal
.getSamplesToHighlightAutosomalTrisomy.bexy <- function(object, ix, sampleNamesToHighlight){
  if (length(sampleNamesToHighlight) == 0){
    return(which(object$posterior_z[ix$ix_z,] > 0.9))
  }
  # get index of match
  return(which(colnames(object$posterior_z) %in% sampleNamesToHighlight))
}

#' Function to create a legend
#' @param x_coord The x-coordinate for legend.
#' @param y_coord The y-coordinate for legend.
#' @param colors A vector of length three with the colors for aneuploid, XX and XY samples, respectively.
#' @return No return value, called for side effects.
#' @keywords internal
.addLegend.bexy <- function(x_coord, y_coord, colors){
  legend(x_coord, y_coord,
         legend = c("Aneuploid", "XX", "XY"),
         col = colors,
         pch = rep(19, 3),
         bty = "n")
}

#' Function to create a legend for autosomal trisomy
#' @param x_coord The x-coordinate for legend.
#' @param y_coord The y-coordinate for legend.
#' @param scaffoldName The name of the scaffold to plot (has to be an autosome).
#' @param colors A vector of length two with the colors for triploid and diploid samples, respectively.
#' @return No return value, called for side effects.
#' @keywords internal
.addLegendAutosomalTrisomy.bexy <- function(x_coord, y_coord, scaffoldName, colors){
  legend(x_coord, y_coord,
         legend = c(paste0("Trisomy ", scaffoldName), "Euploid"),
         col = colors,
         pch = rep(19, 2),
         bty = "n")
}

#' Function to get the width of a string on the current plotting device in user coordinates
#' @param name A string.
#' @return The width of the input string on the current plotting device in user coordinates.
#' @keywords internal
.getStringWidth.bexy <- function(name){
  vec <- strsplit(name, split = "")[[1]]
  return(sum(sapply(vec, graphics::strwidth)))
}

#' Function to plot a ternary plot of the probabilities for being aneuploid, XX or XY
#' @param object A bexy object.
#' @param colors A vector of length three with the colors for aneuploid, XX and XY samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param ... More parameters passed to TernaryPlot()
#' @return The upper boundary of the plot.
#' @keywords internal
.plotTernary.bexy <- function(object,  colors, sampleNamesToHighlight, label, ...){
  # order: aneuploid (top), XX (right), XY (left)
  tern <- .getInputTernary.bexy(object)
  Ternary::TernaryPlot(point = "up", atip = "", btip = "", ctip = "",
              atip.rotate = 0, atip.pos = 3,
              grid.col = "lightgrey", grid.minor.col = "white", # no grid lines
              axis.labels = FALSE, axis.tick = FALSE, ...
  )
  # make nice labels for triangle nodes
  x_range <- Ternary::TernaryXRange()
  y_range <- Ternary::TernaryYRange()
  top <- y_range[2]+0.005
  text(x = mean(x_range), y = top, labels = "Aneuploid", font = 2, cex = 1)
  text(x = x_range[1]-0.03, y = y_range[1]+0.03, labels = "XY", font = 2, cex = 1)
  text(x = x_range[2]+0.03, y = y_range[1]+0.03, labels = "XX", font = 2, cex = 1)

  posterior_mode <- apply(tern, 2, which.max)
  col_per_sample <- colors[posterior_mode]

  Ternary::TernaryPoints(tern, pch = 19, col = col_per_sample, cex = 1)

  # Add labels
  highlight <- .getSamplesToHighlight.bexy(object, sampleNamesToHighlight)
  if (label & length(highlight) > 0){
    pointXY <- Ternary::TernaryToXY(as.matrix(tern[,highlight]))

    bottom <- pointXY[2,] < 0.1
    left <- pointXY[1,] < 0 & !bottom
    right <- pointXY[1,] >= 0 & !bottom

    segment_length_left <- 0.08
    segment_length_right <- 0.18

    # spread x
    xx <- pointXY[1,]
    xx[bottom] <- TeachingDemos::spread.labs(xx[bottom], 2*strheight('A'), maxiter=1000, min=0)
    if (any(left)){
      name_width <- sapply(names(object$posterior_s)[highlight][left], .getStringWidth.bexy)
      xx[left] <- xx[left] - name_width - graphics::strwidth('A') - segment_length_left # add a small extra space (A)
    }
    xx[right] <- xx[right] + 0.2

    # spread y
    yy <- pointXY[2,]
    yy[bottom] <- yy[bottom] - 0.1
    yy[left] <- TeachingDemos::spread.labs(yy[left], 2*strheight('A'), maxiter=1000, min=0)
    yy[right] <- TeachingDemos::spread.labs(yy[right], 2*strheight('A'), maxiter=1000, min=0)

    # define segment length
    x_segment <- pointXY[1,]
    x_segment[left] <- x_segment[left] - segment_length_left
    x_segment[right] <- x_segment[right] + segment_length_right

    # add label
    text(xx, yy, names(object$posterior_s)[highlight], adj = 0, xpd = TRUE, font = 3)
    segments(pointXY[1,], pointXY[2,],
             x_segment, yy)
  }
  return(list(top = top))
}

#' Function to plot a single 7-cell plot displaying posterior probabilities of karyotypes
#' @param object A bexy object.
#' @param ix The sample index.
#' @param label The sample name used as a label.
#' @param x_coords Left and right x-coordinates for 7-cell plot.
#' @param y_coords Bottom and top y-coordinates for 7-cell plot.
#' @param last Boolean indicating if this is the last sample to be plotted.
#' @param colors A vector of length three with the colors for aneuploid, XX and XY samples, respectively.
#' @return No return value, called for side effects.
#' @keywords internal
.plotTicTacToeWithCoordinates.bexy <- function(object, ix, label, x_coords, y_coords, last, colors){
  # draw grid
  width <- diff(x_coords)
  height <- diff(y_coords)
  width_third <- width / 3
  height_third <- height / 3

  xx <- c(x_coords[1], x_coords[1] + width_third, x_coords[1] + 2*width_third, x_coords[2])
  yy <- c(y_coords[1], y_coords[1] + height_third, y_coords[1] + 2*height_third, y_coords[2])

  # draw vertical lines
  lines(x = rep(xx[1], 2), y = y_coords)
  lines(x = rep(xx[2], 2), y = y_coords)
  lines(x = rep(xx[3], 2), y = y_coords)
  lines(x = rep(xx[4], 2), y = c(yy[1], yy[2]))

  # draw horizontal lines
  lines(x = x_coords, y = rep(y_coords[1], 2))
  lines(x = x_coords, y = rep(y_coords[1] + height_third, 2))
  lines(x = c(xx[1], xx[3]), y = rep(yy[3], 2))
  lines(x = c(xx[1], xx[3]), y = rep(yy[4], 2))

  # color according to probabilities
  createCol <- function(col, prob) { return(grDevices::adjustcolor(col, alpha.f = prob)) }

  rect(xleft = xx[1], ybottom = yy[1], xright = xx[2], ytop = yy[2], col = createCol(colors[1], object$posterior_s[3,ix])) # X
  rect(xleft = xx[2], ybottom = yy[1], xright = xx[3], ytop = yy[2], col = createCol(colors[2], object$posterior_s[2,ix])) # XX
  rect(xleft = xx[3], ybottom = yy[1], xright = xx[4], ytop = yy[2], col = createCol(colors[1], object$posterior_s[6,ix])) # XXX
  rect(xleft = xx[1], ybottom = yy[2], xright = xx[2], ytop = yy[3], col = createCol(colors[3], object$posterior_s[1,ix])) # XY
  rect(xleft = xx[2], ybottom = yy[2], xright = xx[3], ytop = yy[3], col = createCol(colors[1], object$posterior_s[4,ix])) # XXY
  rect(xleft = xx[1], ybottom = yy[3], xright = xx[2], ytop = yy[4], col = createCol(colors[1], object$posterior_s[5,ix])) # XYY
  rect(xleft = xx[2], ybottom = yy[3], xright = xx[3], ytop = yy[4], col = createCol(colors[1], object$posterior_s[7,ix])) # XXYY

  # add sample name
  text(x = x_coords[2] + width_third, y = yy[3] - height_third/2, labels = label, adj = 0, font = 3)

  # add Y labels
  text(x = x_coords[1] - width_third/2, y = yy[2] - height_third/2, labels = "0", srt = 0, adj = 1, cex = 0.8)
  text(x = x_coords[1] - width_third/2, y = yy[3] - height_third/2, labels = "Y", srt = 0, adj = 1, cex = 0.8)
  text(x = x_coords[1] - width_third/2, y = yy[4] - height_third/2, labels = "YY", srt = 0, adj = 1, cex = 0.8)

  # add X labels
  if (last){
    text(x = xx[1] + width_third/2, y = yy[1] - height_third/2 , labels = "X", srt = 90, adj = 1, cex = 0.8)
    text(x = xx[2] + width_third/2, y = yy[1] - height_third/2 , labels = "XX", srt = 90, adj = 1, cex = 0.8)
    text(x = xx[3] + width_third/2, y = yy[1] - height_third/2 , labels = "XXX", srt = 90, adj = 1, cex = 0.8)
  }
}

#' Function to plot ternary plot with 7-cell plots displaying posterior probabilities of karyotypes
#' @param object A bexy object.
#' @param colors A vector of length three with the colors for aneuploid, XX and XY samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param addSquares Boolean indicating whether 7-cell posterior probability square plots should be drawn.
#' @param ... Other parameters used for plotting
#' @return No return value, called for side effects.
#' @keywords internal
.plotTernaryAndTicTacToe.bexy <- function(object, colors, sampleNamesToHighlight, label, addSquares, ...){
  highlight <- .getSamplesToHighlight.bexy(object, sampleNamesToHighlight)

  if (!addSquares | !label | length(highlight) == 0){
    # plot ternary, no tic-tac-toe
    .plotTernary.bexy(object, colors, sampleNamesToHighlight, label, ...)
    .addLegend.bexy(par("usr")[1], par("usr")[4], colors)
  } else {
    oldPar <- par(no.readonly = TRUE)
    on.exit(par(oldPar))
    par(mar = c(3, 1, 1, 10), xpd = TRUE)

    # plot ternary
    top <- .plotTernary.bexy(object, colors, sampleNamesToHighlight, label, ...)

    # determine coordinates for tictactoes
    x_left <- par("usr")[2] + 0.01
    y_bottom <- par("usr")[3]
    y_top <- top$top

    height <- (y_top - y_bottom) / (length(highlight) + 1)
    height <- min(height, 0.1)
    dist <- height / (length(highlight) - 1)
    width <- height

    # add tictactoes
    cur_y_top <- y_top
    cur_y_bottom <- y_top - height
    names <- names(object$posterior_s)[highlight]
    numTicTacToes <- min(length(highlight), 11) # make sure we only plot as many tictactoe as there is space
    for (i in 1:numTicTacToes){
      x_coords <- c(x_left, x_left + width)
      y_coords <- c(cur_y_bottom, cur_y_top)

      .plotTicTacToeWithCoordinates.bexy(object = object,
                                         ix = highlight[i],
                                         label = names[i],
                                         x_coords = x_coords,
                                         y_coords = y_coords,
                                         last = (i == length(highlight)),
                                         colors = colors)

      cur_y_top <- cur_y_bottom - dist
      cur_y_bottom <- cur_y_top - height
    }
    .addLegend.bexy(par("usr")[1], par("usr")[4], colors)
  }
}

#' Function to plot X vs Y counts and color them according to BeXY posterior mode
#' @param object A bexy object.
#' @param colors A vector of length three with the colors for aneuploid, XX and XY samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param ... Other parameters used for plotting
#' @return No return value, called for side effects.
#' @keywords internal
.plotXVsY.bexy <- function(object, colors, sampleNamesToHighlight, label, ...){
  # find Y-chromosome, allow only a single hit
  Y <- which(object$posterior_mode_t == 1)
  if (length(Y) == 0){ stop("BeXY did not detect any Y-linked scaffolds, can therefore not plot X vs Y!") }
  if (length(Y) > 1){ stop("BeXY detected multiple Y-linked scaffolds, can therefore not plot X vs Y!") }

  # find X-chromosome, allow only a single hit
  X <- which(object$posterior_mode_t == 2)
  if (length(X) == 0){ stop("BeXY did not detect any X-linked scaffolds, can therefore not plot X vs Y!") }
  if (length(X) > 1){ stop("BeXY detected multiple X-linked scaffolds, can therefore not plot X vs Y!") }

  # get colors
  posterior_mode <- apply(.getInputTernary.bexy(object), 2, which.max)
  col_per_sample <- colors[posterior_mode]

  # get point types
  pchs <- rep(19, length(col_per_sample))

  # plot normalized counts
  n <- object$counts
  n$individual <- NULL
  n$sequencing_type <- NULL
  norm_counts <- sweep(n, 1, rowSums(n), "/")
  xx <- norm_counts[,X] * 100
  yy <- norm_counts[,Y] * 100
  plot(xx, yy, col = col_per_sample, pch = pchs,
       xlab = "Reads mapped to X (%)",
       ylab = "Reads mapped to Y (%)")
  box()

  # identify samples to highlight
  highlight <- .getSamplesToHighlight.bexy(object, sampleNamesToHighlight)
  if (label & length(highlight) > 0){
    label_y <- yy[highlight]

    # add sample labels
    breaks <- hist(xx[highlight], plot = F)$breaks
    breaks <- unique(round(breaks, digits = 1))
    tmp_x <- sapply(xx[highlight], function(this) return(breaks[which.min(abs(breaks - this))]))
    for (i in unique(tmp_x)){
      label_y[tmp_x == i] <- TeachingDemos::spread.labs(label_y[tmp_x == i],
                                                        2*strheight('A'),
                                                        maxiter=1000,
                                                        min=0,
                                                        max = par("usr")[4] - 0.1*(par("usr")[4]-par("usr")[3]))
    }
    # overshoot x-axis on the right?
    label_x <- xx[highlight] + 0.2
    segment_x <- xx[highlight] + 0.15
    name_width <- sapply(names(object$posterior_s)[highlight], .getStringWidth.bexy)
    overshoot <- which(label_x + name_width > par("usr")[2])
    if (length(overshoot) > 0){
      label_y[overshoot] <- label_y[overshoot] + 2*strheight('A') # move up
      label_x[overshoot] <- par("usr")[2] - name_width[overshoot] - graphics::strwidth('A') # move left (with a bit extra: A)
      segment_x[overshoot] <- par("usr")[2] - graphics::strwidth('A')
    }

    # label!
    text(label_x, label_y,
         names(object$posterior_s)[highlight], adj = 0, xpd = T, font = 3)
    segments(xx[highlight], yy[highlight],
             segment_x, label_y )
  }

  # add legend
  .addLegend.bexy(par("usr")[1] + 0.75*(par("usr")[2] - par("usr")[1]),
                  par("usr")[4],
                  colors)
}

#' Function to plot counts of a specific autosome and color them according to BeXY posterior mode of autosomal trisomy
#' @param object A bexy object.
#' @param scaffoldName The name of the scaffold to plot (has to be an autosome)
#' @param plotCounts A boolean. If TRUE, the counts are plotted. If FALSE, the posterior probabilities are plotted.
#' @param colors A vector of length two with the colors for diploid and triploid samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param ... Other parameters used for plotting
#' @return No return value, called for side effects.
#' @keywords internal
.plotAutosomalTrisomy.bexy <- function(object, scaffoldName, plotCounts, colors, sampleNamesToHighlight, label, ...){
  ix <- .getIndexScaffoldAutosomalTrisomy.bexy(object, scaffoldName)

  # get colors
  posterior_mode <- object$posterior_mode_z[ix$ix_z,] + 1
  col_per_sample <- colors[posterior_mode]

  # get point types
  pchs <- rep(19, length(col_per_sample))


  if (plotCounts){  # plot normalized counts
    n <- object$counts
    n$individual <- NULL
    n$sequencing_type <- NULL
    norm_counts <- sweep(n, 1, rowSums(n), "/")
    xx <- 1:nrow(n)
    yy <- norm_counts[,ix$ix_counts] * 100
    plot(yy, col = col_per_sample, pch = pchs,
         xlab = "Sample",
         ylab = paste0("Reads mapped to scaffold ", scaffoldName, " (%)"))
    box()
  } else { # plot posterior probabilities
    xx <- 1:ncol(object$posterior_z)
    yy <- object$posterior_z[ix$ix_z,] * 100
    plot(yy, col = col_per_sample, pch = pchs,
         xlab = "Sample",
         ylab = paste0("Posterior probability of trisomy ", scaffoldName, " (%)"))
  }


  # identify samples to highlight
  highlight <- .getSamplesToHighlightAutosomalTrisomy.bexy(object, ix, sampleNamesToHighlight)
  if (label & length(highlight) > 0){
    label_y <- yy[highlight]

    # add sample labels
    breaks <- hist(xx[highlight], plot = FALSE)$breaks
    breaks <- c(breaks[1] - diff(breaks[1:2]), breaks)
    breaks <- c(breaks, breaks[length(breaks)] + diff(breaks[1:2]))
    tmp_x <- sapply(xx[highlight], function(this) return(breaks[this >= breaks & this < (breaks+1)])) # ceiling(xx[highlight]*2) / 2
    for (i in unique(tmp_x)){
      label_y[tmp_x == i] <- TeachingDemos::spread.labs(label_y[tmp_x == i],
                                                        2*strheight('A'),
                                                        maxiter=1000,
                                                        min=0 )
    }
    text(xx[highlight]+length(xx)*0.05, label_y,
         names(object$posterior_s)[highlight], adj = 0, xpd = TRUE)
    segments(xx[highlight], yy[highlight],
             xx[highlight]+length(xx)*0.04, label_y )
  }

  # add legend
  .addLegendAutosomalTrisomy.bexy(par("usr")[1],
                                  par("usr")[3] + 0.2*(par("usr")[4] - par("usr")[3]),
                                  scaffoldName,
                                  colors)
}

#' Function to extract specific columns from a data frame or matrix
#' @param what A matrix or data frame
#' @param name The name of the parameter to extract.
#' @return A subset of the matrix or data frame containing only the columns whose names match the argument 'name'
#' @keywords internal
.getColumns.bexy <- function(what, name){
  return(what[,which(grepl(paste0("^", name), colnames(what)))])
}

#-------------------------------------------------------------------------------
# Constructor for BeXY
#-------------------------------------------------------------------------------

#' Accurate Bayesian inference of sex chromosome karyotypes and sex-linked scaffolds from low-depth sequencing data
#'
#' Visualization of results produced by BeXY
#' @param path The path where all the output files of BeXY are located. If NULL, the example files will be loaded.
#' @param readMCMCTrace If TRUE, read the full trace of the MCMC of all parameters (potentially slow).
#' @return An object of type bexy.
#' @export

bexy <- function(path = NULL, readMCMCTrace = FALSE){
  if (is.null(path)) {
    # use example files
    path <- paste0(system.file("extdata", package = "bexy"), "/")
    files <- dir(path)
  } else {
    if (!is.character(path)){ stop("Path should be of type 'character'!") }
    # get all files present in that directory
    files <- list.files(path)
  }

  if (length(files) == 0){ stop(paste0("Directory ", path, " is empty!")) }

  # read file containing state posterior probabilities of s, t and z
  posterior_s <- .openFile.bexy(path, files, "s_statePosteriors")
  posterior_t <- .openFile.bexy(path, files, "t_statePosteriors")
  posterior_z <- .openFile.bexy(path, files, "z_statePosteriors")

  # trim away s_ and t_ from header
  posterior_s <- .renameHeader.bexy(posterior_s, "s_")
  posterior_t <- .renameHeader.bexy(posterior_t, "t_")

  # calculate posterior mode for s, t and z
  posterior_mode_s <- apply(posterior_s, 2, which.max) - 1
  posterior_mode_t <- apply(posterior_t, 2, which.max) - 1
  posterior_mode_z <- apply(posterior_z, 2, which.max) - 1

  # re-format z as a matrix
  hasZ <- ncol(posterior_z) > 0
  if (hasZ){
    numScaffoldsForZ <- ncol(posterior_z) / length(posterior_s)
    scaffoldNamesZ <- .getScaffoldNames.bexy(colnames(posterior_z), length(posterior_s))
    posterior_z <- matrix(as.numeric(posterior_z[2,]), # just keep P(z = 1)
                          nrow = numScaffoldsForZ,
                          ncol = length(posterior_s),
                          byrow = TRUE,
                          dimnames = list(scaffoldNamesZ, names(posterior_s)))
    posterior_mode_z <- matrix(as.numeric(posterior_mode_z),
                               nrow = numScaffoldsForZ,
                               ncol = length(posterior_s),
                               byrow = TRUE,
                               dimnames = list(scaffoldNamesZ, names(posterior_s)))
  }

  # read file containing posterior mean and variance
  meanVar <- .openFile.bexy(path, files, "meanVar.txt")
  posterior_mean <- meanVar[1,]
  posterior_var <- meanVar[2,]

  # read counts
  counts <- .openFile.bexy(path, files, "_counts_filtered.txt")

  # read MCMC trace (if needed)
  trace <- NULL
  if (readMCMCTrace){
    trace <- .openFile.bexy(path, files, "trace.txt")
  }

  result <- list(posterior_s = posterior_s,
                 posterior_t = posterior_t,
                 hasZ = hasZ,
                 posterior_z = posterior_z,
                 posterior_mode_s = posterior_mode_s,
                 posterior_mode_t = posterior_mode_t,
                 posterior_mode_z = posterior_mode_z,
                 posterior_mean = posterior_mean,
                 posterior_var = posterior_var,
                 counts = counts,
                 trace = trace)
  class(result) <- "bexy"

  return(result)
}

#-------------------------------------------------------------------------------
# Methods for printing
#-------------------------------------------------------------------------------

#' Printing a bexy object
#'
#' @param x A bexy object.
#' @param ... Additional parameters passed to print functions.
#' @return No return value, called for side effects.
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' print(bex)
print.bexy <- function(x, ...){
  cat("BeXY posterior mode estimates:\n");
  # print scaffold types
  cat(" - scaffolds: ", sum(x$posterior_mode_t == 0), " autosomes, ",
                        sum(x$posterior_mode_t == 1), " Y-linked scaffolds, ",
                        sum(x$posterior_mode_t == 2), " X-linked scaffolds, ",
                        sum(x$posterior_mode_t == 3), " different scaffolds. ",
                        "\n", sep="");

  # print karyotypes
  string_karyo <- c()
  for (s in 1:length(.getAllPossibleSexKaryotypes.bexy())){
    string_karyo <- c(string_karyo, paste(sum(x$posterior_mode_s == (s-1)), .getAllPossibleSexKaryotypes.bexy()[s]))
  }
  string_karyo <- paste0(string_karyo, collapse = ", ")
  cat(" - sex karyotypes: ", string_karyo, ".\n", sep="");

  # print autosomal trisomies
  if (x$hasZ){
    for (c in 1:nrow(x$posterior_mode_z)){
      numSamplesWithTrisomy <- sum(x$posterior_mode_z[c,])
      if (numSamplesWithTrisomy > 0){
        cat(" - autosomal trisomy of scaffold ", row.names(x$posterior_mode_z)[c], ": ", numSamplesWithTrisomy, " individuals.\n", sep="");
      }
    }
  }

  invisible(x)
}

#' Summarizing a bexy object
#'
#' @param object A bexy object.
#' @param ... Additional parameters passed to summary functions.
#' @return No return value, called for side effects.
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' print(bex)
summary.bexy <- function(object, ...){
  print.bexy(object, ...);
}

#-------------------------------------------------------------------------------
# Methods for parsing results
#-------------------------------------------------------------------------------

#' Getting the posterior mode for each sex karyotype
#'
#' @param object A bexy object.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All samples that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return A character vector containing the sex karyotype classification for each sample
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' getPosteriorModeSexKaryotypes(bex)
getPosteriorModeSexKaryotypes <- function(object, threshold_certainty = 0.9){
  karyo <- .getAllPossibleSexKaryotypes.bexy()[object$posterior_mode_s + 1]
  certain <- apply(object$posterior_s, 2, function(x) any(x > threshold_certainty))
  karyo[!certain] <- "uncertain"
  names(karyo) <- names(object$posterior_mode_s)
  return(karyo)
}

#' Getting the posterior mode for each scaffold type
#'
#' @param object A bexy object.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All scaffolds that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return A character vector containing the scaffold type classification for each scaffold
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' getPosteriorModeScaffoldTypes(bex)
getPosteriorModeScaffoldTypes <- function(object, threshold_certainty = 0.9){
  type <- .getAllPossibleScaffoldTypes.bexy()[object$posterior_mode_t + 1]
  certain <- apply(object$posterior_t, 2, function(x) any(x > threshold_certainty))
  type[!certain] <- "uncertain"
  names(type) <- names(object$posterior_mode_t)
  return(type)
}

#' Getting the posterior mode for each autosomal trisomy
#'
#' @param object A bexy object.
#' @param scaffoldName The scaffold name, must be an autosome.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All samples that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return A character vector containing the trisomy classification for each sample for the given scaffold
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' getPosteriorModeAutosomalTrisomy(bex, "21")
getPosteriorModeAutosomalTrisomy <- function(object, scaffoldName, threshold_certainty = 0.9){
  ix <- .getIndexScaffoldAutosomalTrisomy.bexy(object, scaffoldName)
  classification <- .getAllPossibleAutosomalTrisomy.bexy(scaffoldName)[object$posterior_mode_z[ix$ix_z,] + 1]
  certain <- sapply(object$posterior_z[ix$ix_z,], function(x) return(x > threshold_certainty | x < (1 - threshold_certainty)))
  classification[!certain] <- "uncertain"
  names(classification) <- colnames(object$posterior_mode_z)
  return(classification)
}

#' Write a file with the posterior mode for each sex karyotype
#'
#' @param object A bexy object.
#' @param file The name of the output file.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All samples that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return No return value, called for side effects.
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' tf <- tempfile(fileext = ".txt")
#' writePosteriorModeSexKaryotypes(bex, tf)
writePosteriorModeSexKaryotypes <- function(object, file, threshold_certainty = 0.9){
  karyo <- getPosteriorModeSexKaryotypes(object, threshold_certainty = threshold_certainty)
  karyo <- cbind(names(karyo), as.character(karyo))

  write.table(x = karyo, file = file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = c("sample", "sex_karyotype"))
}

#' Write a file with the posterior mode for each scaffold type
#'
#' @param object A bexy object.
#' @param file The name of the output file.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All scaffolds that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return No return value, called for side effects.
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' tf <- tempfile(fileext = ".txt")
#' writePosteriorModeScaffoldTypes(bex, tf)
writePosteriorModeScaffoldTypes <- function(object, file, threshold_certainty = 0.9){
  types <- getPosteriorModeScaffoldTypes(object, threshold_certainty = threshold_certainty)
  types <- cbind(names(types), as.character(types))

  write.table(x = types, file = file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = c("scaffold", "type"))
}

#' Write a file with the posterior mode for each autosomal trisomy
#'
#' @param object A bexy object.
#' @param scaffoldName The scaffold name, must be an autosome.
#' @param file The name of the output file.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All samples that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return No return value, called for side effects.
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' tf <- tempfile(fileext = ".txt")
#' writePosteriorModeAutosomalTrisomies(bex, "21", tf)
writePosteriorModeAutosomalTrisomies <- function(object, scaffoldName, file, threshold_certainty = 0.9){
  types <- getPosteriorModeAutosomalTrisomy(object, scaffoldName, threshold_certainty = threshold_certainty)
  types <- cbind(names(types), as.character(types))

  write.table(x = types, file = file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = c("sample", paste0("autosomal_ploidy_", scaffoldName)))
}

#' Getting all samples classified as a specific sex karyotype
#'
#' @param object A bexy object.
#' @param karyotype One specific sex karyotype, can be XY, XX, X0, XXY, XYY, XXX or XXYY.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All samples that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return A character vector containing all sample names that are classified as the sex karyotype given by argument 'karyotype'.
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' getSamplesWithThisKaryotype(bex, "XX")
getSamplesWithThisKaryotype <- function(object, karyotype, threshold_certainty = 0.9){
  if (!(karyotype %in% .getAllPossibleSexKaryotypes.bexy())){
    stop(paste0("Unknown karyotype ", karyotype, "!"))
  }
  karyo <- getPosteriorModeSexKaryotypes(object, threshold_certainty = threshold_certainty)
  these <- karyo[karyo == karyotype]
  return(names(these))
}

#' Getting all scaffolds classified as a specific scaffold type
#'
#' @param object A bexy object.
#' @param type One specific scaffold type, can be Autosome, Y-linked, X-linked or Different.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All scaffolds that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return A character vector containing all scaffold names that are classified as the scaffold type given by argument 'type'.
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' getScaffoldsWithThisType(bex, "X-linked")
getScaffoldsWithThisType <- function(object, type, threshold_certainty = 0.9){
  if (!(type %in% .getAllPossibleScaffoldTypes.bexy())){
    stop(paste0("Unknown type ", type, "! Possible types are: ", paste0(.getAllPossibleScaffoldTypes.bexy(), collapse = ", "), "."))
  }
  types <- getPosteriorModeScaffoldTypes(object, threshold_certainty = threshold_certainty)
  these <- types[types == type]
  return(names(these))
}

#' Getting all samples classified as a certain autosomal trisomy
#'
#' @param object A bexy object.
#' @param scaffoldName The scaffold name, must be an autosome.
#' @param threshold_certainty The threshold for certainty on the posterior probabilities. All samples that have a posterior probability less than this threshold are classified as 'uncertain'.
#' @return A character vector containing all sample names that are classified as having an autosomal trisomy of scaffold 'scaffoldName'
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' getSamplesWithAutosomalTrisomy(bex, "21")
getSamplesWithAutosomalTrisomy <- function(object, scaffoldName, threshold_certainty = 0.9){
  classification <- getPosteriorModeAutosomalTrisomy(object, scaffoldName, threshold_certainty = threshold_certainty)
  these <- classification[classification != "diploid" & classification != "uncertain"]
  return(names(these))
}

#' Getting posterior mean rho (ploidy ratio parameter) for each scaffold
#'
#' @param object A bexy object.
#' @return A numeric vector containing the posterior mean of rho for each scaffold
#'
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' getPosteriorMeanRho(bex)
getPosteriorMeanRho <- function(object){
  rho <- .getColumns.bexy(object$posterior_mean, "rho")
  names(rho) <- names(object$posterior_mode_t)
  return(rho)
}

#-------------------------------------------------------------------------------
# Methods for plotting
#-------------------------------------------------------------------------------

#' Plotting a bexy object
#'
#' @param x A bexy object.
#' @param colors A vector of length three with the colors for aneuploid, XX and XY samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param addSquares Boolean indicating whether 7-cell posterior probability square plots should be drawn.
#' @param ... Other parameters used for plotting.
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' plot(bex)
plot.bexy <- function(x, colors = c("turquoise3", "darkorange", "royalblue4"), sampleNamesToHighlight = c(), label = TRUE, addSquares = TRUE, ...){
  plotTernary(x, colors, sampleNamesToHighlight, label, addSquares, ...)
  plotXY(x, colors, sampleNamesToHighlight, label, ...)
}

#' Plotting the ternary of a bexy object
#'
#' @param x A bexy object.
#' @param colors A vector of length three with the colors for aneuploid, XX and XY samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param addSquares Boolean indicating whether 7-cell posterior probability square plots should be drawn.
#' @param ... Other parameters used for plotting.
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' plotTernary(bex)
plotTernary <- function(x, colors = c("turquoise3", "darkorange", "royalblue4"), sampleNamesToHighlight = c(), label = TRUE, addSquares = TRUE, ...){
  if (length(colors) != 3){
    stop("Invalid colors! Should be of length 3 (aneuploid, XX and XY).")
  }
  .plotTernaryAndTicTacToe.bexy(x, colors, sampleNamesToHighlight, label, addSquares, ...)
}

#' Plotting the counts of a bexy object, colored by sex karyotype
#'
#' @param x A bexy object.
#' @param colors A vector of length three with the colors for aneuploid, XX and XY samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param ... Other parameters used for plotting.
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' plotXY(bex)
plotXY <- function(x, colors = c("turquoise3", "darkorange", "royalblue4"), sampleNamesToHighlight = c(), label = TRUE, ...){
  if (length(colors) != 3){
    stop("Invalid colors! Should be of length 3 (aneuploid, XX and XY).")
  }
  .plotXVsY.bexy(x, colors, sampleNamesToHighlight, label, ...)
}

#' Plotting the autosomal trisomies
#'
#' @param x A bexy object.
#' @param scaffoldName The name of the scaffold to plot (has to be an autosome)
#' @param colors A vector of length two with the colors for diploid and triploid samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param ... Other parameters used for plotting.
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' plotAutosomalTrisomy(bex, "21")
plotAutosomalTrisomy <- function(x, scaffoldName, colors = c("deepskyblue4", "darkred"), sampleNamesToHighlight = c(), label = TRUE, ...){
  plotAutosomalTrisomyPosteriorProbabilities(x, scaffoldName, colors, sampleNamesToHighlight, label, ...)
  plotAutosomalTrisomyCounts(x, scaffoldName, colors, sampleNamesToHighlight, label, ...)
}

#' Plotting the autosomal trisomies: posterior probabilities
#'
#' @param x A bexy object.
#' @param scaffoldName The name of the scaffold to plot (has to be an autosome)
#' @param colors A vector of length two with the colors for diploid and triploid samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param ... Other parameters used for plotting.
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' plotAutosomalTrisomyPosteriorProbabilities(bex, "21")
plotAutosomalTrisomyPosteriorProbabilities <- function(x, scaffoldName, colors = c("deepskyblue4", "darkred"), sampleNamesToHighlight = c(), label = TRUE, ...){
  if (length(colors) != 2){
    stop("Invalid colors! Should be of length 2 (diploid and triploid).")
  }
  # plot posterior probabilities
  .plotAutosomalTrisomy.bexy(object = x,
                             scaffoldName = scaffoldName,
                             plotCounts = FALSE,
                             colors = colors,
                             sampleNamesToHighlight = sampleNamesToHighlight,
                             label = label,
                             ...)
}

#' Plotting the autosomal trisomies: counts on that scaffold
#'
#' @param x A bexy object.
#' @param scaffoldName The name of the scaffold to plot (has to be an autosome)
#' @param colors A vector of length two with the colors for diploid and triploid samples, respectively.
#' @param sampleNamesToHighlight A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilites.
#' @param label Boolean indicating whether samples should be labeled.
#' @param ... Other parameters used for plotting.
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy()
#' plotAutosomalTrisomyCounts(bex, "21")
plotAutosomalTrisomyCounts <- function(x, scaffoldName, colors = c("deepskyblue4", "darkred"), sampleNamesToHighlight = c(), label = TRUE, ...){
  if (length(colors) != 2){
    stop("Invalid colors! Should be of length 2 (diploid and triploid).")
  }
  # plot counts
  .plotAutosomalTrisomy.bexy(object = x,
                             scaffoldName = scaffoldName,
                             plotCounts = TRUE,
                             colors = colors,
                             sampleNamesToHighlight = sampleNamesToHighlight,
                             label = label,
                             ...)
}

#' Plot the MCMC trace for a specific parameter
#'
#' @param object A bexy object.
#' @param parameterName The name of the parameter to plot (e.g. rho, logSigma, epsilon, s, t, a, f, ...)
#' @param maxNumPlots The maximum number of plots to plot, default 20.
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{bexy}}
#' @examples
#' bex <- bexy(readMCMCTrace = TRUE)
#' plotMCMCTrace(bex, "rho", maxNumPlots = 1)
plotMCMCTrace <- function(object, parameterName, maxNumPlots = 20){
  if (is.null(object$trace)){
    stop("Trace was not loaded in bexy object! Set 'readMCMCTrace = TRUE' when calling bexy().")
  }
  tr <- .getColumns.bexy(object$trace, parameterName)
  if (is.numeric(tr)){ # single parameter
    plot(tr, xlab = "MCMC iteration (thinned)", ylab = parameterName, main = parameterName)
  } else {
    numPlots <- min(ncol(tr), maxNumPlots)
    for (c in 1:numPlots){ # multiple parameters
      plot(tr[,c], xlab = "MCMC iteration (thinned)", ylab = names(tr)[c], main = names(tr)[c])
    }
  }
}
