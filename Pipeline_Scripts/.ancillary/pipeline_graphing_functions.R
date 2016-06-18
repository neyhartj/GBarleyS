## pipeline_graphing_function.R
# Dependencies:
## ggplot2
## gridExtra

# # Define a function that implements all of the graphing and processing
fastqc.stats.grapher <- function(input.file = NULL, # Name of the input file from the fastqc_stats.py script.
                                 output.file = NULL, # Name of the output PDF. Defaults to a modification of the input filename.
                                 dir = getwd(), # Directory for file output. Defaults to the current working directory.
                                 fancy = TRUE
                                 ) {
  
  # Check if libraries are present. If not, use base R
  if (length(grep(pattern = "ggplot2", installed.packages())) < 1 | length(grep(pattern = "gridExtra", installed.packages())) < 1) {
    fancy <- FALSE
  } 

  if (fancy) {
    library(gridExtra)
    library(ggplot2)
  }
  
  # Error reporting
  if (is.null(input.file)) { stop("The input.file was not specified.") }
  
  # Read in data
  qc.stats <- read.table(input.file, header = T, check.names = F, as.is = T)
  
  #Defining functions
  # Define the color palette
  gg_colors <- c("red", "green", "orange", "blue")
    
  # Create a function to find the last column of base quality
  find.end <- function(x) {
    for(i in 2:length(x)) {
      if(x[i] < x[i-1]) {
        return(i - 1)
      } else {
        if (i == length(x)) {
          return(i-1)
        }}}}
  
  # Histogram of read counts
  # Penalize for scientific notation
  options(scipen = 10000)
  
##########
##########

  if (fancy) {
  # Create the histogram
    seq.count.plot <- ggplot(data = qc.stats, aes_string(x = "Total.Sequences")) + 
      geom_histogram(position = "identity", fill = gg_colors[4], color = "black") +
      xlab(label = "Read Count") + 
      ylab(label = "Frequency") +
      ggtitle("Read Count Distribution") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
  }
  
##########
##########
  
  # Bargraphs of FastQC tests
  # Calculate percentages for pass/fail/warn
  tests.perc <- t(data.frame(apply(X = qc.stats[,c(3:14)], MARGIN = 2, FUN = function(x) { 
    table(factor(x, levels = c("pass", "warn", "fail")))/length(x)
    }), row.names = c("pass", "warn", "fail")))
  #row.names(tests.perc) <- abbreviate(row.names(tests.perc))
  
  
  if (fancy) {
    # Melt the data.frame so it's easier for ggplot to handle
    tests.stack <- stack(as.data.frame(tests.perc))
    tests.stack$test <- rep(row.names(tests.perc), ncol(tests.perc))
    
    # Create plot
    test.results.plot <- ggplot(data = tests.stack, aes_string(x = "test", y = "values", fill = "ind")) +
      geom_bar(stat = "identity") +
      coord_flip() +
      ylab("Percentage of Data Files With Grade") +
      xlab("") +
      scale_fill_manual(values = gg_colors[c(2,3,1)], labels = c("Pass", "Warning", "Fail"), guide = guide_legend(title = "Grade")) +
      ggtitle("FastQC Test Results") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
  }

##########
##########
  
  # Plot of base quality scores
  start <- which(colnames(qc.stats) == "%GC") + 1 # Find the start column of the base quality scores
  bases <- as.numeric(colnames(qc.stats)[start:ncol(qc.stats)])
  
  end <- find.end(bases) + start - 1 # Find the end position
  # Extract the base quality data
  base.qual <- qc.stats[,c(start:end)]
  n <- ncol(base.qual)
  nf <- nrow(base.qual)
  
  # Find the positions of the base quality stats
  med.pos <- seq(1, n, 5)
  per25.pos <- seq(2, n, 5)
  per75.pos <- seq(3, n, 5)
  per10.pos <- seq(4, n, 5)
  per90.pos <- seq(5, n, 5)
  
  # Separate the base quality stats
  med.qual <- apply(X = base.qual[,med.pos], MARGIN = 2, FUN = mean, na.rm = T)
  per25.qual <- apply(X = base.qual[,per25.pos], MARGIN = 2, FUN = mean, na.rm = T)
  per75.qual <- apply(X = base.qual[,per75.pos], MARGIN = 2, FUN = mean, na.rm = T)
  per10.qual <- apply(X = base.qual[,per10.pos], MARGIN = 2, FUN = mean, na.rm = T)
  per90.qual <- apply(X = base.qual[,per90.pos], MARGIN = 2, FUN = mean, na.rm = T)
  base.pos <- unique(bases[c(1:(end-start+1))]) # Vector of base positions
  # Combine to a single data.frame
  base.qual <- as.data.frame(cbind(base.pos, med.qual, per25.qual, per75.qual, per10.qual, per90.qual))
  # Remove NAs across rows
  base.qual <- base.qual[apply(X = base.qual, MARGIN = 1, FUN = function(x) all(!is.na(x))),]
  
  if (fancy) {
    # Graph
    read.qual.plot <- ggplot(data = base.qual, aes(x = base.pos, y = med.qual)) +
      geom_ribbon(aes(ymin = per10.qual, ymax = per90.qual, fill = "red"), alpha = 0.5) + # Add ribbon for 10th and 90th percentile
      geom_ribbon(aes(ymin = per25.qual, ymax = per75.qual, fill = "blue"), alpha = 0.9) + # Add a ribbon for the quartiles
      geom_point(aes(color = 'black'), size = 2) + # Add scatter points
      geom_line(aes(color = 'black'), size = 0.5) + # Add line
      # geom_abline(intercept = 20, slope = 0, size = 1.5, color = base.qual$cols.r[1]) +
      ylim(c(0, (max(per90.qual) + 5))) + 
      ylab("Phred Quality Score") +
      xlab("Base Position") +
      scale_color_manual(name = "", values = c('black' = 'black'), labels = "Median") +
      scale_fill_manual(values = gg_colors[c(4,1)], labels = c("Inter-Quartile Range", "80th Percentile Range"), guide = guide_legend(title = "Ranges")) +
      ggtitle("Per-Base Sequence Quality") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
  }
  
##########  
##########
  
  # Plot of read length distribution
  start <- which(colnames(qc.stats) == "%GC") + 1 # Find the start column of the base quality scores
  bases <- as.numeric(colnames(qc.stats)[start:ncol(qc.stats)])
  end <- find.end(bases) + start - 1 # Find the end position
  start = end + 1 # Determine the start position
  
  # Create a separate data.frame
  read.lengths <- as.data.frame(qc.stats[,c(start:ncol(qc.stats))]); colnames(read.lengths) <- colnames(qc.stats)[start:ncol(qc.stats)]
  # Calculate the mean and sd of the number of reads of a particular length
  if(ncol(read.lengths) > 1) { # If there is more than one read length, continue regularly
    read.lengths.sum <- as.data.frame(t(as.data.frame(apply(X = read.lengths, MARGIN = 2, FUN = function(x) return(quantile(x, na.rm = T)[2:4])))))
    lengths <- as.numeric(row.names(read.lengths.sum)) # Create a vector of the read lengths
  } else {
    dummy <- rep(0, nrow(read.lengths)) # Set a dummy vector of zeros
    single.length <- as.numeric(colnames(read.lengths)) # Pull out the length of the single read length
    read.lengths <- cbind(dummy, read.lengths, dummy); colnames(read.lengths)[c(1,3)] <- c(single.length - 1, single.length + 1)
    read.lengths.sum <- as.data.frame(t(as.data.frame(apply(X = read.lengths, MARGIN = 2, FUN = function(x) return(quantile(x, na.rm = T)[2:4])))))
  }
  # Rename columns
  colnames(read.lengths.sum) <- c("IQR1", "med", "IQR2")
  # Add additional data
  read.lengths.sum$cols.b <- rep(gg_colors[4], nrow(read.lengths.sum))
  read.lengths.sum$read_len <- as.numeric(row.names(read.lengths.sum)) # Create a vector of the read lengths
  # Remove NAs across rows
  read.lengths.sum <- read.lengths.sum[apply(X = read.lengths.sum, MARGIN = 1, FUN = function(x) all(!is.na(x))),]
  
  
  if (fancy) {
    # Set the plot aesthetics
    plot.aes <- aes(x = read_len, y = med, ymin = IQR1, ymax = IQR2)
    
    read.dist.plot <- ggplot(data = read.lengths.sum, plot.aes) +
      geom_ribbon(aes(fill = "blue"), alpha = 0.5) +
      geom_point(aes(color = 'black'), size = 2) +
      geom_line(aes(color = 'black'), size = 0.5) +
      xlab("Read Length") +
      ylab("Number of Reads") +
      scale_color_manual(name = "", values = c('black' = 'black'), labels = "Median") +
      scale_fill_manual(values = "blue", labels = "Inter-Quartile Range", guide = guide_legend(title = "Ranges")) +
      ggtitle("Read Length Distribution") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
  }
  
#####################
#####################
  
  
  # Pepare the output filename
  if (is.null(output.file)) {
    output.file <- paste(strsplit(x = input.file, split = ".txt")[[1]], ".pdf", sep = "")
  }
  
  # Attach the directory
  output.file <- paste(dir, output.file, sep = "/")
  
  # Decide which plotting method to use
  if (fancy) {
    # Using ggplot2
    # Arrange the plots on the results output
    pdf(output.file, width = 8.5, height = 11, title = strsplit(x = input.file, split = ".txt")[[1]])
    grid.arrange(ncol = 2, nrow = 2,
                 seq.count.plot, 
                 test.results.plot, 
                 read.qual.plot, 
                 read.dist.plot
                 )
    dev.off()
  } else {
    pdf(output.file, width = 8.5, height = 11, title = strsplit(x = input.file, split = ".txt")[[1]])
    # Plots in base R
    # Read Number distribution
    par(mfrow = c(2,2), mar=c(5,7,4,2))
    hist(qc.stats$Total.Sequences,
      xlim = range(pretty(range(qc.stats$Total.Sequences))), 
      xlab = "Number of Reads", 
      main = "Read Number Distribution")
    par(las = 2)
    barplot(height = t(tests.perc), 
            horiz = T, 
            col = c("green", "orange", "red"),
            xlab = "Proportion of Tests With Indicated Grade", 
            main = "FastQC Test Grade Proportion")
    # legend("bottom", legend = colnames(tests.perc), fill = c("green", "orange", "red"))
    par(las = 1)
    plot(x = base.pos, y = base.qual$per25.qual,
         pch = 16,
         col = "blue",
         cex = 0.75,
         xlim = range(pretty(range(base.pos))),
         ylim = c(0, (max(base.qual$per90.qual) + 5)),
         xlab = "Base Position in Read",
         ylab = "Phred Scale Quality Score", 
         main = "Per-Base Aggregate Quality Score")
    points(x = base.pos, y = base.qual$per25.qual, cex = 0.75, pch = 16, col = "blue")
    lines(x = base.pos, y = base.qual$per25.qual, cex = 0.75, pch = 16, col = "blue")
    points(x = base.pos, y = base.qual$per75.qual, cex = 0.75, pch = 16, col = "blue")
    lines(x = base.pos, y = base.qual$per75.qual, cex = 0.75, pch = 16, col = "blue")
    points(x = base.pos, y = base.qual$per10.qual, cex = 0.5, pch = 16, col = "red")
    lines(x = base.pos, y = base.qual$per10.qual, cex = 0.5, pch = 16, col = "red")
    points(x = base.pos, y = base.qual$per90.qual, cex = 0.5, pch = 16, col = "red")
    lines(x = base.pos, y = base.qual$per90.qual, cex = 0.5, pch = 16, col = "red")
    points(x = base.pos, y = base.qual$med.qual, cex = 1.25, pch = 16)
    lines(x = base.pos, y = base.qual$med.qual)
    legend(x = "bottomleft", 
           legend = c("Median", "Inter-Quartile Range", "80th Percentile Range"), 
           col = c("black", "blue", "red"), 
           pch = 16, 
           lty = 1, 
           cex = 0.5)
    # Read length dist
    plot(x = read.lengths.sum$read_len, y = read.lengths.sum$IQR1,
         cex = 0.75, pch = 16, col = "blue",
         ylim = range(pretty(range(read.lengths.sum$IQR1, read.lengths.sum$IQR2))),
         ylab = "",
         xlab = "",
         main = "Aggregate Read Length Distribution")
    lines(x = read.lengths.sum$read_len, y = read.lengths.sum$IQR1, lwd = 0.75, col = "blue", pch = 16)
    points(x = read.lengths.sum$read_len, y = read.lengths.sum$IQR2, cex = 0.75, col = "blue", pch = 16)
    lines(x = read.lengths.sum$read_len, y = read.lengths.sum$IQR2, lwd = 0.75, col = "blue", pch = 16)
    points(x = read.lengths.sum$read_len, y = read.lengths.sum$med, cex = 1.25, pch = 16)
    lines(x = read.lengths.sum$read_len, y = read.lengths.sum$med, lwd = 2, pch = 16)
    legend(x = "topleft", 
           legend = c("Median", "Inter-Quartile Range"), 
           col = c("black", "blue"), 
           pch = 16,
           lty = 1, 
           cex = 0.5)
    title(xlab = "Read Length")
    par(mgp = c(5, 1, 0))
    title(ylab = "Number of Reads")
    # Close the PDF
    dev.off()
  } # Close the if-else statement
  
} # Close the function



##############################################
# Define the VCF stats plotting function

vcf.stats.grapher <- function(input.file = NULL, # Name of the input file from the VCF_processor.py script.
                                 output.file = NULL, # Name of the output PDF. Defaults to a modification of the input filename.
                                 dir = getwd() # Directory for file output. Defaults to the current working directory.
) {
  
  # Load libraries
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  
  # Error reporting
  if (is.null(input.file)) { stop("The input.file was not specified.") }
  
  # Read in data
  vcf.stats <- read.table(input.file, header = T)
  
  # Define the color palette
  gg_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  
  # Graph mean site quality
  mean.qual.plot <- ggplot(data = vcf.stats, aes(x = mean.site.quality)) +
    geom_histogram(position = "identity", fill = gg_colors[2], color = "black", alpha = 0.5) +
    xlab(label = "Phred Scale Average Site Genotype Quality\n(Note: Quality Scores > 99 Are Saved as 99)") + 
    ylab(label = "Frequency") +
    ggtitle("Average Site Genotye Quality Distribution") + 
    theme(plot.title = element_text(lineheight = 1, face = "bold"), 
          legend.position = "bottom",
          text = element_text(size = 10),
          axis.title.y = element_text(vjust = 1.3))
  
  # Print the plot
  # mean.qual.plot
  
  # Graph mean site depth
  mean.depth.plot <- ggplot(data = vcf.stats, aes(x = mean.site.depth)) +
    geom_histogram(position = "identity", fill = gg_colors[3], color = "black", alpha = 0.5) +
    xlab(label = "Average Site Depth\n(Read Depth Per Site / Number of Non-Missing Genotypes)") + 
    ylab(label = "Frequency") +
    xlim(c(0,50)) +
    ggtitle("Average Site Depth Distribution") + 
    theme(plot.title = element_text(lineheight = 1, face = "bold"), 
          legend.position = "bottom",
          text = element_text(size = 10),
          axis.title.y = element_text(vjust = 1.3))
  
  # Print the plot
  # mean.depth.plot
  
  # Graph homozygous/heterozygous site depths
  # Adjust data first
  to.melt <- match(x = c("mean.heterozygote.site.depth", "mean.homozygote.alt.site.depth", "mean.homozygote.ref.site.depth"), table = colnames(vcf.stats))
  vcf.stats.melt <- melt(vcf.stats[,to.melt])
  
  depth.by.geno.plot <- ggplot(data = vcf.stats.melt) +
    geom_histogram(aes(x = value, fill = variable), alpha = 0.5, color = "black") +
    facet_grid(variable ~ ., ) +
    xlim(c(0,50)) +
    xlab(label = "Average Read Number Supporting a Genotype") + 
    ylab(label = "Frequency") +
    scale_fill_manual(values = gg_colors[c(4:6)], labels = c("Heterozygous", "Homozygous\nAlternate", "Homozygous\nReference"), guide = guide_legend(title = "Genotype")) +
    ggtitle("Average Depth Per Genotype Call") + 
    theme(plot.title = element_text(lineheight = 1, face = "bold"),
          legend.position = "bottom",
          text = element_text(size = 10),
          axis.title.y = element_text(vjust = 1.3),
          strip.background = element_blank(),
          strip.text.y = element_blank())
  
  # Print the plot
  # depth.by.geno.plot
  
  # Pepare the output filename
  if (is.null(output.file)) {
    output.file <- paste(strsplit(x = input.file, split = ".txt")[[1]], ".pdf", sep = "")
  }
  
  # Attach the directory
  output.file <- paste(dir, output.file, sep = "/")
  
  # Arrange the plots in a grid
  pdf(output.file, width = 8.5, height = 11)
  grid.arrange(arrangeGrob(mean.qual.plot, mean.depth.plot, ncol = 1),
               depth.by.geno.plot,
               ncol = 2)
  dev.off()

} # Close the function




# # Define a function to use the cutadapt/fastx log files to create graphs
cf.stats.grapher <- function(input.file = NULL, # Name of the input file from the fastqc_stats.py script.
                                 output.file = NULL, # Name of the output PDF. Defaults to a modification of the input filename.
                                 dir = getwd(), # Directory for file output. Defaults to the current working directory.
                                 fancy = TRUE
) {
  
  # Check if libraries are present. If not, use base R
  if (length(grep(pattern = "ggplot2", installed.packages())) < 1 | length(grep(pattern = "gridExtra", installed.packages())) < 1) {
    fancy <- FALSE
  } 
  
  if (fancy) {
    library(gridExtra)
    library(ggplot2)
  }
  
  # Error reporting
  if (is.null(input.file)) { stop("The input.file was not specified.") }
  
  # Read in data
  cf_log <- read.table(input.file, header = T)
  
  # Split into flowcell-lane
  flowcell_lane <- gsub(pattern = "_BC[0-9]{1,}.fastq.gz", replacement = "", x = cf_log$Filename)
  cf_log$flowcell_lane <- flowcell_lane
  cf_log$flowcell <- sapply(strsplit(x = flowcell_lane, split = "_"), FUN = function(x) return(x[1]))
  cf_log$lane <- sapply(strsplit(x = flowcell_lane, split = "_"), FUN = function(x) return(x[2]))
  
  # Define parameters
  parameters <- colnames(cf_log)[-c(1,11,12,13)]
  ANOVA.list <- list()
  # Iterate through parameters and run ANOVAs
  for (p in parameters) {
    # Fit the anova model
    aov.fit <- aov(as.formula(paste(p, " ~ flowcell + lane + flowcell_lane", sep = "")), data = cf_log)
    # Summary
    anova.sum <- summary(aov.fit)
    sig.inputs <- anova.sum[[1]][,"Pr(>F)"] < 0.05
    # If any of the predictor vriables are significant, record that
    if (any(sig.inputs)) {
      # Find the significant predictor variables
      sig.inputs <- cbind(row.names(anova.sum[[1]])[sig.inputs], anova.sum[[1]][sig.inputs,5])
      sig.inputs <- gsub(pattern = "[[:space:]]", replacement = "", x = sig.inputs[!is.na(sig.inputs[,1]),])
      # Perform a TukeyHSD
      HSD.fit <- TukeyHSD(aov.fit)
      HSD.sub <- lapply(X = HSD.fit, FUN = function(x) cbind(x[x[,4] < 0.05,4]) )
      
      # Append to the ANOVA list
      ANOVA.list[[p]] <- list(sig.predictors = as.data.frame(sig.inputs), sig.comparisons = HSD.sub)
    } else {
      ANOVA.list[[p]] <- "No significant predictor variables"
    }
  } # Close the for-loop
  
  
  
  #Defining functions
  # Define the color palette
  gg_colors <- c("red", "green", "orange", "blue")
  
  # Penalize for scientific notation
  options(scipen = 10000)
  
  ##### Graph 1 - Histogram of percentage of reads with adapters
  per.with.adapters <- as.data.frame(cf_log$Reads.with.adapters / cf_log$Reads.processed)
  colnames(per.with.adapters) <- "per.with.adapter"
  
  if (fancy) {
    
    adapter.hist <- ggplot(data = per.with.adapters, aes(x = per.with.adapter)) +
      geom_histogram(position = "identity", fill = gg_colors[1], color = "black") +
      xlim(range(pretty(range(per.with.adapters)))) +
      xlab("Proportion of Reads") +
      ylab("Frequency") +
      ggtitle("Distribution of Proportion\nof Reads with Adapter") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
      
  }
  
  ##### Graph 2 - Histogram of Proportion of basepairs written
  bases.processed <- as.data.frame((cf_log$Basepairs.processed - cf_log$Basepairs.written) / cf_log$Basepairs.processed)
  colnames(bases.processed) <- "per.processed"
  
  if (fancy) {
    
    bases.hist <- ggplot(data = bases.processed, aes(x = per.processed)) +
      geom_histogram(position = "identity", fill = gg_colors[2], color = "black") +
      xlim(range(pretty(range(bases.processed)))) +
      xlab("Proportion of Bases") +
      ylab("Frequency") +
      ggtitle("Distribution of Proportion\nof Removed Bases") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
    
  }
  
  ##### Graph 3 - Histogram of the proportion of removed reads
  reads.processed <- as.data.frame(cf_log$Discarded.reads / cf_log$Input.reads)
  colnames(reads.processed) <- "per.processed"
  
  if (fancy) {
    
    reads.hist <- ggplot(data = reads.processed, aes(x = per.processed)) +
      geom_histogram(position = "identity", fill = gg_colors[3], color = "black") +
      xlim(range(pretty(range(reads.processed)))) +
      xlab("Proportion of Reads") +
      ylab("Frequency") +
      ggtitle("Distribution of Proportion\nof Removed Reads") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
    
  }
  
  
  ##### Output
  
  # Pepare the output filename
  if (is.null(output.file)) {
    output.file <- paste(strsplit(x = input.file, split = ".txt")[[1]], ".pdf", sep = "")
  }
  
  # Attach the directory
  output.file <- paste(dir, output.file, sep = "/")
  
  # Decide which plotting method to use
  if (fancy) {
    # Using ggplot2
    # Arrange the plots on the results output
    pdf(output.file, width = 8.5, height = 11, title = strsplit(x = input.file, split = ".txt")[[1]])
    grid.arrange(ncol = 2, nrow = 2,
                 adapter.hist,
                 bases.hist,
                 reads.hist
    )
    dev.off()
    
  } else { # If not fancy, use base R
    
    pdf(output.file, width = 8.5, height = 11, title = strsplit(x = input.file, split = ".txt")[[1]])
    # Plots in base R
    # Read Number distribution
    par(mfrow = c(2,2), mar=c(5,7,4,2))
    hist(per.with.adapters$per.with.adapter, 
         xlim = range(pretty(range(per.with.adapters$per.with.adapter))),
         col = gg_colors[1],
         xlab = "Proportion of Reads",
         main = "Distribution of Proportion of\nReads with Adapter",
         breaks = 20)
    hist(bases.processed$per.processed,
         xlim = range(pretty(range(bases.processed$per.processed))),
         col = gg_colors[2],
         xlab = "Proportion of Bases",
         main = "Distribution of Proportion of\nRemoved Bases",
         breaks = 20)
    hist(reads.processed$per.processed,
         xlim = range(pretty(range(reads.processed$per.processed))),
         col = gg_colors[3],
         xlab = "Proportion of Reads",
         main = "Distribution of Proportion of\nRemoved Reads",
         breaks = 20)
    
    # Close the PDF
    dev.off()
    
  } # Close the if-else statement
  
} # Close the function


  
# # Define a function to use the read mapping log file to create graphs
mapping.stats.grapher <- function(input.file = NULL, # Name of the input file from the fastqc_stats.py script.
                             output.file = NULL, # Name of the output PDF. Defaults to a modification of the input filename.
                             dir = getwd(), # Directory for file output. Defaults to the current working directory.
                             fancy = TRUE
) {
  
  # Check if libraries are present. If not, use base R
  if (length(grep(pattern = "ggplot2", installed.packages())) < 1 | length(grep(pattern = "gridExtra", installed.packages())) < 1) {
    fancy <- FALSE
  } 
  
  if (fancy) {
    library(gridExtra)
    library(ggplot2)
  }
  
  # Error reporting
  if (is.null(input.file)) { stop("The input.file was not specified.") }
  
  # Read in data
  map_log <- read.table(input.file, header = T)
  
  # Split into flowcell-lane
  flowcell_lane <- gsub(pattern = "_BC[0-9]{1,}", replacement = "", x = map_log$Filename)
  map_log$flowcell_lane <- flowcell_lane
  map_log$flowcell <- sapply(strsplit(x = flowcell_lane, split = "_"), FUN = function(x) return(x[1]))
  map_log$lane <- sapply(strsplit(x = flowcell_lane, split = "_"), FUN = function(x) return(x[2]))
  
  # Manipulation of data
  map_log$p.Reads.not.aligned <- map_log$Reads.not.aligned / map_log$Total.reads
  map_log$p.Reads.once.aligned <- map_log$Reads.once.aligned / map_log$Total.reads
  map_log$p.Reads.multiple.aligned <- map_log$Reads.multiple.aligned / map_log$Total.reads
  
  # ANOVA of relevant predictor variables
  # Define parameters
  parameters <- colnames(map_log)[c(9,10,11)]
  ANOVA.list <- list()
  # Iterate through parameters and run ANOVAs
  for (p in parameters) {
    # Fit the anova model
    aov.fit <- aov(as.formula(paste(p, " ~ flowcell + lane + flowcell_lane", sep = "")), data = map_log)
    # Summary
    anova.sum <- summary(aov.fit)
    sig.inputs <- anova.sum[[1]][,"Pr(>F)"] < 0.05
    # If any of the predictor vriables are significant, record that
    if (any(sig.inputs)) {
      # Find the significant predictor variables
      sig.inputs <- cbind(row.names(anova.sum[[1]])[sig.inputs], anova.sum[[1]][sig.inputs,5])
      sig.inputs <- gsub(pattern = "[[:space:]]", replacement = "", x = sig.inputs[!is.na(sig.inputs[,1]),])
      # Perform a TukeyHSD
      HSD.fit <- TukeyHSD(aov.fit)
      HSD.sub <- lapply(X = HSD.fit, FUN = function(x) cbind(x[x[,4] < 0.05,4]) )
      
      # Append to the ANOVA list
      ANOVA.list[[p]] <- list(sig.predictors = as.data.frame(sig.inputs), sig.comparisons = HSD.sub)
    } else {
      ANOVA.list[[p]] <- "No significant predictor variables"
    }
  } # Close the for-loop
  
  # Define the color palette
  gg_colors <- c("red", "green", "orange", "blue")
  
  # Penalize for scientific notation
  options(scipen = 10000)
  
  ## START GRAPHING ##
  
  ##### Graph 1 - Histogram of the proportion of reads not aligning

  if (fancy) {
    
    zero.hist <- ggplot(data = map_log, aes(x = p.Reads.not.aligned)) +
      geom_histogram(position = "identity", fill = gg_colors[1], color = "black") +
      xlim(range(pretty(range(map_log$p.Reads.not.aligned)))) +
      xlab("Proportion of Reads") +
      ylab("Frequency") +
      ggtitle("Distribution of Proportion of\nUnaligned Reads") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
    
  }
  
  
  ##### Graph 2 - Histogram of the proportion of reads aligning once
  
  if (fancy) {
    
    once.hist <- ggplot(data = map_log, aes(x = p.Reads.once.aligned)) +
      geom_histogram(position = "identity", fill = gg_colors[2], color = "black") +
      xlim(range(pretty(range(map_log$p.Reads.once.aligned)))) +
      xlab("Proportion of Reads") +
      ylab("Frequency") +
      ggtitle("Distribution of Proportion of\nReads Aligning Once") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
    
  }
  
  ##### Graph 3 - Histogram of the proportion of reads aligning more than once
  
  if (fancy) {
    
    multiple.hist <- ggplot(data = map_log, aes(x = p.Reads.multiple.aligned)) +
      geom_histogram(position = "identity", fill = gg_colors[3], color = "black") +
      xlim(range(pretty(range(map_log$p.Reads.multiple.aligned)))) +
      xlab("Proportion of Reads") +
      ylab("Frequency") +
      ggtitle("Distribution of Proportion of\nReads Aligning Multiple Times") + 
      theme(plot.title = element_text(lineheight = 1, face = "bold"), 
            legend.position = "bottom",
            text = element_text(size = 10),
            axis.title.y = element_text(vjust = 1.3))
    
  }
  
  
  ##### Output
  
  # Pepare the output filename
  if (is.null(output.file)) {
    output.file <- paste(strsplit(x = input.file, split = ".txt")[[1]], ".pdf", sep = "")
  }
  
  # Attach the directory
  output.file <- paste(dir, output.file, sep = "/")
  
  # Decide which plotting method to use
  if (fancy) {
    # Using ggplot2
    # Arrange the plots on the results output
    pdf(output.file, width = 8.5, height = 11, title = strsplit(x = input.file, split = ".txt")[[1]])
    grid.arrange(ncol = 2, nrow = 2,
                 zero.hist,
                 once.hist,
                 multiple.hist
    )
    dev.off()
    
  } else { # If not fancy, use base R
    
    pdf(output.file, width = 8.5, height = 11, title = strsplit(x = input.file, split = ".txt")[[1]])
    # Plots in base R
    # Read Number distribution
    par(mfrow = c(2,2), mar=c(5,7,4,2))
    hist(map_log$p.Reads.not.aligned, 
         xlim = range(pretty(range(map_log$p.Reads.not.aligned))),
         col = gg_colors[1],
         xlab = "Proportion of Reads",
         main = "Distribution of Proportion of\nUnaligned Reads",
         breaks = 20)
    hist(map_log$p.Reads.once.aligned,
         xlim = range(pretty(range(map_log$p.Reads.once.aligned))),
         col = gg_colors[2],
         xlab = "Proportion of Reads",
         main = "Distribution of Proportion of\nReads Aligning Once",
         breaks = 20)
    hist(map_log$p.Reads.multiple.aligned,
         xlim = range(pretty(range(map_log$p.Reads.multiple.aligned))),
         col = gg_colors[3],
         xlab = "Proportion of Reads",
         main = "Distribution of Proportion of\nReads Aligning Multiple Times",
         breaks = 20)
    
    # Close the PDF
    dev.off()
  } # Close the if-else statement
  
} # CLose the function

# Put the functions in a list
f.list <- list(fastqc = fastqc.stats.grapher, 
               vcf = vcf.stats.grapher, 
               cutfast = cf.stats.grapher,
               readmap = mapping.stats.grapher)


# If there are any arguments, run the function
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
	args <- commandArgs(trailingOnly = TRUE)
	
	# Use the argument to select the appropriate function
	f.list[[args[2]]](input.file = args[1])
}

