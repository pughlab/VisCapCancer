# VisCap: Visualize normalized capture coverage
# Generates heatmap of exon coverage from a directory of sample_interval_summary files
# By Trevor Pugh
# March 2012 - April 2013, Laboratory for Molecular Medicine, PCPGM
# May 2013 - June 2013, Center for Advanced Molecular Diagnostics, Brigham and Women's Hospital
# Sept 2013, Advanced Molecular Diagnostics Lab, Princess Margaret Cancer Centre

svn.revision <- "$Id$"

#Used if dev_dir is set. Disable for deployment
#config.file  <- "VisCapCancer.cfg" 
#source(config.file)

###########
#Libraries#
###########

#source("http://bioconductor.org/biocLite.R")
#biocLite("DNAcopy")

library("cluster")
library("gplots")
library("zoo")
library("DNAcopy")

###########
#Functions#
###########

#Modified winDialog function to run in non-interactive mode
winDialog.nonint <- function (type = c("ok", "okcancel", "yesno", "yesnocancel"), message) 
    {
        #if (!interactive())
        #    stop("winDialog() cannot be used non-interactively")
        type <- match.arg(type)
        res <- .Internal(winDialog(type, message))
        if (res == 10L) 
            return(invisible(NULL))
        c("NO", "CANCEL", "YES", "OK")[res + 2L]
    }

collect_arguments_from_user <- function(lane_dir.prompt, out_dir.prompt, explorer_file, cov_file_pattern, clobber.output.directory) {
    #Collect input and output information from user

    #Input directory
    lane_dir <- choose.dir(caption = "Select a lane directory (e.g. L001):", default = lane_dir.prompt)
    if(is.na(lane_dir)) {
        try( winDialog.nonint(type="ok", "Run canceled. No input lane directory provided."), silent=TRUE)
        q(save="no")
    }
    
    #Output diretory: Attempt to derive batch information from file name, prompt user if unsuccessful
    file1 <- list.files(lane_dir, full.names=TRUE, pattern=cov_file_pattern, recursive=TRUE)[1]
    batch.regex <- "__(B*[0-9]+)"
    batch.match <- regexec(batch.regex, file1)[[1]]
    batch <- substring(file1, batch.match[2], batch.match[2] + attr(batch.match, "match.length")[2] - 1)
    if((infer.batch.for.sub.out_dir == FALSE) || is.na(batch)) {
        out_dir <- choose.dir(caption = "Select an output directory:", default = out_dir.prompt)
        batch   <- basename(out_dir)
    } else {
        out_dir <- paste(out_dir, batch, sep="/")
    }
    if(is.na(out_dir)) {
        try( winDialog.nonint(type="ok", "Run canceled. No output directory provided."), silent=TRUE)
        q(save="no")
    }

    #If output directory already exists, prompt user to overwrite
    if((clobber.output.directory == FALSE) & (file.exists(out_dir))) {
        overwrite <- try( winDialog.nonint(type="yesno", "Output directory already exists. Overwrite?"), silent=TRUE)
        if(overwrite == "NO") {
            shell(paste(explorer_file, out_dir, sep=" "), wait=FALSE)
            q(save="no")
        }
    }
    return(c(lane_dir, out_dir, batch))
}

get_size_and_parts_from_coord <- function(coord) {
	coord <- gsub("^chr","",coord)
    parts <- unlist(strsplit(coord, ":|-"))
    if(length(parts) < 3) {
        parts <- as.integer(append(parts, parts[2]))
    }
    size  <- as.integer(parts[3]) - as.integer(parts[2]) + 1
    return(c(size,parts))
}

visualize_coverage_values <- function(mat.cov, cov.plot.range) {
    mat.sizes <- as.integer(sapply(rownames(mat.cov), get_size_and_parts_from_coord)[1,])
    mat.base <- sweep(mat.cov, 1, mat.sizes, "/")

    #sort.mat by median coverage of exon
    mat.base.sorted <- mat.base[order(apply(mat.base, 1, median), decreasing = TRUE), order(apply(mat.base > max(cov.plot.range), 2, sum), decreasing = TRUE)]
    mat.base.sorted[mat.base.sorted > max(cov.plot.range)] <- max(cov.plot.range)
    
    #heatmap requires at least two samples in matrix
    if(dim(mat.cov)[2] == 1) {
        mat.base.sorted <- cbind(mat.base.sorted, mat.base.sorted)
        colnames(mat.base.sorted) <- c("", "") #Do not report vertical sample label
    } 
    par(bg = "lightgrey")
    heatmap.2(mat.base.sorted,
              dendrogram="none",
              Rowv=FALSE,
              Colv=FALSE,
              scale="none",
              trace="none",
              breaks=max(cov.plot.range) - min(cov.plot.range),
              col=colorRampPalette(c("black","red", "orange", "yellow", "white")),
              labRow = FALSE,
              labCol = sapply(colnames(mat.base.sorted), substr, start=1, stop=12),
              cexCol = 0.75,
              ylab = paste(dim(mat.base.sorted)[1], "targets", sep=" "),
              main = "Target coverage per sample",
              )
    #Add horizontal sample label
    if(dim(mat.cov)[2] == 1) {
        text(0, 0, "") #Required to enable mtext printing on heatmap
        mtext(paste("              ", colnames(mat.cov)), side=1, line=3)
    }
}

make_matrix_from_cov_files <- function(lane_dir, cov_file_pattern, cov_field, out_dir, cov.plot.range, cov.plot.filename) {
    #Join all cov files into a single matrix
    filenames <- list.files(lane_dir, full.names=TRUE, pattern=cov_file_pattern, recursive=TRUE)

    for(file in filenames) {
        tab      <- read.table(file, header=TRUE, row.names=1)
        rownames(tab) <- gsub("^chr", "", rownames(tab))
        col_name <- colnames(tab)[grep(cov_field, colnames(tab))]

        #Make new matrix labelled with appropriate identifier
        #new_head <- gsub(cov_field, "", col_name)                  #sample identifier only
        new_head  <- basename(gsub(cov_file_pattern, "", file))      #derive name from original filename
        new       <- matrix(tab[,col_name], dimnames=list(rownames(tab), new_head) )
        
        #If this is the first file, make a new matrix. Otherwise, join it to the existing matrix using genome coordinates
        if(file == filenames[1]) {
            mat.cov <- new
        } else {
            mrg               <- merge(new, mat.cov, by = "row.names", all = TRUE)
            mat.cov           <- as.matrix(mrg[-1])
            rownames(mat.cov) <- mrg[,1]
        }
    }

    #Set NA entries to 0
    mat.cov[is.na(mat.cov)] <- 0
    #mat.cov <- mat.cov[complete.cases(mat.cov),] #Deprecated: Remove any rows with NA entries

    #Visualize raw coverage values across the batch
    pdf(paste(out_dir, cov.plot.filename, sep="/"))
        visualize_coverage_values(mat.cov, cov.plot.range)
    dev.off()

    # Extract fraction of total coverage assigned to each exon in each sample
    if(dim(mat.cov)[2] == 1) {
        # Note two possible methods, depending on preference for X-chromosome handling.
        # mat.cov.totals <- sum(mat.cov) #use X-chromosome for normalization
        mat.cov.totals <- sum(mat.cov[grep("X",rownames(mat.cov), invert=TRUE),]) # do not use X-chromosome for normalization
    } else {
        # Extract fraction of total coverage assigned to each exon in each sample
        # Note two possible methods, depending on preference for X-chromosome handling.
        # mat.cov.totals <- colSums(mat.cov) #use X-chromosome for normalization
        mat.cov.totals <- colSums(mat.cov[grep("X",rownames(mat.cov), invert=TRUE),]) # do not use X-chromosome for normalization
    }
    mat.frac_cov <- sweep(mat.cov, 2, mat.cov.totals, "/")

    #Order columns alphanumerically
    mat.frac_cov <- mat.frac_cov[,order(colnames(mat.frac_cov)),drop=FALSE]
    
    return(mat.frac_cov)
}

load_interval_names <- function(interval_list_dir, interval_file_pattern, interval_list_coord_adjust) {
    #Read interval list files and make lookup_table for genome coords and interval names
    lookup <- c()
    filenames <- list.files(interval_list_dir, full.names=TRUE, pattern=interval_file_pattern, recursive=FALSE)
    for(file in filenames) {
        tab <- read.table(file, header=FALSE, comment.char = "#", stringsAsFactors=FALSE)
	       if(ncol(tab) < 4)
	         {stop("The interval file should have at least 4 columns, chromosome, 
                                 start, end, and interval name, all tab separated.")}
		#Support bed files that may not have strand column
		if(grepl(".bed$", file) && dim(tab)[2] == 4) {
			tab <- cbind(tab[,1:3], strand=rep(NA, dim(tab)[1]), tab[,4])
		}
        colnames(tab) <- c("chr", "start", "end", "strand", "interval_name")
		tab$chr <- gsub("^chr", "", tab$chr) #remove chr prefix
		tab$start <- tab$start+as.integer(interval_list_coord_adjust)
		#tab$end <- tab$end+as.integer(interval_list_coord_adjust)		#only adjust start position as UCSC beds are zero-indexed
        rownames(tab) <- paste(tab$chr, ":", tab$start, "-", tab$end, sep = "")
        tab$interval_file <- file
        lookup <- rbind(lookup, tab)
    }
    return(lookup)
}

purity_to_log2_thresholds <- function(purity) {
    gains  <- round((purity/100)*log2(3/2), 2)
    losses <- round((purity/100)*log2(1/2), 2)
    thresholds <- c(losses, gains)
    return(thresholds)
}

annotate_interval_names <- function(coords, interval_lookup) {
    #coords = vector of genome coordinates
    #interval_lookup = matrix of genome coordinates and target names
    
    #Add name column to mat containing interval name from lookup table, if provided
    if(is.null(interval_lookup)) {
        labels <- coords
    } else {
        matches <- match(coords, rownames(interval_lookup)) # Find lines in interval_list that match mat's genome coordinates
        labels <- interval_lookup[matches, "interval_name"] #the column name "interval_name" is created by function "load_interval_names"
    }
    return(labels) #returns a vector
}

divide_by_batch_median <- function(mat) {
    #Fractional coverage normalization by median across each exon
    rmeds <- apply(mat, 1, median)         ## calculate row medians
    mat_norm <- sweep(mat, 1, rmeds, "/")  ## divide each entry by row median
    mat_norm <- mat_norm[complete.cases(mat_norm),] ##remove NA entries
	return(mat_norm)
}

heatmap_by_chrom <- function(nmat, analysis_name, ylimits, out_dir) {
	nmat.capped <- nmat
    #Limit matrix for plotting
    nmat.capped[nmat.capped < min(ylimits)] <- min(ylimits)
    nmat.capped[nmat.capped > max(ylimits)] <- max(ylimits)

	# Plot heatmap for each chromosome, compare each exon's value/median ratio across samples, save to files
	pdf(file=paste(out_dir, "/", analysis_name, ".pdf", sep=""))

    #Per-chromosome plots
	chroms = c("all", 1:22,"X","Y", "MT", "M")
	for(chr in chroms){
        if(chr == "all") {
            main_title <- "All Chromosomes"
            matches <- rownames(nmat.capped)
        } else {
            main_title <- paste("Chromosome",chr)
            matches <- grep(paste("^",chr,":",sep=""), rownames(nmat.capped))
        }
		if(length(matches) > 1) {
			par("cex.main" = 0.5)
            #Set heatmap color scale
            steps <- 100
            color_scale <- bluered(steps)
            color_breaks <- seq(ylimits[1], ylimits[2], by=(ylimits[2] - ylimits[1])/steps)
            #color_breaks = seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out=steps+1)

            #exons by samples
            #heatmap.2(mat[matches,], Rowv=NA, Colv=NA, scale="none", cexCol=0.5, cexRow=0.3, col=color_scale, breaks=color_breaks, dendrogram="none", main=title, symbreaks=TRUE, trace="none")
            #samples vs exons
            xlabels <- gsub("__", "\n", colnames(mat))
            heatmap.2(
                        t(nmat.capped[matches,]),
                        Rowv       = NA,
                        Colv       = NA,
                        scale      = "none",
                        cexCol     = 0.5,
                        cexRow     = 0.5,
                        col        = color_scale,
                        breaks     = color_breaks,
                        dendrogram = "none",
                        symbreaks  = FALSE,
                        trace      = "none",
                        tracecol   = "black",
                        main       = main_title,
                        labRow     = xlabels,
                        labCol     = NA
            )
		}
	}
    dev.off()
}

score_boxplot <- function(bplot, llim= log2(1/2), ulim=log2(3/2)) {
    fail_lower_thresh <- bplot$stats[1,] < llim
    fail_upper_thresh <- bplot$stats[5,] > ulim
    fail <- as.logical(fail_lower_thresh + fail_upper_thresh)
    qc_string <- gsub(TRUE, "FAIL", gsub(FALSE, "PASS", fail))
    return(qc_string)
}

boxplot_cnv_matrix <- function(nmat, bplot_name, out_dir, ylimits, iqr_multiplier) {
    plot_nmat <- nmat
    plot_nmat[nmat < min(ylimits)] <- min(ylimits)
    plot_nmat[nmat > max(ylimits)] <- max(ylimits)
    pdf(file=paste(out_dir, "/", bplot_name, ".pdf", sep=""))
    par(las=2, mar=c(12,4,4,2))
    xlabels = sub("__", "\n", colnames(plot_nmat))
    bplot <- boxplot(plot_nmat, range=iqr_multiplier, ylim=ylimits, srt=90, pch=16, cex.axis=0.6, names=xlabels, ylab="Log2 ratio sample/batch median")
    bplot$names <- sub("\n", "__", bplot$names)
    dev.off()
    
    #Write boxplot information to separate file
    bplot$qc <- score_boxplot(bplot)
    bplot.out <- cbind(round(t(bplot$stats), 3), bplot$qc)
    rownames(bplot.out) <- bplot$names
    colnames(bplot.out) <- c("del_threshold", "Q1", "median", "Q3", "amp_threshold", "qc")
    write.table(bplot.out, paste(out_dir, "/", bplot_name, ".xls", sep=""), quote = FALSE, sep = "\t", col.names = NA)
    return(bplot)
}

filter_cnvs <- function(segs, threshold.min_exons, threshold.cnv_log2_cutoffs) {
    #Apply consecutive exon filter
    fsegs <- segs[(segs[,"CNV"] != 0) & (segs[,"Interval_count"] >= threshold.min_exons), , drop = FALSE]

    #Apply zero in normal range filter
    fsegs <- fsegs[(as.numeric(fsegs[,"Loss_threshold"]) < 0) & (as.numeric(fsegs[,"Gain_threshold"]) > 0),, drop = FALSE]

    #Apply hard CNV threshold
    fsegs <- fsegs[(as.numeric(fsegs[,"Median_log2ratio"]) < min(threshold.cnv_log2_cutoffs)) | (as.numeric(fsegs[,"Median_log2ratio"]) > max(threshold.cnv_log2_cutoffs)),, drop = FALSE]

    return(fsegs)
}

call_cnvs_using_boxplots <- function(nmat, ylimits, interval_lookup, threshold.min_exons, iqr_multiplier, threshold.cnv_log2_cutoffs, out_dir, use.boxplot.whisker.threshold) {
    #Plot boxplots to visualize ranges used to detect copy number variation
    bplot <- boxplot_cnv_matrix(nmat, "QC_cnv_boxplot", out_dir, ylimits, iqr_multiplier)
    batch_size <- dim(bplot$stats)[2]
    
    if(use.boxplot.whisker.threshold == TRUE) {
        #Set outlier thresholds by distribution test
        #Use boxplot$stats so whiskers in pdf accurately represent thresholds used
        lbound <- round(bplot$stats[1,], 3)
        ubound <- round(bplot$stats[5,], 3)

        #Use hard-threshold only if it is greater than boxplot thresholds
        lbound[lbound > min(threshold.cnv_log2_cutoffs)] <- min(threshold.cnv_log2_cutoffs)
        ubound[ubound < max(threshold.cnv_log2_cutoffs)] <- max(threshold.cnv_log2_cutoffs)
    } else {
        lbound <- rep(min(threshold.cnv_log2_cutoffs), batch_size)
        ubound <- rep(max(threshold.cnv_log2_cutoffs), batch_size)
    }

    #Make matrix of thresholds
    lbound_mat <-  matrix(rep(lbound, dim(nmat)[1]), ncol=length(lbound), byrow=TRUE, dimnames=list(rownames(nmat),colnames(nmat)))
    ubound_mat <-  matrix(rep(ubound, dim(nmat)[1]), ncol=length(lbound), byrow=TRUE, dimnames=list(rownames(nmat),colnames(nmat)))

    #Flag outliers
    nmat_loutliers <- (nmat < lbound_mat) + 0 #Adding zero converts TRUE/FALSE to 1/0
    nmat_uoutliers <- (nmat > ubound_mat) + 0 #Adding zero converts TRUE/FALSE to 1/0
    #TODO: Estimate number of copies gained or lost (i.e. support -2 and +n)
    
    #Make tracking matrix of all zero values then subtract copies lost and add copies gained
    nmat_cnvs <- matrix(data = 0, nrow = dim(nmat)[1], ncol = dim(nmat)[2], dimnames = list(rownames(nmat),colnames(nmat)))
    nmat_cnvs <- nmat_cnvs - nmat_loutliers + nmat_uoutliers
    gene_exon <- annotate_interval_names(rownames(nmat_cnvs), interval_lookup)
    write.table(cbind.data.frame(gene_exon, nmat_cnvs), paste(out_dir, "/", "cnv_boxplot_outliers", ".xls", sep=""), quote = FALSE, sep = "\t", col.names = NA)
    
    #Merge consecutive calls and write out to file
    all_fsegs <- c()
    for(id in colnames(nmat_cnvs)) {
        segs <- c()
        segs.header <- c("Sample", "CNV", "Genome_start_interval", "Genome_end_interval", "Start_interval","End_interval", "Interval_count", "Min_log2ratio", "Median_log2ratio", "Max_log2ratio", "Loss_threshold", "Gain_threshold", "Batch_size")
        chroms = c(1:22,"X","Y", "MT", "M")
        for(chr in chroms){
            matches <- grep(paste("^",chr,":",sep=""), rownames(nmat))
            if(length(matches) > 1) {
                #Segmentation: Detect runs of consecutive copy number calls
                rl <- rle(nmat_cnvs[matches, id])
                values <- rl$values
                lengths <- rl$lengths
                starts <- c(rownames(nmat_cnvs[matches,])[1], names(rl$lengths[1:length(rl$lengths)-1]))
                ends <- names(rl$values)

                #Calculate rounded log2 ratios of intervals involved
                log2s <- lapply(1:length(starts), function(x) round(nmat[(match(starts[x], names(nmat[,id])):match(ends[x], names(nmat[,id]))), id], 3))
                log2s_min <- unlist(lapply(log2s, min))
                log2s_med <- unlist(lapply(log2s, median))
                log2s_max <- unlist(lapply(log2s, max))

                #Report genome coordinate range
                coordinates_part1 <- data.frame(strsplit(starts, "-"), stringsAsFactors = FALSE)[1,]
                coordinates_part2 <- data.frame(strsplit(ends, "-"), stringsAsFactors = FALSE)[2,]                
                coordinates <- paste(coordinates_part1, coordinates_part2, sep="-")
                
                #Lookup interval names
                start_names <- annotate_interval_names(starts,  interval_lookup)
                end_names <- annotate_interval_names(ends,  interval_lookup)
                
                #Add threshold columns, ensure consistent header
                segs <- rbind(segs, cbind(rep(id, length(values)), values, starts, ends, start_names, end_names, lengths, log2s_min, log2s_med, log2s_max, lbound_mat[1,id], ubound_mat[1,id], batch_size))
                colnames(segs) <- segs.header
            }
        }

        #Handle case where there are no cnvs on any chromosomes in any samples
        if(is.null(segs)) {
            segs <- matrix(data=rep(NA,length(segs.header)), ncol=length(segs.header), dimnames=list("", segs.header))[-1,,drop=FALSE]
        }
        
        #Filter segments
        fsegs <- filter_cnvs(segs, threshold.min_exons, threshold.cnv_log2_cutoffs)
        
        #Convert CNV type values to text
        fsegs[,"CNV"] <- gsub("^1$", "Gain", fsegs[,"CNV"])
        fsegs[,"CNV"] <- gsub("^-1$", "Loss", fsegs[,"CNV"])

        #Convert infinite values to large, numerical value
        large_value <- "10"
        num_cols <- c("Min_log2ratio", "Median_log2ratio", "Max_log2ratio", "Loss_threshold", "Gain_threshold")
        fsegs[,num_cols] <- gsub("Inf", large_value, fsegs[,num_cols])

        #Number CNVs for later labeling on visual output
        if(dim(fsegs)[1] > 0) {
            CNV_id <- seq(1, dim(fsegs)[1])
        } else {
            CNV_id <- c()
        }
        fsegs <- cbind(fsegs[,1, drop=FALSE],
                       CNV_id,
                       fsegs[,2:dim(fsegs)[2], drop=FALSE])
       
        #Add fsegs to master fsegs tracking matrix
        all_fsegs <- rbind(all_fsegs, fsegs)
    }
    return(list(bplot, all_fsegs))
}

space_label_names <- function(label_names) {
    spacer <- "     "
    if(length(label_names) > 1) {
        alt_positions <- seq(from=2, to=length(label_names), by=2)
        label_names[alt_positions] <- paste(label_names[alt_positions], spacer)
    }
    return(label_names)
}

plot_custom_gridlines <- function(custom_gridlines, nmat_chr, shift) {
    #Parse out nmat genome coordinates
    nmat_chr.parts <- strsplit(rownames(nmat_chr), ":|-")
    nmat_chr.chrs   <- sub("chr", "", sapply(nmat_chr.parts, "[", 1))
    nmat_chr.starts <- sapply(nmat_chr.parts, "[", 2)
    
    #Parse out custom gridlines genome coordinates, assign plot positions
    custom_gridlines.sub <- custom_gridlines[sub("chr", "", custom_gridlines[,1]) %in% nmat_chr.chrs,,drop=FALSE]
    custom_gridlines.sub <- cbind(custom_gridlines.sub, plot.start=
                                  sapply(custom_gridlines.sub[,2], function(x) {
                                    rank(as.numeric(c(x, nmat_chr.starts)),
                                         ties.method="first")[1] - 1 - shift}))
    custom_gridlines.sub <- cbind(custom_gridlines.sub, plot.end=
                                  sapply(custom_gridlines.sub[,3], function(x) {
                                    rank(as.numeric(c(x, nmat_chr.starts)),
                                         ties.method="first")[1] - 1 - shift}))
    custom_gridlines.sub <- cbind(custom_gridlines.sub, plot.text=apply(custom_gridlines.sub[,c("plot.start","plot.end")], 1, mean))
    abline(v=custom_gridlines.sub[,c("plot.start","plot.end")], col="grey")    
    #Label using name from gridlines file
    #mtext(custom_gridlines.sub[,"name"], side=3, at=cutom_gridlines.sub[,"line.text"])
}

plot_viscap_data <- function(nmat_chr, title, name, labels_rle, labels_size, ylimits, draw_grid_lines, custom_gridlines) {
    shift <- -0.5 #shift for plotting gene guide lines
    par(pch=16)
    plot(nmat_chr, xlim=c(1 + shift,length(nmat_chr)), ylim=ylimits, main=title, ylab = "Log2 ratio sample/batch median", xaxt="n", xlab="")
    mtext(name, line=3, cex=0.5, adj=0)

    #Mark data points outside of ylimits and scaled to fit plot
    off_scale <- as.numeric(which((nmat_chr <= ylimits[1]) | (nmat_chr >= ylimits[2])))
    points(off_scale, nmat_chr[off_scale], pch=22)

    #Plot gene or chomosome grid lines
    section_ends <- (unlist(lapply(seq(1:length(labels_rle$lengths)), function(x) sum(labels_rle$lengths[1:x]))))
    section_starts <- c(1, section_ends[-length(section_ends)] + 1)
    grid_marks <- section_starts + shift
    grid_labels <- c(labels_rle$values)
    axis(1, las=3, at=grid_marks, labels=grid_labels, cex.axis=labels_size)
    if(draw_grid_lines == TRUE) {
        #abline(v=grid_marks, col="grey") #Lines for each gene
        plot_custom_gridlines(custom_gridlines, nmat_chr, shift)
    }
}

add_cbs_segments <- function(all.segs, name, nmat_chr) {
    case.segs <- all.segs[(make.names(all.segs[,"Sample"]) == make.names(name)),,drop=FALSE]
    case.segs <- case.segs[(as.character(case.segs$Genome_start_interval) %in% rownames(nmat_chr)),,drop=FALSE] #Filter to chromosome being plotted
    start_indexes <- which(rownames(nmat_chr) %in% case.segs[,"Genome_start_interval"])
    end_indexes   <- which(rownames(nmat_chr) %in% case.segs[,"Genome_end_interval"])
    if(length(start_indexes) > 0) {
        for(i in seq(1:length(start_indexes))) {
            #guideline at CNV median
            segments(start_indexes[i] - 0.25,
                 as.numeric(case.segs[i, "Median_log2ratio"]),
                 end_indexes[i] + 0.25,
                 as.numeric(case.segs[i, "Median_log2ratio"]),
                 col="orange", lwd=3)       
        }                          
    }               
}

replot_cnv_data_points <- function(calls, name, nmat_chr, label_called_CNVs, ylimits) {
    #Extract thresholds and cnvs, determine intervals from genome coordinates
    types     <- matrix(byrow=TRUE, ncol=4,
                        dimnames=list(c(),
                            c("type", "col", "ypos.line", "ypos.label")),
                        data=c(
                              "Loss", "red",    ylimits[1],         ylimits[1]-0.1,
                              "Gain", "blue",	ylimits[2],          ylimits[2]+0.1))
    segs <- calls[(make.names(calls[,"Sample"]) == make.names(name)),,drop=FALSE]
    segs <- segs[(as.character(segs$Genome_start_interval) %in% rownames(nmat_chr)),,drop=FALSE] #Filter to chromosome being plotted
    start_indexes <- which(rownames(nmat_chr) %in% segs[,"Genome_start_interval"])
    end_indexes   <- which(rownames(nmat_chr) %in% segs[,"Genome_end_interval"])
    if(length(start_indexes) > 0) {
        segs <- cbind(segs, types[match(segs[,"CNV"], types[,"type"]),,drop=FALSE], stringsAsFactors=FALSE)
        for(i in seq(1:length(start_indexes))) {
            indexes <- seq(start_indexes[i], end_indexes[i])
            #Plot colored points
            points(indexes, nmat_chr[indexes], col=segs[i, "col"])

            #Label called CNVs with identifier
            if(label_called_CNVs == TRUE) {
                #guideline at CNV median
                segments(start_indexes[i] - 0.25,
                     as.numeric(segs[i, "Median_log2ratio"]),
                     end_indexes[i] + 0.25,
                     as.numeric(segs[i, "Median_log2ratio"]),
                     col="orange", lwd=3)                                 
                 
                #guideline and text at bottom of plot
                segments(start_indexes[i] - 0.25,
                     as.numeric(segs[i, "ypos.line"]),
                     end_indexes[i] + 0.25,
                     as.numeric(segs[i, "ypos.line"]),
                     col="orange", lwd=2)                                 
                text(mean(c(start_indexes[i], end_indexes[i])),
                     as.numeric(segs[i, "ypos.label"]),
                     segs[i, "CNV_id"],
                     col="orange")
            }
        }
    }
}

exon_plot_per_sample <- function(nmat, ylimits, interval_lookup, thresholds, calls, out_dir, cbs.seg, custom_gridlines, name_prefix="") {
    #Reformat cbs seg output for plotting
    if(class(cbs.seg) == "DNAcopy") {
        all.segs <- reformat_cbs_to_viscap(cbs.seg$output, interval_lookup, thresholds, nmat)
    }

    #Limit matrix for plotting
    nmat[nmat < ylimits[1]] <- ylimits[1]
    nmat[nmat > ylimits[2]] <- ylimits[2]

    for(name in colnames(nmat)) {
        #Get exon names from interval lists
        interval_names <- annotate_interval_names(rownames(nmat), interval_lookup)

        #Plot all intervals, then plot by chromosome
        for(output in c("png", "pdf")) {
            #Direct output to appropriate location
            if(output == "png") {
                out_dir.png <- paste(out_dir, "/png/", sep="")
                dir.create(out_dir.png, showWarnings=FALSE)
            } else if(output == "pdf") {
                #Single pdf with pages for each chromosome
                pdf(file=paste(out_dir, "/", name_prefix, name, ".plot.pdf", sep=""))
            }
            chroms <- c("All", 1:22, "X", "Y")
            for(chr in chroms) {
                #Grep rownames with genome coordinates for desired chromosome, support chr1: and 1: formats
                if(chr == "All") {
                    matches <- seq(1, length(rownames(nmat)))
                } else {
                    matches <- grep(paste("(^chr", chr, "|^", chr, "):", sep=""), rownames(nmat))
                }
                if(length(matches) > 1) {
                    #Output separate pngs for each chromosome if requested
                    if(output == "png") {
                        png(file=paste(out_dir.png, name_prefix, name, ".chr", chr, ".plot.png", sep=""),width=960,height=960,pointsize=24)
                    }
                    nmat_chr <- nmat[matches, name, drop=FALSE]
                    #Different title and axes labeling for All vs individual chromosomes
                    if(chr == "All") {
                        title <- "All chromosomes"
                        chrom_names <- sapply(rownames(nmat_chr), function(x) strsplit(x, ":")[[1]][1])
                        labels_rle <- rle(chrom_names)
                        #TODO: Find better method for spacing chromosome labels
                        labels_rle$values <- space_label_names(labels_rle$values)
                        #labels_rle$values <- rep("", length(labels_rle$values-1))
                        labels_size <- 0.6
                        draw_grid_lines <- FALSE
                        label_called_CNVs <- FALSE
                        draw_all_segments <- FALSE
                    } else {
                        title <- paste("Chromosome", chr)
                        #Convert genome coordinates to unique gene names
                        exon_names <- annotate_interval_names(rownames(nmat_chr), interval_lookup)
                        gene_names <- sapply(exon_names, function(x) strsplit(as.character(x), "_|\\.")[[1]][1])
                        
                        #Replace NA gene_names with nearest non-NA value
                        gene_names <- na.locf(gene_names)
                        labels_rle <- rle(gene_names)
                        labels_size <- 1
                        draw_grid_lines <- TRUE
                        label_called_CNVs <- TRUE
                        draw_all_segments <- TRUE
                    }

                    #Plot log2ratios, guidelines, and thresholds
                    plot_viscap_data(nmat_chr, title, name, labels_rle, labels_size, ylimits, draw_grid_lines, custom_gridlines)

                    #Add guidelines and thresholds
                    abline(h=0, col="black")                          #zero line
                    abline(h=c(log2(1/2), log2(3/2)), col="black")     #expected log2 ratios for single copy loss and gain
                    abline(h=thresholds[thresholds[,1] == name, 2:3], lty=2, col="grey")
#                    abline(h=threshold.cnv_log2_cutoffs, col="grey") #hard cnv log2 ratio thresholds
            
                    #Add copy number segments from CBS
                    if(class(cbs.seg) == "DNAcopy" & (draw_all_segments == TRUE)) {
                        add_cbs_segments(all.segs, name, nmat_chr)
                    }
                                                           
                    #Color CNV data points, if calls have been made
                    if(dim(calls)[1] > 0) {
                        replot_cnv_data_points(calls, name, nmat_chr, label_called_CNVs, ylimits)
                    }
                    #Close chromosome-specific png
                    if(output == "png") { dev.off() }
                }
            }
            #Close single pdf with pages for each chromosome
            if(output == "pdf") { dev.off() }
        }
    }
}

rescale_chrX_by_clustering <- function(nmat_badX, ylimits, iqr_multiplier, out_dir) {
    #Isolate chromosome X, remove outlier probes    
    nmatX <- nmat_badX[grep("^(chr)*X:",rownames(nmat_badX)),,drop=FALSE]

    #Remove probes with Inf values
    nmatX <- nmatX[!(rowSums(!is.finite(nmatX)) > 0),,drop=FALSE]
    
    # return original matrix if no probes are found on X chromosome
    if(dim(nmatX)[1] == 0) {
        return(nmat_badX)
    }

    #Assign cases to clusters and scale medians to zero
    #    DEPRECATED METHOD: kmeans is mislead by outliers. multiple solutions possible.    
    #    clusters <- kmeans(t(nmatX),2, nstart=1)
    #    cluster_med <- apply(clusters$centers, 1, median)
    #    scaling_factors <- cluster_med[clusters$cluster]

    #CURRENT METHOD: "Partitioning Around Medoids" followed by calculation of cluster medians
    clusters <- pam(t(nmatX), 2, cluster.only=TRUE) #t() = transpose
    cluster_med <- c()
    for(i in 1:max(clusters)) {
        cluster_med[i] <- median(nmatX[,clusters == i])
    }
    scaling_factors <- cluster_med[clusters]
    nmatX_scaled <- sweep(nmatX, 2, scaling_factors, "-")  ## subtract cluster median

    #Plot pre- and post-scaled values, write to file
    bplotX        <- boxplot_cnv_matrix(nmatX, "QC_chrX_pre-scale", out_dir, ylimits, iqr_multiplier)
    bplotX_scaled <- boxplot_cnv_matrix(nmatX_scaled, "QC_chrX_post-scale", out_dir, ylimits, iqr_multiplier)

    #Infer sexes from two clusters of fractional coverage
    farthest_from_zero <- abs(cluster_med) == max(abs(cluster_med))
    if(cluster_med[farthest_from_zero] > 0) {
        cluster_female <- farthest_from_zero
    } else {
        cluster_female <- !farthest_from_zero
    }
    sexes <- sapply(cluster_female[clusters], function(x) if(x){"Female"} else {"Male"})
    names(sexes) <- names(clusters)
    write.table(sexes, paste(out_dir, "sexes.xls", sep="/") , quote = FALSE, sep = "\t", col.names = FALSE)

    #Overwrite unscaled X-chromosome values with new, scaled values
    nmat <- nmat_badX
    nmat[rownames(nmatX_scaled),] <- nmatX_scaled #Match interval names from nmatX_scaled with original nmat
    
    return(nmat)
}

rescale_chrX_by_known_sex <- function(nmat_badX, ylimits, iqr_multiplier, out_dir, sex.table) {
    #TODO: Write function to scale X-chromosomes +/- 1 copy, depending on sex table
    return(nmat_badX)
}

remove_failed_samples <- function(mat, cnv_bplot_and_calls) {
    passes <- score_boxplot(cnv_bplot_and_calls)
    mat.trimmed <- mat[,passes]
    return(mat.trimmed)
}

call_segments_using_cbs <- function(nmat, out_dir, threshold.min_exons, cbs.sd.undo, cbs.alpha, sam="sample") {
    require("DNAcopy")

    #Setup CNA object using matrix of log2 ratios
    size_parts <- sapply(rownames(nmat), get_size_and_parts_from_coord) #split out chr, start, end, and length
    chrom <- ordered(size_parts[2,], levels=c(1:22,"X","Y")) #ordering the chromosomes
    maploc <- as.integer(size_parts[3,]) #convert object from string to integer

    #calls to DNAcopy library
    cna <- CNA(nmat, chrom, maploc, data.type=c("logratio"), sampleid=colnames(nmat)) 
    #cna.smooth <- smooth.CNA(cna) #smooth CNA data
    seg <- segment(cna, alpha = cbs.alpha, min.width=threshold.min_exons, undo.splits="sdundo", undo.SD=cbs.sd.undo)     #Segment data

    #Write out all segments
    write.table(seg$output, paste(out_dir, "/", sam, ".segments.seg", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

    return(seg)
}

make_genome_coord_names <- function(mat.all) {
    parts <- strsplit(rownames(mat.all), ":|-")
    parts.incomplete <- !(lapply(parts, length) == 3)

}

filter_and_reformat_cbs <- function(seg, interval_lookup, out_dir, thresholds) {
    seg.filtered <- filter_seg(seg, thresholds)
    calls        <- reformat_cbs_to_viscap(seg, interval_lookup, out_dir, thresholds)
    return(calls)
}

number_consecutive <- function(samples) {
    consec <- c()
    runs <- rle(samples)
    for(i in 1:length(runs$values)) {
        consec <- append(consec, seq(1, runs$lengths[i]))    
    }
    return(consec)
}

match_interval <- function(interval_lookup, info, chr, pos) {
    match <- grep(paste("^", chr, ":", as.integer(pos), "-[0-9]+$", sep=""), rownames(interval_lookup))
    if(length(match) == 0) {
        found <- NA
    } else if(info %in% colnames(interval_lookup)) {
        found <- interval_lookup[match, info]
    } else {
        found <- rownames(interval_lookup)[match]
    }
    return(found)
}

get_range_of_values <- function(call, nmat, func=function(x) {x}) {
    index.start <- match(call["Genome_start_interval"], rownames(nmat))
    index.end <- match(call["Genome_end_interval"], rownames(nmat))
    if(NA %in% c(index.start, index.end)) {
        warning("Coverage files contain intervals not found in interval_list.")
        values <- NA
    } else {
        range <- seq(index.start, index.end)
        values <- as.vector(nmat[range, as.character(call["Sample"])])
    }
    value <- func(values)
    return(value)
}

get_genes_between_intervals <- function(interval_lookup, name.field, start.interval, end.interval) {
    start.pos <- match(start.interval, rownames(interval_lookup))
    end.pos   <- match(end.interval, rownames(interval_lookup))
    range <- seq(start.pos, end.pos)
    interval.names <- interval_lookup[range, name.field]
    interval.genes <- unique(sapply(strsplit(as.character(interval.names), "_|\\."), "[", 1))
    interval.genes.string <- paste(interval.genes, collapse=",")
    return(interval.genes.string)
}

reformat_cbs_to_viscap <- function(seg.filtered, interval_lookup, thresholds, nmat) {
    #Reformat cbs data
    calls.required.header <- c("Sample", "CNV_id", "CNV", "Genes",
                               "Genome_start_interval", "Genome_end_interval",
                               "Start_interval","End_interval", "Interval_count",
                               "Min_log2ratio", "Median_log2ratio", "Max_log2ratio",
                               "Loss_threshold", "Gain_threshold", "Batch_size")
    if(dim(seg.filtered)[1] == 0) {
        calls.ordered <- as.data.frame(matrix(ncol=length(calls.required.header), dimnames=list(c(),calls.required.header))[-1,])
    } else {
        #Rename sample names in calls matrix to be consistent with nmat
        calls <- seg.filtered
        calls$ID <- colnames(nmat)[match(calls$ID, make.names(colnames(nmat)))]  

        #Replace "ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean" #assumes median and mean are same
        colnames(calls) <- c("Sample", "Chromosome", "Start", "End", "Interval_count", "Median_log2ratio")

        #Fill in CNV numbers for review and CNV type
        calls <- cbind(calls, CNV_id=number_consecutive(calls$Sample))
        calls <- cbind(calls, CNV=sapply(calls$Median_log2ratio, function (x) if(x<0) {"Loss"} else {"Gain"}))
    
        #Look up interval information
        calls <- cbind(calls, Genome_start_interval   = unlist(apply(calls, 1, function(x) match_interval(interval_lookup, "coordinates",  x[2], x[3]))), stringsAsFactors=FALSE)
        calls <- cbind(calls, Genome_end_interval     = unlist(apply(calls, 1, function(x) match_interval(interval_lookup, "coordinates", x[2], x[4]))), stringsAsFactors=FALSE)
        calls <- cbind(calls, Start_interval          = unlist(apply(calls, 1, function(x) match_interval(interval_lookup, "interval_name", x[2], x[3]))), stringsAsFactors=FALSE)
        calls <- cbind(calls, End_interval            = unlist(apply(calls, 1, function(x) match_interval(interval_lookup, "interval_name", x[2], x[4]))), stringsAsFactors=FALSE)
        calls <- cbind(calls, Genes                   = unlist(apply(calls, 1, function(x) get_genes_between_intervals(interval_lookup, "interval_name", x["Genome_start_interval"], x["Genome_end_interval"]))), stringsAsFactors=FALSE)

        #Get min and max values from nmat
        calls <- cbind(calls, Min_log2ratio=apply(calls, 1, get_range_of_values, nmat, min))
        calls <- cbind(calls, Max_log2ratio=apply(calls, 1, get_range_of_values, nmat, max))             

        #Attach threshold information
        calls <- cbind(calls, Loss_threshold=thresholds[match(calls[,1], thresholds[,1]),2])
        calls <- cbind(calls, Gain_threshold=thresholds[match(calls[,1], thresholds[,1]),3])
        calls <- cbind(calls, Batch_size=dim(nmat)[2])
      
        #Reorder header and restrict to required fields
        calls.ordered <- calls[calls.required.header]
      }
    return(calls.ordered)
}

call_cnvs_from_segments <- function(seg, thresholds) {
    #setup filters as vector of thresholds for each segment
    thresholds[,1] <- make.names(thresholds[,1])
    loss_thresholds <- as.numeric(thresholds[match(seg$output$ID,thresholds[,1]), 2]) #2=lower threshold
    gain_thresholds <- as.numeric(thresholds[match(seg$output$ID,thresholds[,1]), 3]) #3=upper threshold

    #main filtering step that removes all segments inside thresholds
    seg.filtered <- seg$output[(seg$output$seg.mean <= loss_thresholds) |
                                 (seg$output$seg.mean >= gain_thresholds),]
    return(seg.filtered)
}

run_VisCap_algorithm <- function(mat, ylimits, iqr_multiplier, out_dir.iteration,
                                 call.method, threshold.cnv_log2_cutoffs, interval_lookup,
                                 threshold.min_exons, use.boxplot.whisker.threshold,
                                 samples.of.interest, pct.purities.to.threshold, custom_gridlines) {
    # Normalize exon coverage by exon
    nmat_badX <- log2(divide_by_batch_median(mat))
    #nmat <- nmat_badX #TODO: Use provided sexes to normalize X chromosome
    nmat <- rescale_chrX_by_clustering(nmat_badX, ylimits, iqr_multiplier, out_dir.iteration)

    # Write out matrix of log2 ratios to file
    outfile <- paste(out_dir.iteration, "/", "log2_ratio_table", ".xls", sep="")
    gene_exon <- annotate_interval_names(rownames(nmat), interval_lookup)
    nmat.with.gene_exon <- cbind(gene_exon, nmat) #not used for algorithm, output for end-users only
    write.table(nmat.with.gene_exon, outfile, , quote = FALSE, sep = "\t", col.names = NA)

    # Call cnvs then plot heatmaps by chromosome and per-sample exon coverage
    heatmap_by_chrom(nmat, "QC_cnv_heatmap", ylimits, out_dir.iteration)

    # Run segmentation algorithm before iterating on thresholds
    if(call.method == "cbs") {
        #Circular Binary Segmentation
        seg <- call_segments_using_cbs(nmat[,samples.of.interest,drop=FALSE], out_dir.iteration, threshold.min_exons, cbs.sd.undo, cbs.alpha, sam)
    }
    
    #Iteratively apply log2ratio thresholds by purity
    for(purity in pct.purities.to.threshold) { #purities read from config file
        if(purity == as.character("FALSE")) {
            log2.thresholds          <- threshold.cnv_log2_cutoffs
            name_prefix              <- ""
        } else {
            # This applies the purity bins
            log2.thresholds          <- purity_to_log2_thresholds(purity)
            name_prefix              <- paste("purity", purity, "_", sep="_")
        }

        #Use user-defined segmentation and calling method
        if(call.method == "cbs") {
            thresholds <- matrix(ncol=3, data=c(
                            colnames(nmat),
                            rep(min(log2.thresholds), length(colnames(nmat))),
                            rep(max(log2.thresholds), length(colnames(nmat)))
                            ))
            seg.filtered <- call_cnvs_from_segments(seg, thresholds) #reduces to segments that pass sample-specific thresholds
            calls <- reformat_cbs_to_viscap(seg.filtered, interval_lookup, thresholds, nmat)
            passes <- colnames(nmat)
        } else if (call.method == "boxplot") {
            #Distribution thresholding
            cnv_bplot_and_calls <- call_cnvs_using_boxplots(nmat, ylimits, interval_lookup, threshold.min_exons, iqr_multiplier, log2.thresholds, out_dir.iteration, use.boxplot.whisker.threshold)
            thresholds <- t(rbind(cnv_bplot_and_calls[[1]]$names, cnv_bplot_and_calls[[1]]$stats))[,c(1,2,6)]
            calls <- cnv_bplot_and_calls[[2]]
            seg <- FALSE
            passes <- cnv_bplot_and_calls[[1]]$names[cnv_bplot_and_calls[[1]]$qc == "PASS"]        
        } else {
            stop("Config file setting call.method must equal 'cbs' or 'boxplot")
        }

        #Limit to samples of interest for output (used when running against panel of normals)
        nmat.sam  <- nmat[,colnames(nmat) %in% samples.of.interest,drop=FALSE]

        #Write cnv calls for each sample
        for(sam in samples.of.interest) {
            if( (dim(calls)[1] > 0) )  {
                calls.sam <- calls[calls$Sample == sam,,drop=FALSE]
                if( (dim(calls.sam)[1] > 0) ) {
                    #Write table for individual sample
                    write.table(calls.sam[calls.sam$Sample == sam,], paste(out_dir.iteration, "/", name_prefix, sam, ".cnvs.xls", sep=""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
                } else {
                    #If no calls were made in this sample, write-out empty table            
                    write.table(paste(colnames(calls),collapse="\t"),
                                paste(out_dir.iteration, "/", name_prefix, sam, ".cnvs.xls", sep=""),
                                quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)                  
                }
            } else {
                #If no calls were made in any sample, write-out an empty table for this sample
                write.table(paste(colnames(calls),collapse="\t"),
                            paste(out_dir.iteration, "/", name_prefix, sam, ".cnvs.xls", sep=""),
                            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
            }
        }
        
        #Plot data and cnv calls
        if(print_pdf_plots != FALSE) {
            exon_plot_per_sample(nmat.sam, ylimits, interval_lookup, thresholds, calls, out_dir.iteration, seg, custom_gridlines, name_prefix)
        }
    }
    return(passes)
}

write_VisCap_run_info <- function(svn.revision, lane_dir, cov_file_pattern, cov_field, out_dir.iteration, interval_list_dir, interval_file_pattern, nmat, mat.all, ylimits, threshold.min_exons, iqr_multiplier, threshold.cnv_log2_cutoffs, iteration) {
    #Write out VisCap run information
    run_info_table <- matrix(ncol = 2, byrow=TRUE, data = c(
            "Date",                                         date(),
            "VisCap command",                               paste(commandArgs(), collapse=" "),
            "Subversion revision information",              svn.revision,
            "Batch directory",                              lane_dir,      
            "Coverage file pattern",                        cov_file_pattern,
            "Field within coverage file",                   cov_field,              
            "Output directory",                             out_dir.iteration,                
            "Interval name lookup files",                   interval_list_dir,
            "Interval name lookup file pattern",            interval_file_pattern,
            "Exons used for CNV detection",                 dim(nmat)[1],
            "Samples used for CNV detection",               dim(nmat)[2],
            "Samples not used for CNV detection",           paste(c(colnames(mat.all)[!(colnames(mat.all) %in% colnames(nmat))], ""), sep=",", collapse=","),
            "Plot y-axis limits",                           paste(ylimits, collapse=","),
            "Minimum consecutive exons to call CNV",        threshold.min_exons,
            "IQR multiplier used for boxplots",             iqr_multiplier,
            "Static log2 ratio thresholds to call CNVs",    paste(threshold.cnv_log2_cutoffs, collapse=","),
            "Iteration",                                    iteration
            ))
    run_info_outfile <- paste(out_dir.iteration, "/", "VisCap_run_info", ".xls", sep="")
    write.table(run_info_table, run_info_outfile, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    #save.image(file=paste(out_dir.iteration, "/", "session", ".Rdata", sep=""))
}

###########
#Arguments#
###########
usage <- "Usage: VisCap.R config_file lane_directory output_directory\n"

#Argument collection and parsing
arguments <- commandArgs(trailingOnly = TRUE)

if(length(arguments) == 0) {  
    if(dev_dir == FALSE) {
        windows.attempt <- try(arguments <- collect_arguments_from_user(lane_dir.prompt, out_dir.prompt, explorer_file, cov_file_pattern, clobber.output.directory), silent=TRUE)
        if(attr(windows.attempt, "class") == "try-error") {
            #Usage statement
            stop(usage)
        } else {
            user.arguments <- collect_arguments_from_user(lane_dir.prompt, out_dir.prompt, explorer_file, cov_file_pattern, clobber.output.directory)
            config.file <- sub("\\.R$", "\\.cfg", sub("--file=","", grep("--file", commandArgs(), value=TRUE))) #use config file next to VisCap.R file
            lane_dir <- user.arguments[1]
            out_dir  <- user.arguments[2]
            batch    <- user.arguments[3]
        }
    } else {
        #Skips prompts if dev_dir is set
        lane_dir <- dev_dir
        out_dir  <- dev_dir
        if(!exists("config.file")) {
            config.file <- sub("\\.R$", "\\.cfg", sub("--file=","", grep("--file", commandArgs(), value=TRUE))) #use config file next to VisCap.R file
        }
        batch    <- "dev"
    }   
} else if(length(arguments) == 3) {
    #Uses provided command line arguments
    config.file <- arguments[1]    
    lane_dir    <- arguments[2]
    out_dir     <- arguments[3]
    batch       <- basename(out_dir)
} else {
    #Usage statement
    cat(prompt=usage)
    q(save="no")
}

#Load config file
if(file.exists(config.file)) {
    source(config.file)
} else {
    stop(paste("Config file not found:", config.file))
}

######
#Main#
######

# Make output directory
dir.create(out_dir, showWarnings = FALSE, recursive=TRUE)

# Read coverage tables
mat.batch <- make_matrix_from_cov_files(lane_dir, cov_file_pattern, cov_field, out_dir, cov.plot.range, "QC_coverage_project_samples.pdf")
if(panel.of.normals != FALSE) {
    #Read Panel of Normals and add suffix to normals with names identical to a sample in mat.batch
    mat.normals <- make_matrix_from_cov_files(panel.of.normals, cov_file_pattern, cov_field, out_dir, cov.plot.range, "QC_coverage_panel_of_normals.pdf")    
    colnames(mat.normals)[colnames(mat.normals) %in% colnames(mat.batch)] <- paste(colnames(mat.normals)[colnames(mat.normals) %in% colnames(mat.batch)], "present_in_batch", sep="__")

    #Merge matricies    
    mat.all <- merge(mat.batch, mat.normals, by = "row.names")
    rownames(mat.all) <- mat.all[,"Row.names"]
    mat.all <- as.matrix(mat.all[,2:dim(mat.all)[2]])
} else {
    mat.normals <- FALSE
    mat.all <- mat.batch
}

# Remove intervals without consistent genome format (e.g. 1:2345-6789, chr1:2345-6789)
mat.all <- mat.all[grep(".+:[0-9]+-[0-9]+", rownames(mat.all)),]

# Sort matrix by genome coordinates found in rownames
chroms <- c(1:22,"X","Y", "MT", "M")
chroms <- factor(chroms, levels=chroms, labels=chroms, ordered=TRUE)
coords <- matrix(unlist(strsplit(rownames(mat.all), ":|-")), ncol=3, byrow=TRUE, dimnames=list(rownames(mat.all)))
coords <- coords[order(match(coords[,1], chroms), as.numeric(coords[,2]), as.numeric(coords[,3])),]
mat <- mat.all[rownames(coords),]


# Read interval name files & check that coverage matrix and interval lookup file are consistent
interval_lookup <- load_interval_names(interval_list_dir, interval_file_pattern, interval_list_coord_adjust)

# Bug Fix: If interval_lookup table is not sorted, the genes in CNA regions are identified incorrectly.
# Sort interval_lookup table by chromosome then start position
interval_lookup <- interval_lookup[order(interval_lookup$chr,interval_lookup$start),]

if(length(rownames(mat)) != sum(rownames(mat) %in% rownames(interval_lookup))) {
	warning("Coverage files contain coordinate ranges not found in interval_list. There may be overlapping intervals that were merged by GATK.")
}


# Read custom gridlines file (table downloaded from UCSC Genome Browser website)
custom_gridlines <- read.table(custom_gridlines_file, header=TRUE, comment.char="", stringsAsFactors=FALSE)

#Iteratively run VisCap algorithm, removing bad samples after each run
if(iterative.calling.limit == 0) { iterative.calling.limit <- dim(mat)[2] }
for(iteration in 1:iterative.calling.limit) {
    if(iterative.calling.limit == 1) {
        out_dir.iteration <- out_dir
    } else {
        out_dir.iteration <- paste(out_dir, "/", batch, "_run", iteration, sep="")
    }
   dir.create(out_dir.iteration, showWarnings = FALSE, recursive=TRUE)
        
    #Two modes: Batch-level normalization or Samples-against-panel-of-normals
    if(run.against.normals.only == FALSE) {
        #Batch-level normalization
        passes <- run_VisCap_algorithm(mat, ylimits, iqr_multiplier, out_dir.iteration,
                                       call.method, threshold.cnv_log2_cutoffs,
                                       interval_lookup, threshold.min_exons,
                                       use.boxplot.whisker.threshold, colnames(mat),
                                       pct.purities.to.threshold, custom_gridlines)
        write_VisCap_run_info(svn.revision, lane_dir, cov_file_pattern, cov_field,
                              out_dir.iteration, interval_list_dir, interval_file_pattern,
                              mat, mat.all, ylimits, threshold.min_exons, iqr_multiplier,
                              threshold.cnv_log2_cutoffs, iteration)
    } else {
        #Samples-against-panel-of-normals
        for(sam in colnames(mat.batch)) {
            mat.sam <- mat[,c(sam, colnames(mat.normals))]
            out_dir.iteration.sam <- paste(out_dir.iteration, sam, sep="/")
            dir.create(out_dir.iteration.sam, showWarnings = FALSE, recursive=TRUE)
            passes.sam <- run_VisCap_algorithm(mat.sam, ylimits, iqr_multiplier, out_dir.iteration.sam,
                                           call.method, threshold.cnv_log2_cutoffs,
                                           interval_lookup, threshold.min_exons,
                                           use.boxplot.whisker.threshold, sam, 
                                           pct.purities.to.threshold, custom_gridlines)
            write_VisCap_run_info(svn.revision, lane_dir, cov_file_pattern, cov_field,
                                  out_dir.iteration.sam, interval_list_dir, interval_file_pattern,
                                  mat.sam, mat.all, ylimits, threshold.min_exons, iqr_multiplier,
                                  threshold.cnv_log2_cutoffs, iteration)
        #TODO: Write copy number segments for all samples to a file
        }
        #In this mode, pass every sample
        passes <- colnames(mat)
    }
    #Remove failed samples from matrix for next run
    if(length(passes) == dim(mat)[2]) {
        break
    } else {
        #Restrict mat only to samples that pass boxplot qc
        mat <- mat[,passes]
    }
}

# Open output directory and quit R
if(dev_dir == FALSE) {
    if(length(arguments) == 0) {
        shell(paste(explorer_file, out_dir, sep=" "), wait=FALSE)
    }
    quit(save="no")
}
