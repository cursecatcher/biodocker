# Welcome to CircHunter
# This script allows to obtain backsplicing junction sequences from a circRNA list

source("backsplicing_seq.R")

backsplicing_junction_sequences <- function(data.folder, circrna.file) {
	circ <- circrna.file

	con <- file(circ, "r")
	first_line <- read.table(con, nrows=1, sep='\t')
	close(con)

	circtoseq_file <- file.path(data.folder, "circtoseq")

	# Cheching if strand conversion is required and copying file
	if (first_line$V5 == '1' | first_line$V5 == '-1'){

		cat_command <- paste("cat", circ)
		awk_command <- "awk -v OFS='\t' '{if ($5 == 1) {$5 = \"+\"} else if ($5 == -1) {$5 = \"-\"}; print }'"
		cut_command <- paste("cut -f 4,5 > ",  circtoseq_file)

		conversion <- paste(cat_command, awk_command, cut_command, sep=" | ")
		system(conversion, wait=TRUE)

	} else {
		copy_cf <- paste("cut -f 4,5", circ, ">", circtoseq_file)
		system(copy_cf, wait=TRUE)
	}

	# Obtaining backsplicing junction sequences
	backsplicing_coord <- file.path(data.folder, "backsplicing_coord")

	backsplicing_command <- paste("python 70bp_junction_coordinates.py", circtoseq_file, ">", backsplicing_coord)
	system(backsplicing_command, wait=TRUE)

	cleanup <- paste("rm", circtoseq_file)
	system(cleanup, wait=TRUE)

	# Header
	create_header <- paste("circRNA_NAME", "chromosome", "strand", "ex5_start", "ex5_end", "ex3_start", "ex3_end", sep='\t')
	header_file <- file.path(data.folder, "header_backsplicing")
	out_file <- file(header_file)

	writeLines(create_header, out_file)
	close(out_file)


	# Creating csv file containing backsplicing junction sequences
	backsplicing_csv <- file.path(data.folder, "backsplicing_coord.csv")
	cat_command <- paste("cat", header_file, backsplicing_coord)
	cut_command <- paste("cut -f 1,2,3,4,5,6,7", backsplicing_csv, sep=" > ")
	backsplicing_coordinates <- paste(cat_command, cut_command, sep=" | ")

	system(backsplicing_coordinates, wait=TRUE)

	cleanup <- paste("rm", header_file, backsplicing_coord)
	system(cleanup, wait=TRUE)

	# Obtaining backsplicing junction sequences
	backsplicing_seq(data.folder, assembly)

	cleanup <- paste("rm", backsplicing_csv)
	system(cleanup, wait=TRUE)
}
