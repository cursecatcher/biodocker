# FUNCTION: converter
## Converts FASTA files into FASTQ files
## Renames the FASTQ files and places them in the working directory
converter <- function(element) {
	# Variables for file renaming
	tokens <- strsplit(element, "/")[[1]]
	newname <- paste("CHTEMP", tokens[length(tokens)], ".fastq", sep="")

	# Check first character of file
	check <- file(element, "r")
	first_line <- readLines(check, n=1)

	#check if conversion is needed
	if (substr(first_line,0,1) == ">"){
		# Mounted FASTA files are converted into FASTQ files and renamed
		cat ("FASTA file", element, "will be converted to FASTQ\n", sep=" ")
		command <- paste("python /scripts/circhunter/fasta_to_fastq.py", element, ">", newname)

	} else if ((substr(first_line,0,1) == "@")) {
		# Mounted FASTQ files are renamed creating a symbolic link in the working directory
		cat ("FASTQ file", element, "will be only renamed\n", sep=" ")
		command <- paste("ln -s", element, newname)

	} else {
		cat ("ERROR: Please provide a valid FASTA or FASTQ file\n")
		quit(save="no")
	}
	system(command, wait=TRUE)
}

# FUNCTION: splitter
## Uses fastqsplitter to split FASTQ files
## Renames the resulting FASTQ files for hashcirc
splitter <- function(data.folder, num_parts) {
	cat("FASTQ files will be split into", num_parts, "parts.\n")

	target_files <- file.path(data.folder, "CHTEMP*")
	splitcommand <- paste("/scripts/circhunter/fastq-splitter.pl -n", num_parts, target_files)
	system(splitcommand, wait=TRUE)
}

# FUNCTION: renamer
## Renames FASTQ files in a suitable way for HashCirc execution
renamer <- function(data.folder, oldname) {
	tokens <- strsplit(oldname, "/")[[1]]
	name <- tokens[length(tokens)]				#obtains file name (without path)
	nameparts <- strsplit(name, "\\.")[[1]]		#splits file name into components

	index <- as.integer(strsplit(nameparts[2], "-")[[1]][2]) - 1
	newname <- paste(data.folder, "/", nameparts[1], index, ".", nameparts[3], sep="")

	renamingcommand <- paste("mv", oldname, newname)
	system(renamingcommand, wait=TRUE)
}



r_hashcirc <- function(data.folder, rnaseq, junctionseqs, hcargs, samplename) {
	# =======================================
	# PRELIMINARY OPERATIONS
	# =======================================

	#prefix temporary files
	temp_prefix <- "CHTEMP"
	all_temps <- paste(temp_prefix, "*")
	prev_wd <- getwd()

	cat("Inside r_hashcirc.R\n")
	setwd(data.folder)

	# Make variable with both filenames
	filecheck <- c(rnaseq, junctionseqs)
	#Purges temporary files
	system("rm -f CHTEMP* outputhashcirc* readcount*", wait=TRUE)

	# =======================================
	# MAIN PROGRAM
	# =======================================

	# Conversion step
	for (providedfile in filecheck){
		converter(providedfile)
	}
	# Backup of unsplitted backsplicing junction sequences file
	backup_bksj <- "ALLbksj.fastq"
	temp_bjsk <- paste(temp_prefix, "bksj.fastq", sep="")
	bkp <- paste("cp", temp_bjsk, backup_bksj)
	system(bkp, wait=TRUE)

	# Acquiring data for split and renaming step
	convertedfiles <- system("ls CHTEMP*", intern=TRUE)	# Get list of FASTQ files
	threadnumber <- as.integer(hcargs[2])				# Get thread number

	if (threadnumber < 1){
		cat("ERROR: the 'thread number' argument must be greater than 1.\nQuitting program.\n")
		quit(save="no")
	}
	# Split FASTQ files to perform parallel execution
	splitter(data.folder, threadnumber)

	# Remove unsplitted files
	cleanconv <- paste("rm", paste(convertedfiles, collapse=" "))
	system(cleanconv, wait=TRUE)

	# Renaming splitted files
	for (filename in system("ls CHTEMP*", intern=TRUE)) {
		renamer(data.folder, filename)
	}

	# HashCirc execution step 1
	hashcircexecution <- paste(
		"HashCheckerFilter",		# HashCirc program (first step)
		"CHTEMPbksj fastq",					# Backsplicing junction sequences files
		threadnumber,						# File parts (supplied thread nmumber)
		"CHTEMPrnaseq fastq",				# RNA-Seq files
		threadnumber,						# File parts (supplied thread number)
		"outputhashcirc",					# HashCirc output file
		paste(hcargs[1:5], collapse=" ")	# HashCirc arguments
	)
	print(hashcircexecution)
	system(hashcircexecution, wait=TRUE)

	# HashCirc execution step 2
	for (i in 0:(threadnumber-1)) {
		step1_output <- paste("outputhashcircSim", i, sep="")
		current_output <- paste("readcount", i, sep="")

		hashcircexecution_second <- paste(
			"AlignmentCirRNA",	# HashCirc program (second step)
			"ALLbksj.fastq",				# Backsplicing junction sequences file
			step1_output,					# Output file of step 1
			hcargs[6],						# Matches argument
			current_output					# Output file
		)
		cat("HashCirc second step will run with the following command:\n", hashcircexecution_second, "\n")
		system(hashcircexecution_second, wait=TRUE)
	}

	# Merge output files
	countfile_command <- paste(
		"cd", data.folder, ";",
		"for i in readcount* ; do cut -f 2 $i > $i.num ; done ;",
		"paste *.num | awk '{c=0;for(i=1;i<=NF;++i){c+=$i}; print c}' > count ;"
	)
    output.hashcirc <- paste0(samplename, ".hashcirc")
	merge_all <- paste(
		"cut -f 1 readcount1 > id ;",
        "paste id count >", output.hashcirc, ";",
        "sed -i '1iID\tCounts'", output.hashcirc, ";",
		"rm -f *.num *temp readcount* id count ;",
		"chmod 666 *"
	)
	merging <- paste(countfile_command, merge_all)
	system(merging, wait=TRUE)
	# Cleanup of temporary files
	final_cleanup <- "rm -f CHTEMP* outputhashcirc* ALL*"
	system(final_cleanup, wait=TRUE)
	#set previous working directory
	setwd(prev_wd)
}
