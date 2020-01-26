setwd("/scripts/circhunter")

source("backsplicing_junction_sequences.R")
source("circRNA_classification.R")
source("get_exon_data.R")
source("isoform.R")
source("isoform_info.R")
source("r_hashcirc.R")

main_scratch <- "/scratch"

#create a scratch folder for the current execution
tmp.id <- gsub("[:\ ]", "-", date()) 
my_scratch <- file.path(main_scratch, tmp.id)
writeLines(paste("Current scratch folder ID: ", tmp.id))
dir.create(file.path(my_scratch))

## Useful paths
data.folder <- "/data"
path_exondata <- file.path(data.folder, "genome")
path_isoformdata <- file.path(data.folder, "isoformdata")
path_circrna <- file.path(data.folder, "circRNA")
path_bksjunctions <- file.path(data.folder, "bksj")
path_rnaseq <- file.path(data.folder, "rnaseq")
#path_output_folder <- data.folder
path_output_folder <- file.path(data.folder, "out")

#my_files <- c(path_exondata, path_isoformdata, path_isoformdata, path_circrna, path_bksjunctions, path_rnaseq)

## Now main() begin
writeLines("* Inside R")
start <- Sys.time()
args <- commandArgs(TRUE)


# Checking required genome assembly and assessing validity
assembly <- "hg19"
if ("-as" %in% (args) | "--assembly" %in% (args)) {
    argtouse <- ifelse("-as" %in% args, match("-as", args), match("--assembly", args))
	assembly <- args[argtouse + 1]
}
if (assembly %in% c("hg18", "hg19", "hg38", "mm9", "mm10", "ce11", "dm6", "rn6") == FALSE) {
    writeLines("Please enter a valid genome assembly (hg18, hg19, hg38, mm9, mm10, ce11, dm6, rn6).")
	quit(save="no")
}

# Checking required Ensembl version and assessing validity
version <- "NULL"
if ("-v" %in% (args) | "--version" %in% (args)) {
    argtouse <- ifelse("-v" %in% args, match("-v", args), match("--version", args))
	version <- args[argtouse + 1]
}
if (version %in% c(67,75,77, 79:98) || version > 98 == FALSE) {
    writeLines("Please enter a valid Ensembl version.")
	quit(save="no")
}

## Download isoform and exon data
if ("-f" %in% (args) | "--full" %in% (args) |
    "-p" %in% (args) | "--preparedata" %in% (args)) {
    writeLines("Preparing isoform and exon data.")

    writeLines("Downloading exon data...")
    get_exon_data(my_scratch, assembly, version)

    writeLines("Downloading isoform data...")
    get_isoform_data(my_scratch, assembly, version)

    writeLines("Downloads completed.")
}


## CircRNA classification
if ("-f" %in% (args) | "--full" %in% (args) |
    "-c" %in% (args) | "--classification" %in% (args)) {

    writeLines("circRNA classification in progress...")

    # Checks if a biomart export was supplied with -sg
    if (("-sg" %in% (args) | "--suppliedgenome" %in% (args)) == FALSE) {
        writeLines("Downloading exon data...")
        get_exon_data(my_scratch, assembly, version)
    }
    circRNA_classification(my_scratch, path_exondata, path_circrna)
	writeLines("-- Classification completed --")

    #=================================================
    # CODE FOR MAIN ISOFORM
    #=================================================

    writeLines("Main isoform circRNA classification in progress")
    if(("-id" %in% (args) | "--isoformdata" %in% (args)) == FALSE) {
        writeLines("Downloading isoform data\n")
        get_isoform_data(my_scratch, assembly, version)
    }
    isoform(my_scratch, path_isoformdata)
    writeLines("-- Main isoform classification obtained --")
}

# Obtains CircRNA backsplicing junctions
if ("-f" %in% (args) | "--full" %in%(args) |
    "-s" %in% (args) | "--sequences" %in% (args)) {

	writeLines("Searching backsplicing junction sequences coordinates")

    backsplicing_junction_sequences(my_scratch, path_circrna)
	writeLines("-- Sequences obtained --")
}

if ("-r" %in% (args) | "--readcount" %in% (args)) { # -f was excluded on purpose
    argtouse <- ifelse("-hc" %in% args, match("-hc", args), match("--hashcirc", args))

    hcargs <- c(20, 8, 1000001, 1000001, 17, 1) #default values: very arbitrary ones
    if (!is.na(argtouse)) {
        hcargs <- args[argtouse + 1:6]
    }

	writeLines("Running HashCirc analysis")

    argtouse <- match("--samplename", args)
    samplename <- ifelse(is.na(argtouse), "Readcount", args[argtouse+1])

    r_hashcirc(my_scratch, path_rnaseq, path_bksjunctions, hcargs, samplename)
	writeLines("-- HashCirc analysis completed --")
}

if ("-q" %in% args) {
    writeLines("Converting FASTA to FASTQ")

    fasta_file <- paste(my_scratch, "my_fasta.fastq", sep="/")
    conversion <- paste("python fasta_to_fastq.py", path_circrna, ">", fasta_file)
    system(conversion, wait=TRUE)

    writeLines("-- Conversion completed --")
}


# Output folder step
if ("-of" %in% (args) | "--outputfolder" %in% (args)) {

	writeLines("Placing output files in provided folder")
    data_in_scratch <- file.path(my_scratch, "*")
    move_output <- paste("rsync -avhP", data_in_scratch, path_output_folder)
    system(move_output, wait=TRUE)
    #chmodding everything with the number of the beast
    system(paste("chmod 666", file.path(path_output_folder, "*")), wait=TRUE)
    writeLines("-- Output moved --")

	system(paste("rm -f", data_in_scratch))

    writeLines(paste("Removing temporary directory", my_scratch))
    system(paste("rmdir", my_scratch))  #remove directory if it is empty
}

writeLines("Exiting CircHunter successfully.")
writeLines(paste("Execution time:", Sys.time()-start))
quit(save="no")
