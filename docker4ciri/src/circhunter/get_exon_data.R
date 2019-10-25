# Welcome to CircHunter
# This script extracts the required exon data from ensembl

library(biomaRt)

get_genome <- function(assembly, version) {
	# Choosing the right dataset
	if (assembly == "hg18"){
		my_genome <- useMart(
			"ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
			host="may2009.archive.ensembl.org", path="/biomart/martservice",
			archive=FALSE
		)
	} else if (assembly == "hg19") {
		my_genome <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37, version=version)
	} else if (assembly == "hg38") {
		my_genome <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=version)
	} else if (assembly == "mm9") {
		my_genome <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=67)		
	} else if (assembly == "mm10") {
		my_genome <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=version)
	} else if (assembly == "ce11") {
		my_genome <- useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl", version=version)
	} else if (assembly == "dm6") {
		my_genome <- useEnsembl(biomart="ensembl", dataset="dmelanogaster_gene_ensembl", version=version)
	} else if (assembly == "rn6") {
		my_genome <- useEnsembl(biomart="ensembl", dataset="rnorvegicus_gene_ensembl", version=version)
	}

	return (my_genome)
}

get_exon_data <- function(data.folder, assembly, version) {

	my_genome <- get_genome(assembly, version)

	my_attributes <- c(
		"ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id",
		"chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "rank",
		"start_position", "end_position", "transcript_start", "transcript_end")
	# Extracting data
	exons <- getBM(attributes=my_attributes, mart=my_genome)

	# Correcting chromosome names
	exons$chromosome_name <- paste0("chr", exons$chromosome_name)

	destination_path = file.path(data.folder, "genome")

	# Saving results to exons file
	write.table(exons, file=destination_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

	rm(exons)
}
