# Welcome to CircHunter
# This script extracts the required exon data from ensembl

library(biomaRt)

get_genome <- function(assembly) {
	# Choosing the right dataset
	if (assembly == "hg18"){
		my_genome <- useMart(
			"ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
			host="may2009.archive.ensembl.org", path="/biomart/martservice",
			archive=FALSE
		)
	} else if (assembly == "hg19") {
		my_genome <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
	} else if (assembly == "hg38"){
		my_genome <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	}

	return (my_genome)
}

get_exon_data <- function(data.folder, assembly) {

	my_genome <- get_genome(assembly)

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
