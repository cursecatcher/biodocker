isoform <- function(data.folder, isoformdata_file) {
	cat("Obtaining circRNA main isoform\n")

	circRNA_classification <- file.path(data.folder, "circRNA_classification")
	#prefix temporary files
	temp_prefix <- file.path(data.folder, "CHMAINISO")

	cleanup <- paste0("rm -f ", temp_prefix, "*")
	system(cleanup, wait=TRUE)

	# ===============================
	# COMMANDS SECTION
	# ===============================

	# Filtering for canonical chromosomes
	canonical_isoforms <- paste0(temp_prefix, "_genes_isoform_canonical")

	filter <- paste("awk -v OFS='\\t' '{if(",
		" $3 == \"chr1\" || $3 == \"chr2\" || $3 == \"chr3\" ||",
		" $3 == \"chr4\" || $3 == \"chr5\" || $3 == \"chr6\" ||",
		" $3 == \"chr7\" || $3 == \"chr8\" || $3 == \"chr9\" ||",
		" $3 == \"chr10\" || $3 == \"chr11\" || $3 == \"chr12\" ||",
		" $3 == \"chr13\" || $3 == \"chr14\" || $3 == \"chr15\" ||",
		" $3 == \"chr16\" || $3 == \"chr17\" || $3 == \"chr18\" ||",
		" $3 == \"chr19\" || $3 == \"chr20\" || $3 == \"chr21\" ||",
		" $3 == \"chr22\" || $3 == \"chrX\" || $3 == \"chrY\")",
		"{print $0}}' ",
		isoformdata_file, ">", 			#input file
		canonical_isoforms				#output file
	)

	# Obtaining isoform value
	isoform_value <- paste0(temp_prefix, "_isoform_value")
	getIsoformValue <- paste(
		"cut -f 4", canonical_isoforms, 	#input file
		"| rev | cut -c1-3 | rev", 			#obtaining isoform value
		">", isoform_value					#output file with isoform values
	)

	# Merge gene isoform and isoform values files
	unsorted_isoform_values <- paste0(temp_prefix, "_genes_isoforms_value_unsorted")
	assignGeneIsoform <- paste("paste", canonical_isoforms, isoform_value, ">", unsorted_isoform_values)

	sorted_genes <- paste0(temp_prefix, "_01_genes_isoforms_value")
	sortGenes <- paste("sort -k 4,4 -k 5,5", unsorted_isoform_values, ">", sorted_genes)

	# Obtaining transcripts from classification
	transcripts_temp <- paste0(temp_prefix, "_ENST")

	getTranscripts <- paste(
		"cut -f 1", circRNA_classification, 	#input data
		"| awk 'BEGIN {FS=\"_\"} {print $4}'",	# command
		">", transcripts_temp					#output file
	)

	# Obtaining isoform data from ENST
	transcripts_def <- paste0(temp_prefix, "_ENST_ISO")

	getIsoformData <- paste(
		"python data_extractor.py",  # used program
		sorted_genes, 				#input data file
		transcripts_temp, 			#input query file
		"1 0", 						#column parameters
		"| cut -f 4,5",				#extracted data
		">", transcripts_def		#output file
	)

	# Pasting results
	circ_isoforms <- paste0(temp_prefix, "_circ_iso")
	pasteIsoCircData <- paste("paste",
		transcripts_def, 			# input enst-iso file
		circRNA_classification, 	# input classification file
		">", circ_isoforms			#output file
	)

	# Obtaining unique classification
	univocal_temp <- paste0(temp_prefix, "_circRNA_univocal_classification")
	uniqueClassification <- paste("python univocal_classifier.py", circ_isoforms, ">", univocal_temp)

	univocal_unsorted <- paste0(temp_prefix, "_unsorted_circRNA_univocal_classification")
	postClassification <- paste("python post_univocal.py", univocal_temp, ">", univocal_unsorted)


	univocal_def <- file.path(data.folder, "circRNA_univocal_classification")
	resultSort <- paste("sort", univocal_unsorted, ">", univocal_def)

	# ===============================
	# PIPELINE SECTION
	# ===============================

	# Building pipeline
	command_list <- c(
		filter,   				# Filtering for canonical chromosomes
		getIsoformValue,        # Obtaining isoform value
		assignGeneIsoform,      # Obtaining isoform value
		sortGenes,              # Obtaining isoform value
		getTranscripts,         # Obtaining transcripts from classification
		getIsoformData,         # Obtaining isoform data from ENST
		pasteIsoCircData,       # Pasting results
		uniqueClassification,   # Obtaining unique classification
		postClassification,     # Post-classification final step
		resultSort,             # Result sorting
		cleanup                # Final CHRMAINISO temp file cleanup
	)
	# ===============================
	# EXECUTION SECTION
	# ===============================

	for (command in command_list) {
		system(command, wait=TRUE)
	}

    #set headers
    header <- paste("ID", "Class", "Isoform_ID", "Isoform_Rank", "Ensembl_Gene_ID_1",
        "Ensembl_Trnascript_ID_1", "Ensembl_Exon_ID_1", "Exon_rank_1", "strand_1",
        "Ensembl_Gene_ID_2", "Ensembl_Trnascript_ID_2", "Ensembl_Exon_ID_2",
        "Exon_rank_2",	"strand_2", sep="\t")

    for (file_res in c(univocal_def, circRNA_classification)) {
        header.command <- paste("sed -i", paste0("'1i", header, "'"), file_res)
        system(header.command, wait=TRUE)
    }
}
