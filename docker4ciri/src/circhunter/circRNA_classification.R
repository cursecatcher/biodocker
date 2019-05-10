# Welcome to CircHunter

circRNA_classification <- function(data.folder, genome, circrna_file) {
    #paths where exon data and circRNA files are mounted
    exon_data <- file.path(data.folder, "exon_data")

    ### Filtering exon data in order to obtain only exons from canonical chromosomes
    awk_condition <- "awk -v OFS='\t' '{if($4 == \"chr1\" || $4 == \"chr2\" || $4 == \"chr3\" || $4 == \"chr4\" || $4 == \"chr5\" || $4 == \"chr6\" || $4 == \"chr7\" || $4 == \"chr8\" || $4 == \"chr9\" || $4 == \"chr10\" || $4 == \"chr11\" || $4 == \"chr12\" || $4 == \"chr13\" || $4 == \"chr14\" || $4 == \"chr15\" || $4 == \"chr16\" || $4 == \"chr17\" || $4 == \"chr18\" || $4 == \"chr19\" || $4 == \"chr20\" || $4 == \"chr21\" || $4 == \"chr22\" || $4 == \"chrX\" || $4 == \"chrY\") {print $0}}' "
    filtering <- paste(awk_condition, genome, ">", exon_data)
    system(filtering, wait=TRUE)

    ### Obtaining intron data from the available exon data
    # Sorting exons
    exon_sorted <- file.path(data.folder, "exon_sorted")
    sorting <- paste("sort -k 4,4 -k 1,1 -k 2,2 -k8,8 -V", exon_data, ">", exon_sorted)
    system(sorting, wait=TRUE)

    # Obtaining intron information
    full_genome <- file.path(data.folder, "full_genome")
    introns <- paste("python intron_exporter.py", exon_sorted, ">", full_genome)
    system(introns, wait=TRUE)

    ### Creating BED files from circRNA and exon data files
    # Genome
    genome_tofix <- file.path(data.folder, "genome_tofix.bed")
    awk_command <- paste("awk -v OFS='\t' '{print $4,$5,$6,$2\"_\"$3,\".\",$7}'", full_genome)
    extract_columns <- paste(awk_command, genome_tofix, sep=" > ")
    system(extract_columns, wait=TRUE)


    genome_bed <- file.path(data.folder, "genome.bed")
    cat_command <- paste("cat", genome_tofix)
    awk_command <- paste("awk -v OFS='\t' '{if ($6 == 1) {$6 = \"+\"} else if ($6 == -1) {$6 = \"-\"}; print }'", genome_bed, sep=" > ")
    fix_strand = paste(cat_command, awk_command, sep=" | ")
    system(fix_strand, wait=TRUE)

    cleanup <- paste("rm", genome_tofix)
    system(cleanup, wait=TRUE)

    # circRNAs
    circrna_tofix <- file.path(data.folder, "circRNA_tofix.bed")
    extract_circRNA_cols <- paste(
        "awk -v OFS='\t' '{print $1,$2,$3,$4,\".\",$5}'",
        circrna_file, ">", circrna_tofix
    )
    system(extract_circRNA_cols, wait=TRUE)

    circrna_bed <- file.path(data.folder, "circRNA.bed")
    cat_command <- paste("cat", circrna_tofix)
    awk_command <- paste(
        "awk -v OFS='\t' '{if ($6 == 1) {$6 = \"+\"} else if ($6 == -1) {$6 = \"-\"}; print }'",
        ">", circrna_bed
    )
    fix_circRNA_strand <- paste(cat_command, awk_command, sep=" | ")
    system(fix_circRNA_strand, wait=TRUE)

    cleanup <- paste("rm", circrna_tofix)
    system(cleanup, wait=TRUE)

    ### Overlapping between circRNAs and genomic features using bedtools
    overlap_tofix <- file.path(data.folder, "overlap_tofix")
    overlap <- paste("bedtools intersect -a", genome_bed, "-b", circrna_bed, "-wa -wb -s >", overlap_tofix)
    system(overlap, wait=TRUE)

    #Removing score columns
    overlap_file <- file.path(data.folder, "overlap")
    fix <- paste(
        "cat", overlap_tofix, "|",
        "awk -v OFS='\t' '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$12}'", ">", overlap_file
    )
    system(fix, wait=TRUE)

    cleanup <- paste("rm", overlap_tofix)
    system(cleanup, wait=TRUE)

    ### CircRNA classification using circ_classifier.py
    classification_file <- file.path(data.folder, "circRNA_classification")
    classification_command <- paste("python circ_classifier.py", full_genome, overlap_file, ">", classification_file)
    system(classification_command, wait=TRUE)


    # circRNAs from classification
    names <- file.path(data.folder, "names")
    fullnames <- file.path(data.folder, "fullnames")
    missed <- file.path(data.folder, "missed")
    total <- file.path(data.folder, "total")

    addmissing <- paste(
        "cat", classification_file,
        "| awk 'BEGIN {FS=\"_|\t\"} {print $1\"_\"$2\"_\"$3}' | sort -u >", names
    )
    system(addmissing, wait=TRUE)

    # all circRNAs
    addmissing <- paste("cut -f 4", circrna_file, "| sort -u >", fullnames)
    system(addmissing, wait=TRUE)

    # Identification of circRNAs missed by bedtools
    addmissing <- paste(
        "comm -3", fullnames, names, "| awk 'BEGIN {FS=\"\t\"} {print $0\"_NA\tintergenic\"}' >", missed
    )
    system(addmissing, wait=TRUE)

    # Adding circRNAs to classification
    addmissing <- paste("cat", classification_file, missed, ">", total)
    system(addmissing, wait=TRUE)

    # Cleaning files
    addmissing <- paste("rm", fullnames, names, classification_file)
    system(addmissing, wait=TRUE)

    # Renaming output file
    addmissing <- paste("mv", total, classification_file)
    system(addmissing, wait=TRUE)

    # Final cleanup
    finalcleanup <- paste("rm", exon_data, exon_sorted, full_genome, missed, overlap_file, circrna_bed, genome_bed)
    system(finalcleanup, wait=TRUE)
}
