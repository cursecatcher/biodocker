#library(plotly)
library(heatmaply)


getHeatmap <- function(csv.file, html.filename, scratch.folder, output.folder, range.value) {
    
    html.filename.scratch = file.path(scratch.folder, html.filename)
    html.filename.output = file.path(output.folder, html.filename)

    x <- read.table(file=csv.file, header=TRUE, row.names=1)
    h <- heatmaply(x, plot_method="plotly", file=html.filename.scratch, limits=range.value)

    cat("COPYING FROM scratch TO output...\n")

    file.copy(html.filename.scratch, html.filename.output, overwrite = TRUE)
}


args <- commandArgs(TRUE)

myscratch <- args[1]
filename <- args[2]
lower.range <- as.integer(args[3])
upper.range <- as.integer(args[4])

#define heatmap ranges 
limits <- c(lower.range, upper.range)
if (is.na(lower.range) || is.na(upper.range) || lower.range == upper.range) {
    limits <- NULL
}



#define input file, produced in preprocessing step 
myfile = file.path(myscratch, "myfile.csv")

cat("Obtaining the heatmap...\n")

getHeatmap(
    csv.file = myfile,
    html.filename = paste0(tools::file_path_sans_ext(basename(filename)), ".html"), 
    scratch.folder = myscratch,
    output.folder = "/data/out", 
    range.value = limits
)