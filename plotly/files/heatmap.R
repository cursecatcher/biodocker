#library(plotly)
library(heatmaply)


getHeatmap <- function(csv.file, html.filename, scratch.folder, output.folder, color.palette) {

    html.filename.scratch = file.path(scratch.folder, html.filename)
    html.filename.output = file.path(output.folder, html.filename)

    x <- read.table(file=csv.file, header=TRUE, row.names=1)
    h <- heatmaply(x, plot_method="plotly", file=html.filename.scratch, colors=color.palette)

    if (file.exists(html.filename.scratch)) {
      cat("COPYING FROM scratch TO output...\n")
      file.copy(html.filename.scratch, html.filename.output, overwrite = TRUE)
      system(paste("chmod 777", html.filename.output), wait=TRUE)
    } else {
      quit(save="no", status=1)
    }
}


args <- commandArgs(TRUE)

myscratch <- args[1]
palette.string <- args[2]
#filename <- args[2]

#define input file, produced in preprocessing step
myfile = file.path(myscratch, "myfile.csv")

palette <- switch(tolower(palette.string), 
  "viridis" = viridis::viridis, 
  "brbg" = BrBG, 
  "magma" = viridis::magma, 
  "plasma" = viridis::plasma, 
  "cividis" = viridis::cividis
)

if (identical(palette, NULL)) {
  palette <- viridis::viridis
}

cat("Obtaining the heatmap...\n")

getHeatmap(
    csv.file = myfile,
    html.filename = "heatmaply.html",
    scratch.folder = myscratch,
    output.folder = "/data/out", 
    color.palette = palette 
)
