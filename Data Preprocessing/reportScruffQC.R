require(Biostrings)
require(optparse)
require(readxl)
require(rmarkdown)
require(parallel)
require(stringr)


### parse arguments
option_list <- list(optparse::make_option(c("-P","--Path"),
                                          type="character",
                                          default="/restricted/projectnb/pcga/Conor",
                                          help="Where to find --File and also output data"),
                    optparse::make_option(c("-B","--Batch"),
                                          type="integer",
                                          help="The batch of plates"),
                    optparse::make_option(c("-F", "--File"),
                                          type="character",
                                          help="Name of the file with experiment and read paths"),
                    optparse::make_option(c("-v", "--verbose"), 
                                          action="store_true", default=TRUE,
                                          help="Print extra output [default]"))

arguments <- optparse::parse_args(optparse::OptionParser(option_list=option_list), positional_arguments=TRUE)
opt <- arguments$options
path <- opt[["Path"]]
#print(path)
batch <- opt[["Batch"]]
f <- opt[["File"]]
v <- opt[["verbose"]]

setwd(path)
cat("Present working directory is",getwd(),"\n")
file <- read_excel(f, col_names = TRUE, col_types = "text")

project <- "PCGA2"
experiment <- file$experiment
cores <- min(4, detectCores() - 2)
lane <- paste0("L",rep(x = sprintf(fmt = "%03d", ... = 1:length(experiment)), length.out = length(experiment)))
read1Path <- str_trim(string = file$R1, side = "both")
read2Path <- str_trim(string = file$R2, side = "both")
showSetup <- showSession <- v
output_file <- paste0("Scruff_ResultReport_Batch",batch)
output_dir <- path

render(input = "/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/plotScruffQCReport.Rmd", 
       params = list(project = project, experiment = experiment,
                     lane = lane, read1Path = read1Path, 
                     read2Path = read2Path, cores = cores, 
                     showSetup = showSetup, showSession = showSession), 
                     output_file = output_file, output_dir = output_dir, 
                     intermediates_dir = output_dir, knit_root_dir = output_dir)
 
