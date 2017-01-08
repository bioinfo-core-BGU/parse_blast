library(optparse)
library(tools)

# This script 

# paste("Rscript",
#       "compare_blast_parsed_reports.R",
#       "--blastind",    "Amir_Sagi/parse_perSample_CLC/parsed_BLAST_files_index.txt",
#       "--variable",  "bitscore",
#       "--name",    "qseqid",
#       "--full_txt_output"
#       ) %>% system
# paste("Rscript","compare_blast_parsed_reports.R","-h") %>% system


args    = commandArgs(trailingOnly = F) 
filepos = args %>% grep("--file=",x = .)
curfile = sub(pattern     = "--file=",
              replacement = "",
              x           = args[filepos]) %>% file_path_as_absolute
print(curfile)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-b", "--blastind"),      type="character", default=NULL, 
                help="Path to index file of blast results",   metavar="character"),
    make_option(c("-v", "--variable"),      type="character", default="bitscore", 
                help="Variable to keep in table (bitscore).", metavar="character"),
    make_option(c("-n", "--name"),      type="character", default="name", 
                help="Column to use as gene names, i.e. row identifier.", metavar="character"),
    make_option(c("-o", "--output"),        type="character", default=NULL, 
                help="Path to output file (default=<blastind>.output)", metavar="character"),
    make_option(c("-u", "--full_txt_output"), type="logical", action = "store_true", 
                help="Path to output file", metavar="character")
); 

opt_parser = optparse::OptionParser(usage       = "
This script takes an index file of form 'ID\tfilename' where filename is a BLAST tabular output WITH headers. 
The files are read and merged and then reshaped such that columns are sample IDs, rows are genes (--name) and the values of the tables are those in column --variable.


usage: %prog [options]", 
                                    option_list = option_list,
                                    epilogue    = "NOTE: The BLAST tables must have a header line and --variable and --name must be column titles.\n\nAuthor: Menachem Sklarz");
opt = optparse::parse_args(opt_parser);


if (is.null(opt$blastind)){
    paste("Rscript",curfile,"-h") %>% system
    stop(cat("No --blastind was passed!\n"))
}

if (is.null(opt$output)){
    opt$output <- sprintf("%s.output",opt$blastind)
    sprintf("No --output passed. Using %s as output\n", opt$output) %>% cat
}

# Libraries are loaded here so that the opt parsing is done before failure due to uninstalled packages.
library(magrittr)
library(plyr )
library(reshape2 )


# For debugging:
# opt$output = "Amir_Sagi/parse_perSample_CLC/output"
# opt$blastind = "Amir_Sagi/parse_perSample_CLC/parsed_BLAST_files_index.txt"
# User defined parameters (future):
#   blast : blast_report
#   blast_names : names of blast columns


# Read table of files
file_index = read.delim(opt$blastind,
                        he=F,
                        stringsAsFactors = F,
                        comment.char = "#")

# For each file in table:

#   Create list of tables, one for each sample:
blast_df <- apply(  X      = file_index,
                    MARGIN = 1,
                    FUN    = function(x) {
                                            # Reading parsed results file:
                                            blast_table <- try_default(expr = read.delim(x[2],
                                                                                         he = T),
                                                                       default = NULL) 
                                            # No file read:
                                            if(blast_table %>% is.null) {
                                                cat(sprintf("Could not read file %s\n",x[2]))   
                                                return(NULL)
                                            }
                                            # Solving issue of empty files:
                                            if(dim(blast_table)[1]==0) blast_table[1,1] <- NA
                                            # Setting sample name column:
                                            blast_table$sample_name <- x[1];
                                            return(blast_table)
                                          })
#   concatenate by row (rbind) = merge all tables into one long table
blast_df  <- do.call("rbind", blast_df)

# Find samples with no table (i.e. samples in file_index with no entry in blast_df)
empty_samples <- setdiff(file_index$V1,
                         blast_df$sample_name)
# If there were missiong tables:
if(length(empty_samples) > 0)  {
    # Add one row for each missing table
    blast_df[dim(blast_df)[1]+length(empty_samples),
             dim(blast_df)[2]] <- NA
    # Set sample names for nes rows 
    blast_df$sample_name[(dim(blast_df)[1]-length(empty_samples)+1):(dim(blast_df)[1])] <- empty_samples
    
}
        
# Chaeck that both --name and --variable exist in blast_df header line
if(!all(c(opt$name, opt$variable) %in% names(blast_df))) {
    stop(sprintf("--name (%s) or --variable (%s) do not exist in table header (%s)",
                 opt$name, 
                 opt$variable,
                 paste(names(blast_df),
                       collapse = " ")
                 )
         )
    }


# Convert to tabular form (reshape or similar) using a particular variable (bitscore by default)
blast_tbl <- dcast(data       = blast_df,
                    formula   = sprintf("%s ~ sample_name",opt$name),
                    value.var = opt$variable)
# Remove rows where the name is NA:
# These rows are introduced above for samples with non-existant results tables.
# This is a good thing because now these samples have columns with all-NA
blast_tbl <- blast_tbl[not(is.na(blast_tbl[opt$name])),]

cat("Producing outputs:\n------------------\n")
cat(sprintf("Writing summary table to: %s\n", opt$output))
write.table(x         = blast_tbl, 
            file      = opt$output,
            append    = F,
            quote     = F,
            sep       = "\t",
            row.names = F)

# Save either as text flat file or as R structure, as per full_txt_output parameter.
if(exists(x = "full_txt_output", where = opt)) {
    full_name <- paste(dirname(file_path_as_absolute(opt$output)),
                       .Platform$file.sep,
                       "flat_file_table.txt",
                       sep="")
    cat(sprintf("Writing full flat table to: %s\n", full_name))
    write.table(x         = blast_df, 
                file      = full_name,
                append    = F,
                quote     = F,
                sep       = "\t",
                row.names = F)
} else {
    full_name <- paste(dirname(file_path_as_absolute(opt$output)),
                       .Platform$file.sep,
                       "flat_file_table.txt",
                       sep="")
    cat(sprintf("Writing full flat table to: %s\n", full_name))
    save(blast_df,
         file = full_name
    )
    
}
