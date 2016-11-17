# AUTHOR: Menachem Sklarz & Michal Gordon

library(magrittr)
library(plyr )
library(optparse)
library(tools)



# paste("Rscript",
#       "parse_blast_general.R",
#       "--blast",           "MLST/SAH1503/SAH1503.blast.out",
#       "--dbtable",         "MLST/MLST_scheme.tab",
#       "--names",          "\"qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe\"",
#       "--merge_blast",    "qseqid",
#       "--merge_metadata", "Allele",
#       "--group_dif_name", "Gene",
#       "--min_pident",     "98",
#       "--min_coverage",   "98",
#       "--max_evalue",     "1e-7",
#       "--num_hits",       "1",
#       "--sort_str",       "bitscore,d",
#       "--fasta2extract",       "MLST/SAH1503.nucl.merge.fa",
#  #     "--fasta_allsubject",
#       "--fasta_header",       "sallseqid"
# ) %>% system

# paste("Rscript","parse_blast_general.R","-h") %>% system


args    = commandArgs(trailingOnly = F) 
filepos = args %>% grep("--file=",x = .)
curfile = sub(pattern     = "--file=",
              replacement = "",
              x           = args[filepos]) %>% file_path_as_absolute
print(curfile)


args = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-b", "--blast"), type="character", default=NULL, 
                help="Path to blast results", metavar="character"),
    make_option(c("-t", "--dbtable"), type="character", default=NULL, 
                help="Path to table of database sequence (metadata file).", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Path to output file (default=<blast_input>.parsed)", metavar="character"),
    make_option(c("-n", "--names"), type="character", default="qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe", 
                help="List of names in blast table (default: \"qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe\").", metavar="character"),
    make_option(c("-m", "--merge_blast"), type="character", default="sallseqid", 
                help="Column in blast file to merge with metadata (default: sallseqid)", metavar="character"),
    make_option(c("-r", "--merge_metadata"), type="character", default=NULL, 
                help="Column in metadata file to merge with blast results", metavar="character"),			
    make_option("--group_dif_name", type="character", default="name", 
                help="Name of field in metadata file that distinguishes between groups (default: \"name\").", metavar="character"),
    make_option(c("-s", "--min_bitscore"), type="double", default=0, 
                help="Minimum bitscore to permit. (default: 0)", metavar="character"),
    make_option(c("-c", "--min_coverage"), type="double", default=0, 
                help="Minimum coverage to permit. NOTE: Coverage is defined as fraction: alignment_length/query_length. (default: 0)", metavar="character"),
    make_option(c("-i", "--min_pident"), type="double", default=0, 
                help="Minimum %identity to permit. (default: 0)", metavar="numeric"),
    make_option(c("-p", "--min_align_len"), type="double", default=0, 
                help="Minimum protein length to permit. (default: 0)", metavar="numeric"),
    make_option(c("-e", "--max_evalue"), type="double", default=1, 
                help="Maximum evalue to permit. (default: 1)", metavar="numeric"),
    make_option(c("--num_hits"), type="character", default="1", 
                help="Number of hits to report (numeric or 'all') (default: 1)", metavar="numeric"),
    make_option(c("--sort_str"), type="character", default="bitscore,d,evalue,i,pident,d", 
                help="Method to select best HSPs. A comma-separated string of column names.
                Each column name should be followed by 'd' for decreasing (higher is better) 
                or 'i' for increasing (lower is better) (default: 'bitscore,d,evalue,i,pident,d')
                NOTE: You can also use 2 special columns: 'coverage' and 'align_len', e.g. 'align_len,d'", metavar="numeric"),
    make_option(c("-k","--columns2keep"), type="character", default=NULL, 
                help="Columns from blast and metadata files to keep. Must be a subset of --names and metadata title line (default: all columns)", metavar="numeric"),
    make_option(c("-f","--fasta2extract"), type="character", default=NULL, 
                help="fasta file to extract sequences from. (default: NULL - no extraction) ", metavar="numeric"),
    make_option(c("--fasta_header"), type="character", default="sallseqid", 
                help="Field in BLAST equivalent to FASTA headers (default: 'sallseqid') ", metavar="numeric"),
    make_option(c("--fasta_allsubject"), type="logical", action="store_false", 
                help="If set, the whole target sequence will be extracted rather than just the alignment.", metavar="numeric"),
    make_option(c("-g", "--single"), type="logical", action = "store_true",
                help="Results represent single gene (default: FALSE)"),
    make_option(c("--store_opts"), type="character", default=NULL, 
                help="Target to store the opts R structure", metavar="numeric")
); 



opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor: Menachem Sklarz");
opt = optparse::parse_args(opt_parser);

if (is.null(opt$blast)){
  opt_parser %>% print_help
  stop(cat("No --blast was passed!"))
}
if (is.null(opt$dbtable)){
  # opt_parser %>% print_help
  cat(sprintf("No --dbtable was passed! Using raw blast results...\n"))
}
if (is.null(opt$merge_metadata) & !is.null(opt$dbtable)){
    opt_parser %>% print_help
    stop(cat("--dbtable was passed but no --merge_metadata was passed!\n"))
}
if (is.null(opt$merge_blast) & !is.null(opt$dbtable)){
    opt_parser %>% print_help
    stop(cat("--dbtable was passed but no --merge_blast was passed!\n"))
}
if (is.null(opt$output)){
    opt$output <- sprintf("%s.parsed",opt$blast)
    sprintf("No --output passed. Using %s as output\n", opt$output) %>% cat
}

# opt$dbtable        = "MLST/MLST_scheme.tab"
# opt$blast          =  "MLST/SAH1503/SAH1503.blast.out"
# opt$output         = paste(opt$blast,"parse.out",sep=".")
# opt$merge_blast    = "qseqid"
# opt$merge_metadata = "Allele"
# opt$group_dif_name = "Gene"
# opt$names = "qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe"
# opt$columns2keep = "qseqid sallseqid qlen"
# opt$sort_str       = "bitscore,d"
# opt$min_pident     = 98
# opt$min_coverage   = 98
# opt$fasta_header   = "sallseqid"
# opt$fasta2extract  = "MLST/SAH1503.nucl.merge.fa"
# opt$num_hits       = "1"

if(!is.null(opt$store_opts)) {
    save(opt, file = opt$store_opts)
}

# Read blast results:
blast <- read.delim(opt$blast, 
                    he = F, 
                    stringsAsFactors = F)
# print(opt$merge_metadata)
my.colnames = strsplit(opt$names,
                       split = " ") %>% unlist
# print(opt$merge_blast)
# 
# Check matche between names and blast table:
if(length(colnames(blast)) != length(my.colnames)){
  opt_parser %>% print_help
  stop(cat("length of names and blast file don't match!"))
}
# Set column names 
names(blast) = strsplit(opt$names,
                        split = " ") %>% unlist

# The following section merges in the metadata. Do only if metadata is passed!
if(!is.null(opt$merge_metadata)) {
    
    # Make sure the column names used for merging with metadata exists:
    stopifnot(opt$merge_blast %in% colnames(blast) )
    
    # Read gene metadata:
    metadata = read.delim(file             = opt$dbtable,  
                          header           = T,  
                          comment.char     = "#",
                          quote            = "",
                          stringsAsFactors = F)
    # Test passes merge field exists in metadata file:
    if(!(opt$merge_metadata %in% names(metadata))){
        opt_parser %>% print_help
        stop(sprintf("--merge_metadata ('%s') does not exist in metadata columns! (%s)\n",
                     opt$merge_metadata,
                     paste(names(metadata),collapse = ", ")))
    }
    

    # Merge metadata into main blast results:
    blast <- merge(blast,
                   metadata,
                   by.y = opt$merge_metadata,
                   by.x = opt$merge_blast)
    
}

# Check that opt$group_dif_name exists in blast df:
if(!(opt$group_dif_name %in% names(blast))){
    opt_parser %>% print_help
    stop(sprintf("--group_dif_name ('%s') does not exist in blast or metadata columns! (%s)\n",
                 opt$group_dif_name,
                 paste(names(blast),collapse = ", ")))
}

# Calculating coverage as length of hit (hsp) / total length
# This is what Eliad used (possibly inherited from Nir Gilad)
# See /fastspace/bioinfo_projects/SPICE/01.Legionella/12.Pathogenicity_VFDB_Acinetobacter/02e.Pipe.VFDB_updated/virulanceAllGenes.pl
# Line 96
# Not exacly the same - he parsed a full blast output and here I'm parsing a tabular output...
if (all("qstart" %in% names(blast), 
        "qend"   %in% names(blast), 
        "qlen"   %in% names(blast))) {
    
    blast$align_len = abs(blast$qstart - blast$qend) + 1
    blast$coverage  = blast$align_len / blast$qlen * 100

} else {
    if(any(opt$min_coverage > 0, opt$min_align_len > 0)) {
        opt_parser %>% print_help
        stop(cat("Fields in blast report don't enable filtering by coverage and protein length
                 The BLAST report must have 'qstart', 'qend' and 'qlen' fields."))
    } 
}



# If no columns2keep was passed, set to all columns,
# otherwise, check that all columns2keep are included in the blast column names
if (is.null(opt$columns2keep)){
    opt$columns2keep <- names(blast)
    sprintf("No --columns2keep passed. Outputting all columns\n") %>% cat
} else {
    opt$columns2keep <- strsplit(x     = opt$columns2keep,
                                 split = "\\s+",
                                 perl  = T) %>% unlist
    if(!all(opt$columns2keep %in% names(blast))) {
        stop(sprintf("--columns2keep contains unrecognised columns... (%s).\n",
                     paste(setdiff(opt$columns2keep,names(blast)),
                           collapse = ", ")))
    }
}

# Convert sort_str from input string into valid string for sorting below:
opt$sort_str <- 
    opt$sort_str                           %>%
    strsplit(.,",")                    %>%     # Split by comma
    unlist                             %>%     # Convert to character vector
    matrix(.,
           ncol=2,
           byrow = T)                  %>%     # Convert to 2-column matrix
    data.frame(.,stringsAsFactors = F) %>%     # Convert to data.frame
    setNames(.,c("column","direction"))        # Set names of df

# Checking that second column is only 'd' or 'i'
if(length(setdiff(opt$sort_str$direction,c("d","i"))) > 0) {
    opt_parser %>% print_help
    stop(cat("--sort_str is badly formatted. Make sure every field has a direction tag ('d' or 'i') following it\n"))
}
# Create new column with '-' for decreasing order:
opt$sort_str$direction_tag <- ifelse(test = opt$sort_str$direction=="d","-","")
# Storing column names:
sort_cols <- opt$sort_str$column
# Paste to create string:
opt$sort_str <- paste(paste(opt$sort_str$direction_tag,
                            opt$sort_str$column,
                            sep=""),
                      collapse = ",")


if(!all(sort_cols %in% names(blast))) {
    stop(sprintf("--sort_str contains unrecognised columns... (%s).\n",
                 paste(setdiff(sort_cols,
                               names(blast)),
                       collapse = ", ")))
}

# ------------------------------------
# Retrieve best ortholog (based on bitscore) for each VFG query:

# Helper function. Sorts a sub-data.frame by bitscore, then by evalue and returns the first row.
get_best_hit <- function(x) {
    
    # Step 1 - filter out lines that don't match the requierments (give in opts)
    if(any(exists("coverage", where=x), 
           exists("align_len", where=x))) {
        x <- subset(x,subset = bitscore  >= opt$min_bitscore  &
                        evalue    <= opt$max_evalue    &
                        coverage  >= opt$min_coverage  &
                        align_len >= opt$min_align_len &
                        pident    >= opt$min_pident)
        
    } else {
        x <- subset(x,subset = bitscore  >= opt$min_bitscore  &
                        evalue    <= opt$max_evalue    &
                        pident    >= opt$min_pident)
        
    }


    # If all lines are required or more lines than exist in the result group, set num hits to number of rows.    
    if((opt$num_hits == "all") | (as.numeric(opt$num_hits) > dim(x)[1])) {
        opt$num_hits <- dim(x)[1]
    }

    # Step 2 - take 'opt$num_hits' results of those that pass filtration
    retval <- eval(parse(text = sprintf("x[with(x,order(%s)),][1:as.numeric(opt$num_hits),]",
                                        opt$sort_str)))
    
    # If all results for a gene were filtered out, retval will be all NAs
    # return NULL to remove from final filtered blast table
    if(all(is.na(retval))) return(NULL)
    
    return(retval)
}


if(exists("single", where = opt)) {
    # Use entire BLAST report as result for one gene. Return best line in entire table:
    blast <- get_best_hit(blast)
} else {
    # Use 
    cat("Running ddply. Takes a while...\n")
    # Cuts the 'blast' df into pieces by 'name', and runs get_best_hit on each slice.
    # The lines returned by get_best_hit are concatenated into a new df, 'blast'...
    blast_df_names <- names(blast)
    blast <- ddply(.data      = blast, 
                   .variables = as.quoted(opt$group_dif_name),
                   .fun       = get_best_hit)
    # If there are no results, set blast to an empty df with the same col names as before:
    if(dim(blast)[1]==0) {
        blast <- data.frame(matrix(vector(), 0, length(blast_df_names),
                                   dimnames=list(c(), blast_df_names)),
                            stringsAsFactors=F)
    }

}

# Print output to file (if there are results at all):
write.table(blast[,opt$columns2keep],
            file      = opt$output,
            quote     = F,
            row.names = F,
            sep       = "\t")

cat("Done...\n")
cat(sprintf("See output file at: %s\n", opt$output))

# If user requested extraction of sequences:
if(exists("fasta2extract",opt)) {
    # stop("This option is not implemented yet. Sorry...")
    # Load Biostrings library:
    library(Biostrings)
    # Check fasta file exists:
    if(!file.exists(opt$fasta2extract)) {
        stop(sprintf("FASTA file %s does not exist!\n",opt$fasta2extract))
    }
    # Create fasta file index
    fastaindex <- fasta.index(filepath =  opt$fasta2extract)
    # For each blast result:
    for(i in 1:dim(blast)[1]) {
        # If it is in reverse: (see documentation below for non-reverse)
        newseq <- readBStringSet(fastaindex[grep(pattern = blast[i,opt$fasta_header], 
                                                 x       = fastaindex$desc),])
        if(blast[i,"sstart"] > blast[i,"send"]) {
            if(!exists("fasta_allsubject",opt)) {
                newseq <- subseq(newseq, 
                               start = blast[i,"send"], 
                               end = blast[i,"sstart"])    %>%
                            DNAStringSet                   %>%
                            reverseComplement                
                names(newseq) <- paste(names(newseq),blast[i,"send"],blast[i,"sstart"],"revcomp",
                                       sep="_")
                if(exists("dbtable",opt)) {
                    names(newseq) <- paste(names(newseq), blast[i,opt$merge_blast], sep="|")
                }
            }
            writeXStringSet(newseq, 
                            filepath = sprintf("%s.fasta",opt$output), 
                            append = TRUE,
                            compress=FALSE, 
                            format="fasta")
        } else {
            # Get the full sequence from the fasta file 
            # get subsequence of it based on blast table (sstart and send)
            # newseq <- fastaindex[which(fastaindex$desc==blast[i,opt$fasta_header])]  
            if(!exists("fasta_allsubject",opt)) {
                newseq <- subseq(newseq, 
                               start = blast[i,"sstart"], 
                               end = blast[i,"send"])
                names(newseq) <- paste(names(newseq),blast[i,"sstart"],blast[i,"send"],
                                       sep="_")
                if(exists("dbtable",opt)) {
                    names(newseq) <- paste(names(newseq), blast[i,opt$merge_blast], sep="|")
                }
                
            }
            # Add coordinates to the seq name:
            # Write the sequence to file
            writeXStringSet(newseq, 
                            filepath = sprintf("%s.fasta",opt$output), 
                            append = TRUE,
                            compress=FALSE, 
                            format="fasta")
            
            
        }
        
    }
}

