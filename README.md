# Scripts for parsing and merging BLAST reports

## **`parse_blast`**: A script for parsing tabular BLAST output.

### Brief Description

`parse_blast.R` performs the following tasks:

1. It adds annotation to raw tabular BLAST output files,
2. filters the BLAST results by several possible fields,
3. selects the best hit for a group when passed a grouping field and
4. extracts the sequences equivalent to the alignments.

The program receives two input files:

1. `--blast`: A tabular BLAST report, as created with `blast --outfmt 6` (For non-default output, you have to pass the list of columns the --names parameter). 
2. `--dbtable`: A table containing the annotation for the query sequences.

`parse_blast.R`:

- merges the two tables by the fields specified in parameters `--merge_blast` and `--merge_metadata`, respectively; 
- filters the table by the values supplied in the `--min_...` and `--max_...` parameters;
- sorts the table in the order specified in `--sort_str`;
- groups the hits by the field defined in `--group_dif_name`;
- optionally, adds a field with the fasta sequence of the alignment (or the whole target) and
- exports the resulting table (optionally with only a subset of the original columns)

If no `--output` is specified, will use the `<--blast>.parsed` as output destination.

### Program Usage:


    Usage: bin/parse_blast.R [options]
    
    
    Options:
        -b CHARACTER, --blast=CHARACTER
            Path to blast results
    
        -t CHARACTER, --dbtable=CHARACTER
            Path to table of database sequence (metadata file).
    
        -o CHARACTER, --output=CHARACTER
            Path to output file (default=<blast_input>.parsed)
    
        -n CHARACTER, --names=CHARACTER
            List of names in blast table (default: "qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe").
    
        -m CHARACTER, --merge_blast=CHARACTER
            Column in blast file to merge with metadata (default: sallseqid)
    
        -r CHARACTER, --merge_metadata=CHARACTER
            Column in metadata file to merge with blast results
    
        --group_dif_name=CHARACTER
            Name of field in metadata file that distinguishes between groups (default: "name").
    
        -s CHARACTER, --min_bitscore=CHARACTER
            Minimum bitscore to permit. (default: 0)
    
        -c CHARACTER, --min_coverage=CHARACTER
            Minimum coverage to permit. NOTE: Coverage is defined as fraction: alignment_length/query_length. (default: 0)
    
        -i NUMERIC, --min_pident=NUMERIC
            Minimum %identity to permit. (default: 0)
    
        -p NUMERIC, --min_align_len=NUMERIC
            Minimum protein length to permit. (default: 0)
    
        -e NUMERIC, --max_evalue=NUMERIC
            Maximum evalue to permit. (default: 1)
    
        --num_hits=NUMERIC
            Number of hits to report (numeric or 'all') (default: 1)
    
        --sort_str=NUMERIC
            Method to select best HSPs. A comma-separated string of column names.
                    Each column name should be followed by 'd' for decreasing (higher is better) 
                    or 'i' for increasing (lower is better) (default: 'bitscore,d,evalue,i,pident,d')
                    NOTE: You can also use 2 special columns: 'coverage' and 'align_len', e.g. 'align_len,d'
    
        -k NUMERIC, --columns2keep=NUMERIC
            Columns from blast and metadata files to keep. Must be a subset of --names and metadata title line (default: all columns)
    
        -f NUMERIC, --fasta2extract=NUMERIC
            fasta file to extract sequences from. (default: NULL - no extraction) 
    
        --fasta_header=NUMERIC
            Field in BLAST equivalent to FASTA headers (default: 'sallseqid') 
    
        --fasta_allsubject
            If set, the whole target sequence will be extracted rather than just the alignment.
    
        -g, --single
            Results represent single gene (default: FALSE)
    
        --store_opts=NUMERIC
            Target to store the opts R structure
    
        -h, --help
            Show this help message and exit
    

## **`compare_blast_parsed_reports.R`**: A script for comparing tabular BLAST outputs.

### Brief Description

This script takes an index file of form 'ID<TAB>filename' where filename is a BLAST tabular output WITH headers (as `parse_blast.R` produces).

The files are read and merged and then reshaped such that columns are sample IDs, rows are genes (or the column in the blast reports specified with the `--name` parameter) and the values of the cells in the table are those in column `--variable`.
    

### Program Usage:

    
    usage: bin/compare_blast_parsed_reports.R [options]
    
    
    Options:
            -b CHARACTER, --blastind=CHARACTER
                    Path to index file of blast results
    
            -v CHARACTER, --variable=CHARACTER
                    Variable to keep in table (bitscore).
    
            -n CHARACTER, --name=CHARACTER
                    Column to use as gene names, i.e. row identifier.
    
            -o CHARACTER, --output=CHARACTER
                    Path to output file (default=<blastind>.output)
    
            -u, --full_txt_output
                    Path to output file
    
            -h, --help
                    Show this help message and exit





Contact
---------

Please contact Menachem Sklarz at: [sklarz@bgu.ac.il](mailto:sklarz@bgu.ac.il)

