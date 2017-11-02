
paste("Rscript",
      "./src/parse_blast.R",
      "--blast",           "MLST/SAH1503/SAH1503.blast.out",
      "--dbtable",         "MLST/MLST_scheme.tab",
      "--names",          "\"qseqid sallseqid qlen slen qstart qend sstart send length evalue bitscore score pident qframe\"",
      "--merge_blast",    "qseqid",
      "--merge_metadata", "Allele",
      "--group_dif_name", "Gene",
      "--min_pident",     "98",
      "--min_coverage",   "98",
      "--max_evalue",     "1e-7",
      "--num_hits",       "1",
      "--sort_str",       "bitscore,d",
      "--fasta2extract",       "MLST/SAH1503.nucl.merge.fa",
      "--fasta_allsubject",
      "--fasta_header",       "sallseqid"
) %>% system
#) %>% paste

paste("Rscript","src/parse_blast.R","-h") %>% system



# This script

paste("Rscript",
      "compare_blast_parsed_reports.R",
      "--blastind",    "",
      "--variable",  "bitscore",
      "--name",    "qseqid",
      "--full_txt_output"
      ) %>% system
paste("Rscript","compare_blast_parsed_reports.R","-h") %>% system
