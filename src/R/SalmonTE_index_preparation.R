# Aim: Prepare files to make an index with spruce whole-length LTR-TEs for SalmonTE  
# Prepare a list of whole-length LTR-TEs found by Andrea to extend a SalmonTE reference database and add family names to the fasta identifiers  
  
# Load libraries
library(Biostrings)

# 1) Add our set of LTR-TEs to clades_extended.csv  

# Import fasta file with LTR-TEs
TE_fa <- readDNAStringSet("/mnt/picea/projects/spruce/nstreet/sRNA/LTR-TE/fasta/LTR-TE_CopiaGypsy_whole-length.fa")

# Get the names of fasta sequences
names_TE <- names(TE_fa)
# Make a list with name of TE, family and a superfamily, which will be appended to clades_extended.csv
names_family <- str_extract(names_TE, "[C,G][0-9]+[A-E]?(_Alisei)?")
names_superfamily <- str_extract(names_TE, "C|G")
names_superfamily <- sub("C", "Copia", names_superfamily)
names_superfamily <- sub("G", "Gypsy", names_superfamily)
names_LTRtransposon <- rep("LTR Retrotransposon", length(names_TE))
names_element <- rep("Transposable Element", length(names_TE))
# Construct a table
hierarchy_table <- cbind(names_family, names_superfamily, names_LTRtransposon, names_element)
hierarchy_table <- unique(hierarchy_table)

# Append table to existing file and export
clades_extended <- read.csv("/mnt/picea/projects/spruce/nstreet/sRNA/LTR-TE/SalmonTE/clades_extended_v0.4.csv", header = FALSE)
colnames(hierarchy_table) <- NULL
clades_extended_updated <- rbind(clades_extended, hierarchy_table)
# Export new file
write.table(clades_extended_updated, "/mnt/picea/projects/spruce/nstreet/sRNA/LTR-TE/SalmonTE/clades_extended.csv", 
            sep = ",",
            quote = FALSE, 
            row.names = FALSE,
            col.names = FALSE)

# 2) Add a family name to the fasta sequence identifier, separated from a name by a tab
names(TE_fa) <- paste(names(TE_fa), names_family, sep ="\t")
writeXStringSet(TE_fa, "/mnt/picea/projects/spruce/nstreet/sRNA/LTR-TE/SalmonTE/LTR-TE_CopiaGypsy_whole-length_familyName.fa")