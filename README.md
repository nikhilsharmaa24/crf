## This work is about finding conserved region in sequence data. Clustal omega software is used to do multiple sequence alignment using msa package from bioconductor.

### Installing the packge use devtools
devtools::install_github('nikhilsharmaa24/crf')

library(crf)

### You may not need to call these libraries, but if some functions are not working trying calling and running code
#library(msa)

#library(Biostrings)

#library(plyr)

#library(tidyverse)

### Example usage of the functions
dat <- analyze_fasta_alignment('path_to_your_txt_file.txt',  alignment_method = 'ClustalOmega')

### Analyse for finding conserved region 
res <- analyze_fragments(dat, 12)

### Arrange in desending order  
res$Final_table %>% arrange(-Count) -> res$Final_table

### Saving conserved region information 
write.csv(res$Final_table, 'results.csv', row.names = FALSE)
