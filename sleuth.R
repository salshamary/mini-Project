library(sleuth)
library(data.table)
library(dplyr)
#tale for sleuth
stab <- read.table("sample_covariates.txt", header=TRUE, stringsAsFactors=FALSE, sep = '\t')
so <- sleuth_prep(stab)
so
#fits model compares the conditions
so <- sleuth_fit(so, ~condition, 'full')
so
so <- sleuth_fit(so, ~1, 'reduced')
so<- sleuth_lrt(so, 'reduced', 'full')
#gets reslts
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#significant results
sleuth_significant<- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)
sleuth_sig <- data.frame(sleuth_significant$target_id, sleuth_significant$test_stat, sleuth_significant$pval, sleuth_significant$qval)
colnames(sleuth_sig) <- c("target_id", "test_stat", "pval","qval")
sleuth_sig
#writes table
write.table(sleuth_significant, file="sleuth_output.txt",quote = FALSE,row.names = FALSE)
