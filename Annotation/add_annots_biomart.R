
# using the R package biomaRt to add gene descriptions to peptide matches from Trinotate annotation

library("biomaRt")
ensembl=useMart("ensembl")
listDatasets(ensembl)

## for zebrafish

ensemblZ = useDataset("drerio_gene_ensembl",mart=ensembl)
la <- listAttributes(ensemblZ)
View(la)

getBM(attributes = c("ensembl_peptide_id_version", "ensembl_gene_id","description"),
      filters = 'ensembl_peptide_id_version',
      values = "ENSDARP00000103076.3",
      mart = ensemblZ)

ENSONIP00000022996.1


## get list from file
pep_list <- scan("~/WRASSE/Anns_zebrafish_peptides.txt", what = "", quiet = TRUE)
length(pep_list)

## now run full biomart on all

pepInfo <- getBM(attributes = c("ensembl_peptide_id_version", "ensembl_gene_id","description"),
      filters = 'ensembl_peptide_id_version',
      values = pep_list,
      mart = ensemblZ)


View(pepInfo)

## now output to file

write.table(pepInfo, file = "Anns_zebrafish_peptide_gene_descriptions.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)









