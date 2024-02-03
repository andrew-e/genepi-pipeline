kis <- vroom::vroom("scratch/data/targets/ki_database.tsv")

kis_drug_name <- c("ARIPIPRAZOLE", "CLOZAPINE", "OLANZAPINE", "Paliperidone", "Quetiapine", "RISPERIDONE")

ki_values <- kis[!is.na(kis$unigene) & !is.na(kis$ki_val) & kis$ligand_name %in% kis_drug_name & kis$species == "HUMAN", c("unigene", "ligand_name", "ki_val")]

genes <- unique(ki_values$unigene)
metabrain_top_hits <- "/mnt/storage/private/mrcieu/data/qtl_top_hits/metabrain/top_hits/"
regions <- c("basalganglia", "cerebellum", "cortex", "hippocampus", "spinalcord")

all_hits <- lapply(regions, function(region) {
	file <- paste0(metabrain_top_hits, region, "_eur.tsv.gz")
	top_hits <- vroom::vroom(file)
	top_hits$region <- region
	return(top_hits)
}) |> dplyr::bind_rows()

for (gene in genes) {

}