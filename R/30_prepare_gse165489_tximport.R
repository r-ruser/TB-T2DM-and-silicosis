source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "tximport", "AnnotationDbi", "org.Hs.eg.db"))

acc <- "GSE165489"
raw_tar <- file.path(path_data, "02_GEO_bulk", acc, "raw", paste0(acc, "_RAW.tar"))
processed_dir <- file.path(path_data, "02_GEO_bulk", acc, "processed")
nested_dir <- file.path(processed_dir, "salmon_nested_archives")
quant_dir <- file.path(processed_dir, "salmon_quant")
dir.create(nested_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(quant_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(raw_tar)) stop("GSE165489 RAW archive missing")

entries <- utils::untar(raw_tar, list = TRUE)
entries <- entries[grepl("_quant[.]tar[.]gz$", entries)]
if (length(entries) != 86L) stop("Expected 86 nested Salmon archives, found ", length(entries))
missing_nested <- entries[!file.exists(file.path(nested_dir, entries))]
if (length(missing_nested)) {
  write_log("Extracting ", length(missing_nested), " GSE165489 nested Salmon archives")
  utils::untar(raw_tar, files = missing_nested, exdir = nested_dir)
}

quant_files <- character(length(entries))
sample_ids <- sub("_.*$", "", basename(entries))
for (i in seq_along(entries)) {
  inner <- file.path(nested_dir, entries[i])
  q_entries <- utils::untar(inner, list = TRUE)
  q_entry <- q_entries[grepl("/quant[.]sf$", q_entries)]
  if (length(q_entry) != 1L) stop("quant.sf not uniquely found in ", basename(inner))
  sample_dir <- file.path(quant_dir, sample_ids[i])
  final_q <- file.path(sample_dir, "quant.sf")
  if (!file.exists(final_q)) {
    td <- tempfile(paste0(sample_ids[i], "_")); dir.create(td)
    utils::untar(inner, files = q_entry, exdir = td)
    extracted <- file.path(td, q_entry)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(extracted, final_q, overwrite = FALSE)) stop("Could not stage quant.sf for ", sample_ids[i])
  }
  quant_files[i] <- final_q
}
names(quant_files) <- sample_ids
if (anyDuplicated(sample_ids)) stop("Duplicate GSM IDs in GSE165489 archive")

meta_file <- file.path(processed_dir, paste0(acc, "_metadata_curated.csv"))
meta <- data.table::fread(meta_file, encoding = "UTF-8")
if (!setequal(meta$sample_id, sample_ids)) stop("GSE165489 metadata and Salmon archive GSM IDs differ")
meta <- meta[match(sample_ids, sample_id)]

tx <- data.table::fread(quant_files[1], select = "Name")$Name
tx_stripped <- sub("[.][0-9]+$", "", tx)
entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(tx_stripped),
                                keytype = "REFSEQ", column = "ENTREZID", multiVals = "first")
tx2gene <- data.table::data.table(tx = names(entrez), gene = unname(entrez))
tx2gene <- unique(tx2gene[!is.na(gene) & nzchar(gene)])
if (nrow(tx2gene) < 15000L) stop("Insufficient RefSeq-to-gene mappings: ", nrow(tx2gene))

write_log("Importing GSE165489 Salmon estimates with tximport")
txi <- tximport::tximport(quant_files, type = "salmon", tx2gene = tx2gene,
                         ignoreTxVersion = TRUE, countsFromAbundance = "no",
                         dropInfReps = TRUE)
txi_scaled <- tximport::tximport(quant_files, type = "salmon", tx2gene = tx2gene,
                                ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM",
                                dropInfReps = TRUE)
if (!identical(colnames(txi$counts), meta$sample_id)) stop("tximport sample ordering failed")

symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(txi$counts),
                                keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
gene_map <- data.table::data.table(entrez_id = rownames(txi$counts), gene_symbol = unname(symbol))
safe_write_csv(gene_map, file.path(processed_dir, paste0(acc, "_gene_map.csv")))
safe_write_csv(meta, file.path(processed_dir, paste0(acc, "_metadata_tximport_order.csv")))
saveRDS(list(raw = txi, length_scaled = txi_scaled, tx2gene = tx2gene, metadata = meta),
        file.path(processed_dir, paste0(acc, "_tximport.rds")), compress = FALSE)

manifest <- data.table::data.table(
  sample_id = sample_ids, quant_file = normalizePath(quant_files, winslash = "/"),
  quant_bytes = file.info(quant_files)$size
)
manifest[, total_salmon_numreads := vapply(quant_files, function(f) {
  sum(data.table::fread(f, select = "NumReads")$NumReads)
}, numeric(1))]
manifest[, imported_gene_numreads := colSums(txi$counts)[sample_id]]
manifest[, imported_count_fraction := imported_gene_numreads / total_salmon_numreads]
safe_write_csv(manifest, file.path(path_result, "02_bulk", acc, "A_GSE165489_salmon_quant_manifest.csv"))
write_log("GSE165489 tximport preparation completed: samples=", ncol(txi$counts),
          "; genes=", nrow(txi$counts), "; mapped transcripts=", nrow(tx2gene))
