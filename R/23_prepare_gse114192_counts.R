source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table"))

acc <- "GSE114192"
tar_file <- file.path(path_data, "02_GEO_bulk", acc, "raw", paste0(acc, "_RAW.tar"))
meta_file <- file.path(path_data, "02_GEO_bulk", acc, "processed", paste0(acc, "_metadata_curated.csv"))
if (!file.exists(tar_file)) stop("Verified GSE114192 RAW archive not found")
if (!file.exists(meta_file)) stop("Run R/21_curate_primary_metadata.R first")

members <- utils::untar(tar_file, list = TRUE)
members <- members[grepl("[.]txt[.]gz$", members, ignore.case = TRUE)]
if (length(members) != 249L) warning("Expected 249 count files; found ", length(members))
extract_dir <- tempfile("GSE114192_counts_")
dir.create(extract_dir)
on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
status <- utils::untar(tar_file, files = members, exdir = extract_dir)
if (!identical(status, 0L)) stop("Failed to extract GSE114192 count archive")

read_count <- function(member) {
  f <- file.path(extract_dir, member)
  x <- data.table::fread(f, header = FALSE, col.names = c("gene_id", "count"))
  if (ncol(x) != 2L || anyDuplicated(x$gene_id)) stop("Invalid count file: ", member)
  x
}
first <- read_count(members[1])
gene_id <- first$gene_id
counts <- matrix(0, nrow = nrow(first), ncol = length(members),
                 dimnames = list(gene_id, sub("_.*$", "", basename(members))))
counts[, 1] <- first$count
for (i in 2:length(members)) {
  x <- read_count(members[i])
  if (!identical(x$gene_id, gene_id)) stop("Gene order mismatch in ", members[i])
  counts[, i] <- x$count
}
storage.mode(counts) <- "integer"

meta <- data.table::fread(meta_file, encoding = "UTF-8")
if (!setequal(colnames(counts), meta$sample_id)) stop("Count sample IDs do not match curated metadata")
counts <- counts[, meta$sample_id, drop = FALSE]

out <- data.table::data.table(gene_id = rownames(counts))
out <- cbind(out, data.table::as.data.table(counts))
out_file <- file.path(path_data, "02_GEO_bulk", acc, "processed", paste0(acc, "_counts.tsv.gz"))
data.table::fwrite(out, out_file, sep = "\t", compress = "gzip")

qc <- data.table::data.table(
  sample_id = colnames(counts), library_size = colSums(counts),
  zero_fraction = colMeans(counts == 0), expressed_genes = colSums(counts > 0)
)
qc <- merge(qc, meta[, .(sample_id, group, site, include_primary)], by = "sample_id", all.x = TRUE)
safe_write_csv(qc, file.path(path_result, "02_bulk", acc, "T_GSE114192_count_QC.csv"))

audit <- data.frame(accession = acc, n_genes = nrow(counts), n_samples = ncol(counts),
                    count_min = min(counts), count_max = max(counts),
                    all_integer = all(counts == round(counts)),
                    md5 = unname(tools::md5sum(out_file)))
safe_write_csv(audit, file.path(path_result, "00_audit", "A09_GSE114192_count_matrix_audit.csv"))
write_log("GSE114192 counts prepared: ", nrow(counts), " genes x ", ncol(counts), " samples")

