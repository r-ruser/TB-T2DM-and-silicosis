source(file.path("R", "00_config.R"), encoding = "UTF-8")
assert_packages(c("data.table", "readxl"))

acc <- "GSE326212"
xlsx <- file.path(path_data, "03_scRNA", acc, "raw", paste0(acc, "_Branchett_et_al_Table_1.xlsx"))
result_dir <- file.path(path_result, "04_scRNA", acc)
if (!file.exists(xlsx)) stop("Local Table 1 workbook is missing")
sheets <- readxl::excel_sheets(xlsx)
sheet_audit <- data.table::rbindlist(lapply(seq_along(sheets), function(i) {
  x <- readxl::read_excel(xlsx, sheet = sheets[i])
  safe_name <- gsub("[^A-Za-z0-9]+", "_", sheets[i])
  safe_write_csv(data.table::as.data.table(x), file.path(result_dir, "source_data",
    paste0("SD03_GSE326212_Table1_sheet", sprintf("%02d", i), "_", safe_name, ".csv")))
  data.table::data.table(sheet_index = i, sheet_name = sheets[i], rows = nrow(x), columns = ncol(x),
    column_names = paste(names(x), collapse = " | "))
}))
safe_write_csv(sheet_audit, file.path(result_dir, "tables", "T04_GSE326212_Table1_sheet_audit.csv"))
write_log("GSE326212 Table 1 audit completed: ", paste(sheets, collapse = ", "))
