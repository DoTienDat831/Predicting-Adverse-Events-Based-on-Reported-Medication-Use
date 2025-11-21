library(data.table)
library(readxl)
library(stringr)

# ---------------------------------------------------------
# PATHS
# ---------------------------------------------------------
input_dir  <- "/home/juan/Work/Midterm project/"
output_dir <- "/home/juan/Work/Midterm project/splited/"
dir.create(output_dir, FALSE)

# ---------------------------------------------------------
# READ + MERGE EXCEL
# ---------------------------------------------------------
excel_files <- list.files(input_dir, "List of Reports_.*\\.xlsx$", full.names = TRUE)
dt_all <- rbindlist(lapply(excel_files, \(f) as.data.table(read_excel(f))), fill = TRUE)
dt_all[, c("Case number","Report entry date") := NULL]

# ---------------------------------------------------------
# SPLIT INTO CHUNKS
# ---------------------------------------------------------
chunk_size <- 5000
n_chunks <- ceiling(nrow(dt_all) / chunk_size)

for (i in 1:n_chunks) {
  fwrite(
    dt_all[((i-1)*chunk_size+1):min(i*chunk_size, nrow(dt_all))],
    file.path(output_dir, paste0("part_", i, ".csv"))
  )
}

# ---------------------------------------------------------
# BUILD GLOBAL DICTIONARY
# ---------------------------------------------------------
files <- list.files(output_dir, "part_.*\\.csv$", full.names = TRUE)

extract_terms <- function(path) {
  df <- fread(path, select = c("Medicines reported as being taken","MedDRA reaction terms"))
  ing <- unique(str_trim(unlist(str_split(gsub("[()]", "", unlist(str_extract_all(df[[1]], "\\((.*?)\\)"))), "[;,/]"))))
  react <- unique(str_trim(unlist(str_split(df[[2]], "•"))))
  list(ing = ing[ing != ""], react = react[react != ""]))
}

dict <- lapply(files, extract_terms)
all_ingredients <- sort(unique(unlist(lapply(dict, `[[`, "ing"))))
all_reactions   <- sort(unique(unlist(lapply(dict, `[[`, "react"))))

# ---------------------------------------------------------
# PROCESS EACH FILE
# ---------------------------------------------------------
for (i in seq_along(files)) {
  df <- fread(files[i])
  
  # ---------------- AGE CLEAN ----------------
  df[, Age := as.numeric(`Age (years)`)]
  df[Age < 0 | is.na(Age), Age := mean(Age, na.rm = TRUE)]

  # ---------------- AGE GROUPS ----------------
  df[, AgeGroup := fcase(
        Age <= 17, "G1",
        Age <= 35, "G2",
        Age <= 50, "G3",
        Age <= 65, "G4",
        Age > 65, "G5"
      )]

  # Row ID for merge later
  df[, RowID := .I]

  # ---------------- GENDER ENCODING ----------------
  df[, Gender := fcase(
                     tolower(Gender)=="female", 0,
                     tolower(Gender)=="male", 1,
                     default = mean(
                       fifelse(tolower(Gender)=="female",0,
                               fifelse(tolower(Gender)=="male",1,NA_real_)),
                       na.rm=TRUE)
                   )]
  
  # ---------------- INGREDIENT MATRIX ----------------
  df[, temp_ing := str_extract_all(`Medicines reported as being taken`, "\\((.*?)\\)")]
  df[, susp_val := ifelse(grepl("Suspected", `Medicines reported as being taken`, TRUE), 2, 1)]
  
  dt_ing <- df[, .(raw = unlist(temp_ing)), .(RowID, susp_val)]
  dt_ing[, raw := gsub("[()]", "", raw)]
  dt_ing <- dt_ing[, .(ing = str_trim(unlist(str_split(raw, "[;,/]")))), .(RowID, susp_val)]
  dt_ing <- dt_ing[ing != ""]
  
  matrix_ing <- if (nrow(dt_ing))
    dcast(dt_ing, RowID ~ ing, value.var = "susp_val", fun.aggregate = max, fill = 0)
  else data.table(RowID = df$RowID)
  
  setnames(matrix_ing, names(matrix_ing)[-1], paste0("x ", names(matrix_ing)[-1]))

  # ---------------- REACTION MATRIX ----------------
  dt_react <- df[, .(react = str_trim(unlist(str_split(`MedDRA reaction terms`, "•")))), RowID]
  dt_react <- dt_react[react != ""]
  
  matrix_react <- if (nrow(dt_react))
    dcast(dt_react, RowID ~ react, fun.aggregate = length, fill = 0)
  else data.table(RowID = df$RowID)
  
  setnames(matrix_react, names(matrix_react)[-1], paste0("y ", names(matrix_react)[-1]))
  react_cols <- names(matrix_react)[-1]
  matrix_react[, (react_cols) := lapply(.SD, \(x) as.integer(x > 0)), .SDcols = react_cols]

  # ---------------------------------------------------------
  # AGE-RISK ENCODING (weight = adverse event rate)
  # ---------------------------------------------------------
  
  # Any AE?
  df[, AnyReact := as.integer(rowSums(matrix_react[.SD, .SDcols = react_cols]) > 0)]
  
  # Compute AE frequency for each age group
  age_weight_table <- df[, .(Risk = mean(AnyReact, na.rm = TRUE)), by = AgeGroup]

  # Add weight back to df
  df <- merge(df, age_weight_table, by = "AgeGroup", all.x = TRUE)
  
  # Final numeric encoded feature
  df[, `x AgeRisk` := Risk]

  # Drop helper columns
  df[, c("Risk", "AnyReact") := NULL]

  # ---------------- MERGE EVERYTHING ----------------
  final <- data.table(RowID = df$RowID,
                      `x Age` = df$Age,
                      `x Gender` = df$Gender,
                      `x AgeRisk` = df$`x AgeRisk`)

  final <- merge(final, matrix_ing, "RowID", TRUE)
  final <- merge(final, matrix_react, "RowID", TRUE)

  final[is.na(final)] <- 0
  final[, RowID := NULL]
  
  # ---------------- ALIGN COLUMNS (GLOBAL DICTIONARY) ----------------
  miss_ing   <- setdiff(paste0("x ", all_ingredients), names(final))
  miss_react <- setdiff(paste0("y ", all_reactions), names(final))
  
  if (length(miss_ing))   final[, (miss_ing) := 0]
  if (length(miss_react)) final[, (miss_react) := 0]
  
  # ---------------- WRITE ----------------
  fwrite(final, file.path(output_dir, paste0("processed_", i, ".csv")))
}
