library(data.table)
library(readxl)
library(stringr)

# ---------------------------------------------------------
# PATHS
# ---------------------------------------------------------
input_dir  <- "/home/juan/Work/Midterm project/"
output_dir <- "/home/juan/Work/Midterm project/splited/"
if(!dir.exists(output_dir)) dir.create(output_dir, FALSE)

# ---------------------------------------------------------
# READ + MERGE EXCEL
# ---------------------------------------------------------
excel_files <- list.files(input_dir, "List of Reports_.*\\.xlsx$", full.names = TRUE)
dt_all <- rbindlist(lapply(excel_files, \(f) as.data.table(read_excel(f))), fill = TRUE)

# Remove unnecessary columns
dt_all[, c("Case number","Report entry date") := NULL]

# ---------------------------------------------------------
# PRE-CALCULATE GLOBAL AGE BREAKS
# ---------------------------------------------------------
temp_age <- suppressWarnings(as.numeric(dt_all$`Age (years)`))
temp_age[is.na(temp_age) | temp_age < 0] <- mean(temp_age[temp_age >= 0], na.rm = TRUE)

age_breaks <- unique(quantile(temp_age, probs = seq(0, 1, 0.2), na.rm = TRUE))
age_breaks[1] <- -Inf
age_breaks[length(age_breaks)] <- Inf

rm(temp_age)

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

rm(dt_all)
gc()

# ---------------------------------------------------------
# BUILD GLOBAL DICTIONARY
# ---------------------------------------------------------
files <- list.files(output_dir, "part_.*\\.csv$", full.names = TRUE)

extract_terms <- function(path) {
  df <- fread(path, select = c("Medicines reported as being taken","MedDRA reaction terms"))
  
  ing_raw <- unlist(str_extract_all(df[[1]], "\\((.*?)\\)"))
  ing_clean <- gsub("[()]", "", ing_raw)
  ing <- unique(str_trim(unlist(str_split(ing_clean, "[;,/]"))))
  
  react <- unique(str_trim(unlist(str_split(df[[2]], "•"))))
  
  list(ing = ing[ing != "" & !is.na(ing)], react = react[react != "" & !is.na(react)])
}

dict <- lapply(files, extract_terms)
all_ingredients <- sort(unique(unlist(lapply(dict, `[[`, "ing"))))
all_reactions   <- sort(unique(unlist(lapply(dict, `[[`, "react"))))

# ---------------------------------------------------------
# PROCESS EACH FILE
# ---------------------------------------------------------
for (i in seq_along(files)) {
  cat(sprintf("Processing file %d of %d...\n", i, length(files)))
  df <- fread(files[i])
  
  # --- 1. FIX AGE & CATEGORIZE ---
  df[, age_num := suppressWarnings(as.numeric(`Age (years)`))]
  
  global_mean_age <- mean(df$age_num, na.rm = TRUE)
  if(is.nan(global_mean_age)) global_mean_age <- 0
  
  # Impute missing age with Mean
  df[is.na(age_num) | age_num < 0, age_num := global_mean_age]
  
  # Create Groups (For Linear Models)
  df[, `x AgeGroup` := cut(age_num, breaks = age_breaks, labels = FALSE, include.lowest = TRUE)]
  
  # --- 2. FIX GENDER (PROPORTIONAL RANDOM FILL) ---
  df[, Gender := fcase(
    tolower(Gender) == "female", 0,
    tolower(Gender) == "male", 1,
    default = NA_real_
  )]
  
  # Calculate probability of being Male in this file
  prob_male <- mean(df$Gender, na.rm = TRUE)
  if(is.nan(prob_male)) prob_male <- 0.5
  
  # Identify missing rows
  missing_idx <- which(is.na(df$Gender))
  n_missing <- length(missing_idx)
  
  # Fill randomly respecting the ratio (Preserves distribution for NN/Lasso)
  if(n_missing > 0) {
    random_fills <- sample(c(0, 1), size = n_missing, replace = TRUE, prob = c(1 - prob_male, prob_male))
    df[missing_idx, Gender := random_fills]
  }
  
  # --- 3. PROCESS INGREDIENTS ---
  df[, RowID := .I]
  df[, temp_ing := str_extract_all(`Medicines reported as being taken`, "\\((.*?)\\)")]
  
  df[, susp_val := ifelse(grepl("Suspected", `Medicines reported as being taken`, ignore.case = TRUE), 2, 1)]
  
  dt_ing <- df[, .(raw = unlist(temp_ing)), by = .(RowID, susp_val)]
  dt_ing[, raw := gsub("[()]", "", raw)]
  dt_ing <- dt_ing[, .(ing = str_trim(unlist(str_split(raw, "[;,/]")))), by = .(RowID, susp_val)]
  dt_ing <- dt_ing[ing != ""]
  
  if (nrow(dt_ing) > 0) {
    matrix_ing <- dcast(dt_ing, RowID ~ ing, value.var = "susp_val", fun.aggregate = max, fill = 0)
  } else {
    matrix_ing <- data.table(RowID = df$RowID)
  }
  
  if(ncol(matrix_ing) > 1) {
    setnames(matrix_ing, names(matrix_ing)[-1], paste0("x ", names(matrix_ing)[-1]))
  }
  
  # --- 4. PROCESS REACTIONS ---
  dt_react <- df[, .(react = str_trim(unlist(str_split(`MedDRA reaction terms`, "•")))), by = RowID]
  dt_react <- dt_react[react != ""]
  
  if (nrow(dt_react) > 0) {
    matrix_react <- dcast(dt_react, RowID ~ react, fun.aggregate = length, fill = 0)
  } else {
    matrix_react <- data.table(RowID = df$RowID)
  }
  
  if(ncol(matrix_react) > 1) {
    setnames(matrix_react, names(matrix_react)[-1], paste0("y ", names(matrix_react)[-1]))
    cols_y <- names(matrix_react)[-1]
    matrix_react[, (cols_y) := lapply(.SD, \(x) as.integer(x > 0)), .SDcols = cols_y]
  }
  
  # --- 5. MERGE (UPDATED: KEEP RAW AGE) ---
  # We keep 'x Age' (raw number) for XGBoost/RF
  # We keep 'x AgeGroup' (bins) for analysis/Linear trends
  base_df <- df[, .(RowID, `x Age` = age_num, `x AgeGroup`, `x Gender` = Gender)]
  
  final <- merge(base_df, matrix_ing, by = "RowID", all = TRUE)
  final <- merge(final, matrix_react, by = "RowID", all = TRUE)
  
  final[is.na(final)] <- 0
  final[, RowID := NULL]
  
  # --- 6. FILL MISSING GLOBAL COLUMNS ---
  miss_ing   <- setdiff(paste0("x ", all_ingredients), names(final))
  miss_react <- setdiff(paste0("y ", all_reactions), names(final))
  
  if (length(miss_ing))   final[, (miss_ing) := 0]
  if (length(miss_react)) final[, (miss_react) := 0]
  
  fwrite(final, file.path(output_dir, paste0("processed_", i, ".csv")))
}