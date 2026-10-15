#!/usr/bin/env Rscript
library(dplyr)
library(argparse)

#parse data
parser <- ArgumentParser(description= 'barcode_assignment')
parser$add_argument('--slx_id', help = "SLX_ID")
parser$add_argument('--sample', help = 'sample_index')
parser$add_argument('--cb_table', help= 'cb_table_full.txt')
parser$add_argument('--table_full_filtered', help= 'EnrichPCR_table_full_filtered.txt')
parser$add_argument('--output_file', '-o', help= 'Output file')
xargs<- parser$parse_args()


# Read the input file
barcode.10x <- read.table(xargs$cb_table, header = TRUE)
barcode.10x$CBC <- gsub("CB:Z:", "", barcode.10x$CBC)

# Aggregate by CBC and barcode_name
aggregate <- barcode.10x |>
  dplyr::group_by(CBC, barcode_name) |>
  dplyr::summarize(UMI_count = length(UMI), .groups = 'drop')

# Get total UMIs per CBC
CBC_counts <- aggregate |>
  dplyr::group_by(CBC) |>
  dplyr::summarize(total_UMI = sum(UMI_count), .groups = 'drop')

# Merge counts and calculate fraction
barcode.10x.table <- dplyr::left_join(aggregate, CBC_counts, by = "CBC") |>
    dplyr::mutate(fraction = UMI_count/total_UMI) |>
    dplyr::filter(UMI_count >=2, fraction > 0.5)


# Save the filtered table to disk
filtered_fn=paste0("results/",xargs$slx_id,"/",xargs$sample,"/",xargs$sample, ".10x.table.csv")
write.csv(barcode.10x.table, filtered_fn, row.names = FALSE)



#summarising and filtered the barcode read counts for the enrichment library


enrich <- read.table(xargs$table_full_filtered, header = T)

colnames(enrich) <- c("count", "CBC", "UMI", "barcode_name")

enrich <- enrich |> 
  dplyr::mutate(barcode_number = lengths(strsplit(barcode_name, ","))) |>
  dplyr::filter(barcode_number <= 2) |>
  dplyr::mutate(barcode_score = ifelse(barcode_number == 2, 0.5, 1),
                barcode_name = strsplit(barcode_name, ",")) |>
  tidyr::unnest(barcode_name)



aggregate.enr <- enrich |> dplyr::group_by(CBC,barcode_name) |> 
  dplyr::summarize(barcode_score=sum(barcode_score), .groups = 'drop')

CBC_counts.enr <- aggregate.enr |> dplyr::group_by(CBC) |> 
  dplyr::summarize(total_score=sum(barcode_score), .groups = 'drop')

enrich <- dplyr::left_join(aggregate.enr, CBC_counts.enr, by = "CBC") |>
    dplyr::mutate(fraction = barcode_score/total_score)

enrich.top <- enrich |> dplyr::group_by(CBC) |> dplyr::slice_max(n = 1, order_by = fraction)
enrich.second <- enrich |> dplyr::group_by(CBC) |> dplyr::arrange(desc(fraction)) |> dplyr::slice(2)

enrich <- dplyr::full_join(enrich.top, enrich.second, by = "CBC") |>
    dplyr::mutate(ratio = fraction.x / fraction.y)

filtered <- enrich |> dplyr::filter(is.na(ratio) | ratio >= 2) |>
    dplyr::select(CBC, barcode_name.x, barcode_score.x, total_score.x, fraction.x, fraction.y, ratio) |>
    dplyr::rename(barcode_name=barcode_name.x, barcode_score = barcode_score.x, total_UMI = total_score.x,
                  fraction.first = fraction.x, fraction.second = fraction.y)

#enrich <- sample_metadata[which(sample_metadata$index == sample),2]

merged <- dplyr::full_join(barcode.10x.table,filtered,by="CBC") |>
  dplyr::mutate(clone = dplyr::case_when(
      is.na(barcode_name.x) ~ barcode_name.y,                  # If .x is NA, take .y
      is.na(barcode_name.y) ~ barcode_name.x,                  # If .y is NA, take .x
      barcode_name.x == barcode_name.y ~ barcode_name.x,       # If they match, keep either
      TRUE ~ "non-match")) |>                                 # Otherwise, label as "non-match"
  dplyr::mutate(origin =  dplyr::case_when(
      is.na(barcode_name.x) ~ "pcr-enrichment",
      is.na(barcode_name.y) ~ "10X",
      barcode_name.x == barcode_name.y ~ "both",
      TRUE ~ "non-match")) |>
  dplyr::filter(!(origin == "10X" & UMI_count < 5), !(origin=="pcr-enrichment" & barcode_score <30)) |>
  dplyr::mutate(CBC = paste0(CBC,"-1"))

write.table(merged, xargs$output_file, sep = "\t", row.names = F, quote = F)

