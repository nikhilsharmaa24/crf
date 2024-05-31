#' Analyze FASTA Alignment
#'
#' This function reads sequences from a FASTA file, performs an alignment using the specified method,
#' and returns the alignment as a matrix.
#'
#' @param sequence_file Path to the FASTA file containing the sequences.
#' @param alignment_method The alignment method to use (e.g., "ClustalOmega").
#' @return A matrix representing the aligned sequences.
#' @import msa
#' @importFrom Biostrings readAAStringSet
#' @importFrom Biostrings as.matrix
#' @export
analyze_fasta_alignment <- function(sequence_file, alignment_method) {
  mySequences <- readAAStringSet(sequence_file)
  myAlignment <- msa(mySequences, alignment_method)
  dat_aligned <- Biostrings::as.matrix(myAlignment)
  row_names <- rownames(dat_aligned)
  rownames(dat_aligned) <- sub("^(\\S+).*", "\\1", row_names)
  return(dat_aligned)
}

#' Count Unique Characters
#'
#' This function counts the number of unique characters in a given column, excluding dashes ("-").
#'
#' @param column A vector representing a column from a matrix.
#' @return The count of unique characters excluding dashes.
#' @importFrom plyr ldply
#' @export
count_unique_characters <- function(column) {
  unique_chars <- unique(column)
  unique_chars <- unique_chars[unique_chars != "-"]
  return(length(unique_chars))
}

#' Find Similar Rows
#'
#' This function finds rows in a matrix that have similar values, ignoring dashes.
#'
#' @param matrix_data A matrix of aligned sequences.
#' @return A list of rows with similar values.
#' @importFrom plyr ldply
#' @export
find_similar_rows <- function(matrix_data) {
  row_concatenated <- apply(matrix_data, 1, paste, collapse = "")
  target_strings <- unique(row_concatenated)
  target_strings <- target_strings[grep("-", target_strings, fixed = TRUE, invert = TRUE)]
  similar_rows_list <- lapply(target_strings, function(target) {
    rownames(matrix_data)[row_concatenated == target]
  })
  similar_rows_list <- setNames(similar_rows_list, target_strings)
  return(similar_rows_list)
}

#' Analyze Fragments
#'
#' This function analyzes fragments of the aligned data and returns statistics on the fragments.
#'
#' @param dat_aligned A matrix of aligned sequences.
#' @param fragment_size The size of each fragment to analyze.
#' @return A list containing percentage of alphabetic characters and a final table with fragment analysis.
#' @import dplyr
#' @importFrom plyr ldply
#' @importFrom purrr imap_dfr
#' @export
analyze_fragments <- function(dat_aligned, fragment_size) {
  fragment_list <- list()
  for (start_col in 1:(ncol(dat_aligned) - fragment_size + 1)) {
    end_col <- start_col + fragment_size - 1
    fragment <- dat_aligned[, start_col:end_col]
    fragment_name <- paste(start_col, ":", end_col, sep = "")
    fragment_list[[fragment_name]] <- fragment
  }
  perc_alpha_list <- list()
  unique_alpha_list <- list()
  result_list <- list()
  for (fragment_name in names(fragment_list)) {
    fragment <- fragment_list[[fragment_name]]
    alphabet_count <- sum(grepl("-", fragment))
    total_count <- length(fragment)
    percentage_alphabet <- round((alphabet_count / total_count) * 100, 2)
    perc_alpha_list[[fragment_name]] <- percentage_alphabet
    unique_alpha_list[[fragment_name]] <- sum(apply(fragment, 2, count_unique_characters))
    result_list[[fragment_name]] <- find_similar_rows(fragment)
  }
  alpha_list <- ldply(perc_alpha_list)
  names(alpha_list) <- c('fragment', 'perc_dash')
  unique_list <- ldply(unique_alpha_list)
  names(unique_list) <- c('fragment', 'total_unique_character')
  fragment_result <- left_join(alpha_list, unique_list)
  fragment_result %>%
    mutate(fragment = gsub(":", "_", fragment)) %>%
    arrange(perc_dash, total_unique_character) -> fragment_result
  extracted_data <- result_list %>%
    imap_dfr(~ {
      if (length(.x) > 0) {
        sublist <- .x
        if (!is.list(sublist)) {
          sublist <- list(sublist)
        }
        data.frame(
          Fragment_coverage = gsub(":", "_", .y),
          Fragment = names(sublist),
          Count = sapply(sublist, function(x) length(x)),
          Identifier = sapply(sublist, function(x) paste(x, collapse = " ; "))
        )
      } else {
        NULL
      }
    }) %>%
    filter(!is.null(Fragment))
  rownames(extracted_data) <- NULL
  extracted_data <- extracted_data %>%
    arrange(match(Fragment_coverage, fragment_result$fragment))
  return(list(PercentAlpha = fragment_result, Final_table = extracted_data))
}
