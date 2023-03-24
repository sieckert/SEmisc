#' Standardize data based on internal standard
#'
#' @description This function standardizes a given data frame based on an internal standard (ISTD).
#' @import 'dplyr' 'tibble' 'stringi'

#' @param data A data frame to be standardized.
#' @param remove_ISTD Logical. Should the column containing the ISTD be removed after standardization? `FALSE` by default.
#' @param show_original Logical. Should the original data be returned additionally? `TRUE` by default.
#' @param save_as_textfile Logical. Should the standardized output be saved in the working directory? `FALSE` by default.
#' @return A list of tidy data frames of standardized data and, optionally, the original data. Additionally, the list returns information on what samples where ISTD was not detected.
#' @importFrom utils write.table
#' @importFrom sjmisc is_empty
#' @importFrom here here
#' @importFrom stringi stri_rand_strings
#' @importFrom purrr set_names
#' @importFrom wakefield name
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # generate example data
#' set.seed(123)
#' df <- data.frame(replicate(10, sample(1:10000, 100, T))) %>%
#' set_names(name(10)) %>%
#' mutate(Sample = stri_rand_strings(100, 5),
#'        Treatment = sample(LETTERS[1:3], 100, T),
#'        compound_ISTD = sample(c(NA, 100:120), 100, T),
#'        .before = everything(.)) %>%
#'        tibble(); df
#'
#' # standardize example data
#' standardised_output <- ISTD(data,
#'                             remove_ISTD = F,
#'                             show_original = T,
#'                             save_as_textfile = F)
#' dframe_ISTD$standardised_data
#' dframe_ISTD$original_data
#' dframe_ISTD$no_ISTD
#' }

ISTD <- function(data,
                 remove_ISTD = F,
                 show_original = T,
                 save_as_textfile = F){
  meta <- data %>%
    select(-where(is.numeric))
  if(remove_ISTD == T){
    df <- data %>%
      select_if(is.numeric) %>%
      # select(-.data$Sample) %>%
      relocate(contains("ISTD"))
    m_ISTD <- as.matrix(df)
    m_standardised <- apply(m_ISTD, 1, function(x) (x)/x[1])
    df_standardised <- as.data.frame(t(m_standardised)) %>%
      mutate(meta, .before=everything(.data)) %>%
      select(-contains("ISTD")) %>%
      as_tibble()
  } else{
    df <- data %>%
      select_if(is.numeric) %>%
      # select(-.data$Sample) %>%
      relocate(contains("ISTD"))
    m_ISTD <- as.matrix(df)
    m_standardised <- apply(m_ISTD, 1, function(x) (x)/x[1])
    df_standardised <- as.data.frame(t(m_standardised)) %>%
      mutate(meta, .before=everything(.data)) %>%
      as_tibble()
  }
  if(save_as_textfile == T){
    utils::write.table(x = df_standardised,
                       file = here::here(paste0(format(Sys.time(), "%Y%m%d"), "_", "standardised_data.txt")),
                       row.names=F, sep="\t")
    cat("Standardised output generated and saved as text file.\n")
  } else {
    cat("Standardised output generated but not saved as text file.\n")
  }
  if(is_empty(data %>%
     select(contains("ISTD")) %>%
     filter(if_any(matches("ISTD"), ~ . == 0) |
            if_any(matches("ISTD"), ~ is.na(.data))) %>%
     pull(),
     all.na.empty = F) == F){
    cat("ISTD not detected in the following samples:",
        paste0(no_ISTD <- data %>%
          select(.data$Sample,
                 contains("ISTD")) %>%
          filter(if_any(matches("ISTD"), ~ . == 0) |
                   if_any(matches("ISTD"), ~ is.na(.data))) %>%
          select(.data$Sample) %>%
          pull()
    ), "\n\n")
  }
  if(show_original == T){
    return(list(standardised_data = df_standardised,
                original_data = data,
                no_ISTD = no_ISTD))
  } else {
    return(list(standardised_data = df_standardised,
                no_ISTD = no_ISTD))
  }
}
