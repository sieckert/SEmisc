#' Tidy batch processing output
#'
#' @description Creates a list of tidy output from batch processing based on post-run analysis in Shimadzu GSCMSsolution software v4.45. The function is a parser for output files from batch processing with the "Compound Quantitative Result" checkbox selected in the settings.
#' @import 'dplyr' 'reshape2' 'tibble'

#' @param file A string containing the name of the output file after batch processing. If not renamed, the default file name suggested by the software will be "ASCIIData.txt".
#' @param samples An integer giving the number of samples analyzed.
#' @param compounds An integer giving the number of compounds scored (as shown in the compound table of the method file).
#' @param results_file A string giving the file name to save the tidy output.
#' @param ref_ions An integer giving the number of reference ions (m/z) to add (default is three, maximum is five).
#' @param save_as_textfiles Logical. Should the cleaned output be saved in the working directory? `FALSE` by default.
#' @return A list of tidy output from the cleaned batch processing. The first element contains a tibble with sample names as rows and the compound areas as columns. The second element contains a tibble with additional information for each samples like base peak, peak height, retention index, and retention time. The third element contains a tibble with the summary of all detected compounds listed with their base peak, and the median of their retention index and retention time.
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @importFrom stats median
#' @export
#'
#' @examples
#' # example output file
#' \dontrun{
#' tidy_output <- parse_GCMSsolution(system.file(package="SEmisc",
#'                                     "extdata/example_batchprocessing_GCMSsolution.txt"),
#'                                     samples = 3,
#'                                     compounds = 80,
#'                                     results_file = "tidy_output",
#'                                     ref_ions = 3,
#'                                     save_as_textfiles = F)
#' tidy_output[[1]] # tidy output in wide format
#' tidy_output[[2]] # tidy output in long format
#' tidy_output[[3]] # tidy compound summary}

parse_GCMSsolution <- function(file,
                               samples,
                               compounds,
                               results_file,
                               ref_ions = 3,
                               save_as_textfiles = F){
  skip_rows = 8
  df_temp <- readr::read_file(file)
  df_temp <- gsub(" ", "_", df_temp)
  df_temp <- gsub("\t\r\n", "\r\n", df_temp)
  temp <- "temporary_helper_file_will_be_deleted.txt"
  write.table(df_temp, temp, row.names=T)
  number_samples = samples
  original_start = skip_rows
  var_start = skip_rows
  var_stop = compounds
  var_sum = var_start + var_stop + 1
  next_sample=1
  original_lab=2
  lab_start=2
  lab_stop=1
  lab_sum=lab_start + lab_stop + 1
  df_list = list()
  while (number_samples > 0){
    df_new <- read.table(temp,
                         comment.char = "",
                         header=T, sep="\t",
                         skip=var_start,
                         nrows=var_stop)
    if(ref_ions > 5){
      stop("Error in ref_ions: Only values between 0 and 5 allowed")
    } else {
      if(ref_ions == 0){
        df_new <- df_new[c("Name",
                           "Mass",
                           "Ret.Time",
                           "Area",
                           "Height",
                           "Ret._Index")]
        names(df_new) <- c("Compound",
                           "Base_Peak",
                           "Retention_Time",
                           "Peak_Area",
                           "Peak_Height",
                           "Retention_Index")
      } else {
      df_new <- df_new[c("Name",
                         "Mass",
                         "Ret.Time",
                         "Area",
                         "Height",
                         "Ret._Index",
                         sprintf("Ref.Ion%s_m.z", seq(1:ref_ions)))]
      names(df_new) <- c("Compound",
                         "Base_Peak",
                         "Retention_Time",
                         "Peak_Area",
                         "Peak_Height",
                         "Retention_Index",
                         sprintf("Reference_Ion_%s",seq(1:ref_ions)))
      }
      }
    label_var <- read.table(temp,
                            comment.char = "",
                            header=F, sep="\t",
                            skip=lab_start,
                            nrows=1)
    list_label <- gsub(".*\\\\", "", label_var[[2]])
    list_label_updated <- gsub("\\..*", "", list_label)
    df_list[[list_label_updated]] <- df_new
    var_start = (next_sample * var_sum) + original_start
    lab_start = (next_sample * var_sum) + original_lab
    next_sample=next_sample + 1
    number_samples = number_samples - 1
  }
  df_unlist <- dplyr::bind_rows(df_list, .id = "Sample")
  df_cast <- reshape2::dcast(df_unlist[c("Sample",
                                         "Compound",
                                         "Peak_Area")],#,,-c(3,4,6,7)],
                             Sample ~ Compound,
                             value.var = "Peak_Area",
                             fun.aggregate=NULL)
  df_cast[is.na(df_cast)] <- 0
  df_cast <- df_cast[, colSums(df_cast != 0) > 0]
  df_final <- df_cast
  if(save_as_textfiles == T){
    utils::write.table(df_final,
                       paste0(results_file, "_WideFormat.txt"), row.names=F, sep="\t")
    utils::write.table(df_unlist,
                       paste0(results_file, "_LongFormat.txt"), row.names=F, sep="\t")
    tidy_output_text <- "Tidy output generated and saved as text files.\n"
  } else {
    tidy_output_text <- "Tidy output generated, but not saved as text files (try save_as_textfiles = T).\n"
  }
  df_summary <- df_unlist %>%
    dplyr::filter(.data$Retention_Time != "Not_Identified") %>%
    dplyr::mutate(Retention_Time = as.numeric(.data$Retention_Time)) %>%
    dplyr::group_by(.data$Compound) %>%
    dplyr::summarise(Base_Peak = median(.data$Base_Peak, na.rm=T),
              Retention_Index = median(.data$Retention_Index, na.rm=T),
              Retention_Time = median(.data$Retention_Time, na.rm=T)) %>%
    dplyr::ungroup()
  cat(paste0(tidy_output_text,
             "In total, ", dim(df_summary)[1], " out of ", compounds, " compounds detected after batch processing of ",
             length(unique(df_unlist$Sample)), " samples."))
  file.remove(temp)
  return(list(tibble::tibble(df_final),
              tibble::tibble(df_unlist),
              tibble::tibble(df_summary)))
}
