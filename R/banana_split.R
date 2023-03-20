#' Tidy a string for splitting in half.
#' @description Cleans up a messy string so that it can be split at the separator of choice. The string may contain different separators. A two-element string is not ordered and returned as is.
#' @import 'stringi' 'stringr'
#'
#' @param x A string containing characters with different separators.
#' @param sep Character element at which the resulting string will be split. The default is a single space. Letters are allowed and will be retained in the original character element. Not working if multiple letters occur as separators.
#' @param pos Integer giving the separator position if it occurs twice in a row. Not working if separator occurs more than twice in a row.
#' @param join Character element joining the two resulting strings. The default is a single space.
#'
#' @return A string of two character elements joined by character element of choice.
#' @export
#'
#' @examples
#' # a string with different separators
#' string1 <- "blue_mono green mixed"; string1
#' banana_split(string1, sep = " ")
#'
#' # a string with a letter as separator that also occurs as part of the character elements
#' string2 <- "bluedmonodtypedgreendmixeddtype"; string2
#' banana_split(string2, sep = "d", join = "_")
#'
#' # first character is separator, keep second occurrence in a row (as part of element)
#' string3 <- "blue_monoggreengmixed"; string3
#' banana_split(string3, sep = "g", pos=1, join = " ")
#'
#' # second character is separator, keep first occurrence in a row (as part of element)
#' string4 <- "mono.yellowwmixed red"; string4
#' banana_split(string4, sep = "w", pos=2, join = " ")
#'
#' # will not change anything in the string (already two elements present)
#' string5 <- "1 2"; string5
#' banana_split(string5, sep = "_") # separator is ignored.

banana_split <- function(x,
                         sep = " ",
                         pos = 1,
                         join = " "){
  ifelse(pos==1,
         char_vec <- stringi::stri_split_regex(str=x,
                                               pattern=paste0(sep,"(?<!",sep,sep,")|(?<!",sep,")",sep,"(?!", sep,")|[[:punct:][:space:]]"),
                                               omit_empty=T)[[1]],
         char_vec <- stringi::stri_split_regex(str=x,
                                               pattern=paste0("(?<!",sep,sep,")",sep,"(?!",sep,"|",sep,"(?=",sep,"))|[[:punct:][:space:]]"),
                                               omit_empty=T)[[1]])
  char_vec <- char_vec[char_vec != ""]
  midpoint <- length(char_vec)/2
  new_char_vec <- ifelse(
    (midpoint == 0.5 | midpoint == 1),
    x,
    new_char_vec <- paste0(stringi::stri_paste(char_vec[1:midpoint], collapse=''),
                           join,
                           stringi::stri_paste(char_vec[(midpoint+1):length(char_vec)], collapse=''))
  )
  return(new_char_vec)
}
