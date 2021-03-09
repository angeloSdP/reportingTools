#' Transpose report table
#' @note At the moment this function only works when there is only one variable in report.
#' @param report_data table to transpose
#'
#' @return transposed table
#' @export report_transpose
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select
#'
#' @examples
#' df <- data.frame(V1=rbinom(500,1,0.4),V2=rbinom(500,1,0.5),
#'    group=gl(5,100))
#' df %>% report_nPctBin01(V1,groupVar = group) %>%
#'    report_transpose()
#'
report_transpose <- function(report_data){
  nc <- ncol(report_data)
  varName <- strsplit(names(report_data)[3],"=")[[1]][1]
  names(report_data)[3:(nc-1)] <- gsub(paste0(varName,"="),"",names(report_data)[3:(nc-1)])
  names(report_data)[3:(nc-1)] <- paste0(gsub("\n"," (",names(report_data)[3:(nc-1)]),")")
  report_data %>%
    select(-c(1,2,nc)) %>%
    pivot_longer(everything(),names_to=varName,
                 values_to = paste(report_data$Variable,"n (%)",sep="\n")) %>%
    mutate(`P-value`=c(report_data$`P-value`,rep("",nc-4)))
}
