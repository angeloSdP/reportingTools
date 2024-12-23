#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Arithmetic mean ± sd
#'
#' @param x A vector of numeric values.
#' @param digits Number of decimal digits for showing the mean and sd of x.
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#'
#' @return A string with the value of mean ± sd of the vector x.
#' @export meanSd
#'
#' @examples
#' meanSd(rnorm(50))
#'
meanSd <- function(x,digits=2,na.rm=TRUE){
  m <- mean(x,na.rm=na.rm)
  s <- sd(x,na.rm=na.rm)
  if (is.na(m)) "-" else{
    format <- paste("%.",digits,"f ± %.",digits,"f",sep="")
    sprintf(format,m,s)
  }
}

#' Median and Interquartile Range
#'
#' @param x A vector of numeric values.
#' @param digits Number of decimal digits for showing the median and Interquartile Range of x.
#' @param probs Pair of quantiles to be computed around the median. Default are q25 and q75.
#' @param roundFrom If median is greater than this value, median and IQR are rounded to zero decimals even if digits>0.
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#'
#' @return A character with the value of the median and Interquartile Range.
#' @export medianIQR
#'
#' @examples
#' medianIQR(rnorm(50))
#' medianIQR(runif(500,0,3))
#'
medianIQR <- function(x,digits=2,probs=c(0.25,0.75),roundFrom=30,na.rm=TRUE){
  if ((!na.rm & any(is.na(x)))|all(is.na(x))) "-" else{
    q=quantile(x,na.rm=TRUE,probs=c(0.5,probs))
    format <- ifelse(q[1] >= roundFrom, "%.0f (%.0f; %.0f)",
                     paste("%.",digits,"f (%.",digits,"f; %.",digits,"f)",sep=""))
    sprintf(format,q[1],q[2],q[3])
  }
}

#' Absolute and relative frequencies of 1's in a binary 0-1 variable.
#'
#' @param x A vector of numeric values containing only 1's and 0's. Function returns an error if there are other values in x.
#' @param digits Number of decimal digits for showing relative frequency.
#'
#' @return Number and percent of 1's in a numeric vector x which contains only 1's and 0's.
#' @export nPctBin01
#'
#' @examples
#' nPctBin01(rbinom(100,1,0.4))
#' nPctBin01(sample(c(0,1),50,replace=TRUE))
nPctBin01 <- function(x,digits=1){
  if (!all(x%in%c(0,1,NA))) stop(paste(deparse(substitute(x)),"is not a binary (0,1) variable"))
  format <- paste("%d (%.",digits,"f%%)",sep="")
  n1 <- sum(x,na.rm=TRUE)
  n <- length(na.omit(x))
  sprintf(format,n1,100*n1/n)
}

#' Testing if a variable follows a normal distribution using Shapiro-Wilk test.
#'
#' @param x Variable whose normality is to be tested.
#' @param alpha Significance level for testing normality.
#'
#' @return TRUE (if normality can be accepted) or FALSE (if normality is rejected).
#' @export isNormal
#'
#' @examples
#' isNormal(rnorm(100))
#' isNormal(runif(100))
isNormal <- function(x, alpha=0.05){
  if (!is.numeric(x)) FALSE else {
    stp <- shapiro.test(x)$p.value
    if (stp<alpha) FALSE else TRUE
  }
}

#' Describe a variable depending on normality
#'
#' @param x variable to be described.
#' @param digits decimal digits in the result.
#' @param probs quantiles to be used in case of non normality.
#' @param alpha significance level for testing normality (0.05 by default).
#' @param na.rm Should missing values be considered?. If TRUE and there are missing values present, the funcion returns NA.
#'
#' @return Mean ± sd of variable if normal or Median and Interquartile range if not normal.
#' @export describe
#'
#' @examples
#' describe(rnorm(100))
#' describe(runif(100))
describe <- function(x, digits=2,probs=c(0.25,0.75),alpha=0.05,na.rm=TRUE){
  if (isNormal(x,alpha))
    meanSd(x,digits,na.rm) else
      medianIQR(x,digits,probs,na.rm)
}

#' Format p-values
#'
#' @param x p-value to be formatted.
#' @param pvdigits number of digits in the result.
#'
#' @return p-value rounded to pvdigits decimal digits, or an expression like "<0.0001" when x<1e-pvdigits.
#' @export formatPval
#'
#' @examples
#' formatPval(0.0000001,pvdigits=3)
#' formatPval(0.345251,pvdigits=4)
formatPval <- function(x, pvdigits=4){
  if (pvdigits<=0) pvdigits=1
  ifelse(is.na(x), "-",
         ifelse(x<10^(-pvdigits),
                paste("<.",paste(rep(0,pvdigits-1),collapse=""),1,sep=""),
                sprintf(paste("%.",pvdigits,"f",sep=""),x)))
}

#' P-value of analysis of Variance
#'
#' @param y Response variable.
#' @param g Grouping variable.
#' @param pvdigits number of digits for the p-value.
#'
#' @return p-value of t-test when comparing two groups or p-value of the analysis
#' of variance test aov(y~g) when comparing more than two groups. For comparing two
#' groups t.test is preferred as it allows for different variances in each group.
#' @export
#'
#' @examples
#' df <- data.frame(g=gl(3,10),y=rnorm(30))
#' aov_t_pval(df$y,df$g,pvdigits=3)
aov_t_pval <- function(y, g, pvdigits=4){
  ng <- length(unique(g))
  if (ng==2)
    p=tryCatch(t.test(y~g)$p.value, error=function(e) NA)
  else
    p=tryCatch(summary(aov(y~g))[[1]][["Pr(>F)"]][1],
             error=function(e) NA)
  formatPval(p, pvdigits)
}

#' P-value of Wilcoxon of Kruskal-Wallis test
#'
#' @param y Response variable.
#' @param g Grouping variable.
#' @param pvdigits number of digits for the p-value.
#'
#' @return p-value of Wilcoxon test (when comparing two groups) or Kruskal-Wallis test
#' (when comparing more than two groups).
#' @export
#'
#' @examples
#' df <- data.frame(g=gl(3,10),y=runif(30))
#' kruskal_Wilcox_pval(df$y,df$g,pvdigits=3)
kruskal_Wilcox_pval <- function(y, g, pvdigits=4){
  ng <- length(unique(g))
  if (ng==2)
    p=tryCatch(wilcox.test(y~g,exact=FALSE)$p.value, error=function(e) NA)
  else
    p=tryCatch(kruskal.test(y~g)$p.value, error=function(e) NA)
  formatPval(p, pvdigits)
}

#' P-value of chi-squared contingency table test
#'
#' @param y A numeric vector (or factor). Row variable in contingency table.
#' @param g A numeric vector (or factor). Column variable in contingency table.
#' @param pvdigits number of digits for the p-value
#'
#' @details If conditions for the validity of chi squared test are not met (expected
#' frequencies less than 5) a Fisher exact test is performed instead.
#'
#' @return p-value of the chi squared contingency table test chisq.test(table(y,g))
#' @export
#'
#' @examples
#' df <- data.frame(g=rbinom(100,4,0.5),y=rnorm(sample(c(1:2),100,replace=TRUE)))
#' chisq_pval(df$y,df$g,pvdigits=3)
chisq_pval <- function(y, g, pvdigits=4){
  tb <- table(y,g)
  if (nrow(tb)==1) p <- NA else{
    me <- suppressWarnings(min(chisq.test(tb)$expected))
    p <- if (me>5) chisq.test(tb,correct=FALSE)$p.value else{
      if (nrow(tb)==2)
        tryCatch(fisher.test(tb)$p.value,
                 error=function(e)
                   fisher.test(tb,simulate.p.value=TRUE,B=1e5)$p.value) else
                     fisher.test(tb,simulate.p.value=TRUE,B=1e5)$p.value
    }
  }
  formatPval(p, pvdigits)
}

#' Report table of means ± sd
#'
#' @description report_meanSd builds a table with the overall means ±
#' standard deviations of one or more variables. A grouping variable can also be
#' specified, in which case mean±sd are also calculated for each group and the
#' p-value of an anova test is displayed to decide if there are significant
#' differences between the means of the different groups.
#' @param data data frame or tibble which contains the data.
#' @param summary_vars Variable or variables whose mean ± sd is to be calculated.
#' @param groupVar Grouping variable.
#' @param digits Number of decimal digits for the results.
#' @param pvdigits Number of decimal digits for the p-value of anova test.
#' @param na.rm Should NA values be removed (possible values are TRUE or FALSE)
#'
#' @return A table with the overall mean±sd of the variables and, if a grouping
#' variable is specified, the means±sd by group and the p-value for comparing the
#' means of the groups, by using t-test when there are two groups (maybe with
#' different variances), and anova test when there are more than two groups (equal
#' variances is assumed).
#' @importFrom dplyr mutate filter group_by count add_row summarize as_label pull full_join select enquo
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom forcats fct_inorder
#' @export report_meanSd
#'
#' @examples
#' df <- data.frame(x=rnorm(100,10,3),y=rnorm(100,50,8), z=runif(100,20,30),g=sample(c("Yes","No"),100,replace=TRUE))
#' df %>%
#' report_meanSd(c(x,y,z))  # Only overall means and sd of variables x, y and z
#' df %>%
#' report_meanSd(c(x,y,z), groupVar=g, digits=1)
report_meanSd <- function(data, summary_vars, groupVar=NULL, digits=2, pvdigits=4, na.rm=TRUE) {
  if(rlang::quo_is_null(enquo(groupVar))){
    data %>%
      pivot_longer({{summary_vars}},names_to="Variable",values_to="value") %>%
      mutate(Variable=fct_inorder(Variable)) %>%
      group_by(Variable) %>%
      summarize(`mean ± sd`=meanSd(value, digits=digits, na.rm=na.rm))
  } else{
    data <- data %>% filter(!is.na({{groupVar}}))
    n <- data %>%
      count({{groupVar}}) %>%
      add_row(summarize(.,n=sum(n)),.before=1) %>%
      pull(n)
    data <- data %>%
      pivot_longer({{summary_vars}},names_to="Variable",values_to="value")%>%
      mutate(Variable=fct_inorder(Variable))
    overall <- data %>%
      group_by(Variable) %>%
      summarize(Overall=meanSd(value, digits=digits, na.rm=na.rm))
    perGroup <- data %>%
      group_by(Variable,{{groupVar}}) %>%
      summarize(value=meanSd(value, digits=digits, na.rm=na.rm)) %>%
      pivot_wider(names_from={{groupVar}},values_from=value)
    gvName <- as_label(enquo(groupVar))
    names(perGroup)[-1]=paste(gvName,names(perGroup)[-1],sep="=")
    P <- data %>%
      group_by(Variable) %>%
      summarize(`P-value`=aov_t_pval(value,{{groupVar}},pvdigits=pvdigits))
    summaryTable <- suppressMessages(overall %>% full_join(perGroup) %>% full_join(P))
    nc <- c(1,ncol(summaryTable))
    names(summaryTable)[-nc] <- paste(names(summaryTable)[-nc],n,sep="\nN =")
    summaryTable
  }
}

#' Report table of medians ± IQR (Interquartile Range)
#'
#' @description report_medianIQR builds a table with the overall medians and (by
#' default) quartiles 25 and 75 of one or more variables. A grouping variable can also be
#' specified, in which case medians and quartiles are also calculated for each group and the
#' p-value of a Kruskal-Wallis (more than 2 groups) or Wilcoxon test (two groups) is displayed
#' to test the null hypothesis that the location parameters of the distribution of the variable
#' are the same in each group.
#' @param data data frame or tibble which contains the data.
#' @param summary_vars Variable or variables whose median and quartiles is to be calculated.
#' @param groupVar Grouping variable.
#' @param digits Number of decimal digits for the results.
#' @param pvdigits Number of decimal digits for the p-value of the test.
#' @param roundFrom If median is greater than this value, median and IQR are rounded to zero decimals even if digits>0.
#' @param probs Pair of quantiles to be computed around the median. Default are q25 and q75.
#' @param na.rm Should NA values be removed (possible values are TRUE or FALSE)
#'
#' @return A table with the overall median and quartiles of the variables and, if a grouping
#' variable is specified, the medians and quartiles by group and the p-value of the wilcoxon (two groups) or kruskal
#' (more than two groups) test for comparing location parameters.
#' @importFrom dplyr mutate filter group_by count add_row summarize as_label pull full_join select enquo
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom forcats fct_inorder
#' @export report_medianIQR
#'
#' @examples
#' df <- data.frame(x=rnorm(100,10,3),y=rnorm(100,50,8), z=runif(100,20,30),g=sample(c("Yes","No"),100,replace=TRUE))
#' df %>%
#' report_medianIQR(c(x,y,z))  # Only overall medians and quartiles of variables x, y and z
#' df %>%
#' report_medianIQR(c(x,y,z), groupVar=g, digits=1)
report_medianIQR <- function(data, summary_vars, groupVar=NULL, digits=2, roundFrom=30,
                             probs=c(0.25,0.75), pvdigits=4, na.rm=TRUE){
  if(rlang::quo_is_null(enquo(groupVar))){
    data %>%
      pivot_longer({{summary_vars}},names_to="Variable",values_to="value") %>%
      mutate(Variable=fct_inorder(Variable)) %>%
      group_by(Variable) %>%
      summarize(`median (IQR)`= medianIQR(value, digits=digits, probs=probs,
                                          roundFrom=roundFrom, na.rm=na.rm))
  } else{
    data <- data %>% filter(!is.na({{groupVar}}))
    n <- data %>%
      count({{groupVar}}) %>%
      add_row(summarize(.,n=sum(n)),.before=1) %>%
      pull(n)
    data <- data %>%
      pivot_longer({{summary_vars}},names_to="Variable",values_to="value") %>%
      mutate(Variable=fct_inorder(Variable))
    overall <- data %>%
      group_by(Variable) %>%
      summarize(Overall=medianIQR(value, digits=digits, probs=probs,
                                  roundFrom=roundFrom, na.rm=na.rm))
    perGroup <- data %>%
      group_by(Variable,{{groupVar}}) %>%
      summarize(value=medianIQR(value, digits=digits, probs=probs,
                                roundFrom=roundFrom, na.rm=na.rm)) %>%
      pivot_wider(names_from={{groupVar}},values_from=value)
    gvName <- as_label(enquo(groupVar))
    names(perGroup)[-1]=paste(gvName,names(perGroup)[-1],sep="=")
    P <- data %>%
      group_by(Variable) %>%
      summarize(`P-value`=kruskal_Wilcox_pval(value,{{groupVar}}, pvdigits=pvdigits))
    summaryTable <- suppressMessages(overall %>% full_join(perGroup) %>% full_join(P))
    nc <- c(1,ncol(summaryTable))
    names(summaryTable)[-nc] <- paste(names(summaryTable)[-nc],n,sep="\nN =")
    summaryTable
  }
}

#' Report summary table of means±sd or medians±IQR (Interquartile Range)
#'
#' @description report_medianIQR builds a summary table with the overall mean and sd or
#' medians and quartiles of one or more variables depending on whether the variable is
#' normal or not. A grouping variable can also be specified, in which case means and sd
#' (or medians and quartiles) are also calculated for each group. When means and sd
#' are reported, the p-value of an anova test for testing the significance of the
#' differences between the means is shown. When the reported values are medians and
#' quartiles, the p-value of a Kruskal-Wallis test is displayed to test the null hypothesis
#' that the location parameters of the distribution of the variable are the same in
#' each group.
#' @param data data frame or tibble which contains the data.
#' @param summary_vars Variable or variables whose summary measures (mean±sd or median±IQR)
#' are to be calculated.
#' @param groupVar Grouping variable.
#' @param digits Number of decimal digits for the results.
#' @param pvdigits Number of decimal digits for the p-value of Kruskal test.
#' @param roundFrom If median is greater than this value, median and IQR are rounded to zero
#' decimals even if digits>0.
#' @param probs Pair of quantiles to be computed around the median. Default are q25 and q75.
#' @param alpha Decision on computing mean±sd or median±IQR is made depending on the result of
#' a Shapiro Wilk test of normality. The value of alpha is the significance level for that test.
#' @param na.rm Should NA values be removed (possible values are TRUE or FALSE)
#'
#' @return A table with the overall median and quartiles of the variables and, if a grouping
#' variable is specified, the medians and quartiles by group and the p-value of the kruskal test for
#' comparing location parameters.
#' @importFrom dplyr filter summarize pull bind_rows across enquo
#' @importFrom tidyr pivot_longer
#' @export report_continuous
#'
#' @examples
#' df <- data.frame(x=rnorm(100,10,3),y=rnorm(100,50,8), z=runif(100,20,30),g=sample(c("Yes","No"),100,replace=TRUE))
#' df %>%
#' report_continuous(c(x,y,z))  # Only overall summary of variables x, y and z
#' df %>%
#' report_continuous(c(x,y,z), groupVar=g, digits=1)
report_continuous <- function(data, summary_vars, groupVar=NULL, digits=2, probs=c(0.25,0.75),
                              pvdigits=4, alpha=0.05, na.rm=TRUE) {
  groupVar <- rlang::enquo(groupVar)
  data <- if (rlang::quo_is_null(groupVar)) data %>% mutate(g=1) else
    data %>% mutate(g=!!groupVar)
  # Testing for normality in each group and variable. Variable is declared normal
  # if it is normal in every group
  normalTest <- data %>% select(g,{{summary_vars}}) %>%
    pivot_longer(-g,names_to = "variable",values_to = "value") %>%
    group_by(g,variable) %>%
    nest() %>%
    mutate(
      shapiro_test = map(data, ~ shapiro.test(.x$value)),
      normal = map_dbl(shapiro_test, ~ .x$p.value)>=0.05
    ) %>%
    select(g,variable,normal) %>%
    group_by(variable) %>%
    nest() %>%
    mutate(normal=map_lgl(data,\(x) all(x[["normal"]]))) %>%
    select(variable,normal) %>%
    ungroup()
  normales <- normalTest %>% filter(normal) %>% pull(variable)
  noNormales <- normalTest %>% filter(!normal) %>% pull(variable)
  normalSummary <- NULL
  if (length(normales)>0)
    normalSummary <- data %>%
    report_meanSd(all_of(normales),{{groupVar}},digits=digits,pvdigits=pvdigits,na.rm=na.rm)
  nonNormalSummary <- NULL
  if (length(noNormales)>0)
    nonNormalSummary <- data %>%
    report_medianIQR(all_of(noNormales),{{groupVar}},digits=digits,probs=probs,
                     pvdigits=pvdigits,na.rm=na.rm)
  if(rlang::quo_is_null(enquo(groupVar))){
    if (!is.null(normalSummary)) names(normalSummary)[2]="Overall"
    if (!is.null(nonNormalSummary)) names(nonNormalSummary)[2]="Overall"
  }
  bind_rows(normalSummary,nonNormalSummary)
}


#' Frequency table for binary 0-1 variables
#'
#' @description report_nPctBin01 builds a summary table with the frequency and
#' percentage of times that one or more binary 0-1 variables takes the value 1. If
#' a grouping variable is specified, absolute frequencies and percentages of 1's
#' are computed in each group, and a chi-squared test is performed to test the
#' difference of proportions of 1's between groups. If conditions for the validity
#' of the chi-squared test are not met, a Fisher exact test is performed instead.
#' @param data data frame or tibble which contains the variables.
#' @param summary_vars Binary variables which frequencies and percentages are to
#' be computed
#' @param groupVar Grouping variables.
#' @param digits Number of decimal digits of the result.
#' @param pvdigits Number of decimal digits in the p-value of chi-squared test.
#'
#' @return A table with the overall frequencies and percentages of 1's in each one
#' of the required variables and, if a grouping variable is specified, the frequencies
#' and percentages of 1's in each group as well as a chi-squared test for comparing
#' proportions of 1's between groups.
#' @importFrom dplyr mutate filter group_by count add_row summarize as_label pull full_join select enquo
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom forcats fct_inorder
#' @export report_nPctBin01
#'
#' @examples
#'  df <- data.frame(x=rbinom(90,1,0.3),y=rbinom(90,1,0.8), z=rbinom(90,1,0.5),
#'  g=sample(c("Yes","No"),90,replace=TRUE))
#' df %>%
#' report_nPctBin01(c(x,y,z))  # Overall summary of variables x, y and z
#' df %>%
#' report_nPctBin01(c(x,y,z), groupVar=g, digits=3)
report_nPctBin01 <- function(data, summary_vars, groupVar=NULL, digits=1, pvdigits=4) {
  if(rlang::quo_is_null(enquo(groupVar))){
    data %>%
      pivot_longer({{summary_vars}},names_to="Variable",values_to="value") %>%
      mutate(Variable=fct_inorder(Variable)) %>%
      group_by(Variable) %>%
      summarize(`n (%)`= nPctBin01(value, digits=digits))
  } else{
    data <- data %>% filter(!is.na({{groupVar}}))
    n <- data %>%
      count({{groupVar}}) %>%
      add_row(summarize(.,n=sum(n)),.before=1) %>%
      pull(n)
    data <- data %>%
      select({{groupVar}},{{summary_vars}}) %>%
      pivot_longer(-{{groupVar}},names_to="Variable",values_to="value") %>%
      mutate(Variable=fct_inorder(Variable))
    overall <- data %>%
      group_by(Variable) %>%
      summarize(Overall=nPctBin01(value, digits=digits))
    perGroup <- data %>%
      group_by(Variable,{{groupVar}}) %>%
      summarize(value=nPctBin01(value, digits=digits)) %>%
      pivot_wider(names_from={{groupVar}},values_from=value)
    gvName <- as_label(enquo(groupVar))
    names(perGroup)[-1]=paste(gvName,names(perGroup)[-1],sep="=")
    P <- data %>%
      group_by(Variable) %>%
      summarize(`P-value`=chisq_pval(value,{{groupVar}},pvdigits=pvdigits))
    summaryTable <- suppressMessages(overall %>% full_join(perGroup) %>% full_join(P))
    nc <- c(1,ncol(summaryTable))
    names(summaryTable)[-nc] <- paste(names(summaryTable)[-nc],n,sep="\nN =")
    summaryTable
  }
}

#' Frequency table for categorical variables
#'
#' @description report_nPct builds a summary table with the frequency and
#' percentage of the different values of one or more categorical variables. If
#' a grouping variable is specified, absolute frequencies and percentages o values
#' are computed in each group, and a chi-squared test is performed to test the
#' association with the grouping variable. If conditions for the validity
#' of the chi-squared test are not met, a Fisher exact test is performed instead.
#' @param data data frame or tibble which contains the variables.
#' @param summary_vars Variables whose frequencies and percentages are to
#' be computed
#' @param groupVar Grouping variables.
#' @param digits Number of decimal digits of the result.
#' @param pvdigits Number of decimal digits in the p-value of chi-squared test.
#' @param na.rm Should missing values be included in the table? Percentages are
#' always computed excluding missing values.
#'
#' @return A table with the overall frequencies and percentages of the values of
#' the categorical variable and, if a grouping variable is specified, the frequencies
#' and percentages of the values in each group as well as a chi-squared test for
#' association between the variable and the grouping variable.
#' @importFrom dplyr mutate filter group_by count add_row summarize as_label pull full_join select arrange rename ungroup
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom forcats fct_inorder
#' @importFrom rlang :=
#' @export report_nPct
#'
#' @examples
#' df <- data.frame(x=rbinom(90,3,0.6),y=rbinom(90,4,0.8), z=rbinom(90,5,0.5),
#' g=sample(c("Yes","No"),90,replace=TRUE))
#' df %>%
#' report_nPct(c(x,y,z))  # Overall summary of variables x, y and z
#' df %>%
#' report_nPct(c(x,y,z), groupVar=g, digits=3)
#' df$x[sample(1:90,15)]=NA  # Some missing values are included in x
#' df$y[sample(1:90,8)]=NA  # Some missing values are included in y
#' df %>% report_nPct(c(x,y,z), groupVar=g, digits=3, na.rm=FALSE)

report_nPct <- function(data, summary_vars, groupVar=NULL, digits=1, pvdigits=4, na.rm=TRUE) {
  report_nPct1 <- function(data,Variable,digits=1){
    format <- paste("%d (%.",digits,"f%%)",sep="")
    data %>%
      count({{Variable}}) %>%
      mutate(nPct=ifelse(is.na({{Variable}}),n,
                         sprintf(format,n,100*n/sum(n[!is.na({{Variable}})]))),
             {{Variable}}:=ifelse(is.na({{Variable}}),"(NA) Missing Values",
                                  as.character({{Variable}})))
  }
  message("\nPercentages are computed over the total number of valid values\n")
  if(rlang::quo_is_null(enquo(groupVar))){
    dt1 <- data %>%
      pivot_longer({{summary_vars}},names_to="Variable",values_to="value") %>%
      mutate(Variable=fct_inorder(Variable)) %>%
      filter(if (na.rm) !is.na(value) else TRUE) %>%
      group_by(Variable) %>%
      report_nPct1(value)
    dt1 %>%
      summarize(nPct=paste("Overall N =",as.character(sum(n))), value=" ") %>%
      full_join(dt1) %>%
      select(-n)%>%
      arrange(Variable) %>%
      mutate(Variable=ifelse(value!=" ",paste("|  ",value), as.character(Variable))) %>%
      select(-value)
  } else{
    data <- data %>% filter(!is.na({{groupVar}}))
    n <- data %>%
      count({{groupVar}}) %>%
      add_row(summarize(.,n=sum(n)),.before=1) %>%
      pull(n)
    data <- data %>%
      select({{groupVar}},{{summary_vars}}) %>%
      pivot_longer(-{{groupVar}}, names_to="Variable") %>%
      mutate(Variable=fct_inorder(Variable)) %>%
      filter(if (na.rm) !is.na(value) else TRUE)
    overall <- data %>%
      group_by(Variable) %>%
      report_nPct1(value) %>%
      select(-n) %>%
      rename(Overall=nPct)
    perGroup <-data %>%
      group_by({{groupVar}},Variable) %>%
      report_nPct1(value) %>%
      select(-n) %>%
      pivot_wider(everything(),names_from={{groupVar}},values_from=nPct)
    gvName <- as_label(enquo(groupVar))
    names(perGroup)[-(1:2)]=paste(gvName,names(perGroup)[-(1:2)],sep="=")
    P <- data %>%
      group_by(Variable) %>%
      summarize(`P-value`=chisq_pval(value,{{groupVar}},pvdigits=pvdigits))
    summaryTable <- suppressMessages(
      full_join(overall,perGroup) %>%
        mutate(`P-value`=""))
    summaryTable <- suppressMessages(
      full_join(P,summaryTable) %>%
        replace(is.na(.), "") %>%
        arrange(Variable) %>%
        select(1,3:ncol(.),2) %>%
        mutate(Variable=ifelse(value!="",paste("|  ",value), as.character(Variable))) %>%
        select(-value)
    )
    nc <- c(1,ncol(summaryTable))
    names(summaryTable)[-nc] <- paste(names(summaryTable)[-nc],n,sep="\nN =")
    summaryTable %>% ungroup()
  }
}
#' Elegant format for report tables
#'
#' @param report_data report_data Report table to be formatted.
#'
#' @return Table formatted with flextable.
#' @importFrom dplyr mutate
#' @importFrom flextable flextable set_header_df merge_h merge_v color fontsize width theme_booktabs padding fix_border_issues set_caption autofit
#' @export report_format
#'
#' @examples
#' df <- data.frame(x=rbinom(90,3,0.6),y=rbinom(90,4,0.8), z=rbinom(90,5,0.5),
#' g=sample(c("Yes","No"),90,replace=TRUE))
#' df %>%
#' report_nPct(c(x,y,z)) %>% report_format()  # Overall summary of variables x, y and z
#' df %>%
#' report_nPct(c(x,y,z), groupVar=g, digits=3) %>% report_format(caption="My report table")
#'

report_format <- function (report_data, caption = NULL, fontsize = NULL)
{
  stopifnot(is.data.frame(report_data), ncol(report_data) > 0)
  typology <- function(df) {
    n <- ncol(df)
    col_keys <- names(df)
    grpVar <- strsplit(col_keys[3], "=")[[1]][1]
    if (has_pValue)
      data.frame(col_keys = names(df),
                 what = c(col_keys[1:2], rep(grpVar, n - 3), col_keys[n]),
                 values = c(col_keys[1:2], gsub(paste0(grpVar, "="), "", col_keys[3:(n - 1)]),
                            col_keys[n]), stringsAsFactors = FALSE) else
                              data.frame(col_keys = names(df),
                                         what = c(col_keys[1:2], rep(grpVar, n - 2)),
                                         values = c(col_keys[1:2], gsub(paste0(grpVar, "="), "", col_keys[3:n])),
                                         stringsAsFactors = FALSE)
  }
  nc <- ncol(report_data)
  has_pValue=names(report_data)[nc]=="P-value"
  report_data[[1]] <- gsub("[|]", "   ", report_data[[1]])
  blueLines <- grepl("[(]NA[])]", report_data[[1]])
  valueLines <- grepl("   ", report_data[[1]])
  out <- report_data %>% flextable() %>% {
    if (nc > 3) {
      set_header_df(., mapping = typology(report_data),
                    key = "col_keys") %>% merge_h(part = "header") %>%
        merge_v(j = c(1, 2, nc), part = "header")
    }
    else .
  } %>% theme_booktabs() %>% {
    if (has_pValue)
      color(., color = "red", j = nc, i = ~`P-value` <
              0.05)
    else .
  } %>% {
    if (!is.null(fontsize))
      fontsize(., size = fontsize)
    else .
  } %>% {
    if (nc > 5)
      width(., j = 1:nc, width = c(2, rep(6/(nc - 1), nc -
                                            1)))
    else autofit(.)
  } %>% color(color = "blue", j = 1:max(2, nc - 1), i = blueLines) %>%
    padding(j = 1, i = valueLines, padding.left = 30) %>%
    fix_border_issues() %>% {
      if (!is.null(caption))
        set_caption(., caption = caption)
      else .
    }
  out
}
