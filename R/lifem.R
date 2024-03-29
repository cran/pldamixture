#' LIFE-M Data
#'
#' The `lifem` data set contains a subset of data from the Life-M project (\url{https://life-m.org/}) on 3,238 individuals born between 1883
#' to 1906. These records were obtained from linking birth certificates and death certificates either of two
#' ways. A fraction of the records (2,159 records) were randomly sampled to be “hand-linked at some level” (HL).
#' These records are high quality and were manually linked at some point by trained research assistants.
#' The remaining records were “purely machine-linked” (ML) based on probabilistic record linkage without clerical
#' review. The Life-M team expects the mismatch rate among these records to be around 5% (Bailey et al.
#' 2022). Of interest is the relationship between age at death and year of birth. The `lifem` demo data set consists of 2,159 hand-linked records
#' and 1,079 records that were randomly sampled from the purely machine-linked records (~2:1 HL-ML ratio).
#'
#' \itemize{
#'   \item yob: year of birth (value from 1883 and 1906)
#'   \item unit_yob: yob re-scaled to the unit interval for analysis (between 0 and 1). If X is the yob, we use the following:
#'   (X – min(X)) / (max(X) – min(X)) = a * X + b, a = 1/(max(X) – min(X)), b = -min(X)*a
#'   \item age_at_death: age at death (in years)
#'   \item hndlnk: whether record was purely machine-linked or hand-linked at some level.
#'   \item commf: commonness score of first name (between 0 and 1). It is based on the 1940 census.
#'   It is a ratio of the log count of the individual’s first name over the log count of the most
#'   commonly occurring first name in the census.
#'   \item comml: commonness score of last name (between 0 and 1). It is based on the 1940 census.
#'   It is a ratio of the log count of the individual’s last name over the log count of the most
#'   commonly occurring last name in the census.
#' }
#'
#' @docType data
#' @name lifem
#' @usage data(lifem)
#' @format a data frame with 3,238 rows and 6 variables
#'
#' @references Bailey, Martha J., Lin, Peter Z., Mohammed, A.R. Shariq, Mohnen, Paul, Murray, Jared,
#' Zhang, Mengying, and Prettyman, Alexa. LIFE-M: The Longitudinal, Intergenerational
#' Family Electronic Micro-Database. Ann Arbor, MI: Inter-university Consortium for
#' Political and Social Research (distributor), 2022-12-21. < \doi{10.3886/E155186V5} >
"lifem"
