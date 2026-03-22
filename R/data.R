#' GAD-7 Anxiety Dataset
#' @description Binary GAD-7 responses from a national survey of Chinese
#'   adults (N = 2,404). Original 4-point responses dichotomised:
#'   "never" = 0; all other responses = 1. Exemplar dataset for the
#'   NIRApost tutorial.
#' @format A \code{data.frame} with 2,404 rows and 7 binary (0/1) columns.
#' @source National survey on anxiety among Chinese adults.
#' @docType data
#' @name single_gds
#' @examples
#' data("single_gds")
#' dim(single_gds)
#' checkMissing(single_gds)
"single_gds"

#' PHQ-9 Depression Dataset
#' @description Binary PHQ-9 responses from the 2023 PBICR database
#'   (N = 45,829 mainland Chinese participants). "Not at all" = 0;
#'   all other responses = 1.
#' @format A \code{data.frame} with 45,829 rows and 9 binary (0/1) columns.
#' @source 2023 PBICR database.
#' @docType data
#' @name PHQ9
#' @examples
#' data("PHQ9"); dim(PHQ9)
"PHQ9"

#' GDS-9 Geriatric Depression Dataset
#' @description GDS-9 (Geriatric Depression Scale, 9-item version) data
#'   from 3,097 Chinese elderly participants. Already binary; no
#'   transformation required before NIRA.
#' @format A \code{data.frame} with 3,097 rows and 9 binary (0/1) columns.
#' @source Clinical survey of geriatric depression in Chinese elderly.
#' @docType data
#' @name GDS9
#' @examples
#' data("GDS9"); dim(GDS9)
"GDS9"

#' ACE (Adverse Childhood Experiences) Dataset
#' @description ACE data from 525 Chinese university students.
#'   Ten binary items; directly compatible with NIRA.
#' @format A \code{data.frame} with 525 rows and 10 binary (0/1) columns.
#' @source Survey of adverse childhood experiences among Chinese students.
#' @docType data
#' @name ACE
#' @examples
#' data("ACE"); dim(ACE)
"ACE"

#' Chronic Disease Dataset
#' @description Disease presence (1) / absence (0) for 8 conditions
#'   across 27,353 Chinese adults: hypertension, diabetes,
#'   hyperlipidemia, stroke, coronary heart disease, digestive disorders,
#'   respiratory diseases, and urinary system diseases.
#' @format A \code{data.frame} with 27,353 rows and 8 binary (0/1) columns.
#' @source Epidemiological survey of chronic diseases among Chinese adults.
#' @docType data
#' @name disease
#' @examples
#' data("disease"); dim(disease)
"disease"

