#' Case-control study on oesophageal cancer in Chinese Singapore men.
#'
#' This dataset is borrowed from "Aetiological factors in oesophageal cancer in Singapore Chinese" by De Jong UW, Breslow N, Hong JG, Sridharan M, Shanmugaratnam K (1974).
#'
#' @docType data
#' @name singapore
#' @usage data(singapore)
#' @format The dataset contains the following variables:
#' \describe{
#' \item{Age}{age of the patient.}
#' \item{Dial}{dialect group where 1 represent "\code{Hokhien/Teochew}" and 0 represent "\code{Cantonese/Other}".}
#' \item{Samsu}{a numeric indicator of whether the patient consumes Samsu wine or not.}
#' \item{Cigs}{number of cigarettes smoked per day.}
#' \item{Bev}{number of beverage at "burning hot" temperatures ranging between 0 to 3 different drinks per day.}
#' \item{Everhotbev}{a numeric indicator of whether the patients ever drinks "burning hot beverage" or not. Recoded from the variable "\code{Bev}".}
#' \item{Set}{matched set identification number.}
#' \item{CC}{a numeric variable where 1 represent if the patient is a case, 2 represent if the patient is a control from the same ward as the case and 3 represent if the patient is control from orthopedic hospital.}
#' \item{Oesophagealcancer}{a numeric indicator variable of whether the patient is a case of oesophageal cancer or not.}
#' }
#'
#' \strong{The following changes have been made to the data from the original data in De Jong UW (1974):}
#'
#' - The variable "\code{Bev}" is recoded into the numeric indicator variable "\code{Everhotbev}":
#'
#' \code{singapore$Everhotbev <- ifelse(singapore$Bev >= 1, 1, 0)}
#'
#' @references 	De Jong UW, Breslow N, Hong JG, Sridharan M, Shanmugaratnam K. (1974). Aetiological factors in oesophageal cancer in Singapore Chinese. \emph{Int J Cancer} Mar 15;13(3), 291-303.
#' @references \url{http://faculty.washington.edu/heagerty/Courses/b513/WEB2002/datasets.html}

NULL


#' Birthweight data clustered on the mother.
#'
#' This dataset is borrowed from "An introduction to Stata for health reserachers" (Juul and Frydenberg, 2010).
#' The dataset contains data on 189 mothers who have given birth to one or several children. In total, the dataset contains data on 487 births.
#'
#' @docType data
#' @name clslowbwt
#' @usage data(clslowbwt)
#' @format The dataset is structured so that each row corresponds to one birth/child. It contains the following variables:
#' \describe{
#'   \item{id}{the identification number of the mother.}
#'   \item{birth}{the number of the birth, i.e. "1" for the mother's first birth, "2" for the mother's second birth etc.}
#'   \item{smoke}{a categorical variable indicating if the mother is a smoker or not with levels "\code{0. No}" and "\code{1. Yes}".}
#'   \item{race}{the race of the mother with levels "\code{1. White}", "\code{2. Black}" or "\code{3. Other}".}
#'   \item{age}{the age of the mother at childbirth.}
#'   \item{lwt}{weight of the mother at last menstruational period (in pounds).}
#'   \item{bwt}{birthweight of the newborn.}
#'   \item{low}{a categorical variable indicating if the newborn is categorized as a low birthweight baby (<2500 grams) or not with levels "\code{0. No}" and "\code{1. Yes}".}
#'   \item{smoker}{a numeric indicator if the mother is a smoker or not. Recoded version of the variable "\code{smoke}" where "\code{0.No}" is recoded as "0" and "\code{1.Yes}" is recoded as "1".}
#'   \item{lbw}{a numeric indicator of whether the newborn is categorized as a low birthweight baby (<2500 grams) or not. Recoded version of the variable "\code{low}" where "\code{0.No}" is recoded as "0" and "\code{1.Yes}" is recoded as "1".}
#' }
#'
#'
#'\strong{The following changes have been made to the original data in Juul & Frydenberg (2010):}
#'
#' - The variable "\code{low}" is recoded into the numeric indicator variable "\code{lbw}":
#'
#' \code{clslowbwt$lbw <- as.numeric(clslowbwt$low == "1. Yes")}
#'
#' - The variable "\code{smoke}" is recoded into the numeric indicator variable "\code{smoker}":
#'
#' \code{clslowbwt$smoker <- as.numeric(clslowbwt$smoke == "1. Yes")}
#'
#' @references Juul, Svend & Frydenberg, Morten (2010). \emph{An introduction to Stata for health researchers}, Texas, Stata press, 2010 (Third edition).
#' @references \url{http://www.stata-press.com/data/ishr3.html}
#'

NULL
