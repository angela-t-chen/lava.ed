#' Apply rule-based algorithms
#'
#' This function applies original ICD code rule-based algorithms to your data. Seven 
#' algorithms available in ICD-10 native coding are provided for use. Diagnosis codes 
#' are coded using ICD-10 and included at the “category” level, corresponding to 
#' the first three characters of the codes. When indicated, this function follows 
#' the conditional logic of algorithms provided. 
#' 
#' @param dat A dataframe with observations to be classified.
#' @param dx1 The primary diagnosis, coded in ICD-10. Include only the first three 
#'   characters of the ICD-10 code.
#' @param algorithm Selection of rule-based algorithms to apply to the data. Choices
#'   include "billings", "ahrq", "weinick", "aihw", "sundmacher", "nicholl", and 
#'   "scheuttig".
#' @param dx2 A secondary diagnosis, if available, also coded in ICD-10. Include
#'   only the first three characters of the code.
#' @param age Patient age.
#' 
#' @return The dataframe `dat` with additional columns appended for each selected 
#'   algorithm. Appended columns consist of TRUE/FALSE indicators, with TRUE indicating 
#'   that the observation is classified as low acuity by the selected algorithm.
#'
#' @examples 
#' nhamcs <- rule_based_algos(nhamcs, "icd10", dx1 = "diag1", 
#'      algorithm = c("billings", "ahrq", "weinick", "aihw", "sundmacher", 
#'      "nicholl", "scheuttig"), dx2 ="diag2", age = "age")
#' 
#' @export

rule_based_algos <- function(dat, dx1, algorithm, ...) {
     # List of algorithms (i.e., different papers)
     algo <- list(
          Billings = function(dat, dx1, dx2, age, ...) {
               icd10 <- c("A50", "A33", "A34", "A35", "A37", "A80", "G00", 
                          "I01", "G40", "R56", "H66", "J02", "J03", "J06", 
                          "J31", "A15", "A17", "A18", "A19", "J20", "J40", 
                          "J41", "J42", "J43", "J44", "J47", "J13", "J14", 
                          "J15", "J16", "J18", "J45", "I50", "I11", "J81", 
                          "I10", "I11", "I16", "I20", "I24", "L03", "L04", 
                          "L08", "L88", "L98", "E10", "E13", "E11", "E16", 
                          "K52", "N10", "N11", "N12", "E86", "D50", "E40", 
                          "E41", "E43", "E55", "E64", "R62", "N70", "N73", 
                          "K02", "K03", "K04", "K05", "K06", "K08", "K12", 
                          "K13", "M27", "A69", "K09")
               codes <- dat[[dx1]]
               logivec <- codes %in% icd10
               # Modifiers - constraints (changing to FALSE)
               tempvec <- which(grepl("^A50", dat[[dx1]]) == TRUE & dat[[age]] < 1 |
                                     grepl("^G00", dat[[dx1]]) == TRUE & dat[[age]] > 5 |
                                     grepl("^H66|^H67", dat[[dx1]]) == TRUE & grepl("^Z96", dat[[dx2]])|
                                     grepl("^J13|^J14|^J15|^J16|^J18", dat[[dx1]]) == TRUE & grepl("^D57", dat[[dx2]]) |
                                     grepl("^J13|^J14|^J15|^J16|^J18", dat[[dx1]]) == TRUE & dat[[age]] < 1 |
                                     grepl("^D50", dat[[dx1]]) == TRUE & dat[[age]] > 5 |
                                     grepl("^R62", dat[[dx1]]) == TRUE & dat[[age]] < 1)
               logivec[tempvec] <- FALSE
               # Modifiers - checking secondary diagnoses (changing to TRUE)
               tempvec <- which(grepl("^E86", dat[[dx2]]) == TRUE|
                                     grepl("^D50", dat[[dx2]]) == TRUE & dat[[age]] <= 5 |
                                     grepl("^E40|^E41|^E43|^E55|^E64", dat[[dx2]]) == TRUE)
               logivec[tempvec] <- TRUE
               return(logivec)
          },
          Ahrq = function(dat, dx1, ...){
               icd10 <- c("E10", "E11", "J41", "J42", "J43", "J44", "J47", "J45", 
                          "E84", "J84", "P27", "Q25", "Q31", "Q32", "Q33", "Q34", 
                          "Q34", "Q39", "Q31", "Q89", "I10", "I11", "I12", "I13", 
                          "I16", "031", "041", "I09", "I11", "I13", "I50", "J13", 
                          "J14", "J15", "J16", "J18", "D57", "N10", "N12", "N15", 
                          "N16", "N28", "N30", "N39", "N11", "N13", "Q60", "Q61", 
                          "Q62", "Q63", "Q64", "P27")
               codes <- dat[[dx1]]
               logivec <- codes %in% icd10
          },
          Weinick = function(dat, dx1, ...){
               icd10 <- c("J30", "J20", "J21", "J40", "A74", "B30", "B58", "H16", 
                          "H10", "H11", "H01", "H02", "J11", "J12", "J10", "J09", 
                          "H60", "H61", "H62", "H65", "H68", "H69", "H66", "H67", 
                          "J02", "J03", "A38", "R07", "J00", "J01", "J04", "J05", 
                          "J06", "J32", "J34", "N30", "N34", "N39", "R30", "B33", 
                          "B97", "B34", "L02", "L03", "K12", "L04", "L01", "L05", 
                          "L08", "L98", "E83", "S00", "S10", "S05", "S20", "S30", 
                          "T14", "S40", "S50", "S60", "S70", "S80", "S90", "K00", 
                          "K01", "K02", "K03", "K04", "K05", "K06", "S02", "M20", 
                          "S43", "S41", "S53", "S51", "S63", "S83", "S46", "S56", 
                          "S66", "S73", "S76", "S86", "S93", "S96", "S33", "S13", 
                          "S23", "S03", "S29", "S39", "S22", "S32", "S42", "S52", 
                          "S62", "S82", "S92", "S99", "M19", "M18", "M16", "M17", 
                          "M25", "M75", "M77", "M70", "M76", "M65", "D48", "M21", 
                          "M71", "M67", "M66", "M60", "M62", "M79", "M54", "S01", 
                          "S08", "S09", "S21", "S31", "S71", "S81", "S91", "N36", 
                          "N13", "R31", "M15", "S61", "H00")
               codes <- dat[[dx1]]
               logivec <- codes %in% icd10 
          },
          Aihw = function(dat, dx1, dx2, age, ...){
               icd10 <- c("J10", "J11", "J13", "J14", "J15", "J16", "J18", "A35", 
                          "A36", "A37", "A80", "B05", "B06", "B16", "B18", "B26", 
                          "G00", "M01", "E10", "E11", "E13", "E14", "E40", "E41", 
                          "E42", "E43", "E55", "E64", "D50", "I10", "I11", "I50", 
                          "J81", "I20", "I24", "J41", "J42", "J43", "J44", "J47", 
                          "J20", "J45", "J46", "E86", "K52", "G40", "G41", "O15", 
                          "R56", "H66", "H67", "J02", "J03", "J06", "J31", "A69", 
                          "K02", "K03", "K04", "K05", "K06", "K08", "K09", "K00", 
                          "K12", "K13", "K25", "K26", "K27", "K28", "K35", "N10", 
                          "N11", "N12", "N13", "N70", "N73", "N74", "L03", "L04", 
                          "L08", "L88", "L98", "R02")
               codes <- dat[[dx1]]
               logivec <- codes %in% icd10
               # Modifiers - constraints (changing to FALSE)
               tempvec <- which(grepl("^J10|^J11|^J13|^J14|^J15|^J16|^J18", dat[[dx1]]) == TRUE & dat[[age]] < 1 |
                                     grepl("^J10|^J11|^J13|^J14|^J15|^J16|^J18", dat[[dx1]]) == TRUE & grepl("^D57", dat[[dx2]]))
               logivec[tempvec] <- FALSE
               return(logivec)
          },
          Sundmacher = function(dat, dx1, ...){
               icd10 <- c("I20", "I25", "I50", "I05", "I06", "I08", "I49", "I67",
                          "I70", "I73", "I78", "I80", "I83", "I86", "I87", "I95", 
                          "R00", "R47", "J20", "J21", "J40", "J41", "J42", "J43", 
                          "J44", "J47", "F10", "F11", "M42", "M47", "M53", "M54", 
                          "I10", "I11", "I12", "I13", "I14", "I15", "K52", "K57", 
                          "K58", "K59", "A01", "A02", "A04", "A05", "A07", "A08", 
                          "A09", "J10", "J11", "J13", "J14", "J15", "J16", "J18", 
                          "H66", "J01", "J02", "J03", "J06", "J31", "J32", "J35", 
                          "F32", "F33", "E10", "E11", "E13", "E14", "M17", "G56", 
                          "M67", "M71", "M75", "M76", "M77", "M79", "F40", "F41", 
                          "F43", "F45", "F50", "F60", "H25", "H40", "N30", "N34", 
                          "N39", "G47", "A46", "L01", "L02", "L04", "L08", "L60", 
                          "L72", "L98", "E03", "E04", "E05", "E86", "E87", "E89", 
                          "C43", "C44", "K21", "K29", "K30", "K31", "G43", "G44", 
                          "R51", "D50", "D51", "D52", "D53", "D56", "E40", "E41", 
                          "E42", "E43", "E44", "E45", "E46", "E47", "E48", "E49", 
                          "E50", "E51", "E52", "E53", "E54", "E55", "E56", "E57", 
                          "E58", "E59", "E60", "E61", "E62", "E63", "E64", "R63", 
                          "K70", "K02", "K04", "K05", "K06", "K08", "K12", "K13", 
                          "N70", "N71", "N72", "N75", "N76", "N84", "N86", "N87", 
                          "F01", "F03", "O23", "O24", "N41", "N45", "N48", "J45", 
                          "G62", "A15", "A16", "A34", "A35", "A36", "A37", "A50", 
                          "A51", "A52", "A53", "A54", "A55", "A56", "A57", "A58", 
                          "A63", "A64", "A80", "B05", "B06", "B07", "B15", "B16", 
                          "B17", "B18", "B20", "B21", "B22", "B23", "B24", "B26", 
                          "B34", "B51", "B52", "B53", "B54", "B77", "B86", "R56", 
                          "L89", "E66", "K25", "K27", "L97", "F80", "Z73")
               codes <- dat[[dx1]]
               logivec <- codes %in% icd10
          },
          Nicholl = function(dat, dx1, ...){
               icd10 <- c("S00", "R07", "I80", "I81", "I82", "I20", "J40", "J41", 
                          "J42", "J43", "J44", "R10", "L03", "R50", "T83", "E10", 
                          "E11", "E12", "E13", "E14", "E15", "E16", "N39", "G40", 
                          "G41")
               codes <- dat[[dx1]]
               logivec <- codes %in% icd10
          },
          Scheuttig = function(dat, dx1, ...){
               icd10 <- c("Z23", "Z24", "Z25", "Z27", "Z29", "L20", "L22", "L23", 
                          "R21", "M99", "J40", "Z72", "Z74", "Z75", "Z76", "Z01", 
                          "Z09", "M62", "M79", "I10", "M10", "M25", "H66", "J00", 
                          "J02", "J03", "J04", "J06", "J20", "H10", "B34", "B85", 
                          "B86", "B99", "M54", "K52", "K59", "E11", "N30", "N39", 
                          "R30", "R32", "Z41", "Z43", "Z47", "Z48", "Z59", "Z63", 
                          "R50", "R51", "R52", "R53", "K11", "K12", "K14", "H91", 
                          "H93", "T67", "T75", "T78", "T83", "A09", "I80", "I83", 
                          "S00", "S30", "S40", "S50", "S60", "S70", "S80", "S90", 
                          "T00", "T14", "F32", "F33", "F34", "F38", "F39", "R04", 
                          "R05", "R06", "R10", "R11", "N41", "N45", "N48", "N12", 
                          "N28", "G35", "S20", "S63", "S83", "S93", "S92", "A46", 
                          "L03")
               codes <- dat[[dx1]]
               logivec <- codes %in% icd10
          }
     )
     algo_used <- algo[stringr::str_to_title(algorithm)]
     logical_list <- lapply(algo_used, function(f) f(dat, dx1, dx2, age, ...))
     logical_df <- as.data.frame(logical_list)
     dat2 <- cbind(dat, logical_df)
     return(dat2)
}
