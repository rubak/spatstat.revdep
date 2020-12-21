#' Import visual field data from Haag-Streit Eyesuite
#'
#' \code{loadvfEyesuite} imports visual field data from the Eyesuite software by Haag-Streit. These data are converted into a vf object.
#'
#' @param filename the filename of the csv-file
#' @param date_format the order of the dates used in the specific locale
#' @return a vf object

loadvfEyesuite <- function(filename, date_format = "%d.%m.%Y") {
  
  eyesuiteNames <- c(
    "id",
    "lastname",
    "firstname",
    "dateofbirth",
    "gender",
    "ethnicity",
    "V7",
    "V8",
    "V9",
    "V10",
    "apparatus",
    "serial_number",
    "V13",
    "V14",
    "V15",
    "V16",
    "V17",
    "eye",
    "pattern",
    "stimulus_size",
    "stimulus_duration",
    "stiumulus_luminance",
    "strategy",
    "tperimetry",
    "V25",
    "testduration",
    "testdate",
    "test_starting_time",
    "reliability_factor",
    "locnum",
    "questions",
    "repetitions",
    "positive_catch_trials",
    "false_positives",
    "negative_catch_trials",
    "false_negatives",
    "notes",
    "sphere",
    "cylinder",
    "axis",
    "bcva"
  )
  
  # Read the csv-file exported by EyeSuite
  vFieldsRaw <-
    read.csv2(
      filename,
      header = F,
      quote = "",
      stringsAsFactors = F,
      fill = T,
      col.names = paste("V", 1:2000, sep = "")
    )
  
  names(vFieldsRaw)[seq(eyesuiteNames)] <- eyesuiteNames
  
  # use factors to replace values
  vFieldsRaw$eye <- factor(vFieldsRaw$eye,
                             levels = c(0, 1, 3),
                             labels = c("OD", "OS", "binocular"))
  vFieldsRaw$strategy <- factor(
    vFieldsRaw$strategy,
    levels = c(0, 1, 2, 3, 4, 6, 11),
    labels = c(
      "normal",
      "dynamic",
      "2LT/normal",
      "low vision",
      "1LT",
      "TOP",
      "GATE"
    )
  )
  vFieldsRaw$pattern <- factor(vFieldsRaw$pattern,
                                 levels = c("G", "BT"),
                                 labels = c("pG1", "BT"))
  vFieldsRaw$tperimetry <- factor(vFieldsRaw$tperimetry,
                                    levels = c(0, 1),
                                    labels = c("sap", "swap"))
  
  vFieldsRaw <- vFieldsRaw[!is.na(vFieldsRaw$strategy), ]
  vFieldsRaw <- vFieldsRaw[!is.na(vFieldsRaw$pattern), ]
  vFieldsRaw <- vFieldsRaw[!is.na(vFieldsRaw$tperimetry), ]
  
  if (nrow(vFieldsRaw) < 1) stop("There are no valid visual fields in this file.")
  
  # add numbers for each visual field
  vFieldsRaw$i <- 1:nrow(vFieldsRaw)
  
  #extract the header information (id, tperimetry, ...)
  header <- cbind(vFieldsRaw[, seq(eyesuiteNames)], i = vFieldsRaw$i)
  
  # create a data table from the pattern
  patternMatrix <-
    visualFields::saplocmap$pG1[, c("loc", "xod", "yod")]
  
  # exclude binocular visual fields
  if (any(vFieldsRaw$eye == "binocular")) {
    warning(
      "Binocular visual fields are not supported! Binocular visual fields have been removed."
    )
    vFieldsRaw <- vFieldsRaw[vFieldsRaw$eye != "binocular", ]
  }
  
  #function to extract sensitivities for the different loci
  extractLocations <- function(tLine) {
    locMatrix <-
      data.frame(matrix(unlist(tLine[, 44:338]) / 10, 59, 5, byrow = T))
    names(locMatrix) <- c("xod", "yod", "sens1", "sens2", "norm")
    
    if (tLine[1, 18] == "OS")
      locMatrix$xod <- -locMatrix$xod
    
    combinedTable <- merge(locMatrix, patternMatrix)
    returnTable <-
      apply(combinedTable[, c("sens1", "sens2")], 1, mean, na.rm = T)
    names(returnTable) <-
      paste("L", as.character(combinedTable$loc), sep = "")
    rValue <- t(as.matrix(returnTable, 1, 59))
    return(rValue)
    
  }
  
  # apply the extractLocations function on each row
  vFieldsLocs <- data.frame()
  for(j in 1:nrow(vFieldsRaw)) {
    vFieldsLocs <- rbind(vFieldsLocs, extractLocations(vFieldsRaw[j, ]))
  }
  vFieldsRaw <- cbind(vFieldsRaw, vFieldsLocs)
  
  # merge sensitivity data with header
  vFieldsRaw <- merge(header, vFieldsRaw)
  vFieldsRaw$i <- NULL
  
  # convert all text columns into the correct class
  vFieldsRaw$tperimetry <- as.character(vFieldsRaw$tperimetry)
  vFieldsRaw$talgorithm <- as.character(vFieldsRaw$strategy)
  vFieldsRaw$tpattern <- as.character(vFieldsRaw$pattern)
  vFieldsRaw$tdate <- as.Date(vFieldsRaw$testdate, date_format)
  vFieldsRaw$ttime <- as.character(vFieldsRaw$test_starting_time)
  vFieldsRaw$stype <- as.character(vFieldsRaw$notes)
  vFieldsRaw$sage <- agecalc(as.Date(vFieldsRaw$dateofbirth, date_format),
                               as.Date(vFieldsRaw$testdate, date_format))
  vFieldsRaw$seye <- as.character(vFieldsRaw$eye)
  vFieldsRaw$sbsx <- 15
  vFieldsRaw$sbsy <- -3
  vFieldsRaw$sfp <- vFieldsRaw$false_positives / vFieldsRaw$positive_catch_trials
  vFieldsRaw$sfn <- vFieldsRaw$false_negatives / vFieldsRaw$negative_catch_trials
  vFieldsRaw$sfl <- vFieldsRaw$repetitions / vFieldsRaw$questions
  vFieldsRaw$sduration <- as.character(vFieldsRaw$testduration)
  vFieldsRaw$spause <- NA
  
  #exract final table
  finalIndex <-
    c(
      "id",
      "tperimetry",
      "talgorithm",
      "tpattern",
      "tdate",
      "ttime",
      "stype",
      "sage",
      "seye",
      "sbsx",
      "sbsy",
      "sfp",
      "sfn",
      "sfl",
      "sduration",
      "spause"
    )
  finalIndex <- c(finalIndex, paste("L", 1:59, sep = ""))
  
  # return vf-object
  return(vFieldsRaw[, finalIndex])
  
}
