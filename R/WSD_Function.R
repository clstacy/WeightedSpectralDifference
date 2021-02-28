#' Prepare Raw Data for Weighted Spectral Difference
#'
#' This function takes mulitple data.frame objects (or a list of data.frame objects) in wide format and prepares data for calculation of WSD.
#' @param ... input in your raw dataframe(s) either in a single list or as seperate inputs (see example).
#' @param pH a numeric or vector corresponding to pH for each data.frame input into the function
#' @param Temp Default NULL if columns of data represent temperatures. a numeric or vector corresponding to Temp condition for each data.frame input into the function
#' @param ID Identifier string of specimen being sampled, or vector of strings if multiple proteins compared simultaneously. Used for grouping for WSD comparison
#' @return A tidy dataframe with 5 colummns corresponding to: Wavelength, Temperature, MolarElipticity, pH, and ID for each measurement.
#' @export
WSD_dataprep <- function(..., pH, Temp=NULL, ID) {

  if(length(list(...)) == 1L) {
    # check if single input is list (maybe list of dataframes)
    if (is.list(...) == TRUE) {
      x <- list(...)[[1]] #have to list (...) to be able to use it, but already was list, so get 1st list back out with [[1]]
      #purpose of this was to get x <- ... which doesn't work in functions. this was best solution I found.
    } else if (unique(summary(list(...))[,2]) %in% c("data.frame", "tibble", "data.table")) {
      x <- list(...)
    } else {
      stop("data variable input is not of type dataframe or list of dataframes")
    }
  } else {
    x <- list(...)
  }

  if (unique(summary(x)[,2]) %in% c("data.frame", "tibble", "data.table")) {

    df_count <- length(x)
  } else {
    stop("not all unnamed input to function are dataframes of same length")
  }

  ID <- rep(ID, times=df_count, each = max(length(pH), length(Temp)))
  ID <- head(ID, df_count)
  pH <- rep(pH, times=df_count)
  pH <- head(pH, df_count)

  if (df_count == length(pH)) {

    df_tidy <- data.frame() #initialize dataframe to store final output
    for (i in 1:df_count) {
      df <- x[i][[1]]
      colnames(df)[1] <- "WL"
      df_tidy_i <- tidyr::gather(df, Temperature, MolarElipticity, -WL)
      df_tidy_i$pH <- pH[i]
      df_tidy_i$Temperature <- as.numeric(gsub("\\D+", "", df_tidy_i$Temperature))
      df_tidy_i$ID <- ID[i]
      df_tidy <- rbind(df_tidy, df_tidy_i)
    }

    df_tidy_WSD <- df_tidy
    return(df_tidy_WSD)

  } else if (df_count == length(Temp)) {

    df_tidy <- data.frame() #initialize dataframe to store final output

    for (i in 1:df_count) {
      df <- x[i][[1]]
      colnames(df)[1] <- "WL"
      df_tidy_i <- tidyr::gather(df, pH, MolarElipticity, -WL)
      df_tidy_i$Temperature <- Temp[i]
      df_tidy_i$pH <- as.numeric(gsub("\\D+", "", df_tidy_i$pH))
      df_tidy_i$ID <- ID[i]
      df_tidy <- rbind(df_tidy, df_tidy_i)
    }

    df_tidy_WSD <- df_tidy
    return(df_tidy_WSD)

    # } else if (pH == NULL && Temp == NULL) {
    # stop("provide etiher pH or Temperature vector corresponding to each dataframe uploaded")
  } else {
    stop("unknown error, unfamiliar conditions loaded in")
  }


}







#' Weighted Spectral Difference
#'
#' Input data must be in tidy format.
#' This function calculates Weighted Spectral Difference for preprocessed data
#' that has been made with WSD_dataprep or otherwise tidy preprocessed data.
#' Output is a tidy data.frame with all WSD (and SD) values for input raw data.
#' @param data A tidy 5 column dataframe containing wavelength, temperature, Molar Elipticity or Radians Measures from Circular Dichroism Experiment, pH, and protein ID columns
#' @param units (unsupported) The units of input data. Defaults to Molar Elipticity. Options: "MolElip", "Radians", etc
#' @param control_temp Temperature of control condition from which to calculate weighted spectral difference
#' @param control_pH pH of control condition from which to calculate weighted spectral difference
#' @param min_wavelength lowest wavelength (from first column of dataframe) to include in weighted spectral difference calculation
#' @param max_wavelength highest wavelength (from first column of dataframe) to include in weighted spectral difference calculation
#' @return A tidy data.frame that contains WSD estimates for each condition tested
#' @export
WSD <- function(data, control_temp, control_pH, units=NULL, min_wavelength, max_wavelength) {

  ID <- data$ID
  data <- data
  for (protein in unique(ID)) {
    #select control data
    control <- dplyr::filter(data, ID == protein)
    control <- dplyr::filter(control, Temperature == control_temp & pH == control_pH)
    control <- dplyr::filter(control, WL >= min_wavelength, WL <= max_wavelength)
    control <- dplyr::arrange(control, WL)
    control <- dplyr::pull(control, MolarElipticity)

    # select sample data (including control)
    samples <- dplyr::filter(data, ID == protein)
    samples <- dplyr::filter(samples, WL >= min_wavelength, WL <= max_wavelength)
    samples <- tidyr::spread(samples, Temperature, MolarElipticity)
    samples <- dplyr::arrange(samples, pH, WL)
    samples <- dplyr::group_by(samples, pH)
    samples <- dplyr::group_split(samples)

    #find weighted squared difference between y_Ai and y_Bi
    if (!exists("WSD_output")) {
      WSD_output <- data.frame()
    }

    for (i in 1:length(samples)) {
      num_conditions <- ncol(samples[[i]])
      for (j in 4:num_conditions) { #note we start at 3 because col 1 is wavelength and col2 is pH, col 3 is protein
        SD <- sqrt(
          sum(
            (control - samples[[i]][[j]])^2  / length(control)
          )
        )
        WSD <- sqrt(
          sum(
            ((abs(control)) / sum(abs(control))) * (control - samples[[i]][[j]])^2
          )
        )
        if (!exists("WSD_output")) {
          WSD_output <- data.frame(colnames(samples[[i]][j]), samples[[i]][[2]][1], WSD, SD, protein)
        } else {
          WSD_output_new <- data.frame(colnames(samples[[i]][j]), samples[[i]][[2]][1], WSD, SD, protein)
          WSD_output <- rbind(WSD_output, WSD_output_new)
        }
      }

    }

  }
  colnames(WSD_output) <- c("Temperature", "pH", "WSD", "SD", "ID")
  return(WSD_output) # return output for samples from all different proteins
}











# Here is how I imagine workflow to go:
# save raw data as .csv file(s)
# upload files into R env as dataframes
# create list of dataframes to run at once?
#>>better idea: have separate command in package to prep and gather/tidy data for WSD() call
# call WSD()
#returns output dataframe in tidy format with all WSD values for each 'group'
# save WSD() object into env [user command]
# plot WSD() output with ggplot [ggplot2 function in package]


#
#
#
#
# WSD <- function(control, samples, ControlTemp, ControlpH, units, ID, MinWavelength, MaxWavelength) {
#   #find squared difference between y_Ai and y_Bi
#   WSD_output <- data.frame()
#   for (i in 1:length(samples)) {
#     for (j in 3:ncol(samples[[i]])) { #note we start at 3 because col 1 is wavelength and col2 is pH
#       SD <- sqrt(
#         sum(
#           (control - samples[[i]][[j]])^2  / length(control)
#         )
#       )
#       WSD <- sqrt(
#         sum(
#           ((abs(control)) / sum(abs(control))) * (control - samples[[i]][[j]])^2
#         )
#       )
#       if (!exists("WSD_output")) {
#         WSD_output <- data.frame(colnames(samples[[i]][j]), samples[[i]][[2]][1], WSD, SD)
#       } else {
#         WSD_output_new <- data.frame(colnames(samples[[i]][j]), samples[[i]][[2]][1], WSD, SD)
#         WSD_output <- rbind(WSD_output, WSD_output_new)
#       }
#     }
#   }
#   colnames(WSD_output) <- c("Temperature", "pH", "WSD", "SD")
#   return(WSD_output)
#
#
#
#
# }
#
#
#
# WSD_manual <- function(control, samples, ControlTemp, ControlpH, units, ID, MinWavelength, MaxWavelength) {
#   #find squared difference between y_Ai and y_Bi
#   WSD_output <- data.frame()
#   for (i in 1:length(samples)) {
#     for (j in 3:ncol(samples[[i]])) { #note we start at 3 because col 1 is wavelength and col2 is pH
#       SD <- sqrt(
#         sum(
#           (control - samples[[i]][[j]])^2  / length(control)
#         )
#       )
#       WSD <- sqrt(
#         sum(
#           ((abs(control)) / sum(abs(control))) * (control - samples[[i]][[j]])^2
#         )
#       )
#       if (!exists("WSD_output")) {
#         WSD_output <- data.frame(colnames(samples[[i]][j]), samples[[i]][[2]][1], WSD, SD)
#       } else {
#         WSD_output_new <- data.frame(colnames(samples[[i]][j]), samples[[i]][[2]][1], WSD, SD)
#         WSD_output <- rbind(WSD_output, WSD_output_new)
#       }
#     }
#   }
#   colnames(WSD_output) <- c("Temperature", "pH", "WSD", "SD")
#   return(WSD_output)




# }


#TODO

# check if control is either a vector or a dataframe/datatable/matrix with only 1 dimension
# attempt to convert control input to a vector
# error if ctl cannot be made into a vector
# determine type of samples input
# provide error if it is not either a dataframe/datatable/tibble/matrix variable with at least 2 columns
# attempt to convert samples to dataframe (if not already)
# NOTE: samples must have first column be wavelength values corresponding to elipticity measurements.
# give error if samples cannot be converted
# give error if length ctl != nrow(samples)

# provide possible pre-built options for units of input data
# default the unit to be MolarElipticity?
# give error if units is not recognized type

# Allow optional ID column for labelling all results from output. for example, name of protein
#allows for combination of results from multiple

# allow ControlTemp and ControlpH be either numeric values or vectors same length as




# new function to convert from raw data to Molar Elipticity units?
#

#
# WSD_preprocessed <- function(control, samples, ControlTemp, ControlpH, units, ID, MinWavelength, MaxWavelength) {
#   #this function is appropriate to be used on data that is already tidy formatted with col 1 re
#   ####NEED TO FINISH THIS DESCRIPTION######
#
#   #find squared difference between y_Ai and y_Bi
#   WSD_output <- data.frame()
#   for (i in 1:length(samples)) {
#     for (j in 3:ncol(samples[[i]])) { #note we start at 3 because col 1 is wavelength and col2 is pH
#       SD <- sqrt(
#         sum(
#           (control - samples[[i]][[j]])^2  / length(control)
#         )
#       )
#       WSD <- sqrt(
#         sum(
#           ((abs(control)) / sum(abs(control))) * (control - samples[[i]][[j]])^2
#         )
#       )
#       if (!exists("WSD_output")) {
#         WSD_output <- data.frame(colnames(samples[[i]][j]), samples[[i]][[2]][1], WSD, SD)
#       } else {
#         WSD_output_new <- data.frame(colnames(samples[[i]][j]), samples[[i]][[2]][1], WSD, SD)
#         WSD_output <- rbind(WSD_output, WSD_output_new)
#       }
#     }
#   }
#   colnames(WSD_output) <- c("Temperature", "pH", "WSD", "SD")
#   return(WSD_output)
# }




# WSD_Combine <- function() {
#   #used to combine multiple WSD calculations from different trials (useful for comparisons, statistical significance, plotting combined, etc)
#   }




# #' Plot WSD values according to Temperature
# #'
# #' This function plots Weighted Spectral Difference between each condition compared to control and plots as a function of temperature
# #' It plots where color is pH, shape is protein ID.
# #' @param data The saved data.frame output of WSD() which includes weighted spectral difference values for curves in data set
# #' @return A ggplot object that shows WSD values as a function of temperature for each pH and substrate tested.
# #' @export
# WSD_plot <- function(data){
  #able to plot WSD outputs using ggplot2
#  ggplot2::ggplot(data = data, aes(x=Temperature, y=WSD, group= interaction(pH, ID), color=as.factor(pH), shape=ID)) + #convert temp to factor so not continuous
#    geom_abline(intercept = 0, size = 1) + #horizontal line
#    geom_point(size = 2.5) +
#    geom_line(size = 1) +
#    theme_bw() +
#    labs(color="pH", title = "Spectral Difference values of Proteins") +
#    xlab("Temperature (°C)")
#}
