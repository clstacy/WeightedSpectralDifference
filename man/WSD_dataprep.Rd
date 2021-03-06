% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WSD_Function.R
\name{WSD_dataprep}
\alias{WSD_dataprep}
\title{Prepare Raw Data for Weighted Spectral Difference}
\usage{
WSD_dataprep(..., pH, Temp = NULL, ID)
}
\arguments{
\item{...}{input in your raw dataframe(s) either in a single list or as seperate inputs (see example).}

\item{pH}{a numeric or vector corresponding to pH for each data.frame input into the function}

\item{Temp}{Default NULL if columns of data represent temperatures. a numeric or vector corresponding to Temp condition for each data.frame input into the function}

\item{ID}{Identifier string of specimen being sampled, or vector of strings if multiple proteins compared simultaneously. Used for grouping for WSD comparison}
}
\value{
A tidy dataframe with 5 colummns corresponding to: Wavelength, Temperature, MolarElipticity, pH, and ID for each measurement.
}
\description{
This function takes mulitple data.frame objects (or a list of data.frame objects) in wide format and prepares data for calculation of WSD.
}
