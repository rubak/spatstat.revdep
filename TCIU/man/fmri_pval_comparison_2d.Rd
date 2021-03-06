% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmri_pval_comparison_2d.R
\name{fmri_pval_comparison_2d}
\alias{fmri_pval_comparison_2d}
\title{2D comparison visualization between the p-values}
\usage{
fmri_pval_comparison_2d(
  pval_ls,
  pval_name_ls,
  axis_i_lses,
  hemody_data = NULL,
  mask,
  p_threshold = 0.05,
  legend_show = TRUE,
  method = "scale_p",
  color_pal = "YlOrRd",
  multi_pranges = TRUE,
  mask_width = 1.5
)
}
\arguments{
\item{pval_ls}{a list. Each element is a 3D array of p-values data}

\item{pval_name_ls}{a list with the element as name for the p-values data provided in 'pval_ls'}

\item{axis_i_lses}{a list with 3 numeric elements or a list of lists. If the elements are numeric, 
they would specify indices of slice for the three direction. 
If any direction of the slice need not to be shown, make it as NULL for that element. 
If elements are lists, each list provides specified cuts for corresponding 3D p-values data.}

\item{hemody_data}{a parameter to have the plot with/without hemodynamic contour. The default is NULL to make the plot 
without hemodynamic contour, otherwise assign a 3D array of the hemodynamic data.}

\item{mask}{a 3D nifti or 3D array of data to show the shell of the brain}

\item{p_threshold}{NULL or a numeric value that can be selected randomly below 0.05 to 
drop all p-values above the threshold. If 'low5_percent' method is used, 
make 'p_threshold' as NULL. The default is 0.05.}

\item{legend_show}{a logical parameter to specify whether the final plot has legend
for all the subplots or the shared legend for all the subplots. The default is TRUE.}

\item{method}{a string that represents method for the plot. 
There are 3 options: 'min_max', 'scale_p' and 'low5_percent'. The default is 'scale_p'. 
'min_max' is to draw plot based on the color scale of the minimum and maximum of the p-values; 
'scale_p' is to draw the plot with fixed color scale for fixed range of p-values; 
'low5_percent' is to draw the plot for the smallest 5 percent of p-values
when all the p-values are not significant}

\item{color_pal}{the name of the color palettes provided by \code{RColorBrewer} The default is "YlOrRd".}

\item{multi_pranges}{an option under 'scale_p' method to decide whether there are at most 9 colors 
in the legend for the ranges of p-values, or at most 4 colors. 
The default is TRUE, choosing the larger number of colors for the plot.}

\item{mask_width}{a numeric value to specify the width of mask contour. The default is 1.5.}
}
\value{
a combination plot arranged by \code{gridExtra}
}
\description{
a plot arrangement method, which uses \code{gridExtra} to combine multiple 
2D plots of the fMRI data together. This can bring convenience for users to compare 
the result of different statistical tests based on the p values they provide
}
\details{
The function \code{fmri_pval_comparison_2d} is used to combine and compare the 2D plots for different 3D arrays of p-values.
The plots in each row are generated by one specific 3D p value data. 
The first column of the integrated plot specifies the name of the 3D p value data (for generation of 
the plots in that row). The rest of the three columns are the plots from sagittal, coronal and 
axial view for each 3D p value data.
}
\examples{
# sample 3D data of mask provided by the package
dim(mask)
# sample 3D p value provided by the package
dim(phase2_pval)
dim(phase3_pval)

\donttest{
fmri_pval_comparison_2d(list(phase2_pval, phase3_pval), 
                        list('phase2_pval', 'phase3_pval'),
                        list(list(40, 26, 33), list(40, 26, 33)), 
                        hemody_data = NULL, 
                        mask = mask, p_threshold = 0.05, 
                        legend_show = FALSE, method = 'scale_p',
                        color_pal = "YlOrRd", multi_pranges=TRUE)
}                  
}
\author{
SOCR team <\url{http://socr.umich.edu/people/}>
}
