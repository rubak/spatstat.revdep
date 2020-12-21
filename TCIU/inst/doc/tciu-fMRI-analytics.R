## ----warning=FALSE, message=FALSE---------------------------------------------
require(TCIU)
require(DT)
require(AnalyzeFMRI)

## -----------------------------------------------------------------------------
fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask, 
								   ons = c(1, 21, 41, 61, 81, 101, 121, 141), 
								   dur = c(10, 10, 10, 10, 10, 10, 10, 10))

# the outputs include simulated fMRI data, its mask, 
# the starting time points of the stimulated period and its duration 
# as well as all the stimulated time points
dim(fmri_generate$fmri_data)

## ----fig.width = 7, fig.align = "center", warning=FALSE, message=FALSE--------
fmri_time_series(sample[[5]], voxel_location = NULL, is.4d = FALSE, ref = sample[[4]])

## ----fig.width = 7, fig.align = "center", warning=FALSE, message=FALSE--------
# a data-frame with 160 rows and 4 columns: time (1:10), phases (8), states (2), and fMRI data (Complex or Real intensity)
datatable(fmri_kimesurface(fmri_generate$fmri_data, c(44,30,33))[[1]])
# ON Kime-Surface
fmri_kimesurface(fmri_generate$fmri_data, c(44,30,33))[[2]]
# User can try themself to plot the off / on&off figure
# OFF Kime-Surface 
# fmri_kimesurface(fmri_generate$fmri_data, c(44,30,33))[[3]]
# ON&OFF Kime-Surface 
# fmri_kimesurface(fmri_generate$fmri_data, c(44,30,33))[[4]]

## ----fig.width = 7, fig.align = "center", warning=FALSE, message=FALSE--------
fmri_image(fmri_generate$fmri_data, option="manually", voxel_location = c(40,22,33), time=4)

## ----fig.width = 7, fig.align = "center", warning=FALSE, message=FALSE--------
smoothmod<-GaussSmoothArray(fmri_generate$fmri_data, sigma = diag(3,3))
fmri_ts_forecast(smoothmod, voxel_location=c(41,44,33))

## -----------------------------------------------------------------------------
p_simulate_t_test = 
fmri_stimulus_detect(fmridata= fmri_generate$fmri_data, 
                     mask = fmri_generate$mask,
                     stimulus_idx = fmri_generate$on_time,
                     method = "t-test" , 
                     ons = fmri_generate$ons, 
                     dur = fmri_generate$dur)

dim(p_simulate_t_test)
summary(p_simulate_t_test)

## ---- eval = FALSE------------------------------------------------------------
#  # do the FDR correction
#  pval_fdr = fmri_post_hoc(phase2_pval , fdr_corr = "fdr",
#  						 spatial_cluster.thr = NULL,
#  						 spatial_cluster.size = NULL,
#  						 show_comparison = FALSE)
#  
#  # do the spatial clustering
#  pval_posthoc = fmri_post_hoc(pval_fdr, fdr_corr = NULL,
#  							 spatial_cluster.thr = 0.05,
#  							 spatial_cluster.size = 5,
#  							 show_comparison = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  # the output figure is hidden
#  for(axis in c("x", "y", "z")){
#    axis_i = switch(axis,
#                    "x" = {35},
#                    "y" = {30},
#                    "z" = {22})
#    print(fmri_2dvisual(p_simulate_t_test, list(axis, axis_i),
#                        hemody_data=NULL, mask=fmri_generate$mask,
#                        p_threshold = 0.05, legend_show = TRUE,
#                        method = "scale_p",
#                        color_pal = "YlOrRd", multi_pranges=TRUE))
#  }
#  			

## ----fig.width = 9, fig.align = "center", warning=FALSE-----------------------
fmri_3dvisual(p_simulate_t_test, fmri_generate$mask, 
							p_threshold = 0.05, method="scale_p",
              multi_pranges=TRUE)$plot

## ----eval = FALSE-------------------------------------------------------------
#  # the two p value are the p value generated based on the simulated fMRI
#  # and the p value saved in the package and finished post hoc test
#  # the output figure is hidden
#  fmri_pval_comparison_3d(list(p_simulate_t_test, phase3_pval), mask,
#  				                list(0.05, 0.05), list("scale_p", "scale_p"),
#  				                multi_pranges=FALSE)
#  

## ----fig.width = 9, fig.align = "center", warning=FALSE-----------------------
fmri_pval_comparison_2d(list(p_simulate_t_test, phase3_pval), 
                        list('pval_simulated', 'pval_posthoc'),
                        list(list(35, 33, 22), list(40, 26, 33)), 
                        hemody_data = NULL, 
                        mask = mask, p_threshold = 0.05, 
                        legend_show = FALSE, method = 'scale_p',
                        color_pal = "YlOrRd", multi_pranges=FALSE)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  ROI_phase1 = fmri_ROI_phase1(fmri_generate$fmri_data, mask_label, mask_dict, stimulus_idx = fmri_generate$on_time)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  ROI_phase2 = fmri_ROI_phase2(fmridata = fmri_generate$fmridata, label_mask = mask_label,
#                               ROI_label_dict = mask_dict, stimulus_idx = fmri_generate$on_time,
#                               stimulus_dur = fmri_generate$dur, rrr_rank = 3,
#                               fmri.design_order = 2, fmri.stimulus_TR = 3,
#                               method = "t_test", parallel_computing = TRUE, max(detectCores()-2,1))

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  # do the FDR correction
#  # do the spatial clustering
#  ROI_phase3 = fmri_post_hoc(ROI_phase2 , fdr_corr = "fdr",
#                             spatial_cluster.thr = 0.05,
#                             spatial_cluster.size = 5,
#                             show_comparison = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  # the output figure is hidden due to the size of vignettes
#  label_index = mask_dict$index
#  label_name = as.character(mask_dict$name)
#  label_mask = mask_label
#  fmri_3dvisual_region(phase1_pval, mask_label, label_index, label_name, title = "phase1 p-values")

## ----eval = FALSE-------------------------------------------------------------
#  # the output figure is hidden due to the size of vignettes
#  fmri_3dvisual_region(list(phase2_pval,phase3_pval), mask_label,
#                       label_index, label_name, title = "phase2&3 p-values")

