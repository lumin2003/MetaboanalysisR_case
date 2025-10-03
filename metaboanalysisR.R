library(MetaboAnalystR)
library("OptiLCMS")
#MS1 data
mSet <- PerformROIExtraction(datapath = "/.mzML files", rt.idx = 0.9, rmConts = TRUE);
best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "Q-TOF"), ncore = 64)
mSet <- ImportRawMSData(path = "/.mzML/", metadata="metadata.txt",plotSettings = SetPlotParam(Plot = T))
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=TRUE))
annParams <- SetAnnotationParam(polarity = 'positive/negative', mz_abs_add = 0.001) #check best_params
mSet <- PerformPeakAnnotation(mSet, annParams)
mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1)

#MS2 data #ppm1/2,rt_tol, sn, window_size are determinated based on best_param 
mSet <- PerformMSnImport(filesPath = c(list.files("/.mzML/",
                                                  pattern = ".mzML",
                                                  full.names = T, recursive = T)),
                                                  targetFeatures = ft_dt, #ft_dt from MS1 data
                                                  acquisitionMode = "DDA")
system.time(mSet <- PerformDDADeconvolution(mSet,
                                ppm1 = 5,
                                ppm2 = 10,
                                sn = 12,
                                filtering = 0,
                                window_size = 1.5,
                                intensity_thresh = 1.6e5, ###usually 5e4
                                database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                ncores = 6L))

mSet <- PerformSpectrumConsenus (mSet,
                                 ppm2 = 15, 
                                 concensus_fraction = 0.2,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)
mSet <- PerformDBSearchingBatch (mSet,
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                 use_rt = FALSE,
                                 enableNL = FALSE,
                                 ncores = 6L)
mSet <- PerformResultsExport (mSet,
                             type = 0L,
                             topN = 10L,
                             ncores = 3L)
dtx <- FormatMSnAnnotation(mSet, 5L, F)
write.csv(dtx, file="dtx.csv")

###R analysis the results from MS2
rm(list = ls())

mSet<-InitDataObjects("conc", "stat", FALSE)
Meta_data <- Read.TextData(mSet, "data_2.csv", "colu", "disc")
mSet<-SanityCheckData(Meta_data)
mSet<-ReplaceMin(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "_norm_0_", format ="png", dpi=72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "_snorm_0_", format = "png", dpi=72, width=NA)
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=TRUE, equal.var=TRUE, "fdr", FALSE)
mSet<-Volcano.Anal(mSet, TRUE, 2.0, 0, F, 0.1, TRUE, "raw")
mSet<-PlotVolcano(mSet, "_volcano_0_", 1, 0, format ="png", dpi=72, width=NA)
mSet <- ANOVA.Anal(mSet, F, 0.05, "fisher")

mSet<-PlotStaticHeatMap(mSet, "heatmap_1_", "png", 300, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, 8, "overview", T, T, NULL, T, F, T, T, T)
