# Load required libraries
libraries <- c(
  'tidyverse',
  'lavaan',
  'semTools',
  'semPlot',
  'semptools',
  'parallel'
)

invisible(lapply(libraries, require, character.only = TRUE))

rm(libraries)

#//----- VARIABLE DICTIONARY -----

# wide0: SPRINT dataset in wide format with participants who completed baseline MRI scan that passed quality control and had complete data on SVD indicators
# wide01: SPRINT dataset in wide format with participants who completed both baseline and follow-up MRI scans that passed quality control and had complete data on SVD indicators
# r_pvwml.0: rescaled periventricular white matter hyperintensity volume at baseline
# r_pvwml.1: rescaled periventricular white matter hyperintensity volume at follow-up
# r_fw.0: rescaled mean white matter free water at baseline
# r_fw.1: rescaled mean white matter free water at follow-up
# r_bgepvsc.0: rescaled basal ganglia perivascular space count at baseline
# r_bgepvsc.1: rescaled basal ganglia perivascular space count at follow-up
# r_pvwmln.0: rescaled ROI volume-normalized periventricular white matter hyperintensity volume at baseline
# r_pvwmln.1: rescaled ROI volume-normalized periventricular white matter hyperintensity volume at follow-up
# r_bgepvscn.0: rescaled ROI volume-normalized basal ganglia perivascular space count at baseline
# r_bgepvscn.1: rescaled ROI volume-normalized basal ganglia perivascular space count at follow-up
# age.c.0: mean-centered baseline age
# time_years.1: time of follow-up (in years)
# treat: intensive (vs standard) BP group assignment binary indicator
# female: female sex binary indicator
# r_icv.0: rescaled total intracranial volume at baseline
# race: dummy variable vector for race/ethnicity with "White" as reference
# edu: dummy variable vector for education with "College degree" as reference
# smk: dummy variable vector for smoking with "Never" as reference
# polyph: dummy variable vector for polypharmacy with "<5 medications" as reference
# sub_cvd: subclinical cardiovascular disease binary indicator
# sbp: systolic blood pressure at baseline visit
# dbp: diastolic blood pressure at baseline visit
# BMI: body mass index
# HDL: fasting high-density lipoprotein cholesterol
# result_CO2: serum bicarbonate
# egfr: estimated glomerular filtration rate
# log2_umalcr: log urine albumin-to-creatinine ratio
# lm_delayed1: logical memory delayed score
# sbp_group_1: dummy variable indicating attained SBP change from baseline of 0-10 mmHg
# sbp_group_2: dummy variable indicating attained SBP change from baseline of 10-20 mmHg
# sbp_group_3: dummy variable indicating attained SBP change from baseline of ≥20 mmHg
# sbp_group_n: numeric attained SBP change group variable for linear trend calculation
# sbp_delta: attained SBP change from baseline (continuous variable)

#//--------------------------------------------------------------------------- START OF ANALYTIC CODE ---------------------------------------------------------------------------//

#//-------------------------------------------- ANALYSES 3.2: SVD MEASUREMENT MODEL AND LONGITUDINAL MEASUREMENT INVARIANCE --------------------------------------------

#//----- LONGITUDINAL INVARIANCE TESTING WITH RAW SVD INDICATORS -----

# Testing for configural, metric (weak), and scalar (strong) measurement invariance using longitudinal CFA models for multiple indicator data

# Configural invariance model (model 1)
model1 <-paste0('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw.0*r_fw.0 + bg.0*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw.1*r_fw.1 + bg.1*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvsc.0 ~ ibg.0*1
r_pvwml.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvsc.1 ~ ibg.1*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1
                ')

# Metric (weak) invariance model (model 2)
model2 <-paste0('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvsc.0 ~ ibg.0*1
r_pvwml.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvsc.1 ~ ibg.1*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1
                ')

# Scalar (strong) invariance model (model 3)
model3 <-('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv*1
r_fw.0 ~ ifw*1
r_bgepvsc.0 ~ ibg*1
r_pvwml.1 ~ ipv*1
r_fw.1 ~ ifw*1
r_bgepvsc.1 ~ ibg*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 1
          ')

# Fit longitudinal CFA models with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit1 <- lavaan::sem(model1, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
fit2 <- lavaan::sem(model2, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
fit3 <- lavaan::sem(model3, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")

# Compare model fit indices
fitm1 <- fitmeasures(fit1, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                            "cfi.robust", "rmsea.robust", "srmr"))

fitm2 <- fitmeasures(fit2, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                            "cfi.robust", "rmsea.robust", "srmr"))

fitm3 <- fitmeasures(fit3, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                            "cfi.robust", "rmsea.robust", "srmr"))

fitm_all <- as.data.frame(rbind(fitm1, fitm2, fitm3))
fitm_all <- fitm_all %>% mutate(across(where(is.numeric), function(x) {round(x, 3)}))

print(fitm_all) # Supplemental Table 4; Models 1, 2, and 3

# Display model summaries
summary(fit1, fit.measures = TRUE, standardized = TRUE)
summary(fit2, fit.measures = TRUE, standardized = TRUE)
summary(fit3, fit.measures = TRUE, standardized = TRUE)


#//----- LONGITUDINAL INVARIANCE TESTING WITH ROI VOLUME-NORMALIZED INDICATORS -----

# Testing for configural, metric (weak), and scalar (strong) measurement invariance using longitudinal CFA models for multiple indicator data

# Configural invariance model (model 4)
model4 <-paste0('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwmln.0 + fw.0*r_fw.0 + bg.0*r_bgepvscn.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwmln.1 + fw.1*r_fw.1 + bg.1*r_bgepvscn.1

# Mean structure specification
r_pvwmln.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvscn.0 ~ ibg.0*1
r_pvwmln.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvscn.1 ~ ibg.1*1

# Residual covariances
r_pvwmln.0 ~~ r_pvwmln.1
r_fw.0 ~~ r_fw.1
r_bgepvscn.0 ~~ r_bgepvscn.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1
                ')

# Metric (weak) invariance model (model 5)
model5 <-paste0('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwmln.0 + fw*r_fw.0 + bg*r_bgepvscn.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwmln.1 + fw*r_fw.1 + bg*r_bgepvscn.1

# Mean structure specification
r_pvwmln.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvscn.0 ~ ibg.0*1
r_pvwmln.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvscn.1 ~ ibg.1*1

# Residual covariances
r_pvwmln.0 ~~ r_pvwmln.1
r_fw.0 ~~ r_fw.1
r_bgepvscn.0 ~~ r_bgepvscn.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1
                ')

# Scalar (strong) invariance model (model 6)
model6 <-('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwmln.0 + fw*r_fw.0 + bg*r_bgepvscn.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwmln.1 + fw*r_fw.1 + bg*r_bgepvscn.1

# Mean structure specification
r_pvwmln.0 ~ ipv*1
r_fw.0 ~ ifw*1
r_bgepvscn.0 ~ ibg*1
r_pvwmln.1 ~ ipv*1
r_fw.1 ~ ifw*1
r_bgepvscn.1 ~ ibg*1

# Residual covariances
r_pvwmln.0 ~~ r_pvwmln.1
r_fw.0 ~~ r_fw.1
r_bgepvscn.0 ~~ r_bgepvscn.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 1
          ')

# Fit longitudinal CFA models with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit4 <- lavaan::sem(model4, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
fit5 <- lavaan::sem(model5, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
fit6 <- lavaan::sem(model6, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")

# Compare model fit indices
fitm4 <- fitmeasures(fit4, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                             "cfi.robust", "rmsea.robust", "srmr"))

fitm5 <- fitmeasures(fit5, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                             "cfi.robust", "rmsea.robust", "srmr"))

fitm6 <- fitmeasures(fit6, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                             "cfi.robust", "rmsea.robust", "srmr"))

fitm_all <- as.data.frame(rbind(fitm4, fitm5, fitm6))
fitm_all <- fitm_all %>% mutate(across(where(is.numeric), function(x) {round(x, 3)}))

print(fitm_all) # Supplemental Table 4; Models 4, 5, and 6

# Display model summaries
summary(fit4, fit.measures = TRUE, standardized = TRUE)
summary(fit5, fit.measures = TRUE, standardized = TRUE)
summary(fit6, fit.measures = TRUE, standardized = TRUE)


#//----- MULTIPLE INDICATORS MULTIPLE CAUSES (MIMIC) MODELS TO EXAMINE CONCEPTUALLY PLAUSIBLE SOURCES OF DIFFERENTIAL ITEM FUNCTIONING (DIF) OVER TIME -----

# MIMIC Model for periventricular white matter hyperintensity volume (pvwml)
# Metric (weak) invariance model with both direct and indirect (through SVD) effects of age on pvwml (modeldif_wml)
modeldif_wml <-paste0('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvsc.0 ~ ibg.0*1
r_pvwml.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvsc.1 ~ ibg.1*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1

# Regressions on age
svd_0 ~ age.c.0
svd_1 ~ age.c.0 + time_years.1
                
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
                      ')

# Fit longitudinal MIMIC model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fitdif_wml <- lavaan::sem(modeldif_wml, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fitdif_wml, fit.measures = TRUE, standardized = TRUE)

# Path diagram for pvwml MIMIC model (Supplementary Figure 1; panel A)
pl_mod <- semPlotModel(fitdif_wml)

pl_par <- pl_mod@Pars

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "int", ]), ]

pl_mod@Pars <- pl_par

plt <- semPaths(pl_mod, whatLabels = "est", residuals = FALSE, intStyle = "multi", sizeInt = 9, edge.label.cex = 1.2)

plt$layout[1, ] <- c(-1, 0.6) #p.0
plt$layout[2, ] <- c(-0.7, 0.6) #f.0
plt$layout[3, ] <- c(-0.4, 0.6) #b.0
plt$layout[4, ] <- c(0.2, 0.6) #p.1
plt$layout[6, ] <- c(0.5, 0.6) #f.1
plt$layout[7, ] <- c(0.8, 0.6) #b.1
plt$layout[5, ] <- c(-1.3, -0.2) #age.c.0
plt$layout[8, ] <- c(-0.1, -0.2) #dufu
plt$layout[9, ] <- c(-0.7, -0.8) #s_0
plt$layout[10, ] <- c(0.5, -0.8) #s_1

plt$graphAttributes$Nodes$labels <- c("Perive-\nticular\nWMH\nt0", "WM\nFree\nWater\nt0", "Basal\nganglia\nPVS\nt0",
                                      "Perive-\nticular\nWMH\nt1", "Base-\nline\nage", "WM\nFree\nWater\nt1",
                                      "Basal\nganglia\nPVS\nt1", "Follow\nup\ntime",
                                      "SVD\nt0", "SVD\nt1")

plt <- semptools::mark_sig(semPaths_plot = plt, object = fitdif_wml)

plt$graphAttributes$Nodes$label.cex <- 2

plt$graphAttributes$Edges$curve[c(7)] <- 2
plt$graphAttributes$Edges$curve[16] <- -1

plt$graphAttributes$Edges$label.margin[c(7:9)] <- -0.035

plt$graphAttributes$Edges$edge.label.position[14] <- 0.8

plt$graphAttributes$Edges$edge.label.position[17] <- 0.7

plt$plotOptions$label.prop <- c(0.8, 0.8, 0.8, 0.8, 0.7, 0.8, 0.8, 0.7, 0.6, 0.6)

plt$graphAttributes$Edges$color[c(1:6)] <- "coral"

plt$graphAttributes$Edges$color[c(7:9)] <- "lightgreen"

plt$graphAttributes$Edges$color[c(10:12)] <- "magenta2"

plt$graphAttributes$Edges$color[c(13:15)] <- "skyblue"

plot(plt)


# MIMIC Model for mean white matter free water (fw)
# Metric (weak) invariance model with both direct and indirect (through SVD) effects of age on fw (modeldif_fw)
modeldif_fw <-paste0('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvsc.0 ~ ibg.0*1
r_pvwml.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvsc.1 ~ ibg.1*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1

# Regressions on age
svd_0 ~ age.c.0
svd_1 ~ age.c.0 + time_years.1
                
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
                     ')

# Fit longitudinal MIMIC model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fitdif_fw <- lavaan::sem(modeldif_fw, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fitdif_fw, fit.measures = TRUE, standardized = TRUE)

# Path diagram for fw MIMIC model (Supplementary Figure 1; panel B)
pl_mod <- semPlotModel(fitdif_fw)

pl_par <- pl_mod@Pars

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "int", ]), ]

pl_mod@Pars <- pl_par

plt <- semPaths(pl_mod, whatLabels = "est", residuals = FALSE, intStyle = "multi", sizeInt = 9, edge.label.cex = 1.2)

plt$layout[1, ] <- c(-1, 0.6) #p.0
plt$layout[2, ] <- c(-0.7, 0.6) #f.0
plt$layout[3, ] <- c(-0.4, 0.6) #b.0
plt$layout[4, ] <- c(-1.3, -0.2) #age.c.0
plt$layout[5, ] <- c(0.5, 0.6) #f.1
plt$layout[6, ] <- c(0.2, 0.6) #p.1
plt$layout[7, ] <- c(0.8, 0.6) #b.1
plt$layout[8, ] <- c(-0.1, -0.2) #dufu
plt$layout[9, ] <- c(-0.7, -0.8) #s_0
plt$layout[10, ] <- c(0.5, -0.8) #s_1

plt$graphAttributes$Nodes$labels <- c("Perive-\nticular\nWMH\nt0", "WM\nFree\nWater\nt0", "Basal\nganglia\nPVS\nt0",
                                      "Base-\nline\nage",
                                      "WM\nFree\nWater\nt1", "Perive-\nticular\nWMH\nt1", "Basal\nganglia\nPVS\nt1",
                                      "Follow\nup\ntime",
                                      "SVD\nt0", "SVD\nt1")

plt <- semptools::mark_sig(semPaths_plot = plt, object = fitdif_fw)

plt$graphAttributes$Nodes$label.cex <- 2

plt$graphAttributes$Edges$curve[c(8)] <- 2
plt$graphAttributes$Edges$curve[16] <- -1

plt$graphAttributes$Edges$label.margin[c(7:9)] <- -0.035

plt$graphAttributes$Edges$edge.label.position[14] <- 0.72

plt$graphAttributes$Edges$edge.label.position[17] <- 0.7

plt$plotOptions$label.prop <- c(0.8, 0.8, 0.8, 0.7, 0.8, 0.8, 0.8, 0.7, 0.6, 0.6)

plt$graphAttributes$Edges$color[c(1:6)] <- "coral"

plt$graphAttributes$Edges$color[c(7:9)] <- "lightgreen"

plt$graphAttributes$Edges$color[c(10:12)] <- "magenta2"

plt$graphAttributes$Edges$color[c(13:15)] <- "skyblue"

plot(plt)


# MIMIC Model for basal ganglia perivascular space count (bgepvsc)
# Metric (weak) invariance model with both direct and indirect (through SVD) effects of age on bgepvsc (modeldif_pvs)
modeldif_pvs <-paste0('
# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvsc.0 ~ ibg.0*1
r_pvwml.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvsc.1 ~ ibg.1*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1

# Regressions on age
svd_0 ~ age.c.0
svd_1 ~ age.c.0 + time_years.1
                
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1
                      ')

# Fit longitudinal MIMIC model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fitdif_pvs <- lavaan::sem(modeldif_pvs, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fitdif_pvs, fit.measures = TRUE, standardized = TRUE)

# Path diagram for bgepvsc MIMIC model (Supplementary Figure 1; panel C)
pl_mod <- semPlotModel(fitdif_pvs)

pl_par <- pl_mod@Pars

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "int", ]), ]

pl_mod@Pars <- pl_par

plt <- semPaths(pl_mod, whatLabels = "est", residuals = FALSE, intStyle = "multi", sizeInt = 9, edge.label.cex = 1.2)

plt$layout[1, ] <- c(-1, 0.6) #p.0
plt$layout[2, ] <- c(-0.7, 0.6) #f.0
plt$layout[3, ] <- c(-0.4, 0.6) #b.0
plt$layout[4, ] <- c(-1.3, -0.2) #age.c.0
plt$layout[5, ] <- c(0.2, 0.6) #p.1
plt$layout[6, ] <- c(0.8, 0.6) #b.1
plt$layout[7, ] <- c(0.5, 0.6) #f.1
plt$layout[8, ] <- c(-0.1, -0.2) #dufu
plt$layout[9, ] <- c(-0.7, -0.8) #s_0
plt$layout[10, ] <- c(0.5, -0.8) #s_1

plt$graphAttributes$Nodes$labels <- c("Perive-\nticular\nWMH\nt0", "WM\nFree\nWater\nt0", "Basal\nganglia\nPVS\nt0",
                                      "Base-\nline\nage",
                                      "Perive-\nticular\nWMH\nt1", "Basal\nganglia\nPVS\nt1", "WM\nFree\nWater\nt1",
                                      "Follow\nup\ntime",
                                      "SVD\nt0", "SVD\nt1")

plt <- semptools::mark_sig(semPaths_plot = plt, object = fitdif_pvs)

plt$graphAttributes$Nodes$label.cex <- 2

plt$graphAttributes$Edges$curve[c(9)] <- 2
plt$graphAttributes$Edges$curve[16] <- -1

plt$graphAttributes$Edges$label.margin[c(7:9)] <- -0.035

plt$graphAttributes$Edges$edge.label.position[14] <- 0.6

plt$graphAttributes$Edges$edge.label.position[17] <- 0.7

plt$plotOptions$label.prop <- c(0.8, 0.8, 0.8, 0.7, 0.8, 0.8, 0.8, 0.7, 0.6, 0.6)

plt$graphAttributes$Edges$color[c(1:6)] <- "coral"

plt$graphAttributes$Edges$color[c(7:9)] <- "lightgreen"

plt$graphAttributes$Edges$color[c(10:12)] <- "magenta2"

plt$graphAttributes$Edges$color[c(13:15)] <- "skyblue"

plot(plt)


#//----- LONGITUDINAL INVARIANCE TESTING ACCOUNTING FOR TIME-VARYING AGE EFFECTS WITH RAW SVD INDICATORS -----

# Testing for configural, metric (weak), and scalar (strong) measurement invariance using longitudinal CFA models for multiple indicator data

# Configural invariance model (model 7)
model7 <-paste0('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw.0*r_fw.0 + bg.0*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw.1*r_fw.1 + bg.1*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvsc.0 ~ ibg.0*1
r_pvwml.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvsc.1 ~ ibg.1*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1
                ')

# Metric (weak) invariance model (model 8)
model8 <-paste0('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvsc.0 ~ ibg.0*1
r_pvwml.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvsc.1 ~ ibg.1*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1
                ')

# Scalar (strong) invariance model (model 9)
model9 <-('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ ipv*1
r_fw.0 ~ ifw*1
r_bgepvsc.0 ~ ibg*1
r_pvwml.1 ~ ipv*1
r_fw.1 ~ ifw*1
r_bgepvsc.1 ~ ibg*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 1
          ')

# Fit longitudinal CFA models with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit7 <- lavaan::sem(model7, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
fit8 <- lavaan::sem(model8, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
fit9 <- lavaan::sem(model9, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")

# Compare model fit indices
fitm7 <- fitmeasures(fit7, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                            "cfi.robust", "rmsea.robust", "srmr"))

fitm8 <- fitmeasures(fit8, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                            "cfi.robust", "rmsea.robust", "srmr"))

fitm9 <- fitmeasures(fit9, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                            "cfi.robust", "rmsea.robust", "srmr"))

fitm_all <- as.data.frame(rbind(fitm7, fitm8, fitm9))
fitm_all <- fitm_all %>% mutate(across(where(is.numeric), function(x) {round(x, 3)}))

print(fitm_all) # Supplemental Table 4; Models 7, 8, and 9

# Display model summaries
summary(fit7, fit.measures = TRUE, standardized = TRUE)
summary(fit8, fit.measures = TRUE, standardized = TRUE)
summary(fit9, fit.measures = TRUE, standardized = TRUE)


#//----- PARAMETER ESTIMATES FOR THE COVARIANCE AND MEAN STRUCTURE OF THE SVD MEASUREMENT MODEL -----

# Supplemental Tables 5 and 6

# Model summary
sum <- summary(fit9, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# Model diagnostics
lavInspect(fit9, "sampstat") # sample covariance matrix
lavInspect(fit9, "implied") # model-implied covariance matrix
lavInspect(fit9, "resid") # difference between observed and model-implied covariance matrix (unstandardized and unscaled model residuals)
lavResiduals(fit9, type = "cor.bentler")$cov # unstandardized model residuals after transformation to correlation matrix and rescaling (by dividing the elements by the square roots of the corresponding variances of the observed covariance matrix)
lavResiduals(fit9, type = "cor.bentler")$cov.z # standardized model residuals after transformation to correlation matrix and rescaling (by dividing the elements by the square roots of the corresponding variances of the observed covariance matrix)
modindices(fit9, sort. = TRUE, standardized = TRUE) # sorted (from largest to smallest) model modification indices
lavTestScore(fit9, cumulative = TRUE) # Score test (or Lagrange Multiplier test) for releasing one or more fixed or constrained parameters in model

# Display baseline and follow-up loadings
load <- sum$pe[sum$pe$op == "=~", c("rhs", "est", "se", "std.all")]
load <- load %>% mutate(across(where(is.numeric), function(x) {formatC(x, format = "f", digits = 3)}))
print(load)

# Display baseline and follow-up age regression coefficients
reg <- sum$pe[sum$pe$op == "~", c("lhs", "rhs", "est", "se", "std.all")]
reg <- reg %>% mutate(across(where(is.numeric), function(x) {formatC(x, format = "f", digits = 3)}))
print(reg)

# Display baseline and follow-up residual (error) variances
er <- sum$pe[((sum$pe$op == "~~") & (sum$pe$lhs == sum$pe$rhs)), c("lhs", "rhs", "est", "se", "std.all")]
er <- er %>% mutate(across(where(is.numeric), function(x) {formatC(x, format = "f", digits = 3)}))
print(er)

# Display factor covariance
cov <- sum$pe[((sum$pe$op == "~~") & (sum$pe$lhs != sum$pe$rhs)), c("lhs", "op", "rhs", "est", "se", "std.all")]
cov <- cov %>% mutate(across(where(is.numeric), function(x) {formatC(x, format = "f", digits = 3)}))
print(cov)

# Display intercepts
int <- sum$pe[(sum$pe$op == "~1"), c("lhs", "op", "rhs", "est", "se", "std.all")]
int <- int %>% mutate(across(where(is.numeric), function(x) {formatC(x, format = "f", digits = 3)}))
print(int)


#//----- LONGITUDINAL INVARIANCE TESTING ACCOUNTING FOR TIME-VARYING AGE EFFECTS WITH ROI VOLUME-NORMALIZED SVD INDICATORS -----

# Testing for configural, metric (weak), and scalar (strong) measurement invariance using longitudinal CFA models for multiple indicator data

# Configural invariance model (model 10)
model10 <-paste0('
r_pvwmln.0 ~ age.c.0
r_pvwmln.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvscn.0 ~ age.c.0
r_bgepvscn.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwmln.0 + fw.0*r_fw.0 + bg.0*r_bgepvscn.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwmln.1 + fw.1*r_fw.1 + bg.1*r_bgepvscn.1

# Mean structure specification
r_pvwmln.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvscn.0 ~ ibg.0*1
r_pvwmln.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvscn.1 ~ ibg.1*1

# Residual covariances
r_pvwmln.0 ~~ r_pvwmln.1
r_fw.0 ~~ r_fw.1
r_bgepvscn.0 ~~ r_bgepvscn.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1
                 ')

# Metric (weak) invariance model (model 11)
model11 <-paste0('
r_pvwmln.0 ~ age.c.0
r_pvwmln.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvscn.0 ~ age.c.0
r_bgepvscn.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwmln.0 + fw*r_fw.0 + bg*r_bgepvscn.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwmln.1 + fw*r_fw.1 + bg*r_bgepvscn.1

# Mean structure specification
r_pvwmln.0 ~ ipv.0*1
r_fw.0 ~ ifw.0*1
r_bgepvscn.0 ~ ibg.0*1
r_pvwmln.1 ~ ipv.1*1
r_fw.1 ~ ifw.1*1
r_bgepvscn.1 ~ ibg.1*1

# Residual covariances
r_pvwmln.0 ~~ r_pvwmln.1
r_fw.0 ~~ r_fw.1
r_bgepvscn.0 ~~ r_bgepvscn.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 0*1
                 ')

# Scalar (strong) invariance model (model 12)
model12 <-('
r_pvwmln.0 ~ age.c.0
r_pvwmln.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvscn.0 ~ age.c.0
r_bgepvscn.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwmln.0 + fw*r_fw.0 + bg*r_bgepvscn.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwmln.1 + fw*r_fw.1 + bg*r_bgepvscn.1

# Mean structure specification
r_pvwmln.0 ~ ipv*1
r_fw.0 ~ ifw*1
r_bgepvscn.0 ~ ibg*1
r_pvwmln.1 ~ ipv*1
r_fw.1 ~ ifw*1
r_bgepvscn.1 ~ ibg*1

# Residual covariances
r_pvwmln.0 ~~ r_pvwmln.1
r_fw.0 ~~ r_fw.1
r_bgepvscn.0 ~~ r_bgepvscn.1

# Latent variable means
svd_0 ~ 0*1
svd_1 ~ 1
           ')

# Fit longitudinal CFA models with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit10 <- lavaan::sem(model10, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
fit11 <- lavaan::sem(model11, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
fit12 <- lavaan::sem(model12, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")

# Compare model fit indices
fitm10 <- fitmeasures(fit10, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                             "cfi.robust", "rmsea.robust", "srmr"))

fitm11 <- fitmeasures(fit11, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                             "cfi.robust", "rmsea.robust", "srmr"))

fitm12 <- fitmeasures(fit12, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                             "cfi.robust", "rmsea.robust", "srmr"))

fitm_all <- as.data.frame(rbind(fitm10, fitm11, fitm12))
fitm_all <- fitm_all %>% mutate(across(where(is.numeric), function(x) {round(x, 3)}))

print(fitm_all) # Supplemental Table 4; Models 10, 11, and 12

# Display model summaries
summary(fit10, fit.measures = TRUE, standardized = TRUE)
summary(fit11, fit.measures = TRUE, standardized = TRUE)
summary(fit12, fit.measures = TRUE, standardized = TRUE)


#//----- MULTIPLE-INDICATOR LATENT CHANGE SCORE (LCS) MODELS -----

# LCS model equivalent to longitudinal CFA model 9
model13 <-paste0('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ 0*1
r_fw.0 ~ ifw*1
r_bgepvsc.0 ~ ibg*1
r_pvwml.1 ~ 0*1
r_fw.1 ~ ifw*1
r_bgepvsc.1 ~ ibg*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent change
svd_1 ~ 1*svd_0
dC =~ 1*svd_1

# Latent variable mean structure
svd_0 ~ 1
svd_1 ~ 0*1
dC ~ 1

# Latent variable covariance structure
svd_0 ~~ svd_0
svd_1 ~~ 0*svd_1
dC ~~ dC
dC ~~ svd_0
                 ')

# Fit LCS model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit13 <- lavaan::sem(model13, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fit13, fit.measures = TRUE, standardized = TRUE)

# Compare models (sanity check)
lavTestLRT(fit9, fit13)

# Compare model fit indices (sanity check)
fitm9 <- fitmeasures(fit9, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                            "cfi.robust", "tli.robust",
                                            "rmsea.robust", "srmr"))

fitm13 <- fitmeasures(fit13, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                              "cfi.robust", "tli.robust",
                                              "rmsea.robust", "srmr"))

fitm_all <- as.data.frame(rbind(fitm9, fitm13))
fitm_all <- fitm_all %>% mutate(across(where(is.numeric), function(x) {round(x, 3)}))

print(fitm_all) 


#//-------------------------------------------- ANALYSES 3.3: EFFECT OF TREATMENT ASSIGNMENT ON SVD BURDEN --------------------------------------------

#//----- COMPUTE UNADJUSTED TREATMENT EFFECT -----
# LCS model with treatment arm as latent change score predictor (unadjusted)
model1 <-paste0('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ 0*1
r_fw.0 ~ ifw*1
r_bgepvsc.0 ~ ibg*1
r_pvwml.1 ~ 0*1
r_fw.1 ~ ifw*1
r_bgepvsc.1 ~ ibg*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent change
svd_1 ~ 1*svd_0
dC =~ 1*svd_1

# Latent variable mean structure
svd_0 ~ a0*1
svd_1 ~ 0*1
dC ~ 1 + a*treat

# Latent variable covariance structure
svd_0 ~~ svd_0
svd_1 ~~ 0*svd_1
dC ~~ b*dC
dC ~~ svd_0

# Standardized mean difference (Cohens d) and % change relative to mean baseline SVD burden
cohen_d := a/sqrt(b)
prop := a/a0
                ')

# Fit LCS model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit1 <- lavaan::sem(model1, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fit1, fit.measures = TRUE, standardized = TRUE)

# Inspect fit indices
fitmeasures(fit1, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                   "cfi.robust", "rmsea.robust", "srmr"))

# Display parameter estimates of interest
parameterEstimates(fit1)[parameterEstimates(fit1)$label %in%
                           c("a", "cohen_d", "prop"),
                         c("lhs", "est", "ci.lower", "ci.upper")]

# Display latent change score variance
parameterEstimates(fit1)[parameterEstimates(fit1)$label %in% "b", "est"]

# Path diagram for the unadjusted treatment effect (Figure 2; panel A)
pl_mod <- semPlotModel(fit1)

pl_par <- pl_mod@Pars

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "int" & !pl_par$rhs %in% c("dC", "svd_0"), ]), ]

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "<->" & pl_par$rhs %in% c("treat"), ]), ]

pl_mod@Pars <- pl_par

plt <- semPaths(pl_mod, whatLabels = "est", residuals = FALSE, intStyle = "multi", sizeInt = 9, edge.label.cex = 1.2)

plt$layout[1, ] <- c(-1, 0.6) #p.0
plt$layout[2, ] <- c(-0.7, 0.6) #f.0
plt$layout[3, ] <- c(-0.4, 0.6) #b.0
plt$layout[4, ] <- c(0.2, 0.6) #p.1
plt$layout[5, ] <- c(0.5, 0.6) #f.1
plt$layout[6, ] <- c(0.8, 0.6) #b.1
plt$layout[7, ] <- c(0.3, -1) #trt
plt$layout[8, ] <- c(-0.7, 1.35) #age.c.0
plt$layout[9, ] <- c(0.5, 1.35) #dufu
plt$layout[10, ] <- c(-0.7, -0.5) #s_0
plt$layout[11, ] <- c(0.5, -0.5) #s_1
plt$layout[12, ] <- c(-0.2, -1) #dC
plt$layout[13, ] <- c(-1.2, -0.5) #int.s_0
plt$layout[14, ] <- c(-0.8, -1) #int.dC

plt$graphAttributes$Nodes$labels <- c("Perive-\nticular\nWMH\nt0", "WM\nFree\nWater\nt0", "Basal\nganglia\nPVS\nt0",
                                      "Perive-\nticular\nWMH\nt1", "WM\nFree\nWater\nt1", "Basal\nganglia\nPVS\nt1",
                                      "Intensive\nBP\nControl",
                                      "Base-\nline\nage", "Follow\nup\ntime",
                                      "SVD\nt0", "SVD\nt1", "SVD\nChange",
                                      "Mean\nBaseline\nSVD", "Mean\nSVD\nChange")

plt <- semptools::mark_sig(semPaths_plot = plt, object = fit1)

plt$graphAttributes$Nodes$label.cex <- 1.6

plt$graphAttributes$Edges$curve[c(24,25)] <- 0.5

plt$graphAttributes$Edges$label.margin[c(16:18)] <- -0.035
plt$graphAttributes$Edges$label.margin[c(21,23)] <- -0.04
plt$graphAttributes$Edges$label.margin[c(22)] <- -0.03

plt$plotOptions$label.prop <- c(1, 1, 1, 1, 1, 1, 0.95, 0.85, 0.85, 0.7, 0.7, 0.85, 0.45, 0.425)

plt$graphAttributes$Edges$color[c(1,2,4,5,7,8)] <- "skyblue"

plt$graphAttributes$Edges$color[c(3,6,9)] <- "navy"

plt$graphAttributes$Edges$color[c(10:15)] <- "coral"

plt$graphAttributes$Edges$color[c(16:18)] <- "lightgreen"

plt$graphAttributes$Edges$color[c(21:22)] <- "blue"

plt$graphAttributes$Edges$color[c(23)] <- "red2"

plt$graphAttributes$Nodes$width <- c(5, 5, 5, 5, 5, 5, 7, 5, 5, 8, 8, 8, 12, 12)

plot(plt)


#//----- COMPUTE ADJUSTED TREATMENT EFFECT -----
# LCS model with treatment arm as latent change score predictor adjusted for baseline age, sex, race, and icv
model2 <-paste0('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ 0*1
r_fw.0 ~ ifw*1
r_bgepvsc.0 ~ ibg*1
r_pvwml.1 ~ 0*1
r_fw.1 ~ ifw*1
r_bgepvsc.1 ~ ibg*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent change
svd_1 ~ 1*svd_0
dC =~ 1*svd_1

# Latent variable mean structure
svd_0 ~ a0*1 + age.c.0 + female + r_icv.0 + ', paste(c(race), collapse = " + "), '
svd_1 ~ 0*1
dC ~ 1 + a*treat + age.c.0 + female + r_icv.0 + ', paste(c(race), collapse = " + "), '

# Latent variable covariance structure
svd_0 ~~ svd_0
svd_1 ~~ 0*svd_1
dC ~~ b*dC
dC ~~ svd_0

# Standardized mean difference (Cohens d)
cohen_d := a/sqrt(b)
                ')

# Fit LCS model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit2 <- lavaan::sem(model2, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fit2, fit.measures = TRUE, standardized = TRUE)

# Inspect fit indices
fitmeasures(fit2, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                   "cfi.robust", "rmsea.robust", "srmr"))

# Display parameter estimates of interest
parameterEstimates(fit2)[parameterEstimates(fit2)$label %in%
                           c("a", "cohen_d"),
                         c("lhs", "est", "ci.lower", "ci.upper")]

# Display latent change score variance
parameterEstimates(fit2)[parameterEstimates(fit2)$label %in% "b", "est"]

# Path diagram for the adjusted treatment effect (Figure 2; panel B)
pl_mod <- semptools::drop_nodes(
  object = semPlotModel(fit2),
  nodes = c(race, "female", "r_icv.0"))

pl_par <- pl_mod@Pars

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "int" & !pl_par$rhs %in% c("dC", "svd_0"), ]), ]

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "<->" & pl_par$rhs %in% c("treat"), ]), ]

pl_mod@Pars <- pl_par

plt <- semPaths(pl_mod, whatLabels = "est", residuals = FALSE, intStyle = "multi", sizeInt = 9, edge.label.cex = 1.2)

plt$layout[1, ] <- c(-1, 0.6) #p.0
plt$layout[2, ] <- c(-0.7, 0.6) #f.0
plt$layout[3, ] <- c(-0.4, 0.6) #b.0
plt$layout[4, ] <- c(0.2, 0.6) #p.1
plt$layout[5, ] <- c(0.5, 0.6) #f.1
plt$layout[6, ] <- c(0.8, 0.6) #b.1
plt$layout[7, ] <- c(-0.7, 1.35) #age.c.0
plt$layout[8, ] <- c(0.3, -1) #trt
plt$layout[9, ] <- c(0.5, 1.35) #dufu
plt$layout[10, ] <- c(-0.7, -0.5) #s_0
plt$layout[11, ] <- c(0.5, -0.5) #s_1
plt$layout[12, ] <- c(-0.2, -1) #dC
plt$layout[13, ] <- c(-1.2, -0.5) #int.s_0
plt$layout[14, ] <- c(-0.8, -1) #int.dC

plt$graphAttributes$Nodes$labels <- c("Perive-\nticular\nWMH\nt0", "WM\nFree\nWater\nt0", "Basal\nganglia\nPVS\nt0",
                                      "Perive-\nticular\nWMH\nt1", "WM\nFree\nWater\nt1", "Basal\nganglia\nPVS\nt1",
                                      "Base-\nline\nage", "Intensive\nBP\nControl", "Follow\nup\ntime",
                                      "SVD\nt0", "SVD\nt1", "SVD\nChange",
                                      "Mean\nBaseline\nSVD", "Mean\nSVD\nChange")

plt <- semptools::mark_sig(semPaths_plot = plt, object = fit2)

plt$graphAttributes$Nodes$label.cex <- 1.6

plt$graphAttributes$Edges$curve[c(22)] <- -6

plt$graphAttributes$Edges$curve[c(26,27)] <- 0.5

plt$graphAttributes$Edges$label.margin[c(16:18)] <- -0.035
plt$graphAttributes$Edges$label.margin[c(21,23,24)] <- -0.04
plt$graphAttributes$Edges$label.margin[c(23)] <- -0.03

plt$graphAttributes$Edges$edge.label.position[25] <- 0.7

plt$plotOptions$label.prop <- c(1, 1, 1, 1, 1, 1, 0.85, 0.95, 0.85, 0.7, 0.7, 0.85, 0.45, 0.425)

plt$graphAttributes$Edges$color[c(1,2,4,5,7,8)] <- "skyblue"

plt$graphAttributes$Edges$color[c(3,6,9)] <- "navy"

plt$graphAttributes$Edges$color[c(10:15)] <- "coral"

plt$graphAttributes$Edges$color[c(16:18)] <- "lightgreen"

plt$graphAttributes$Edges$color[c(21,23)] <- "blue"

plt$graphAttributes$Edges$color[c(24)] <- "red2"

plt$graphAttributes$Edges$color[c(22, 25)] <- "magenta2"

plt$graphAttributes$Nodes$width <- c(5, 5, 5, 5, 5, 5, 5, 7, 5, 8, 8, 8, 12, 12)

plot(plt)


#//----- COMPUTE TREATMENT EFFECT WITH ROI VOLUME-NORMALIZED SVD INDICATORS -----
# LCS model with treatment arm as latent change score predictor (unadjusted) using ROI volume-adjusted MRI indicators
model3 <-paste0('
r_pvwmln.0 ~ age.c.0
r_pvwmln.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvscn.0 ~ age.c.0
r_bgepvscn.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwmln.0 + fw*r_fw.0 + bg*r_bgepvscn.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwmln.1 + fw*r_fw.1 + bg*r_bgepvscn.1

# Mean structure specification
r_pvwmln.0 ~ 0*1
r_fw.0 ~ ifw*1
r_bgepvscn.0 ~ ibg*1
r_pvwmln.1 ~ 0*1
r_fw.1 ~ ifw*1
r_bgepvscn.1 ~ ibg*1

# Residual covariances
r_pvwmln.0 ~~ r_pvwmln.1
r_fw.0 ~~ r_fw.1
r_bgepvscn.0 ~~ r_bgepvscn.1

# Latent change
svd_1 ~ 1*svd_0
dC =~ 1*svd_1

# Latent variable mean structure
svd_0 ~ a0*1
svd_1 ~ 0*1
dC ~ 1 + a*treat

# Latent variable covariance structure
svd_0 ~~ svd_0
svd_1 ~~ 0*svd_1
dC ~~ b*dC
dC ~~ svd_0

# Standardized mean difference (Cohens d) and % change relative to mean baseline SVD burden
cohen_d := a/sqrt(b)
prop := a/a0
                ')

# Fit LCS model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit3 <- lavaan::sem(model3, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fit3, fit.measures = TRUE, standardized = TRUE)

# Inspect fit indices
fitmeasures(fit3, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                   "cfi.robust", "rmsea.robust", "srmr"))

# Display parameter estimates of interest
parameterEstimates(fit3)[parameterEstimates(fit3)$label %in%
                           c("a", "cohen_d", "prop"),
                         c("lhs", "est", "ci.lower", "ci.upper")]

# Display latent change score variance
parameterEstimates(fit3)[parameterEstimates(fit3)$label %in% "b", "est"]


#//----- COMPUTE TREATMENT EFFECT IN ALL PARTICIPANTS WITH AVAILABLE BASELINE MRI WITH FIML -----
# LCS model for the total sample of participants with baseline MRI with treatment arm as latent change score predictor (unadjusted) with fiml for missing follow-up MRI markers + auxiliary variables

# Create dataset for fiml
fiml <- wide0

# Base auxiliary variable set
base_set <- c('female', race, edu, polyph, smk, 'sub_cvd', 'sbp', 'dbp')
# Inspect for missing values
sapply(base_set, function(x){sum(is.na(fiml[ ,x]))})

# Extended auxiliary variable set
ext_set <- c("BMI", "HDL", "result_CO2", "egfr", "log2_umalcr", "lm_delayed1", "r_icv.0")
# Inspect for missing values
sapply(ext_set, function(x){sum(is.na(fiml[ ,x]))})

# Fit LCS model with Full Information Maximum Likelihood estimation method with auxiliary variables using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit.aux <- auxiliary(model1, aux = c(base_set, ext_set), fun = "sem", data = fiml, estimator = "ML", missing = "FIML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fit.aux, fit.measures = TRUE, standardized = TRUE)

fitmeasures(fit.aux, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                      "cfi.robust", "rmsea.robust", "srmr"))

parameterEstimates(fit.aux)[parameterEstimates(fit.aux)$label %in%
                              c("a", "cohen_d", "prop"),
                            c("lhs", "est", "ci.lower", "ci.upper")]


#//-------------------------------------------- ANALYSES 3.4: ATTAINED SBP REDUCTION AND SVD CHANGE --------------------------------------------
# LCS model with SBP delta group as latent change score predictor
model4 <-paste0('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ 0*1
r_fw.0 ~ ifw*1
r_bgepvsc.0 ~ ibg*1
r_pvwml.1 ~ 0*1
r_fw.1 ~ ifw*1
r_bgepvsc.1 ~ ibg*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent change
svd_1 ~ 1*svd_0
dC =~ 1*svd_1

# Latent variable mean structure
svd_0 ~ a0*1
svd_1 ~ 0*1
dC ~ 1 + a1*sbp_group_1 + a2*sbp_group_2 + a3*sbp_group_3

# Latent variable covariance structure
svd_0 ~~ svd_0
svd_1 ~~ 0*svd_1
dC ~~ dC
dC ~~ svd_0

# % change relative to mean baseline SVD burden
prop_1 := a1/a0
prop_2 := a2/a0
prop_3 := a3/a0
                ')

# Fit LCS model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit4 <- lavaan::sem(model4, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fit4, fit.measures = TRUE, standardized = TRUE)

# Inspect fit indices
fitmeasures(fit4, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                   "cfi.robust", "rmsea.robust", "srmr"))

# Display parameter estimates of interest
parameterEstimates(fit4)[parameterEstimates(fit4)$label %in%
                           c("a1", "a2", "a3",
                             "prop_1", "prop_2", "prop_3"),
                         c("lhs", "est", "ci.lower", "ci.upper")]

# Path diagram for different SBP delta groups and SVD burden change (Figure 3)
pl_mod <- semPlotModel(fit4)

pl_par <- pl_mod@Pars

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "int" & !pl_par$rhs %in% c("dC", "svd_0"), ]), ]

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "<->" & pl_par$rhs %in%
                                                          c("sbp_group_1", "sbp_group_2", "sbp_group_3"), ]), ]

pl_mod@Pars <- pl_par

plt <- semPaths(pl_mod, whatLabels = "est", residuals = FALSE, intStyle = "multi", sizeInt = 9, edge.label.cex = 1.2)

plt$layout[1, ] <- c(-1, 0.6) #p.0
plt$layout[2, ] <- c(-0.7, 0.6) #f.0
plt$layout[3, ] <- c(-0.4, 0.6) #b.0
plt$layout[4, ] <- c(0.2, 0.6) #p.1
plt$layout[5, ] <- c(0.5, 0.6) #f.1
plt$layout[6, ] <- c(0.8, 0.6) #b.1

plt$layout[7, ] <- c(0.8, -0.7) #sbp_1
plt$layout[8, ] <- c(0.8, -1.0) #sbp_2
plt$layout[9, ] <- c(0.8, -1.3) #sbp_3

plt$layout[10, ] <- c(-0.7, 1.35) #age.c.0
plt$layout[11, ] <- c(0.5, 1.35) #dufu
plt$layout[12, ] <- c(-0.7, -0.5) #s_0
plt$layout[13, ] <- c(0.5, -0.5) #s_1
plt$layout[14, ] <- c(-0.2, -1) #dC
plt$layout[15, ] <- c(-1.2, -0.5) #int.s_0
plt$layout[16, ] <- c(-0.8, -1) #int.dC

plt$graphAttributes$Nodes$labels <- c("Perive-\nticular\nWMH\nt0", "WM\nFree\nWater\nt0", "Basal\nganglia\nPVS\nt0",
                                      "Perive-\nticular\nWMH\nt1", "WM\nFree\nWater\nt1", "Basal\nganglia\nPVS\nt1",
                                      "SBP\nDrop\n0 to 10", "SBP\nDrop\n10 to 20", "SBP\nDrop\n≥ 20",
                                      "Base-\nline\nage", "Follow\nup\ntime",
                                      "SVD\nt0", "SVD\nt1", "SVD\nChange",
                                      "Mean\nBaseline\nSVD", "Mean\nSVD\nChange")

plt <- semptools::mark_sig(semPaths_plot = plt, object = fit4)

plt$graphAttributes$Nodes$label.cex <- 1.6

plt$graphAttributes$Edges$curve[c(26,27)] <- 0.5

plt$graphAttributes$Edges$label.margin[c(16:18)] <- -0.035
plt$graphAttributes$Edges$label.margin[c(22, 23)] <- -0.02
plt$graphAttributes$Edges$label.margin[c(21, 24, 25)] <- -0.03

plt$plotOptions$label.prop <- c(1, 1, 1,
                                1, 1, 1,
                                0.9, 1, 0.9,
                                0.85, 0.85,
                                0.7, 0.7, 0.85,
                                0.45, 0.425)

plt$graphAttributes$Edges$color[c(1,2,4,5,7,8)] <- "skyblue"

plt$graphAttributes$Edges$color[c(3,6,9)] <- "navy"

plt$graphAttributes$Edges$color[c(10:15)] <- "coral"

plt$graphAttributes$Edges$color[c(16:18)] <- "lightgreen"

plt$graphAttributes$Edges$color[c(21:22)] <- "blue"

plt$graphAttributes$Edges$color[c(23:25)] <- "red2"

plt$graphAttributes$Nodes$width <- c(5, 5, 5,
                                     5, 5, 5,
                                     5, 5, 5,
                                     5, 5,
                                     8, 8, 8,
                                     12, 12)

plot(plt)


#//----- TREND TEST -----
# LCS model with BP delta group as latent change score predictor - trend test
model5 <-paste0('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ 0*1
r_fw.0 ~ ifw*1
r_bgepvsc.0 ~ ibg*1
r_pvwml.1 ~ 0*1
r_fw.1 ~ ifw*1
r_bgepvsc.1 ~ ibg*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent change
svd_1 ~ 1*svd_0
dC =~ 1*svd_1

# Latent variable mean structure
svd_0 ~ 1
svd_1 ~ 0*1
dC ~ 1 + sbp_group_n

# Latent variable covariance structure
svd_0 ~~ svd_0
svd_1 ~~ 0*svd_1
dC ~~ dC
dC ~~ svd_0
                ')

# Fit LCS model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit5 <- lavaan::sem(model5, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fit5, fit.measures = TRUE, standardized = TRUE)

# Inspect fit indices
fitmeasures(fit5, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                   "cfi.robust", "rmsea.robust", "srmr"))

# Extract parameter estimates of interest
par <- parameterEstimates(fit5)
par[par$lhs == "dC" & grepl("sbp_group_n", par$rhs),]


#//-------------------------------------------- ANALYSES 3.4: BP MEDIATION OF THE TREATMENT ASSIGNMENT EFFECT --------------------------------------------
# LCS model with treatment group as latent change score predictor and delta SBP as mediator
model6 <-paste0('
r_pvwml.0 ~ age.c.0
r_pvwml.1 ~ age.c.0 + time_years.1
r_fw.0 ~ age.c.0
r_fw.1 ~ age.c.0 + time_years.1
r_bgepvsc.0 ~ age.c.0
r_bgepvsc.1 ~ age.c.0 + time_years.1

# Measurement model for svd in time 1
 svd_0 =~ 1*r_pvwml.0 + fw*r_fw.0 + bg*r_bgepvsc.0

# Measurement model for svd in time 2
 svd_1 =~ 1*r_pvwml.1 + fw*r_fw.1 + bg*r_bgepvsc.1

# Mean structure specification
r_pvwml.0 ~ 0*1
r_fw.0 ~ ifw*1
r_bgepvsc.0 ~ ibg*1
r_pvwml.1 ~ 0*1
r_fw.1 ~ ifw*1
r_bgepvsc.1 ~ ibg*1

# Residual covariances
r_pvwml.0 ~~ r_pvwml.1
r_fw.0 ~~ r_fw.1
r_bgepvsc.0 ~~ r_bgepvsc.1

# Latent change
svd_1 ~ 1*svd_0
dC =~ 1*svd_1

# Latent variable mean structure
svd_0 ~ 1
svd_1 ~ 0*1
sbp_delta ~ a*treat
dC ~ 1 + b*sbp_delta + c*treat

# Latent variable covariance structure
svd_0 ~~ svd_0
svd_1 ~~ 0*svd_1
dC ~~ dC
dC ~~ svd_0
                
Indirect := a*b
Direct := c
Total := (a*b) + c
                ')

# Fit LCS model with Maximum Likelihood estimation method using robust (Huber-White) standard errors and a scaled (Yuan-Bentler) test statistic
fit6 <- lavaan::sem(model6, data = wide01, estimator = "ML", se = "robust.huber.white", test = "yuan.bentler.mplus")
summary(fit6, fit.measures = TRUE, standardized = TRUE)

# Inspect fit indices
fitmeasures(fit6, fit.measures = c("chisq.scaled", "pvalue.scaled",
                                   "cfi.robust", "rmsea.robust", "srmr"))

# Extract parameter estimates of interest
parameterEstimates(fit6)[parameterEstimates(fit6)$op %in% ":=",
                         c("lhs", "est", "ci.lower", "ci.upper")]

# Path diagram of indirect (through delta SBP) and direct treatment effect (Supplemental Figure 3)
pl_mod <- semPlotModel(fit6)

pl_par <- pl_mod@Pars

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "int" & !pl_par$rhs %in% c("dC", "svd_0"), ]), ]

pl_par <- pl_par[!rownames(pl_par) %in% rownames(pl_par[pl_par$edge == "<->" & pl_par$rhs %in%
                                                          c("sbp_delta", "treat"), ]), ]

pl_mod@Pars <- pl_par

plt <- semPaths(pl_mod, whatLabels = "est", residuals = FALSE, intStyle = "multi", sizeInt = 9, edge.label.cex = 1.2)

plt$layout[1, ] <- c(-1, 0.6) #p.0
plt$layout[2, ] <- c(-0.7, 0.6) #f.0
plt$layout[3, ] <- c(-0.4, 0.6) #b.0
plt$layout[4, ] <- c(0.2, 0.6) #p.1
plt$layout[5, ] <- c(0.5, 0.6) #f.1
plt$layout[6, ] <- c(0.8, 0.6) #b.1

plt$layout[7, ] <- c(0.8, -0.7) #delta_sbp
plt$layout[8, ] <- c(0.8, -1.3) #treat

plt$layout[9, ] <- c(-0.7, 1.35) #age.c.0
plt$layout[10, ] <- c(0.5, 1.35) #dufu
plt$layout[11, ] <- c(-0.7, -0.5) #s_0
plt$layout[12, ] <- c(0.5, -0.5) #s_1
plt$layout[13, ] <- c(-0.2, -1) #dC
plt$layout[14, ] <- c(-1.2, -0.5) #int.s_0
plt$layout[15, ] <- c(-0.8, -1) #int.dC

plt$graphAttributes$Nodes$labels <- c("Perive-\nticular\nWMH\nt0", "WM\nFree\nWater\nt0", "Basal\nganglia\nPVS\nt0",
                                      "Perive-\nticular\nWMH\nt1", "WM\nFree\nWater\nt1", "Basal\nganglia\nPVS\nt1",
                                      "\u0394SBP", "Intensive\nBP\nControl",
                                      "Base-\nline\nage", "Follow\nup\ntime",
                                      "SVD\nt0", "SVD\nt1", "SVD\nChange",
                                      "Mean\nBaseline\nSVD", "Mean\nSVD\nChange")

plt <- semptools::mark_sig(semPaths_plot = plt, object = fit6)

plt$graphAttributes$Nodes$label.cex <- 1.6

plt$graphAttributes$Edges$curve[c(26,27)] <- 0.5

plt$graphAttributes$Edges$label.margin[c(16:18)] <- -0.035
plt$graphAttributes$Edges$label.margin[c(21)] <- -0.03
plt$graphAttributes$Edges$label.margin[c(23)] <- -0.02
plt$graphAttributes$Edges$label.margin[c(24)] <- -0.02

plt$plotOptions$label.prop <- c(1, 1, 1,
                                1, 1, 1,
                                0.9, 0.95,
                                0.85, 0.85,
                                0.7, 0.7, 0.85,
                                0.45, 0.425)

plt$graphAttributes$Edges$color[c(1,2,4,5,7,8)] <- "skyblue"

plt$graphAttributes$Edges$color[c(3,6,9)] <- "navy"

plt$graphAttributes$Edges$color[c(10:15)] <- "coral"

plt$graphAttributes$Edges$color[c(16:18)] <- "lightgreen"

plt$graphAttributes$Edges$color[c(21,23)] <- "blue"

plt$graphAttributes$Edges$color[c(22,24)] <- "red2"

plt$graphAttributes$Nodes$width <- c(5, 5, 5,
                                     5, 5, 5,
                                     5, 5,
                                     5, 5,
                                     8, 8, 8,
                                     12, 12)

plot(plt)

# Bootstrapped solution
# Detect cpu number
ncpus = max(1, parallel::detectCores() - 1)

# Fit LCS model with Maximum Likelihood estimation method using bootstrapped standard errors and a scaled (Yuan-Bentler) test statistic
fit6_boot <- lavaan::sem(model6, data = wide01, estimator = "ML", test = "yuan.bentler.mplus",
                         se = "bootstrap", bootstrap = 5000, parallel = "snow",
                         ncpus = ncpus, iseed = 2025)

# Extract parameter estimates of interest from bootstrapped solution
sum <- parameterEstimates(fit6_boot)
sum <- sum[sum$op == ":=" , c("label", "est", "ci.lower", "ci.upper", "pvalue")]
rownames(sum) <- sum[,"label"]
sum <- sum[ ,-1]
sum


#//--------------------------------------------------------------------------- END OF ANALYTIC CODE ---------------------------------------------------------------------------//



