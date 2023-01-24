library(tidyverse)
library(mgcv)
library(lme4)
library(gratia)

# # # LOAD XY data 
# (Processed, Curated, and Combined in Python)
region = 'Kruger'
datadir = '/n/home02/pbb/scripts/halo-metadata-server/GrassBiomass/data'
XYdata = read_csv(file=paste0(datadir, '/out/', region,'/XYdata.csv'))
colnames(XYdata)


# # # Clean & Scale data, and break into training/testing
# how to replace column names all at once: 
# https://dnidzgorski.wordpress.com/2017/06/09/r-fix-column-names-spaces/
# Take out spaces in column names
colnames(XYdata)<-str_replace_all(colnames(XYdata), c(" " = "_", "#" = "number"))
colnames(XYdata)<-str_replace_all(colnames(XYdata), c("25" = "perc25"))
colnames(XYdata)<-str_replace_all(colnames(XYdata), c("50" = "perc50"))
colnames(XYdata)<-str_replace_all(colnames(XYdata), c("75" = "perc75"))
colnames(XYdata)<-str_replace_all(colnames(XYdata), c("98" = "perc98"))
colnames(XYdata)<-str_replace_all(colnames(XYdata), c("100" = "perc100"))

# Note, have to use \\ with this one, otherwise it thinks "(" is a regex
colnames(XYdata)<-str_replace_all(colnames(XYdata), c("\\(" = "", "\\)" = ""))
colnames(XYdata)


# Add variables, PAIherb
# Note: You have to do it as > 0 and < 1 to avoid NaNs
XYdata = XYdata %>% mutate(PAIherb_5cm = if_else((coverherb_5cm<1) & (coverherb_5cm>0),
                                                 log(1 - coverherb_5cm)/0.5,
                                                 NaN),
                           PAIherb_0m = if_else((coverherb_0m<1) & (coverherb_0m>0),
                                                log(1 - coverherb_0m)/0.5,
                                                NaN),
                           PAIherb_5cm_norm = if_else((coverherb_5cm<1) & (coverherb_5cm>0),
                                                      (log(1 - coverherb_5cm)/0.5)/herbh,
                                                      NaN),
                           PAIherb_0m_norm = if_else((coverherb_0m<1) & (coverherb_0m>0),
                                                     (log(1 - coverherb_0m)/0.5)/herbh,
                                                     NaN))


# Add a log mean height
XYdata = XYdata %>% mutate(logmean = if_else(mean>0, log(mean), mean))

# Break into X and Y for scaling
X_lidar = subset(XYdata, select = c(coverherb_0m, coverherb_5cm,
                                    cover_0m, cover_5cm,
                                    perc25, perc50, perc75,
                                    perc98, perc100, mean, std,
                                    nlayers, gapsize, maxpeakh,
                                    ptoh, cscore, FHD, VDR,
                                    meanpeakh, stdpeakh, cvpeakh,
                                    herbh, PAIherb_5cm, PAIherb_0m,
                                    PAIherb_5cm_norm, PAIherb_0m_norm, logmean,
                                    PAIsum_0m, PAIsum_5cm, PAIsum_0mto1p5m,
                                    PAIsum_5cmto1p5m, PAImean_0mto1p5m, PAImean_5cmto1p5m))

# Scale numeric X vars, and combine back with Y vars
XY_scale = X_lidar %>%
  scale() %>%
  data.frame() %>%
  bind_cols(select(XYdata, c('id', 'Dry_Weight_g', 'logDryWeight')))

colnames(XY_scale)


# Split into a training and testing dataset

#make this example reproducible
set.seed(42)

# Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE),
                 nrow(XY_scale),
                 replace=TRUE,
                 prob=c(0.7,0.3))

# From: https://www.statology.org/train-test-split-r/
train  <- XY_scale[sample, ]
test   <- XY_scale[!sample, ]

nrow(train)
nrow(test)



# # # Fit GAMs

# KRUGER - top pairwise corr - 1/23
# non-logged dry weight (g)
# cover_5cm             0.865090
# herbh                 0.864227
# coverherb_5cm         0.863829
# PAIsum_5cm           0.818907
# 50                    0.812874
# 25                    0.811870
# mean                  0.799406
# PAImean_5cmto1p5m    0.797188
# PAIsum_5cmto1p5m     0.797188
# 75                    0.791930
# cover_0m              0.739724
# 98                    0.738793
# coverherb_0m          0.679265
kruger_m1 = gam(Dry_Weight_g ~ s(mean) +
                 s(PAIsum_5cm),
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")

draw(kruger_m1, residuals = TRUE)
summary(kruger_m1)
appraise(kruger_m1)
gam.check(kruger_m1, rep=500)
lines(c(0, 370), c(0, 370))


pred = predict(kruger_m1, data=XY_scale)
RMSE = sqrt(mean((pred - XY_scale$Dry_Weight_g)^2))
rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
RMSE
rRMSE

kruger_m2 = gam(Dry_Weight_g ~ s(herbh) +
                  s(cover_5cm),
                data=XY_scale,
                family=gaussian,
                select=TRUE,
                method="REML")

draw(kruger_m2, residuals = TRUE)
summary(kruger_m2)
appraise(kruger_m2)
gam.check(kruger_m2, rep=500)
lines(c(0, 370), c(0, 370))


kruger_m2 = gam(Dry_Weight_g ~ s(herbh) +
                  s(PAIsum_5cm),
                data=XY_scale,
                family=gaussian,
                select=TRUE,
                method="REML")

draw(kruger_m2, residuals = TRUE)
summary(kruger_m2)
appraise(kruger_m2)
gam.check(kruger_m2, rep=500)
lines(c(0, 370), c(0, 370))



kruger_m3 = gam(Dry_Weight_g ~ s(herbh) +
                  s(PAIherb_5cm),
                data=XY_scale,
                family=gaussian,
                select=TRUE,
                method="REML")

draw(kruger_m3, residuals = TRUE)
summary(kruger_m3)
appraise(kruger_m3)
gam.check(kruger_m3, rep=500)
lines(c(0, 370), c(0, 370))


kruger_m4 = gam(Dry_Weight_g ~ s(herbh) +
                  s(PAImean_5cmto1p5m),
                data=XY_scale,
                family=gaussian,
                select=TRUE,
                method="REML")

draw(kruger_m4, residuals = TRUE)
summary(kruger_m4)
appraise(kruger_m4)
gam.check(kruger_m4, rep=500)
lines(c(0, 370), c(0, 370))

pred = predict(kruger_m4, data=XY_scale)
RMSE = sqrt(mean((pred - XY_scale$Dry_Weight_g)^2))
rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
RMSE
rRMSE
