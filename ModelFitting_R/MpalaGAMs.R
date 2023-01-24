library(tidyverse)
library(mgcv)
library(lme4)
library(gratia)

# # # LOAD XY data 
# (Processed, Curated, and Combined in Python)

# Set region (kruger or mpala)
region = 'Mpala'

# Set dataset label (initial, spicy)
label = 'spicy'

datadir = '/n/home02/pbb/scripts/halo-metadata-server/GrassBiomass/data'
XYdata = read_csv(file=paste0(datadir, '/out/',
                              label, '/',
                              region,'/XYdata.csv'))
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

# # Add variables, PAIherb
# XYdata = XYdata %>% mutate(PAIherb_5cm = -log(1 - coverherb_5cm)/0.5,
#                            PAIherb_0m = -log(1 - coverherb_0m)/0.5,
#                            PAIherb_5cm_norm = (-log(1 - coverherb_5cm)/0.5)/herbh,
#                            PAIherb_0m_norm = (-log(1 - coverherb_0m)/0.5)/herbh)

# Add variables, PAIherb
# Note: You have to do it as > 0 and < 1 to avoid NaNs
XYdata = XYdata %>% mutate(PAIherb_5cm = if_else((coverherb_5cm<1) & (coverherb_5cm>0),
                                                 -log(1 - coverherb_5cm)/0.5,
                                                 NaN),
                           PAIherb_0m = if_else((coverherb_0m<1) & (coverherb_0m>0),
                                                -log(1 - coverherb_0m)/0.5,
                                                NaN),
                           PAIherb_5cm_norm = if_else((coverherb_5cm<1) & (coverherb_5cm>0),
                                                      (-log(1 - coverherb_5cm)/0.5)/herbh,
                                                      NaN),
                           PAIherb_0m_norm = if_else((coverherb_0m<1) & (coverherb_0m>0),
                                                     (-log(1 - coverherb_0m)/0.5)/herbh,
                                                     NaN))

# Add a log mean height
XYdata = XYdata %>% mutate(logmean = if_else(mean>0, log(mean), mean))

# Add Coefficient of variation
XYdata = XYdata %>% mutate(CV = mean/std,
                           prat = perc100 - mean / perc100,
                           prat50 = perc100 - perc50 / perc100)

# Remove na rows where biomass is na
XYdata = XYdata %>% drop_na(Dry_Weight_g)


#Add a coverXheight
XYdata = XYdata %>% mutate(CCxH = cover_0m*mean,
                           CCxH100 = cover_0m*perc100)

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
                                    PAIsum_5cmto1p5m, PAImean_0mto1p5m, PAImean_5cmto1p5m,
                                    prat, prat50, CV, CCxH, CCxH100))

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

# # # SPICY
# Mpala - Vars with pairwise correlation with Non-Logged DryWeight
# For SPICY dataset 1/24
# 100                           0.804155
# 98                            0.784037
# 75                            0.760976
# std                           0.750670
# mean                          0.744222
# coverherb_5cm                 0.729503
# cover_5cm                     0.726991
# 50                            0.715745
# PAIsum_5cm                    0.701315
# PAImean_5cmto1p5m             0.701315
# PAIsum_5cmto1p5m              0.701315
# meanpeakh                     0.701067

mpala_spicy1 = gam(Dry_Weight_g ~ s(perc100) +
                         s(coverherb_5cm),
                       data=XY_scale,
                       family=gaussian,
                       select=TRUE,
                       method="REML")


draw(mpala_spicy1, residuals = TRUE)
summary(mpala_spicy1)
appraise(mpala_spicy1)
gam.check(mpala_spicy1, rep=500)
lines(c(0, 370), c(0, 370))

pred = predict(mpala_spicy1, data=XY_scale)
RMSE = sqrt(mean((mpala_spicy1$residuals)^2))
rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
# sqrt(mean((mpala_spicy1$residuals)^2))
RMSE
rRMSE
# Nice - looks like coverherb_5cm didn't do as well though
# perc100 underneath herbh works v well though

# Try perc100 and std of points
mpala_spicy2 = gam(Dry_Weight_g ~ s(perc100) +
                     s(meanpeakh),
                   data=XY_scale,
                   family=gaussian,
                   select=TRUE,
                   method="REML")


draw(mpala_spicy2, residuals = TRUE)
summary(mpala_spicy2)
appraise(mpala_spicy2)
gam.check(mpala_spicy2, rep=500)
lines(c(0, 370), c(0, 370))

pred = predict(mpala_spicy2, data=XY_scale)
RMSE = sqrt(mean((mpala_spicy2$residuals)^2))
rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
# sqrt(mean((mpala_spicy1$residuals)^2))
RMSE
rRMSE
# Slightly better, but definitely weird with those large std and mean values


# Try Lots of Vars
# Complex - but decent model
mpala_spicy3 = gam(Dry_Weight_g ~ s(perc100) +
                     s(prat) + s(std),
                   data=XY_scale,
                   family=gaussian,
                   select=TRUE,
                   method="REML")


draw(mpala_spicy3, residuals = TRUE)
summary(mpala_spicy3)
appraise(mpala_spicy3)
gam.check(mpala_spicy3, rep=500)
lines(c(0, 370), c(0, 370))

pred = predict(mpala_spicy2, data=XY_scale)
RMSE = sqrt(mean((mpala_spicy2$residuals)^2))
rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
# sqrt(mean((mpala_spicy1$residuals)^2))
RMSE
rRMSE



# Try CCxH
mpala_spicy4 = gam(Dry_Weight_g ~ s(CCxH100),
                   data=XY_scale,
                   family=gaussian,
                   select=TRUE,
                   method="REML")


draw(mpala_spicy4, residuals = TRUE)
summary(mpala_spicy4)
appraise(mpala_spicy4)
gam.check(mpala_spicy4, rep=500)
lines(c(0, 370), c(0, 370))

pred = predict(mpala_spicy4, data=XY_scale)
RMSE = sqrt(mean((mpala_spicy4$residuals)^2))
rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
# sqrt(mean((mpala_spicy1$residuals)^2))
RMSE
rRMSE


p1 = ggplot(XYdata, aes(x = perc75,
                        y = Dry_Weight_g,
                        colour=Tree_Species)) +
  geom_point(size=2.5) +
  xlab("75th percentile height [m]") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p1

p1 = ggplot(XYdata, aes(x = CCxH,
                        y = Dry_Weight_g,
                        colour=Tree_Species)) +
  geom_point(size=2.5) +
  xlab("Cover x Mean Height") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p1

p1 = ggplot(XYdata, aes(x = mean,
                        y = Dry_Weight_g,
                        colour=Tree_Species)) +
  geom_point(size=2.5) +
  xlab("Mean Height") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p1

p1 = ggplot(XYdata, aes(x =cover_5cm,
                        y = Dry_Weight_g,
                        colour=Tree_Species)) +
  geom_point(size=2.5) +
  xlab("Canopy Density") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p1


# # # INITIAL Dataset 1/23
# Looks like the best (for non-logged dry-weight data) are:
# cover_5cm             0.865090
# herbh                 0.864227
# coverherb_5cm         0.863829
# PAIsum_25cm           0.818907
# 50                    0.812874
# 25                    0.811870
# mean                  0.799406
# PAImean_25cmto1p5m    0.797188
# PAIsum_25cmto1p5m     0.797188
# 75                    0.791930
# cover_0m              0.739724
# 98                    0.738793

# initial dataset best model
# mpala_m1 = gam(Dry_Weight_g ~ s(logmean) +
#                          s(PAImean_5cmto1p5m),
#                        data=XY_scale,
#                        family=gaussian,
#                        select=TRUE,
#                        method="REML")
# 
# 
# draw(mpala_m1, residuals = TRUE)
# summary(mpala_m1)
# appraise(mpala_m1)
# gam.check(mpala_m1, rep=500)
# lines(c(0, 370), c(0, 370))

#pred = predict(mpala_m1, data=XY_scale)
# RMSE = sqrt(mean((mpala_m1$residuals)^2))
# rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
# RMSE
# rRMSE

# ----
#  # Another initial dataset good model 
# mpala_m2 = gam(Dry_Weight_g ~ s(herbh) +
#                  s(cover_5cm),
#                data=XY_scale,
#                family=gaussian,
#                select=TRUE,
#                method="REML")
# 
# 
# draw(mpala_m2, residuals = TRUE)
# summary(mpala_m2)
# appraise(mpala_m2)
# gam.check(mpala_m2, rep=500)
# lines(c(0, 370), c(0, 370))
# 
# 
# #pred = predict(mpala_m1, data=XY_scale)
# RMSE = sqrt(mean((mpala_m2$residuals)^2))
# rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
# RMSE
# rRMSE
# 
# 
# mpala_m3 = gam(Dry_Weight_g ~ s(perc100) +
#                  s(cover_5cm),
#                data=XY_scale,
#                family=gaussian,
#                select=TRUE,
#                method="REML")
# 
# 
# draw(mpala_m3, residuals = TRUE)
# summary(mpala_m3)
# appraise(mpala_m3)
# gam.check(mpala_m3, rep=500)
# lines(c(0, 370), c(0, 370))

## The 2 models above also are not bad 1/23
# Looks like an issue with the smaller grass biomass values
# Which has to do with the mean height saturdat
# Let's try looking at some of our variables 
# ------ 


p1 = ggplot(XYdata, aes(x = PAImean_5cmto1p5m, y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("PAI") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p1


p2 = ggplot(XYdata, aes(x = perc100, y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("100th percentile height") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p2

p3 = ggplot(XYdata, aes(x = perc50, y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("50th percentile height") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p3

p3 = ggplot(XYdata, aes(x = perc98, y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("75th percentile height") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p3

p3 = ggplot(XYdata, aes(x = mean, y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Mean height") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p3

p3 = ggplot(XYdata, aes(x = logmean, y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Log Mean height") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p3


p4 = ggplot(XYdata, aes(x = logmean, y = logDryWeight)) +
  geom_point(size=2.5) +
  xlab("PAI") +
  ylab("log(Biomass (g))") + scale_color_hue(direction = -1)

p4


p5 = ggplot(XYdata, aes(x = cover_0m,
                        y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Canopy Density 0 m") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p5


p6 = ggplot(XYdata, aes(x = cover_5cm,
                        y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Canopy Density 5 cm") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p6


p7 = ggplot(XYdata, aes(x = coverherb_5cm,
                        y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Canopy Density Herbaceous 5 cm") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p7

p8 = ggplot(XYdata, aes(x = coverherb_0m,
                        y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Canopy Density Herbaceous 0 m") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p8

p9 = ggplot(XYdata, aes(x = PAIherb_0m,
                        y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("PAI Herbaceous 0 m") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p9

p10 = ggplot(XYdata, aes(x = PAIherb_5cm_norm,
                        y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("PAI Herbaceous 5 cm Normalized") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p10

p11 = ggplot(XYdata, aes(x = PAIherb_0m_norm,
                        y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("PAI Herbaceous 0 m Normalized") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p11

p12 = ggplot(XYdata, aes(x = PAIherb_5cm,
                        y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("PAI Herbaceous 5 cm") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p12

p12 = ggplot(XYdata, aes(x = herbh,
                         y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Herbaceous Height [m]") +
  ylab("Biomass (g)") + scale_color_hue(direction = -1)

p12

p13 = ggplot(XYdata, aes(x = herbh,
                         y = Average_Height)) +
  geom_point(size=2.5) +
  xlab("Herbaceous Height [m]") +
  ylab("Mean Field Height") + scale_color_hue(direction = -1)

p13

p13 = ggplot(XYdata, aes(x = herbh,
                         y = Max_Height_m)) +
  geom_point(size=2.5) +
  xlab("Herbaceous Height [m]") +
  ylab("Max Field Height") + scale_color_hue(direction = -1)

p13

p14 = ggplot(XYdata, aes(x = mean,
                         y = Average_Height)) +
  geom_point(size=2.5) +
  xlab("Mean Lidar Height") +
  ylab("Mean Field Height") + scale_color_hue(direction = -1)

p14

p14 = ggplot(XYdata, aes(x = perc98,
                         y = Max_Height_m)) +
  geom_point(size=2.5) +
  xlab("Mean Lidar Height") +
  ylab("Mean Field Height") + scale_color_hue(direction = -1)

p14

p14 = ggplot(XYdata, aes(x = perc98,
                         y = Max_Height_m)) +
  geom_point(size=2.5) +
  xlab("Mean Lidar Height") +
  ylab("Mean Field Height") + scale_color_hue(direction = -1)

p14

p14 = ggplot(XYdata, aes(x = prat,
                         y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Penetration Ratio") +
  ylab("Biomass [g]") + scale_color_hue(direction = -1)

p14


p14 = ggplot(XYdata, aes(x = CV,
                         y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Coefficient of Variation") +
  ylab("Biomass [g]") + scale_color_hue(direction = -1)

p14

p14 = ggplot(XYdata, aes(x = prat50,
                         y = Dry_Weight_g)) +
  geom_point(size=2.5) +
  xlab("Penetration Ratio (using 50th Percentile)") +
  ylab("Biomass [g]") + scale_color_hue(direction = -1)

p14

# # # With 100 height instead
mpala_m2 = gam(Dry_Weight_g ~ s(PAIsum),
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")


draw(mpala_m2, residuals = TRUE)
summary(mpala_m2)
appraise(mpala_m2)
gam.check(mpala_m2, rep=500)
lines(c(0, 370), c(0, 370))


mpala_m3 = gam(Dry_Weight_g ~ s(cover_0m),
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")


draw(mpala_m3, residuals = TRUE)
summary(mpala_m3)
appraise(mpala_m3)
gam.check(mpala_m3, rep=500)
lines(c(0, 370), c(0, 370))


mpala_m4 = gam(Dry_Weight_g ~ s(cover_5cm),
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")


draw(mpala_m4, residuals = TRUE)
summary(mpala_m4)
appraise(mpala_m4)
gam.check(mpala_m4, rep=500)
lines(c(0, 370), c(0, 370))

mpala_m4 = gam(Dry_Weight_g ~ s(coverherb_5cm),
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")


draw(mpala_m4, residuals = TRUE)
summary(mpala_m4)
appraise(mpala_m4)
gam.check(mpala_m4, rep=500)
lines(c(0, 370), c(0, 370))


### Investigate duplicate rows- why you end up with all these 50 g values

XY_dupes = XYdata[duplicated(mpala_m1$fitted.values), ]
XY_nodupes = XYdata[!duplicated(mpala_m1$fitted.values), ]

XY_dupes$mean
XY_dupes$perc50
XY_dupes$perc100
XY_dupes$coverherb_0m
XY_dupes$coverherb_5cm

# Ok, so they're plots that have all their cover below 5 cm.

boxplot(XY_dupes$Average_Height)
boxplot(XY_dupes$Max_Height_cm)

boxplot(XY_nodupes$Max_Height_cm)
boxplot(XY_nodupes$Average_Height)

XYdata = XYdata %>% mutate(Average_Height_m = Average_Height/100) 

XYdata_nobareground = XYdata %>% subset(Average_Height>0)

XYdata_lidarbareground = XYdata %>% subset(mean==0)

XYdata_lidaraboveground = XYdata %>% subset(mean>0)

XYdata_bareground = XYdata %>% subset(Average_Height==0)

p14 = ggplot(XYdata_lidaraboveground, aes(x = mean,
                         y = Average_Height_m)) +
  geom_point(size=2.5) +
  xlab("Mean Lidar Height") +
  ylab("Mean Field Height") + scale_color_hue(direction = -1)

p14

p15 = ggplot(XYdata_lidarbareground, aes(x = mean,
                                         y = Average_Height_m)) +
  geom_boxplot(size=2.5) +
  xlab("Mean Lidar Height") +
  ylab("Mean Field Height") + scale_color_hue(direction = -1)

p15

p16 = ggplot(data=XYdata_lidarbareground) +
  geom_point(aes(x = herbh,
                 y = Average_Height_m),
             size=2.5) +
  xlab("Herb Lidar Height") +
  ylab("Mean Field Height") + 
  scale_color_hue(direction = -1)

p16

p16 = ggplot(data=XYdata_lidaraboveground) +
  geom_point(aes(x = herbh,
                 y = Average_Height_m),
             size=2.5) +
  xlab("Herb Lidar Height") +
  ylab("Mean Field Height") + 
  scale_color_hue(direction = -1)

p16

p16 = ggplot(data=XYdata_lidaraboveground) +
  geom_point(aes(x = perc50,
                 y = Average_Height_m),
             size=2.5) +
  xlab("Herb Lidar Height") +
  ylab("Mean Field Height") + 
  scale_color_hue(direction = -1)

p16

