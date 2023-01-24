library(tidyverse)
library(mgcv)
library(lme4)
library(gratia)

# # # LOAD XY data 
# (Processed, Curated, and Combined in Python)
region = 'Mpala'
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
XYdata = XYdata %>% mutate(PAIherb_5cm = log(1 - coverherb_5cm)/0.5,
                           PAIherb_0m = log(1 - coverherb_0m)/0.5,
                           PAIherb_5cm_norm = (log(1 - coverherb_5cm)/0.5)/herbh,
                           PAIherb_0m_norm = (log(1 - coverherb_0m)/0.5)/herbh)

# Add a log mean height
XYdata = XYdata %>% mutate(logmean = if_else(mean>0, log(mean), mean))

# Remove na rows where biomass is na
XYdata = XYdata %>% drop_na(Dry_Weight_g)

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

# Non-Logged DryWeight
# Mpala - top pairwise corr - 1/23/23
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
# coverherb_0m          0.679265
# meanpeakh             0.666096
# maxpeakh              0.654558
mpala_m1 = gam(Dry_Weight_g ~ s(logmean) +
                         s(PAImean_5cmto1p5m),
                       data=XY_scale,
                       family=gaussian,
                       select=TRUE,
                       method="REML")


draw(mpala_m1, residuals = TRUE)
summary(mpala_m1)
appraise(mpala_m1)
gam.check(mpala_m1, rep=500)
lines(c(0, 370), c(0, 370))


#pred = predict(mpala_m1, data=XY_scale)
RMSE = sqrt(mean((mpala_m1$residuals)^2))
rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
RMSE
rRMSE

# Above model is prob the best so far
# ----

mpala_m2 = gam(Dry_Weight_g ~ s(logmean) +
                 s(cover_5cm),
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")


draw(mpala_m2, residuals = TRUE)
summary(mpala_m2)
appraise(mpala_m2)
gam.check(mpala_m2, rep=500)
lines(c(0, 370), c(0, 370))


#pred = predict(mpala_m1, data=XY_scale)
RMSE = sqrt(mean((mpala_m2$residuals)^2))
rRMSE = RMSE/mean(XY_scale$Dry_Weight_g)
RMSE
rRMSE


mpala_m3 = gam(Dry_Weight_g ~ s(perc100) +
                 s(cover_5cm),
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")


draw(mpala_m3, residuals = TRUE)
summary(mpala_m3)
appraise(mpala_m3)
gam.check(mpala_m3, rep=500)
lines(c(0, 370), c(0, 370))

## The 2 models above also are not bad


# Looks like an issue with the smaller grass biomass values
# Which has to do with the mean height saturdat
# Let's try looking at some of our variables 
p1 = ggplot(XYdata, aes(x = PAIsum, y = Dry_Weight_g)) +
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


p4 = ggplot(XYdata, aes(x = PAIsum, y = logDryWeight)) +
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


