################################################################################
# The one-stage design for countrywide studies on temperature and health.
# Aurelio Tobias (aurelio.tobias@idaea.csic.es)
# Carmen Iniguez (Carmen.Iniguez@uv.es )
# version 2025.07.31
################################################################################
# Comparison of one-stage and two-stage designs fitting temperature in 
# absolute scale (°C)
################################################################################

# Remove objects 
rm(list=ls())

# Load packages and functions
library(dplyr)
library(dlnm)
library(splines)
library(tsModel)
library(mixmeta)
library(ggplot2)
library(RColorBrewer)

source("findmin.R")  

# Import openaccess England and Wales dataset from Antonio Gasparrini GitHub
url <- "https://raw.githubusercontent.com/gasparrini/2012_gasparrini_StatMed_Rcodedata/master/regEngWales.csv"
data <- read.csv(url)
head(data)

################################################################################
# One-Stage design
################################################################################

# Aggregate data across all regions by date
data_agg <- data %>%
  group_by(date) %>%
  summarise(
    death = sum(death, na.rm = TRUE),
    tmean = mean(tmean, na.rm = TRUE),
    dow = first(dow),
    time = first(time),
    year = first(year)
  ) %>%
  ungroup()

# Define cross-basis for temperature
nlag <- 21
klag <- logknots(nlag, 3)
varknots_agg <- quantile(data_agg$tmean, c(0.1, 0.75, 0.9), na.rm = TRUE)
cb <- crossbasis(data_agg$tmean, lag = nlag,
                 argvar = list(fun = "ns", knots = varknots_agg),
                 arglag = list(fun = "ns", knots = klag))
# Fit model
df <- 10
ny <- length(unique(data_agg$year))
mod1 <- glm(death ~ cb + ns(time, df*ny) + factor(dow), family = quasipoisson(), data = data_agg)

# Predict using temperature in absolute scale (oC)  
temp.seq <- seq(min(data_agg$tmean, na.rm = TRUE), max(data_agg$tmean, na.rm = TRUE), by = 0.1)
mmt1 <- findmin(cb, mod1)
pred1 <- crosspred(cb, mod1, cen = mmt1, by = 0.1, at = temp.seq)

################################################################################
# Two-Stage design
################################################################################

regions <- unique(data$region)
coef <- matrix(NA, nrow = length(regions), ncol = 4) # ncol depens on the num. coefs  
vcov <- vector("list", length(regions))
names(coef) <- regions
names(vcov) <- regions
mmt <- numeric(length(regions))

# Loop for city-specific analysis
for (i in seq_along(regions)) {
  c <- regions[i]
  data_c <- subset(data, region == c)
  varknots_c <- quantile(data_c$tmean, c(0.1, 0.75, 0.9), na.rm = TRUE)
  cb <- crossbasis(data_c$tmean, lag = nlag,
                   argvar = list(fun = "ns", knots = varknots_c),
                   arglag = list(fun = "ns", knots = klag))
  ny_c <- length(unique(data_c$year))
  mod <- glm(death ~ cb + ns(time, df*ny) + factor(dow), family = quasipoisson(), data = data_c)
  mmt[i] <- findmin(cb, mod)
  pred <- crossreduce(cb, mod, cen = mmt[i], by = 0.1)
  coef[i, ] <- pred$coef
  vcov[[i]] <- pred$vcov
}

# Mixed effects meta-analysis
mv <- mixmeta(coef ~ 1, vcov, method = "reml")

# Predict using temperature in absolute scale (oC)
temp.seq   <- seq(min(data$tmean, na.rm = TRUE), max(data$tmean, na.rm = TRUE), by = 0.1)
basis.pred <- onebasis(temp.seq, fun = "ns", knots = quantile(data$tmean, c(0.10, 0.75, 0.90), na.rm = TRUE))
mmt2 <- findmin(basis.pred, coef = coef(mv), vcov = vcov(mv))
pred2 <- crosspred(basis.pred, coef = coef(mv), vcov = vcov(mv), model.link = "log", by = 0.1, cen = mmt2, at = temp.seq)

################################################################################
# Plot
################################################################################

# Colors and theme
colors2 <- brewer.pal(n = 10, name = "RdBu")
colors2 <- colors2[ c(2,9, 3,8)]

theme_cc <- function(){  
  theme_classic() %+replace%  
    theme(
      plot.title = element_text(hjust = 0.5, colour = "black", size = 12, margin = margin(b = 20)),
      panel.border = element_rect(colour= "black", size= 0.5,fill = NA, linetype=1),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = 8),
    )  
}

# Combined data.frame with with predicted values to plot
dfpred <- rbind(
  data.frame(x = pred1$predvar, y = pred1$allRRfit, lo = pred1$allRRlow, hi = pred1$allRRhigh, cen = mmt1, 
             type = "One-stage", coun = "England & Wales"),
  data.frame(x = pred2$predvar, y = pred2$allRRfit, lo = pred2$allRRlow, hi = pred2$allRRhigh, cen = mmt2, 
             type = "Two-stage" , coun = "England & Wales") )

# Overall plot for comparison
ggplot(dfpred, aes(x, y, color=type)) + 
  geom_hline(yintercept = 1, size = 0.25) +
  geom_vline(aes(xintercept = cen, color = type), linetype = "dashed", size = 0.25, show.legend = FALSE) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = .2, size=0.05, show.legend = FALSE) +
  facet_wrap(vars(coun))+
  scale_x_continuous(breaks = seq(-5, 30, by = 5)) +
  scale_y_continuous(breaks = seq(1, 1.4, by = 0.1), labels = seq(1, 1.4, by = 0.1)) +
  scale_color_manual(values = colors2[c(3, 4)], labels = c("One-stage", "Two-stage")) +
  scale_fill_manual(values = colors2[c(3, 4)], labels = c("One-stage", "Two-stage")) +
  coord_cartesian(xlim = c(-5, 30), ylim = c(0.99, 1.4))+
  labs(x = "Temperature, ºC", y = "Relative Risk", color = NULL, fill = NULL) +
  theme_cc()

################################################################################
# Plot in percentile scale 
################################################################################

# One-stage predictions temperature in percentile scale (%)
pct.seq <- seq(1, 99, by = 1)
temp.at <- quantile(data_agg$tmean, probs = pct.seq / 100, na.rm = TRUE)
mmt1 <- findmin(cb, mod1)
pred1 <- crosspred(cb, mod1, cen = mmt1, by = 1, at = temp.at)

# Two-stage predictions temperature in percentile scale (%)
basis.pred <- onebasis(temp.at, fun = "ns", 
                       knots = quantile(data$tmean, c(0.10, 0.75, 0.90), na.rm = TRUE),
                       Boundary.knots = range(data$tmean))
mmt2 <- findmin(basis.pred, coef = coef(mv), vcov = vcov(mv))
pred2 <- crosspred(basis.pred, coef = coef(mv), vcov = vcov(mv),
                   model.link = "log", cen = mmt2, at = temp.at)

# Combined data.frame with with predicted values to plot
dfpred <- rbind(
  data.frame(x = pct.seq, y = pred1$allRRfit, lo = pred1$allRRlow, hi = pred1$allRRhigh, 
             cen = ecdf(data_agg$tmean)(mmt1)*100,
             type = "One-stage", coun = "England & Wales"),
  data.frame(x = pct.seq, y = pred2$allRRfit, lo = pred2$allRRlow, hi = pred2$allRRhigh,
             cen = ecdf(data$tmean)(mmt2)*100,
             type = "Two-stage", coun = "England & Wales")
)

# Overall plot for comparison labeling percentiles (%)
pct.seq <- seq(0, 100, by = 10)
temp.at <- quantile(data$tmean, probs = pct.seq / 100, na.rm = TRUE)
x_labels <- paste0(pct.seq, "\n(", round(temp.at, 1), ")")

ggplot(dfpred, aes(x, y, color = type)) + 
  geom_hline(yintercept = 1, size = 0.25) +
  geom_vline(aes(xintercept = cen, color = type), linetype = "dashed", size = 0.25, show.legend = FALSE) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = .2, size = 0.05, show.legend = FALSE) +
  facet_wrap(vars(coun))+
  scale_x_continuous(breaks = pct.seq, labels = x_labels + limits = c(-1, 101)) +
  scale_y_continuous(breaks = seq(1, 1.4, by = 0.1), labels = seq(1, 1.4, by = 0.1) + limits=c(0.99, 1.4)) +
  scale_color_manual(values = colors2[c(3, 4)], labels = c("One-stage", "Two-stage")) +
  scale_fill_manual(values = colors2[c(3, 4)], labels = c("One-stage", "Two-stage")) +
  labs(x = "Percentile, % (Temperature, °C)", y = "Relative Risk", color = NULL, fill = NULL) +
  theme_cc()

################################################################################
# End of script file
################################################################################
