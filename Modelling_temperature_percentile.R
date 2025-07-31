################################################################################
# The one-stage design for countrywide studies on temperature and health.
# Aurelio Tobias (aurelio.tobias@idaea.csic.es)
# Carmen Iniguez (Carmen.Iniguez@uv.es )
# version 2025.07.31
################################################################################
# Comparison of one-stage and two-stage designs fitting temperature in 
# percentile scale (%)
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

# Import open access England and Wales dataset from Antonio Gasparrini GitHub
url <- "https://raw.githubusercontent.com/gasparrini/2012_gasparrini_StatMed_Rcodedata/master/regEngWales.csv"
data <- read.csv(url)
head(data)

################################################################################
# Transform temperature to percentiles
################################################################################

data$pct <- ave(data$tmean, data$region, FUN = function(x) ecdf(x)(x)*100)
head(data)

################################################################################
# One-Stage Design
################################################################################

# Aggregate data across all regions by date
data_agg <- data %>%
  group_by(date) %>%
  summarise(death = sum(death, na.rm = TRUE),
            pct = mean(pct, na.rm = TRUE), # Temperature to percentiles
            dow = first(dow),
            time = first(time),
            year = first(year)) %>%
  ungroup()

# Define cross-basis for temperature
nlag <- 21
klag <- logknots(nlag, 3)
varknots <- c(10, 75, 90)
cb1 <- crossbasis(data_agg$pct, lag = nlag,
                  argvar = list(fun = "ns", knots = varknots),
                  arglag = list(fun = "ns", knots = klag))

# Fit model
df <- 10
ny <- length(unique(data_agg$year))
mod1 <- glm(death ~ cb1 + ns(time, df*ny) + factor(dow), family = quasipoisson(), data = data_agg)

# Predict using temperature in in percentile scale (%)
pct.seq <- seq(0, 100, by = 1)
mmt1 <- findmin(cb1, mod1)
pred1 <- crosspred(cb1, mod1, cen = mmt1, by = 1, at = pct.seq)

################################################################################
# Two-Stage Design
################################################################################

regions <- unique(data$region)
coef <- matrix(NA, nrow = length(regions), ncol = 4) # ncol depens on the num. coefs  
vcov <- vector("list", length(regions))
mmt <- numeric(length(regions))
names(coef) <- names(vcov) <- regions

# Loop for city-specific analysis
for (i in seq_along(regions)) {
  c <- regions[i]
  data_c <- subset(data, region == c)
  data_c$pct <- ecdf(data_c$tmean)(data_c$tmean)*100 # Temperature to percentiles
  cb <- crossbasis(data_c$pct, lag = nlag,
                   argvar = list(fun = "ns", knots = varknots),
                   arglag = list(fun = "ns", knots = klag))
  ny_c <- length(unique(data_c$year))
  mod <- glm(death ~ cb + ns(time, df*ny_c) + factor(dow), family = quasipoisson(), data = data_c)
  mmt[i] <- findmin(cb, mod)
  red <- crossreduce(cb, mod, cen = mmt[i], by = 1)
  coef[i, ] <- red$coef
  vcov[[i]] <- red$vcov
}

# Mixed effects meta-analysis
mv <- mixmeta(coef ~ 1, vcov, method = "reml")

# Predict using temperature in percentile scale (%)
pct.seq <- seq(0, 100, by = 1)
basis.pred <- onebasis(pct.seq, fun = "ns", knots = varknots)
mmt2 <- findmin(basis.pred, coef = coef(mv), vcov = vcov(mv))
pred2 <- crosspred(basis.pred, coef = coef(mv), vcov = vcov(mv), model.link = "log", cen = mmt2, by = 1, at = pct.seq)

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
             type = "Two-stage", coun = "England & Wales")
)

# Overall plot for comparison labeling percentiles (%)
pct.seq <- seq(0, 100, by = 10)
temp.at <- quantile(data$tmean, probs = pct.seq / 100, na.rm = TRUE)
x_labels <- paste0(pct.seq, "\n(", round(temp.at, 1), ")")

ggplot(dfpred, aes(x, y, color=type)) + 
  geom_hline(yintercept = 1, size = 0.25) +
  geom_vline(aes(xintercept = cen, color = type), linetype = "dashed", size = 0.25, show.legend = FALSE) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = .2, size=0.05, show.legend = FALSE) +
  facet_wrap(vars(coun))+
  scale_x_continuous(breaks = pct.seq, labels = x_labels) +
  scale_y_continuous(breaks = seq(1, 1.4, by = 0.1), labels = seq(1, 1.4, by = 0.1), limits=c(0.99, 1.4)) +
  scale_color_manual(values = colors2[c(3, 4)], labels = c("One-stage", "Two-stage")) +
  scale_fill_manual(values = colors2[c(3, 4)], labels = c("One-stage", "Two-stage")) +
  labs(x = "Percentile, % (Temperature, °C)", y = "Relative Risk", color = NULL, fill = NULL) +
  theme_cc()

################################################################################
# End of script file
################################################################################
