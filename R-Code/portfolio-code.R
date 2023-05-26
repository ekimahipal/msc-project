# R Script - Exploring the Power of Clustering for Portfolio Selection Using Spline Coefficients
# Eki Mahipal 
# Student ID: 20377383
# psxem8@nottingham.ac.uk

# Install and Load Packages
install.packages("pacman")
install.packages('devtools')
install.packages('mgcv')
library(mgcv)
library(devtools)
library(pacman)
p_load(dplyr, tidyverse, stats,
       mgcv, cluster, magrittr, tidyquant,
       quantmod, PerformanceAnalytics,
       zoo, lubridate, factoextra, kableExtra,
       patchwork, gridExtra, ggplot2, xtable)


# 1. DATA COLLECTIOM AMD PREPROCESSING ------------------------------------

# 1a. Load data from csv
datalink <- ""
kompasdata <- read.csv('data/kompas100list.csv', sep=';')
kompasdata$Kode <- paste(kompasdata$Kode, '.JK', sep='') # All of Indonesian stocks symbols end with ".JK"

# 1b. Get the stock prices from Yahoo! Finance
startdate <- as.Date('2012-01-01')
enddate <- as.Date('2022-12-31')
kompas <- as.data.frame(tq_get(kompasdata$Kode,
                 from = startdate,
                 to = enddate))

kompas$date <- as.character(kompas$date) 
print(xtable(head(kompas)), type = "latex") # export R Table to Latex
print(xtable(tail(kompas)), type = "latex") # export R Table to Latex
kompas$date <- as.Date(kompas$date)

kompas <- kompas %>%
  select(date, symbol, close)


# 1c. Check starting date 2012-01-02 
## Function to check for NA 
find_NA <- function (x) {
  x[which(is.na(x$close)),]
}

checkstartdate <- find_NA(kompas100 <- kompas %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = close) %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "close") %>%
  arrange(symbol, date) %>%
  as.data.frame()) %>%
  filter(date == startdate + 1)

## Remove symbols included in checkstardate
kompas100 <- kompas100[!(kompas100$symbol %in% checkstartdate$symbol),]

# 1d. Find NA and treat with interpolation
find_NA(kompas100)
kompas100$close <- na.approx(kompas100$close)

# 1e. Visualisation of Preprocessed KOMPAS100 Index daily prices
k100plot1 <- kompas100 %>%
  ggplot(aes(x = date, y = close)) +
  geom_line(aes(color = symbol), show.legend = FALSE) +
  xlab("Year") +
  ylab("Price") +
  ggtitle("KOMPAS100 Index Price") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        # panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank()) 

ggsave(k100plot1, file = "image/Plot 1 - KOMPAS100 Index Prices 2012-2022.png", 
       width = 2000, height = 1500, units = "px", dpi = 220)

# 1f. Data Normalisation
## Raw Data
kompas100p <- kompas100 %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = close)

## Normalised Data (Z-Score)
kompas100sc <- kompas100p %>%
  mutate_at(c(names(kompas100p[2:ncol(kompas100p)])), ~(scale(.)[,1])) %>%
  as.data.frame()

## Visualisation of Normalised Data
k100plot2 <- kompas100sc %>%
  as_tibble() %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot(aes(x = date, y = price)) +
  geom_line(aes(color = symbol), show.legend = FALSE) +
  xlab("Year") +
  ylab("Normalised Price") +
  ggtitle("KOMPAS100 Index Price Normalised") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        # panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank()) 

k100plot2
ggsave(k100plot2, file = "image/Plot 2 - KOMPAS100 Index Normalised Prices 2012-2022.png", 
       width = 2000, height = 1500, units = "px", dpi = 220)

  

# 2. P-SPLINE Modelling with Normalized Data ---------------------------------------------------

## Test using different parameter
## k=10 on BMRI.JK
bmri <- kompas100sc$BMRI.JK %>%
  t() %>%
  as.matrix()

k10 <- gam(BMRI.JK ~ s(as.numeric(date), k = 10, bs = "ps", m = c(2,3)), 
           data = kompas100sc)
summary(k10)
k10coef <- coef(k10)
k10knots <- cbind('date' = k10$smooth[[1]]$knots)
k10iknots <- k10knots %>%
  as_tibble() %>%
  filter(date >= as.Date("2012-01-01") & date <= as.Date("2022-12-31")) %>%
  as.matrix()

## Visualisation: Comparison between raw and spline k= 10
plot(kompas100sc$date, kompas100sc$BMRI.JK, col = 1, 
     main="k=10", type = "l",
     xlab = "Date", ylab = "Price")
lines(kompas100sc$date, fitted(k10), col = 2, lwd = 2)
abline(v = as.Date(k10knots), lty = 2, col = "grey")

## ggplot version
k100plot3 <- kompas100sc %>%
  as_tibble() %>%
  pivot_longer(cols = BMRI.JK, 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = date, y = price), color = "grey", show.legend = FALSE) +
  xlab("Year") +
  ylab("Price") +
  ggtitle("k = 10") +
  geom_line(aes(x = date, y = fitted(k10)), color = "red") +
  geom_vline(xintercept = k10iknots, 
             linetype = 'dashed', color = "#555555") + 
  geom_vline(xintercept = k10knots[!k10knots[,1] %in% k10iknots[,1]], 
             linetype = 'dashed', color = "#BBBBBB") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank()) 


## Apply to all stocks
spline10 <- data.frame()
spline10list <- list()
for (n in 2:ncol(kompas100sc)) {
  value <- as_tibble(kompas100sc)[,n] %>%
    pull()
  ps10 <- gam(value ~ s(as.numeric(kompas100sc$date), k=10, bs = "ps", m= c(2,3)), method = "REML")
  if (!"index" %in% colnames(spline10)) {
    spline10 <- data.frame(index = 1:length(fitted(ps10)), fitted(ps10))
    names(spline10)[n] <- colnames(kompas100sc[n])
  } else {
    spline10 <- cbind(spline10, fitted(ps10))
    names(spline10)[n] <- colnames(kompas100sc[n])
  }
  spline10list[[n]] <- ps10
}

## Get the spline coefficients
spline10coef <- data.frame(index = 1:length(k10$coefficients))
for (l in 2:length(spline10list)) {
    spline10coef <- cbind(spline10coef, coef(spline10list[[l]]))
    names(spline10coef)[l] <- colnames(kompas100p[,l])
}

## k=20 on BMRI.JK
# bmri <- kompas100sc$BMRI.JK %>%
#   t() %>%
#   as.matrix()

k20 <- gam(BMRI.JK ~ s(as.numeric(date), k = 20, bs = "ps", m = c(2,3)), 
           data = kompas100sc)
summary(k20)
k20coef <- coef(k20)
k20knots <- k20$smooth[[1]]$knots

par(mfrow = c(1,2))
plot(k20coef)
plot(k10)

## Visualisation: Comparison between raw and spline
plot(kompas100sc$date, kompas100sc$BMRI.JK, col = 1, 
     main="k=20", type = "l",
     xlab = "Date", ylab = "Price")
lines(kompas100sc$date, fitted(k20), col = 2, lwd = 2)
abline(v = as.Date(k20knots), lty = 2, col = "grey")


par(mfrow=c(1,1))

## Apply k = 20 to all stocks
spline20 <- data.frame()
spline20list <- list()
for (n in 2:ncol(kompas100sc)) {
  value <- as_tibble(kompas100sc)[,n] %>%
    pull()
  ps20 <- gam(value ~ s(as.numeric(kompas100sc$date), k=20, bs = "ps", m= c(2,3)), method = "REML")
  if (!"index" %in% colnames(spline20)) {
    spline20 <- data.frame(index = 1:length(fitted(ps20)), fitted(ps20))
    names(spline20)[n] <- colnames(kompas100sc[n])
  } else {
    spline20 <- cbind(spline20, fitted(ps20))
    names(spline20)[n] <- colnames(kompas100sc[n])
  }
  spline20list[[n]] <- ps20
}


## Get the spline coefficients
spline20coef <- data.frame(index = 1:length(k20$coefficients))
for (l in 2:length(spline20list)) {
  spline20coef <- cbind(spline20coef, coef(spline20list[[l]]))
  names(spline20coef)[l] <- colnames(kompas100p[,l])
}


## k=40 on BMRI.JK
bmri <- kompas100sc$BMRI.JK %>%
  t() %>%
  as.matrix()
k40 <- gam(BMRI.JK ~ s(as.numeric(date), k = 40, bs = "ps", m = c(2,3)), 
           data = kompas100sc)
summary(k40)
k40coef <- coef(k40)
k40knots <- k40$smooth[[1]]$knots

gam.check(k40)
summary.gam(k40)
lpred <- predict.gam(k40, type = "link")
mat %*% k40coef

print(data.frame(lpred) == fitted(k40))

## Visualisation: Comparison between raw and spline
plot(kompas100sc$date, kompas100sc$BMRI.JK, col = 1, 
     main="k=40", type = "l",
     xlab = "Date", ylab = "Price")
lines(kompas100sc$date, fitted(k40), col = 2, lwd = 2)
lines(k40coef)
abline(v = as.Date(k40knots), lty = 2, col = "grey")

install.packages("scales")
coefindex <- rescale(index(k40coef), 
                     to = c(15341, 
                     19356))

k40coef <- cbind('date' = coefindex, k40coef)
  
## Apply k =40 to all stocks
spline40 <- data.frame()
spline40list <- list()
for (n in 2:ncol(kompas100sc)) {
  value <- as_tibble(kompas100sc)[,n] %>%
    pull()
  ps40 <- gam(value ~ s(as.numeric(kompas100sc$date), k=40, bs = "ps", m= c(2,3)), method = "REML")
  if (!"date" %in% colnames(spline40)) {
    spline40 <- data.frame(date = kompas100sc$date, fitted(ps40))
    names(spline40)[n] <- colnames(kompas100sc[n])
  } else {
    spline40 <- cbind(spline40, fitted(ps40))
    names(spline40)[n] <- colnames(kompas100sc[n])
  }
  spline40list[[n]] <- ps40
}


## Get the spline coefficients
spline40coef <- data.frame(index = 1:length(k40$coefficients))
for (l in 2:length(spline40list)) {
  spline40coef <- cbind(spline40coef, coef(spline40list[[l]]))
  names(spline40coef)[l] <- colnames(kompas100p[,l])
}

# Visualise the spline model
k100plot6 <- spline40 %>%
  as_tibble() %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot(aes(x = date, y = price)) +
  geom_line(aes(color = symbol), show.legend = FALSE) +
  xlab("Year") +
  ylab("Normalised Price") +
  ggtitle("KOMPAS100 Index Price Normalised") +
  theme_light()

# 3. CLUSTERING -----------------------------------------------------------

# 3. Clustering the KOMPAS100 Index
## Raw Data
kompas100raw <- kompas100sc %>%
  select(-date) %>%
  t() %>%
  as.data.frame()

## Spline Coefficients
## K = 10
spline10v <- spline10coef %>%
  select(-index) %>%
  t() %>% 
  as.data.frame()

## K = 20
spline20v <- spline20coef %>%
  select(-index) %>%
  t() %>% 
  as.data.frame()

## K = 40 - This model is selected for dissertation
spline40v <- spline40coef %>%
  select(-index) %>%
  t() %>% 
  as.data.frame()

# 3a. Find optimal number of clusters
## Estimate optimum number of clusters (Elbow Method)
## Create function to find no of clusters using wss plot 
## (https://www.kaggle.com/code/biswam87/exploring-winequality-with-clustering#Clustering-methods-and-plots)
wssplot <- function (data, nc=10) {
  set.seed(300)
  wsstot <- c()
  for (i in 1:nc) {
    cl <- kmeans(data, centers = i, nstart = 50, iter.max = 50)
    wsstot[i] <- cl$tot.withinss
  }
  plot(x = 1:nc, y =  wsstot,
       type = "b", xlab = "Number of Clusters", ylab = "Within Groups Sum of Squares")
}
## Run elbow method on the KOMPAS100 data
wssplot(spline10v)
wssplot(spline20v)
wssplot(spline40v)
wssplot(kompas100raw)

## Estimate number of clusters with factoextra
## (https://uc-r.github.io/kmeans_clustering)
## (https://rstudio-pubs-static.s3.amazonaws.com/375287_5021917f670c435bb0458af333716136.html)
## Elbow method
set.seed(300)
rplot1 <- fviz_nbclust(kompas100raw, kmeans, method = "wss")
cplot1 <- fviz_nbclust(spline10v, kmeans, method = "wss")
cplot12 <- fviz_nbclust(spline20v, kmeans, method = "wss")
cplot14 <- fviz_nbclust(spline40v, kmeans, method = "wss")

## Silhouette method
set.seed(300)
rplot2 <- fviz_nbclust(kompas100raw, kmeans, method = "silhouette")
cplot2 <- fviz_nbclust(spline10v, kmeans, method = "silhouette")
cplot22 <- fviz_nbclust(spline20v, kmeans, method = "silhouette")
cplot24 <- fviz_nbclust(spline40v, kmeans, method = "silhouette")

## Gap statistic method
set.seed(300)
gapstat <- clusGap(spline40v, FUN = kmeans, d.power = 2, nstart = 200, K.max = 10, B=200)
fviz_gap_stat(gapstat)
print(gapstat, method = "firstmax")

rplot3 <- fviz_nbclust(kompas100raw, kmeans, method = "gap_stat", nboot = 200)
cplot3 <- fviz_nbclust(spline10v, kmeans, method = "gap_stat", nboot = 200)
cplot32 <- fviz_nbclust(spline20v, kmeans, method = "gap_stat", nboot = 200)
cplot34 <- fviz_nbclust(spline40v, kmeans, method = "gap_stat", nboot = 200)

grid.arrange(rplot1, rplot2, rplot3, nrow=1) # Shows optimum number of cluster is 2
grid.arrange(cplot1, cplot2, cplot3, nrow=1) # Shows optimum number of cluster is 2
grid.arrange(cplot12, cplot22, cplot32, nrow=1) 
grid.arrange(cplot14, cplot24, cplot34, nrow=1)

## Run cluster algorithm (K-Means)
## 2 Clusters
set.seed(300)
km_k100r <- kmeans(kompas100raw, 2, nstart = 50, iter.max = 50)
fviz_cluster(km_k100r, data = kompas100raw)
km_k100r$cluster

set.seed(300)
km_k100c <- kmeans(spline10v, 2, nstart = 50, iter.max = 50)
fviz_cluster(km_k100c, data = spline10v)
km_k100c$cluster

set.seed(300)
km_k100c2 <- kmeans(spline20v, 3, nstart = 50, iter.max = 50)
fviz_cluster(km_k100c2, data = spline20v)
km_k100c2$cluster

set.seed(300)
km_k100c4 <- kmeans(spline40v, 3, nstart = 50, iter.max = 50)
fviz_cluster(km_k100c4, data = spline40v)
km_k100c4$cluster

## Create dataframe for stock clusters
## 2 Clusters
## Raw Prices
k100r_cl <- km_k100r$cluster %>%
  as.data.frame() %>%
  rename(cluster = 1) %>%
  mutate(symbol = rownames(.)) %>%
  arrange(cluster)
## Fitted Spline
k100f_cl <- km_k100f$cluster %>%
  as.data.frame() %>%
  rename(cluster = 1) %>%
  mutate(symbol = rownames(.)) %>%
  arrange(cluster)