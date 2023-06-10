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
       patchwork, gridExtra, ggplot2, xtable, ggpubr,
       tseries, grid)


# #### 1. DATA COLLECTION AND PREPROCESSING #### #------------------------------------

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
# KOMPAS100 Index Prices for evaluation purposes
k100index <- read.csv('data/IDXKOMPAS100.csv', sep = ",")
k100index$date <- strptime(k100index$Date, "%m/%d/%Y")
k100index$date <- as.Date(k100index$date)
k100index <- k100index %>%
  select(date, Price) %>%
  arrange(date, Price)

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
  scale_y_continuous(breaks = seq(0, 100000, 25000)) +
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

# 1f. Data Normalisation/Scaled
## Raw Data
kompas100p <- kompas100 %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = close)

## Normalised Data (Z-Score)
kompas100sc <- kompas100p %>%
  mutate_at(c(names(kompas100p[2:ncol(kompas100p)])), ~(scale(.)[,1])) %>%
  as.data.frame()

## Visualisation of Scaled Data
k100plot2 <- kompas100sc %>%
  as_tibble() %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot(aes(x = date, y = price)) +
  geom_line(aes(color = symbol), show.legend = FALSE) +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4, 6, 2)) +
  xlab("Year") +
  ylab("Normalised Price") +
  ggtitle("KOMPAS100 Index Price (Scaled)") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        # panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank()) 

k100plot2
ggsave(k100plot2, file = "image/Plot 2 - KOMPAS100 Index Normalised Prices 2012-2022.png", 
       width = 2000, height = 1500, units = "px", dpi = 220)


# #### 2A. P-SPLINE Modelling with Raw Data #### #------------------------------------
bmri <- kompas100sc$BMRI.JK %>%
  t() %>%
  as.matrix()
str(kompas100sc$BMRI.JK)

## Test using different parameter
## k=14 (10 knots) on BMRI.JK
k14raw <- gam(BMRI.JK ~ s(as.numeric(date), k = 14, bs = "ps", m = c(2,3)), 
           data = kompas100p)
summary(k14raw)
k14rawcoef <- coef(k14raw)
k14rawknots <- cbind('date' = k14raw$smooth[[1]]$knots)
k14rawiknots <- k14rawknots %>%
  as_tibble() %>%
  filter(date >= as.Date("2012-01-01") & date <= as.Date("2022-12-31")) %>%
  as.matrix()

## Visualisation
k100plot3 <- kompas100p %>%
  as_tibble() %>%
  pivot_longer(cols = BMRI.JK, 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = date, y = price), color = "black", show.legend = FALSE) +
  xlab("Year") +
  ylab("Price") +
  ggtitle("knots = 10") +
  geom_line(aes(x = date, y = fitted(k14raw)), color = "red") +
  geom_vline(xintercept = k14rawiknots, 
             linetype = 'dashed', color = "#555555") + 
  geom_vline(xintercept = k14rawknots[!k14rawknots[,1] %in% k14rawiknots[,1]], 
             linetype = 'dashed', color = "#BBBBBB") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank()) 

## k=24 (20 knots) on BMRI.JK
k24raw <- gam(BMRI.JK ~ s(as.numeric(date), k = 24, bs = "ps", m = c(2,3)), 
              data = kompas100p)
summary(k24raw)
k24rawcoef <- coef(k24raw)
k24rawknots <- cbind('date' = k24raw$smooth[[1]]$knots)
k24rawiknots <- k24rawknots %>%
  as_tibble() %>%
  filter(date >= as.Date("2012-01-01") & date <= as.Date("2022-12-31")) %>%
  as.matrix()

## Visualisation
k100plot4 <- kompas100p %>%
  as_tibble() %>%
  pivot_longer(cols = BMRI.JK, 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = date, y = price), color = "black", show.legend = FALSE) +
  xlab("Year") +
  ylab("Price") +
  ggtitle("knots = 20") +
  geom_line(aes(x = date, y = fitted(k24raw)), color = "red") +
  geom_vline(xintercept = k24rawiknots, 
             linetype = 'dashed', color = "#555555") + 
  geom_vline(xintercept = k24rawknots[!k24rawknots[,1] %in% k24rawiknots[,1]], 
             linetype = 'dashed', color = "#BBBBBB") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank()) 

## k=44 (40 knots) on BMRI.JK
k44raw <- gam(BMRI.JK ~ s(as.numeric(date), k = 44, bs = "ps", m = c(2,3)), 
              data = kompas100p)
summary(k44raw)
k44rawcoef <- coef(k44raw)
k44rawknots <- cbind('date' = k44raw$smooth[[1]]$knots)
k44rawiknots <- k44rawknots %>%
  as_tibble() %>%
  filter(date >= as.Date("2012-01-01") & date <= as.Date("2022-12-31")) %>%
  as.matrix()


## Visualisation
k100plot5 <- kompas100p %>%
  as_tibble() %>%
  pivot_longer(cols = BMRI.JK, 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = date, y = price), color = "black", show.legend = FALSE) +
  xlab("Year") +
  ylab("Price") +
  ggtitle("knots = 40") +
  geom_line(aes(x = date, y = fitted(k44raw)), color = "red") +
  geom_vline(xintercept = k44rawiknots, 
             linetype = 'dashed', color = "#555555") + 
  geom_vline(xintercept = k44rawknots[!k44rawknots[,1] %in% k44rawiknots[,1]], 
             linetype = 'dashed', color = "#BBBBBB") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank()) 

figcol <- ggarrange(k100plot3,k100plot4, k100plot5, labels = c("a", "b", "c"),
                    nrow = 2, ncol = 2)

figcol
ggsave(figcol, file = "image/Plot 3 - P-Spline Model Compariosn.png", 
       width = 1600, height = 1600, units = "px", dpi = 150)

## Apply to all stocks
spline44raw <- data.frame()
spline44rawlist <- list()
for (n in 2:ncol(kompas100p)) {
  value <- as_tibble(kompas100p)[,n] %>%
    pull()
  ps44raw <- gam(value ~ s(as.numeric(kompas100p$date), k=44, bs = "ps", m= c(2,3)))
  if (!"date" %in% colnames(spline44raw)) {
    spline44raw <- data.frame(date = 1:length(fitted(ps44raw)), fitted(ps44raw))
    names(spline44raw)[n] <- colnames(kompas100p[n])
  } else {
    spline44raw <- cbind(spline44raw, fitted(ps44raw))
    names(spline44raw)[n] <- colnames(kompas100p[n])
  }
  spline44rawlist[[n]] <- ps44raw
}


## Get the spline coefficients
spline44rawcoef <- data.frame(index = 1:length(k44raw$coefficients))
for (l in 2:length(spline44rawlist)) {
  spline44rawcoef <- cbind(spline44rawcoef, coef(spline44rawlist[[l]]))
  names(spline44rawcoef)[l] <- colnames(kompas100p[,l])
}

# Visualise the spline model
k100plot6 <- spline44raw %>%
  as_tibble() %>%
  mutate(date = kompas100p$date) %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot(aes(x = date, y = price)) +
  geom_line(aes(color = symbol), show.legend = FALSE) +
  scale_y_continuous(breaks = seq(0, 100000, 25000)) +
  xlab("Year") +
  ylab("Price") +
  ggtitle("P-Spline Model of KOMPAS100 Index Price") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        # panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank())

figcol2 <- ggarrange(k100plot1, k100plot6, nrow = 1, ncol = 2)

ggsave(figcol2, file = "image/Plot 4 - P-Spline Model All Stocks.png", 
       width = 1600, height = 800, units = "px", dpi = 150)



# #### 2B. P-SPLINE Modelling with Scaled Data ####---------------------------------------------------

## k=44 on BMRI.JK
bmri_sc <- kompas100sc$BMRI.JK %>%
  t() %>%
  as.matrix()

k44sc <- gam(BMRI.JK ~ s(as.numeric(date), k = 44, bs = "ps", m = c(2,3)), 
           data = kompas100sc)
summary(k44sc)
k44sccoef <- coef(k44sc)
k44scknots <- cbind('date' = k44sc$smooth[[1]]$knots)
k44sciknots <- k44scknots %>%
  as_tibble() %>%
  filter(date >= as.Date("2012-01-01") & date <= as.Date("2022-12-31")) %>%
  as.matrix()

  
## Apply k =44 to all stocks
spline44sc <- data.frame()
spline44sclist <- list()
for (n in 2:ncol(kompas100sc)) {
  value <- as_tibble(kompas100sc)[,n] %>%
    pull()
  ps44sc <- gam(value ~ s(as.numeric(kompas100sc$date), k=44, bs = "ps", m= c(2,3)))
  if (!"date" %in% colnames(spline44sc)) {
    spline44sc <- data.frame(date = kompas100sc$date, fitted(ps44sc))
    names(spline44sc)[n] <- colnames(kompas100sc[n])
  } else {
    spline44sc <- cbind(spline44sc, fitted(ps44sc))
    names(spline44sc)[n] <- colnames(kompas100sc[n])
  }
  spline44sclist[[n]] <- ps44sc
}


## Get the spline coefficients
spline44sccoef <- data.frame(index = 1:length(k44sc$coefficients))
for (l in 2:length(spline44sclist)) {
  spline44sccoef <- cbind(spline44sccoef, coef(spline44sclist[[l]]))
  names(spline44sccoef)[l] <- colnames(kompas100p[,l])
}

# Visualise the spline model
k100plot7 <- spline44sc %>%
  as_tibble() %>%
  mutate(date = kompas100sc$date) %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% as.data.frame() %>%
  ggplot(aes(x = date, y = price)) +
  geom_line(aes(color = symbol), show.legend = FALSE) + 
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4, 6, 2)) +
  xlab("Year") +
  ylab("Price") +
  ggtitle("P-Spline Model of KOMPAS100 Index Price (Scaled)") +
  theme_light() + theme(axis.line = element_line(colour = "black"),
                        # panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black"),
                        panel.background = element_blank())

figcol3 <- ggarrange(k100plot2, k100plot7, nrow = 1, ncol = 2)

ggsave(figcol3, file = "image/Plot 8 - P-Spline Model Scaled).png", 
       width = 1600, height = 800, units = "px", dpi = 160)

# #### 3A. CLUSTERING (K-MEANS) #### -----------------------------------------------------------

## Raw Data
data_raw <- kompas100p %>%
  select(-date) %>%
  t() %>%
  as.data.frame()

data_raw_sc <- kompas100sc %>%
  select(-date) %>%
  t() %>%
  as.data.frame()

## Spline Coefficients
data_spline <- spline44rawcoef %>%
  select(-index) %>%
  t() %>% 
  as.data.frame()

## Spline Coefficients Scaled
data_spline_sc <- spline44sccoef %>%
  select(-index) %>%
  t() %>%
  as.data.frame()

## Find Optimal Number of Clusters ----------------------------
## (https://rstudio-pubs-static.s3.amazonaws.com/375287_5021917f670c435bb0458af333716136.html)
## Elbow Method ---------------------------
set.seed(26)
elbow_raw <- fviz_nbclust(data_raw, kmeans, method = "wss")
set.seed(26)
elbow_spline <- fviz_nbclust(data_spline, kmeans, method = "wss")
## Calculate the computation time
set.seed(26)
system.time({fviz_nbclust(data_raw, kmeans, method = "wss")})
set.seed(26)
system.time({fviz_nbclust(data_spline, kmeans, method = "wss")})

## Silhouette Method -----------------------
set.seed(26)
sil_raw <- fviz_nbclust(data_raw, kmeans, method = "silhouette")
set.seed(26)
sil_spline <- fviz_nbclust(data_spline, kmeans, method = "silhouette")
## Calculate computation time
set.seed(26)
system.time({fviz_nbclust(data_raw, kmeans, method = "silhouette")})
set.seed(26)
system.time({fviz_nbclust(data_spline, kmeans, method = "silhouette")})
## Print the result as latex table
sil_raw_score <- t(sil_raw$data) %>% as.data.frame()
sil_spline_score <- t(sil_spline$data) %>% as.data.frame()
print(xtable(sil_raw_score), type = "latex")
print(xtable(sil_spline_score), type = "latex")

## Visualisation of elbow and silhouette --
grid.arrange(elbow_raw, sil_raw, nrow =1)
grid.arrange(elbow_spline, sil_spline, nrow =1)

## Elbow and Silhouette Method on Scaled Data ---------------------
set.seed(26)
elbow_raw_sc <- fviz_nbclust(data_raw_sc, kmeans, method = "wss")
elbow_spline_sc <- fviz_nbclust(data_spline_sc, kmeans, method = "wss")
system.time({fviz_nbclust(data_spline_sc, kmeans, method = "wss")})
set.seed(26)
sil_raw_sc <- fviz_nbclust(data_raw_sc, kmeans, method = "silhouette")
set.seed(26)
sil_spline_sc <- fviz_nbclust(data_spline_sc, kmeans, method = "silhouette")
system.time({fviz_nbclust(data_spline_sc, kmeans, method = "silhouette")})
## Print the result as latex table
sil_spline_sc_score <- t(sil_spline_sc$data) %>% as.data.frame()
print(xtable(sil_spline_sc_score), type = "latex")
## Visualisation of elbow and silhouette
grid.arrange(elbow_raw_sc, sil_raw_sc, nrow = 1)
grid.arrange(sil_raw_sc, sil_spline_sc, nrow =1)

## Run K-Means Algorithm ---------------------------------------
## Raw Data
set.seed(26)
km_raw <- kmeans(data_raw, centers = 2, nstart = 50, iter.max = 50)
system.time({kmeans(data_raw, centers = 2, nstart = 50, iter.max = 50)})
viz_raw <- fviz_cluster(km_raw, data = data_raw, labelsize = 9)

## Raw Data Scaled
set.seed(26)
km_raw_sc <- kmeans(data_raw_sc, centers = 2, nstart = 50, iter.max = 50)
system.time({kmeans(data_raw_sc, centers = 2, nstart = 50, iter.max = 50)})
viz_sc <- fviz_cluster(km_raw_sc, data = data_raw_sc, ggtheme = theme_minimal())
km_raw_sc$cluster

## Spline Coefficients
set.seed(26)
km_spline <- kmeans(data_spline, centers = 2, nstart = 50, iter.max = 50)
system.time({kmeans(data_spline, centers = 2, nstart = 50, iter.max = 50)})
viz_spline <- fviz_cluster(km_spline, data = data_spline, labelsize = 9)
km_spline$cluster

grid.arrange(viz_spline, viz_raw, nrow = 1)

## Spline Coefficients Scaled
set.seed(26)
km_spline_sc <- kmeans(data_spline_sc, centers = 3, nstart = 50, iter.max = 50)
system.time({kmeans(data_spline_sc, centers = 3, nstart = 50, iter.max = 50)})
viz_spline_sc <- fviz_cluster(km_spline_sc, data = data_spline_sc, ggtheme = theme_minimal())
km_spline_sc$cluster

set.seed(26)
km_spline_sc <- kmeans(data_spline_sc, centers = 2, nstart = 50, iter.max = 50)
fviz_cluster(km_spline_sc, data = data_spline_sc, ggtheme = theme_minimal())

figcol5 <- ggarrange(viz_sc, viz_spline_sc, nrow = 1)
ggsave(figcol5, file = "image/Plot 10 - Recognised Cluster (Scaled Data).jpg", 
       width = 2000, height = 1000, units = "px", dpi = 150)

## Create dataframe for visualising stock clusters ------------------------------
## RAW Prices (not scaled) -----------------------------------------
## Raw Data
cl_raw <- km_raw$cluster %>%
  as.data.frame() %>%
  rename(cluster = 1) %>%
  mutate(symbol = rownames(.)) %>%
  arrange(cluster)
## Spline Coefficients
cl_spline <- km_spline$cluster %>%
  as.data.frame() %>%
  rename(cluster = 1) %>%
  mutate(symbol = rownames(.)) %>%
  arrange(cluster)
## Raw stock prices (2 clusters)
cl_raw_plot1 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_raw$symbol[cl_raw$cluster==1]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "red") +
  coord_cartesian(ylim=c(0,100000)) +
  scale_y_continuous(breaks = seq(0, 100000, 25000)) +
  ggtitle("Cluster 1") + 
  theme_light()  

cl_raw_plot2 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_raw$symbol[cl_raw$cluster==2]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "green") +
  coord_cartesian(ylim=c(0,100000)) +
  scale_y_continuous(breaks = seq(0, 100000, 25000)) +
  ggtitle("Cluster 2") + 
  theme_light() 
## Spline Coefficients (2 clusters)
cl_spline_plot1 <- spline44raw %>%
  mutate(date = kompas100p$date) %>%
  pivot_longer(cols = 2:ncol(spline44raw),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_spline$symbol[cl_spline$cluster==1]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "red") +
  coord_cartesian(ylim=c(0,100000)) +
  scale_y_continuous(breaks = seq(0, 100000, 25000)) +
  ggtitle("Cluster 1 Spline Coef") + 
  theme_light()  

cl_spline_plot2 <- spline44raw %>%
  mutate(date = kompas100p$date) %>%
  pivot_longer(cols = 2:ncol(spline44raw),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_spline$symbol[cl_spline$cluster==2]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "green") +
  coord_cartesian(ylim=c(0,100000)) +
  scale_y_continuous(breaks = seq(0, 100000, 25000)) +
  ggtitle("Cluster 2 Spline Coef") + 
  theme_light() 

figcol4 <- ggarrange(cl_raw_plot1, cl_raw_plot2, cl_spline_plot1, cl_spline_plot2, nrow = 2, ncol=2)

ggsave(figcol4, file = "image/Plot 7 - Recognised Cluster (Raw Data).png", 
       width = 1600, height = 800, units = "px", dpi = 180)


## Scaled Data -----------------------------------------
## Raw Data
cl_raw_sc <- km_raw_sc$cluster %>%
  as.data.frame() %>%
  rename(cluster = 1) %>%
  mutate(symbol = rownames(.)) %>%
  arrange(cluster)
## Spline Coefficients
cl_spline_sc <- km_spline_sc$cluster %>%
  as.data.frame() %>%
  rename(cluster = 1) %>%
  mutate(symbol = rownames(.)) %>%
  arrange(cluster)
## Raw stock prices (2 clusters)
cl_raw_sc_plot1 <- kompas100sc %>%
  pivot_longer(cols = 2:ncol(kompas100sc),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_raw_sc$symbol[cl_raw_sc$cluster==1]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "red") +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4,6,2)) +
  ggtitle("Cluster 1") + 
  theme_light()  

cl_raw_sc_plot2 <- kompas100sc %>%
  pivot_longer(cols = 2:ncol(kompas100sc),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_raw_sc$symbol[cl_raw_sc$cluster==2]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "green") +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4,6,2)) +
  ggtitle("Cluster 2") + 
  theme_light() 

## Spline Coefficients Scaled (2 clusters)
cl_spline_sc_plot1 <- spline44sc %>%
  mutate(date = kompas100sc$date) %>%
  pivot_longer(cols = 2:ncol(spline44sc),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_spline_sc$symbol[cl_spline_sc$cluster==2]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "red") +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4,6,2)) +
  ggtitle("Cluster 1 Spline Coef") + 
  theme_light()  

cl_spline_sc_plot2 <- spline44sc %>%
  mutate(date = kompas100sc$date) %>%
  pivot_longer(cols = 2:ncol(spline44sc),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_spline_sc$symbol[cl_spline_sc$cluster==1]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "green") +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4,6,2)) +
  ggtitle("Cluster 2 Spline Coef") + 
  theme_light() 

figcol6 <- ggarrange(cl_raw_sc_plot1, cl_raw_sc_plot2, cl_spline_sc_plot1, cl_spline_sc_plot2, nrow = 2, ncol=2)

ggsave(figcol6, file = "image/Plot 11 - Recognised Cluster (Scaled Data).png", 
       width = 1600, height = 800, units = "px", dpi = 180)




# #### 3B. CLUSTERING (HIERARCHICAL) ##### -------------------------------------------

# Set Data
data_raw_sc
data_spline_sc

# Set Distance Matrix
dist_raw_sc <- dist(data_raw_sc, method = "euclidean")
dist_spline_sc <- dist(data_spline_sc, method = "euclidean")

# Define Linkage 
# (https://www.statology.org/hierarchical-clustering-in-r/)
linkage <- c( "average", "single", "complete", "ward")
names(linkage) <- c( "average", "single", "complete", "ward")

# Function to Compute Agglomerative Coefficient
agcoef_raw_sc <- function(l) {
  agnes(dist_raw_sc, method = l)$ac
}
agcoef_spline_sc <- function(l) {
  agnes(dist_spline_sc, method = l)$ac
}

# Compute and Compare Agglomerative Coefficient of All Linkage Method
agcoef1 <- as.data.frame(sapply(linkage, agcoef_raw_sc)) # ward method is the highest
agcoef2 <- as.data.frame(sapply(linkage, agcoef_spline_sc)) # ward method is the highest
print(xtable(agcoef1), type = "latex")
print(xtable(agcoef2), type = "latex")

# Run Hierarchical with Ward's Linkage Method
# Scaled Data
set.seed(12)
hc_raw_sc <- agnes(data_raw_sc, 
                   metric = "euclidean",
                   method = "complete") 
fviz_nbclust(as.matrix(dist_raw_sc), hcut, method = "silhouette") # find optimum number of cluster
hctree_raw_sc <- fviz_dend(hc_raw_sc, cex = 0.6, lwd = 0.6, 
                           k = 2, main = "Cluster Dendogram - Scaled Data",
                           k_colors = c("jco"),
                           rect = TRUE, 
                           rect_border = "jco", 
                           rect_fill = TRUE) # visualisation
## Create dataframe with cluster result
hcl_raw_sc <- data.frame(cutree(hc_raw_sc, k = 2))
hcl_raw_sc <- hcl_raw_sc %>%
  rename(cluster = 1) %>%
  mutate(symbol = rownames(data_raw_sc)) %>%
  arrange(symbol)


# Spline Coefficients
set.seed(12)
hc_spline_sc <- agnes(data_spline_sc,
                      metric = "euclidean",
                      method = "ward")
fviz_nbclust(as.matrix(dist_spline_sc), hcut, method = "silhouette") # find optimum number of cluster
hctree_spline_sc <- fviz_dend(hc_spline_sc, cex = 0.6, lwd = 0.6, 
                              k = 2, main = "Cluster Dendogram - P-Spline Coefficients",
                              k_colors = c("jco"),
                              rect = TRUE, 
                              rect_border = "jco", 
                              rect_fill = TRUE) # visualisation
## Create dataframe with cluster result
hcl_spline_sc <- data.frame(cutree(hc_spline_sc, k = 2))
hcl_spline_sc <- hcl_spline_sc %>%
  rename(cluster = 1) %>%
  mutate(symbol = rownames(data_spline_sc)) %>%
  arrange(symbol)

grid.arrange(hctree_raw_sc, hctree_spline_sc, nrow = 2)


# Time series Visualisation
## Scaled data (2 clusters)
hcl_raw_sc_plot1 <- kompas100sc %>%
  pivot_longer(cols = 2:ncol(kompas100sc),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_raw_sc$symbol[hcl_raw_sc$cluster==1]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "#0066ff") +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4,6,2)) +
  ggtitle("Cluster 1") + 
  theme_light()  

hcl_raw_sc_plot2 <- kompas100sc %>%
  pivot_longer(cols = 2:ncol(kompas100sc),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_raw_sc$symbol[hcl_raw_sc$cluster==2]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "#ffee00") +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4,6,2)) +
  ggtitle("Cluster 2") + 
  theme_light() 

## Spline Coefficients Scaled (2 clusters)
hcl_spline_sc_plot1 <- spline44sc %>%
  mutate(date = kompas100sc$date) %>%
  pivot_longer(cols = 2:ncol(spline44sc),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_spline_sc$symbol[hcl_spline_sc$cluster==1]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "#0066ff") +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4,6,2)) +
  ggtitle("Cluster 1 Spline Coef") + 
  theme_light()  

hcl_spline_sc_plot2 <- spline44sc %>%
  mutate(date = kompas100sc$date) %>%
  pivot_longer(cols = 2:ncol(spline44sc),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_spline_sc$symbol[hcl_spline_sc$cluster==2]) %>%
  ggplot(aes(x = date, y = price, group = symbol)) +
  geom_line( color = "#ffee00") +
  coord_cartesian(ylim=c(-4,6)) +
  scale_y_continuous(breaks = seq(-4,6,2)) +
  ggtitle("Cluster 2 Spline Coef") + 
  theme_light() 

figcol7 <- ggarrange(hcl_raw_sc_plot1, hcl_raw_sc_plot2, hcl_spline_sc_plot1, hcl_spline_sc_plot2, nrow = 2, ncol=2)

ggsave(figcol7, file = "image/Plot 13 - Recognised Cluster (Scaled Data).png", 
       width = 1600, height = 800, units = "px", dpi = 180)


# ###### 4. PORTFOLIO SELECTION  #######--------------------------------------------------
# Calculate Return from each cluster
km.sc1 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_raw_sc$symbol[cl_raw_sc$cluster==1]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date") %>%
  Return.calculate() %>%
  slice(-1)

km.sc2 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_raw_sc$symbol[cl_raw_sc$cluster==2]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date") %>%
  Return.calculate() %>%
  slice(-1)

km.spline1 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_spline_sc$symbol[cl_spline_sc$cluster==2]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date") %>%
  Return.calculate() %>%
  slice(-1)

km.spline2 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_spline_sc$symbol[cl_spline_sc$cluster==1]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date") %>%
  Return.calculate() %>%
  slice(-1)

hc.sc1 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_raw_sc$symbol[hcl_raw_sc$cluster==1]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date") %>%
  Return.calculate() %>%
  slice(-1)

hc.sc2 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_raw_sc$symbol[hcl_raw_sc$cluster==2]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date") %>%
  Return.calculate() %>%
  slice(-1)

hc.spline1 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_spline_sc$symbol[hcl_spline_sc$cluster==1]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date") %>%
  Return.calculate() %>%
  slice(-1)

hc.spline2 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_spline_sc$symbol[hcl_spline_sc$cluster==2]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date") %>%
  Return.calculate() %>%
  slice(-1)

# ###### Select 10 stocks from each portfolios ###### #
# (https://www.youtube.com/watch?v=2EpIdnA0pPo) #
# Set KOMPAS100 Index as Benchmark
k100index_xts <- xts(k100index[,-1], order.by = k100index$date)
k100index_ret <- Return.calculate(k100index_xts)
k100index_pf <- Return.portfolio(k100index_ret[-1,])
k100index_pf_ <- Return.portfolio(k100index_ret[-1,], rebalance_on = "years",wealth.index = TRUE)

# #### Portfolio 1 #### #
sh1 <- SharpeRatio.annualized(xts(km.sc1, order.by = kompas100p$date[-1])) %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "symbol",
               values_to = "sharpe") %>%
  arrange(desc(sharpe)) %>% head(10)

pf1 <- km.sc1 %>%
  rownames_to_column(var = "date") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "return") %>%
  arrange(symbol) %>%
  filter(symbol %in% sh1$symbol) %>%
  pivot_wider(names_from = symbol, values_from = return)

# --- PF1 MV weighting --- #
pf1 <- xts(pf1[-1], order.by = as.Date(kompas100p$date[-1]))
pf1_mv <- portfolio.optim(pf1)$pw
names(pf1_mv) <- colnames(as.data.frame(pf1))
# Calculate MV Portfolio Return
pf1_mv_ret <- Return.portfolio(pf1, weights = pf1_mv)
pf1_mv_ret_ <- Return.portfolio(pf1, weights = pf1_mv, 
                                wealth.index = TRUE)

# --- PF1 EW Weighting --- #
pf1_ew <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
names(pf1_ew) <- colnames(as.data.frame(pf1))
# Calculate MV Portfolio Return
pf1_ew_ret <- Return.portfolio(pf1, weights = pf1_ew)
pf1_ew_ret_ <- Return.portfolio(pf1, weights = pf1_ew, 
                                wealth.index = TRUE)

# Plot - PF1 Return
plot(pf1_mv_ret, col = "#A89F68", main = "")
lines(pf1_ew_ret, col = "#F5C396")
# Plot - PF1 Wealth Index
plot(pf1_mv_ret_, col = "#A89F68", main = "Portfolio Wealth Index - PF1", lwd = 2)
lines(pf1_ew_ret_, col = "#F5C396" , lwd= 2)
lines(k100index_pf_, col = "#000000", lwd = 2)
addLegend("topleft", c("PF1 MV", "PF1 EW", "KOMPAS100"),
       col = c("#A89F68", "#F5C396", "#000000"), lty = 1, bty="o", lwd = 2)

# Latex Table
print(xtable(as.data.frame(cbind(pf1_mv, pf1_ew))), type = "latex")


# #### Portfolio 2 #### #
sh2 <- SharpeRatio.annualized(xts(km.sc2, order.by = kompas100p$date[-1])) %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "symbol",
               values_to = "sharpe") %>%
  arrange(desc(sharpe)) %>% head(10)

pf2 <- km.sc2 %>%
  rownames_to_column(var = "date") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "return") %>%
  arrange(symbol) %>%
  filter(symbol %in% sh2$symbol) %>%
  pivot_wider(names_from = symbol, values_from = return)

# --- PF2 MV weighting --- #
pf2 <- xts(pf2[-1], order.by = as.Date(kompas100p$date[-1]))
pf2_mv <- portfolio.optim(pf2)$pw
names(pf2_mv) <- colnames(as.data.frame(pf2))

# Calculate MV Portfolio Return
pf2_mv_ret <- Return.portfolio(pf2, weights = pf2_mv)
pf2_mv_ret_ <- Return.portfolio(pf2, weights = pf2_mv, 
                                wealth.index = TRUE)



# --- PF2 EW Weighting --- #
pf2_ew <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
names(pf2_ew) <- colnames(as.data.frame(pf2))
print(xtable(as.data.frame(t(pf2_ew))), type = "latex")
# Calculate MV Portfolio Return
pf2_ew_ret <- Return.portfolio(pf2, weights = pf2_ew)
pf2_ew_ret_ <- Return.portfolio(pf2, weights = pf2_ew, 
                                wealth.index = TRUE)
# Plot - PF2 Wealth Index
plot(pf2_mv_ret_, col = "#EE4266", main = "Portfolio Wealth Index - PF2", lwd = 2)
lines(pf2_ew_ret_, col = "#3CBBB1" , lwd= 2)
lines(k100index_pf_, col = "#000000", lwd = 2)
addLegend("topright", c("PF2 MV", "PF2 EW", "KOMPAS100"),
          col = c("#EE4266", "#3CBBB1", "#000000"), lty = 1, bty="o", lwd = 2)

# Latex Table
print(xtable(as.data.frame(cbind(pf2_mv, pf2_ew))), type = "latex")

# #### Portfolio 3 #### #
sh3 <- SharpeRatio.annualized(xts(km.spline1, order.by = kompas100p$date[-1])) %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "symbol",
               values_to = "sharpe") %>%
  arrange(desc(sharpe)) %>% head(10)


pf3 <- km.spline1 %>%
  rownames_to_column(var = "date") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "return") %>%
  arrange(symbol) %>%
  filter(symbol %in% sh3$symbol) %>%
  pivot_wider(names_from = symbol, values_from = return)

# --- PF3 MV weighting --- #
pf3 <- xts(pf3[-1], order.by = as.Date(kompas100p$date[-1]))
pf3_mv <- portfolio.optim(pf3)$pw
names(pf3_mv) <- colnames(as.data.frame(pf3))
print(xtable(as.data.frame(t(pf3_mv))), type = "latex")
# Calculate MV Portfolio Return
pf3_mv_ret <- Return.portfolio(pf3, weights = pf3_mv)
pf3_mv_ret_ <- Return.portfolio(pf3, weights = pf3_mv, 
                                wealth.index = TRUE)
# Compare MV with KOMPAS100 Index 
lines(pf3_mv_ret_, col = "#BB55DD", legend.loc = "topleft")
lines(k100index_pf_, col = "#555500", main = "KOMPAS100 Index")

# --- PF3 EW Weighting --- #
pf3_ew <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
names(pf3_ew) <- colnames(as.data.frame(pf3))
print(xtable(as.data.frame(t(pf3_ew))), type = "latex")
# Calculate MV Portfolio Return
pf3_ew_ret <- Return.portfolio(pf3, weights = pf3_ew)
pf3_ew_ret_ <- Return.portfolio(pf3, weights = pf3_ew, 
                                wealth.index = TRUE)
# Plot - PF3 Wealth Index
plot(pf3_mv_ret_, col = "#881600", main = "Portfolio Wealth Index - PF3", lwd = 2)
lines(pf3_ew_ret_, col = "#E3E36A" , lwd= 2)
lines(k100index_pf_, col = "#000000", lwd = 2)
addLegend("topleft", c("PF3 MV", "PF3 EW", "KOMPAS100"),
          col = c("#881600", "#E3E36A", "#000000"), lty = 1, bty="o", lwd = 2)

# Latex Table
print(xtable(as.data.frame(cbind(pf3_mv, pf3_ew))), type = "latex")

# #### Portfolio 4 #### #
sh4 <- SharpeRatio.annualized(xts(km.spline2, order.by = kompas100p$date[-1])) %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "symbol",
               values_to = "sharpe") %>%
  arrange(desc(sharpe)) %>% head(10)

pf4 <- km.spline2 %>%
  rownames_to_column(var = "date") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "return") %>%
  arrange(symbol) %>%
  filter(symbol %in% sh4$symbol) %>%
  pivot_wider(names_from = symbol, values_from = return)

# --- PF4 MV weighting --- #
pf4 <- xts(pf4[-1], order.by = as.Date(kompas100p$date[-1]))
pf4_mv <- portfolio.optim(pf4)$pw
names(pf4_mv) <- colnames(as.data.frame(pf4))
print(xtable(as.data.frame(t(pf4_mv))), type = "latex")
# Calculate MV Portfolio Return
pf4_mv_ret <- Return.portfolio(pf4, weights = pf4_mv)
pf4_mv_ret_ <- Return.portfolio(pf4, weights = pf4_mv, 
                                wealth.index = TRUE)

# --- PF4 EW Weighting --- #
pf4_ew <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
names(pf4_ew) <- colnames(as.data.frame(pf4))
print(xtable(as.data.frame(t(pf4_ew))), type = "latex")
# Calculate MV Portfolio Return
pf4_ew_ret <- Return.portfolio(pf4, weights = pf4_ew)
pf4_ew_ret_ <- Return.portfolio(pf4, weights = pf4_ew, 
                                wealth.index = TRUE)

# Plot - PF4 Wealth Index
plot(pf4_mv_ret_, col = "#32936F", main = "Portfolio Wealth Index - PF4", lwd = 2)
lines(pf4_ew_ret_, col = "#E83F6F" , lwd= 2)
lines(k100index_pf_, col = "#000000", lwd = 2)
addLegend("topleft", c("PF4 MV", "PF4 EW", "KOMPAS100"),
          col = c("#32936F", "#E83F6F", "#000000"), lty = 1, bty="o", lwd = 2)

# Latex Table
print(xtable(as.data.frame(cbind(pf4_mv, pf4_ew))), type = "latex")


# #### Portfolio 5 #### #
sh5 <- SharpeRatio.annualized(xts(hc.sc1, order.by = kompas100p$date[-1])) %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "symbol",
               values_to = "sharpe") %>%
  arrange(desc(sharpe)) %>% head(10)

pf5 <- hc.sc1 %>%
  rownames_to_column(var = "date") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "return") %>%
  arrange(symbol) %>%
  filter(symbol %in% sh5$symbol) %>%
  pivot_wider(names_from = symbol, values_from = return)

# --- PF5 MV weighting --- #
pf5 <- xts(pf5[-1], order.by = as.Date(kompas100p$date[-1]))
pf5_mv <- portfolio.optim(pf5)$pw
names(pf5_mv) <- colnames(as.data.frame(pf5))
print(xtable(as.data.frame(t(pf5_mv))), type = "latex")
# Calculate MV Portfolio Return
pf5_mv_ret <- Return.portfolio(pf5, weights = pf5_mv)
pf5_mv_ret_ <- Return.portfolio(pf5, weights = pf5_mv, 
                                wealth.index = TRUE)

# --- PF5 EW Weighting --- #
pf5_ew <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
names(pf5_ew) <- colnames(as.data.frame(pf5))
print(xtable(as.data.frame(t(pf5_ew))), type = "latex")
# Calculate MV Portfolio Return
pf5_ew_ret <- Return.portfolio(pf5, weights = pf5_ew)
pf5_ew_ret_ <- Return.portfolio(pf5, weights = pf5_ew, 
                                wealth.index = TRUE)

# Plot - PF5 Wealth Index
plot(pf5_mv_ret_, col = "#EF233C", main = "Portfolio Wealth Index - PF5", lwd = 2)
lines(pf5_ew_ret_, col = "#8D99AE" , lwd= 2)
lines(k100index_pf_, col = "#000000", lwd = 2)
addLegend("topright", c("PF5 MV", "PF5 EW", "KOMPAS100"),
          col = c("#EF233C", "#8D99AE", "#000000"), lty = 1, bty="o", lwd = 2)

# Latex Table
print(xtable(as.data.frame(cbind(pf5_mv, pf5_ew))), type = "latex")

# #### Portfolio 6 #### # (Same as Portfolio 1)
sh6 <- SharpeRatio.annualized(xts(hc.sc2, order.by = kompas100p$date[-1])) %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "symbol",
               values_to = "sharpe") %>%
  arrange(desc(sharpe)) %>% head(10)

pf6 <- hc.sc2 %>%
  rownames_to_column(var = "date") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "return") %>%
  arrange(symbol) %>%
  filter(symbol %in% sh6$symbol) %>%
  pivot_wider(names_from = symbol, values_from = return)

# --- PF6 MV weighting --- #
pf6 <- xts(pf6[-1], order.by = as.Date(kompas100p$date[-1]))
pf6_mv <- portfolio.optim(pf6)$pw
names(pf6_mv) <- colnames(as.data.frame(pf6))
print(xtable(as.data.frame(t(pf6_mv))), type = "latex")
# Calculate MV Portfolio Return
pf6_mv_ret <- Return.portfolio(pf6, weights = pf6_mv)
pf6_mv_ret_ <- Return.portfolio(pf6, weights = pf6_mv, 
                                wealth.index = TRUE)
# Compare MV with KOMPAS100 Index 
lines(pf6_mv_ret_, col = "#FF22BB", legend.loc = "topleft")
lines(k100index_pf_, col = "#555500", main = "KOMPAS100 Index")

# --- PF6 EW Weighting --- #
pf6_ew <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
names(pf6_ew) <- colnames(as.data.frame(pf6))
print(xtable(as.data.frame(t(pf6_ew))), type = "latex")
# Calculate MV Portfolio Return
pf6_ew_ret <- Return.portfolio(pf6, weights = pf6_ew)
pf6_ew_ret_ <- Return.portfolio(pf6, weights = pf6_ew, 
                                wealth.index = TRUE)
# Plot - PF6 Wealth Index
plot(pf6_mv_ret_, col = "#FF22BB", main = "Portfolio Wealth Index - PF6", lwd = 2)
lines(pf6_ew_ret_, col = "#555500" , lwd= 2)
lines(k100index_pf_, col = "#000000", lwd = 2)
addLegend("topleft", c("PF6 MV", "PF6 EW", "KOMPAS100"),
          col = c("#FF22BB", "#555500", "#000000"), lty = 1, bty="o", lwd = 2)

# Latex Table
print(xtable(as.data.frame(cbind(pf5_mv, pf5_ew))), type = "latex")


# #### Portfolio 7 #### #
sh7 <- SharpeRatio.annualized(xts(hc.spline1, order.by = kompas100p$date[-1])) %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "symbol",
               values_to = "sharpe") %>%
  arrange(desc(sharpe)) %>% head(10)

pf7 <- hc.spline1 %>%
  rownames_to_column(var = "date") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "return") %>%
  arrange(symbol) %>%
  filter(symbol %in% sh7$symbol) %>%
  pivot_wider(names_from = symbol, values_from = return)

# --- PF7 MV weighting --- #
pf7 <- xts(pf7[-1], order.by = as.Date(kompas100p$date[-1]))
pf7_mv <- portfolio.optim(pf7)$pw
names(pf7_mv) <- colnames(as.data.frame(pf7))
print(xtable(as.data.frame(t(pf7_mv))), type = "latex")
# Calculate MV Portfolio Return
pf7_mv_ret <- Return.portfolio(pf7, weights = pf7_mv)
pf7_mv_ret_ <- Return.portfolio(pf7, weights = pf7_mv, 
                                wealth.index = TRUE)
# Compare MV with KOMPAS100 Index 
lines(pf7_mv_ret_, col = "#222200", legend.loc = "topleft")
lines(k100index_pf_, col = "#555500", main = "KOMPAS100 Index")

# --- PF7 EW Weighting --- #
pf7_ew <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
names(pf7_ew) <- colnames(as.data.frame(pf7))
print(xtable(as.data.frame(t(pf7_ew))), type = "latex")
# Calculate MV Portfolio Return
pf7_ew_ret <- Return.portfolio(pf7, weights = pf7_ew)
pf7_ew_ret_ <- Return.portfolio(pf7, weights = pf7_ew, 
                                wealth.index = TRUE)
# Plot - PF7 Wealth Index
plot(pf7_mv_ret_, col = "#87BBA2", main = "Portfolio Wealth Index - PF7", lwd = 2)
lines(pf7_ew_ret_, col = "#F0A868" , lwd= 2)
lines(k100index_pf_, col = "#000000", lwd = 2)
addLegend("topleft", c("PF7 MV", "PF7 EW", "KOMPAS100"),
          col = c("#87BBA2", "#F0A868", "#000000"), lty = 1, bty="o", lwd = 2)

# Latex Table
print(xtable(as.data.frame(cbind(pf7_mv, pf7_ew))), type = "latex")

# #### Portfolio 8 #### #
sh8 <- SharpeRatio.annualized(xts(hc.spline2, order.by = kompas100p$date[-1])) %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "symbol",
               values_to = "sharpe") %>%
  arrange(desc(sharpe)) %>% head(10)

pf8 <- hc.spline2 %>%
  rownames_to_column(var = "date") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "symbol",
               values_to = "return") %>%
  arrange(symbol) %>%
  filter(symbol %in% sh8$symbol) %>%
  pivot_wider(names_from = symbol, values_from = return)

# --- PF8 MV weighting --- #
pf8 <- xts(pf8[-1], order.by = as.Date(kompas100p$date[-1]))
pf8_mv <- portfolio.optim(pf8)$pw
names(pf8_mv) <- colnames(as.data.frame(pf8))
print(xtable(as.data.frame(t(pf8_mv))), type = "latex")
# Calculate MV Portfolio Return
pf8_mv_ret <- Return.portfolio(pf8, weights = pf8_mv)
pf8_mv_ret_ <- Return.portfolio(pf8, weights = pf8_mv, 
                                wealth.index = TRUE)
# Compare MV with KOMPAS100 Index 
lines(pf8_mv_ret_, col = "#00BcFF", legend.loc = "topleft")
lines(k100index_pf_, col = "#555500", main = "KOMPAS100 Index")

# --- PF8 EW Weighting --- #
pf8_ew <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
names(pf8_ew) <- colnames(as.data.frame(pf8))
print(xtable(as.data.frame(t(pf8_ew))), type = "latex")
# Calculate MV Portfolio Return
pf8_ew_ret <- Return.portfolio(pf8, weights = pf8_ew)
pf8_ew_ret_ <- Return.portfolio(pf8, weights = pf8_ew, 
                                wealth.index = TRUE)
# Plot - PF7 Wealth Index
plot(pf8_mv_ret_, col = "#F34213", main = "Portfolio Wealth Index - PF8", lwd = 2)
lines(pf8_ew_ret_, col = "#48A9A6" , lwd= 2)
lines(k100index_pf_, col = "#000000", lwd = 2)
addLegend("topleft", c("PF8 MV", "PF8 EW", "KOMPAS100"),
          col = c("#F34213", "#48A9A6", "#000000"), lty = 1, bty="o", lwd = 2)

# Latex Table
print(xtable(as.data.frame(cbind(pf8_mv, pf8_ew))), type = "latex")

## ### Performance Comparison ### ##
# Create function for performance table #
# (https://www.youtube.com/watch?v=2EpIdnA0pPo) #
perf_table <- function(pfreturn, indexreturn, rfrate) {
  perf_table <- Return.cumulative(pfreturn)
  perf_table <- rbind(perf_table, Return.annualized(pfreturn, scale = 252))
  perf_table <- rbind(perf_table, StdDev.annualized(pfreturn, scale = 252))
  perf_table <- rbind(perf_table, SharpeRatio.annualized(pfreturn, Rf = rfrate))
  perf_table <- rbind(perf_table, maxDrawdown(pfreturn))
  perf_table <- rbind(perf_table, CAPM.beta(pfreturn, indexreturn, Rf = rfrate))
  rownames(perf_table) <- c("TotalReturn", "Return_Ann", "Volatility_Ann", "SharpeRatio_Ann", "MaxDD", "Beta")
  perf_table <- t(perf_table)
  return(perf_table)
}

# Performance Table of All Portfolios
pf1_mv_perf <- perf_table(pf1_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF1 MV"))

pf1_ew_perf <- perf_table(pf1_ew_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF1 EW"))

pf2_mv_perf <- perf_table(pf2_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF2 MV"))

pf2_ew_perf <- perf_table(pf2_ew_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF2 EW"))

pf3_mv_perf <- perf_table(pf3_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF3 MV"))

pf3_ew_perf <- perf_table(pf3_ew_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF3 EW"))

pf4_mv_perf <- perf_table(pf4_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF4 MV"))

pf4_ew_perf <- perf_table(pf4_ew_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF4 EW"))

pf5_mv_perf <- perf_table(pf5_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF5 MV"))

pf5_ew_perf <- perf_table(pf5_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF5 EW"))

pf6_mv_perf <- perf_table(pf6_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF6 MV"))

pf6_ew_perf <- perf_table(pf6_ew_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF6 EW"))

pf7_mv_perf <- perf_table(pf7_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF7 MV"))

pf7_ew_perf <- perf_table(pf7_ew_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF7 EW"))

pf8_mv_perf <- perf_table(pf8_mv_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF8 MV"))

pf8_ew_perf <- perf_table(pf8_ew_ret, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "PF8 EW"))

k100index_perf <- perf_table(k100index_pf, k100index_pf, rfrate = 0) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Portfolio") %>%
  mutate(Portfolio = str_replace(Portfolio, "portfolio.returns", "KOMPAS100"))

AllPerf <- rbind(pf1_mv_perf,pf1_ew_perf,pf2_mv_perf,pf2_ew_perf,pf3_mv_perf,
                 pf3_ew_perf,pf4_mv_perf,pf4_ew_perf,pf5_mv_perf,pf5_ew_perf,
                 pf6_mv_perf,pf6_ew_perf,pf7_mv_perf,pf7_ew_perf,pf8_mv_perf,
                 pf8_ew_perf,k100index_perf)
print(xtable(AllPerf), type = "latex")

pf1perf <- rbind(pf1_mv_perf, pf1_ew_perf, k100index_perf)
print(xtable(pf1perf), type = "latex")

pf2perf <- rbind(pf2_mv_perf, pf2_ew_perf, k100index_perf)
print(xtable(pf2perf), type = "latex")

pf3perf <- rbind(pf3_mv_perf, pf3_ew_perf, k100index_perf)
print(xtable(pf3perf), type = "latex")

pf4perf <- rbind(pf4_mv_perf, pf4_ew_perf, k100index_perf)
print(xtable(pf4perf), type = "latex")

pf5perf <- rbind(pf5_mv_perf, pf5_ew_perf, k100index_perf)
print(xtable(pf5perf), type = "latex")

pf6perf <- rbind(pf6_mv_perf, pf6_ew_perf, k100index_perf)
print(xtable(pf6perf), type = "latex")

pf7perf <- rbind(pf7_mv_perf, pf7_ew_perf, k100index_perf)
print(xtable(pf7perf), type = "latex")

pf8perf <- rbind(pf8_mv_perf, pf8_ew_perf, k100index_perf)
print(xtable(pf8perf), type = "latex")
