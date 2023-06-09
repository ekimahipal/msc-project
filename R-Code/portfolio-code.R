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
       patchwork, gridExtra, ggplot2, xtable, ggpubr)


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
  column_to_rownames(var = "date")

km.spline1 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_spline_sc$symbol[cl_spline_sc$cluster==2]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date")

km.spline2 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% cl_spline_sc$symbol[cl_spline_sc$cluster==1]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date")

hc.sc1 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_raw_sc$symbol[hcl_raw_sc$cluster==1]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date")

hc.sc2 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_raw_sc$symbol[hcl_raw_sc$cluster==2]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date")

hc.spline1 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_spline_sc$symbol[hcl_spline_sc$cluster==1]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date")

hc.spline2 <- kompas100p %>%
  pivot_longer(cols = 2:ncol(kompas100p),
               names_to = "symbol",
               values_to = "price") %>%
  arrange(symbol, date) %>% 
  filter(symbol %in% hcl_spline_sc$symbol[hcl_spline_sc$cluster==2]) %>%
  pivot_wider(id_cols = date, names_from = symbol, values_from = price) %>%
  column_to_rownames(var = "date")

# Calculate Return from each cluster
ew <- 
Return.portfolio(R = km.sc1, weights = )