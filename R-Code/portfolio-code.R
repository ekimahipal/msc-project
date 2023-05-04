# R Script - Exploring the Power of Clustering for Portfolio Selection Using Spline Coefficients
# Eki Mahipal 
# Student ID: 20377383
# psxem8@nottingham.ac.uk

# Install and Load Packages
install.packages("pacman")
library(pacman)


# 1. DATA COLLECTIOM AMD PREPROCESSING ------------------------------------
# 1a. Load Packages
p_load(tidyverse, tidyquant)

# 1b. Load data from csv
datalink <- ""
kompasdata <- read.csv('data/kompas100list.csv', sep=';')
kompasdata$Kode <- paste(kompasdata$Kode, '.JK', sep='') # All of Indonesian stocks symbols end with ".JK"

# 1c. Get the stock prices from Yahoo! Finance
startdate <- as.Date('2012-01-01')
enddate <- as.Date('2022-12-31')
kompas <- as.data.frame(tq_get(kompasdata$Kode,
                 from = startdate,
                 to = enddate))
kompas <- kompas %>%
  select(date, symbol, close)

# 1d. Check starting date 2012-01-02 
## Function to check for NA 
find_NA <- function (x) {
  x[which(is.na(x$close)),]
}

## Find stock that has NA on 2012-01-02
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

## Find NA and treat with interpolation
find_NA(kompas100)
kompas100$close <- na.approx(kompas100$close)

# 1e. Visualisation of Preprocessed KOMPAS100 Index daily prices
k100plot1 <- kompas100 %>%
  ggplot(aes(x = date, y = close)) +
  geom_line(aes(color = symbol), show.legend = FALSE) +
  xlab("Year") +
  ylab("Price") +
  ggtitle("KOMPAS100 Index Price") +
  theme_light()

ggsave(k100plot1, file = "image/Plot 1 - KOMPAS100 Index Prices 2012-2022.png", 
       width = 1600, height = 1200, units = "px", dpi = 150)
