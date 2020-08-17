rm(list = ls())
library(rgdal) # For reading in ESRI shapefile.
library(ggplot2)
library(raster)
library(gridExtra) # For grid.arrange.
library(lemon) # For grid_arrange_shared_legend.
library(MASS) # Ordinal regression.
library(knitr)
# library(mlogit) # Multinomial Probit regression.
library(dplyr)
library(viridis)
library(LaplacesDemon)
library(tidyr)

PGA_list <- c(0.00000001, seq(0.0025, 1, by = 0.0025))
plot_pga <- data.frame("PGA" = PGA_list)
plot_pga$logPGA <- log(plot_pga$PGA)

excp_to_prob <- function(x){
  temp <- rep(NA, length(x) +1)
  temp[1] <- 1 - x$Mean[x$damage_grade == "Grade 0"]
  temp[2] <- x$Mean[x$damage_grade == "Grade 0"] - x$Mean[x$damage_grade == "Grade 1"]
  temp[3] <- x$Mean[x$damage_grade == "Grade 1"] - x$Mean[x$damage_grade == "Grade 2"]
  temp[4] <- x$Mean[x$damage_grade == "Grade 2"] - x$Mean[x$damage_grade == "Grade 3"]
  temp[5] <- x$Mean[x$damage_grade == "Grade 3"] - x$Mean[x$damage_grade == "Grade 4"]
  temp[6] <- x$Mean[x$damage_grade == "Grade 4"]
  return(temp)
}

## Dataset 1: Training/Base Data ##

# Read in Nepal buildings category data:

ordinal_data <- read.csv(file = "D:/Documents/Ordinal_Fragility_Curves/Data/ordinal_data.csv", 
                         stringsAsFactors = FALSE)

# Only use timber superstructure:

timber_data <- ordinal_data[ordinal_data$superstructure == "timber", ]

# Fit ordinal regression:
timber_data$damage_grade <- ordered(timber_data$damage_grade, 
                                    levels = c("Grade 0", "Grade 1", 
                                               "Grade 2", "Grade 3",
                                               "Grade 4", "Grade 5"))
timber_ordinal <- polr(damage_grade ~ logPGA, data = timber_data, 
                       method = "probit", Hess = TRUE)

# Generate Dataset 1 via predicted probabilities:

no.buildings <- 100 
# (Maximum) number of buildings per PGA bin; can delete rows later if want diff. no. per PGA bin.
pga_range <- data.frame("PGA" = seq(0.1, 0.9, by = 0.3))
pga_range$logPGA <- log(pga_range$PGA)
ds1_prob <- predict(timber_ordinal, newdata = pga_range, type = "prob")

set.seed(1)
ds1_count <- apply(ds1_prob, MARGIN = 1, FUN = function(x){sample(c("Grade 0", "Grade 1", 
                                                                    "Grade 2", "Grade 3",
                                                                    "Grade 4", "Grade 5"), 
                                                                  size = no.buildings, 
                                                                  prob = x, replace = TRUE)})

dataset_1 <- data.frame("PGA" = rep(pga_range$PGA, each = no.buildings), 
                        "logPGA" = rep(pga_range$logPGA, each = no.buildings), 
                        "damage_grade" = c(ds1_count[, 1], ds1_count[, 2], ds1_count[, 3]))
table(dataset_1$damage_grade)
dataset_1$Data <- "Dataset 1"

mega_KL_table <-data.frame("Model" = NA, "PGA_0.1" = NA, "PGA_0.4" = NA,
                          "PGA_0.7" = NA, "Average" = NA, "x" = NA)

x_vec <- seq(0, 0.5, by = 0.05)

for (xi in 1:length(x_vec)){

## Dataset 3 (Training set): For illustrating overfitting ##
  
# Generate Dataset 3 as Dataset 1 but x% chosen randomly from a 
# uniform distribution across the damage states:
  
  
temp_dataset_3 <- dataset_1
x <- x_vec[xi]
n <- nrow(temp_dataset_3)
temp_dataset_3$damage_grade <- as.character(temp_dataset_3$damage_grade)

# Choose x% of rows to change damage states for:
set.seed(2)
change_id <- sample(n, x*n, replace = FALSE)
set.seed(3)
# Exclude Grade 0 for now.
change_val <- sample(c("Grade 1", "Grade 2", "Grade 3", "Grade 4", "Grade 5"),
                     size = x*n, replace = TRUE)

temp_dataset_3[change_id, "damage_grade"] <- change_val
dataset_3 <- temp_dataset_3

table(dataset_3$damage_grade)

# Fit fragility curves to Dataset 3 using separate probit regressions:

ds3_ex0 <- dataset_3
ds3_ex0$Damage <- 1
ds3_ex0$Damage[ds3_ex0$damage_grade %in% c("Grade 0")] <- 0
ds3_probit0 <- glm(Damage ~ logPGA, family = binomial(link = "probit"), data = ds3_ex0)
# summary(ds3_probit0)

ds3_ex1 <- dataset_3
ds3_ex1$Damage <- 1
ds3_ex1$Damage[ds3_ex1$damage_grade %in% c("Grade 0", "Grade 1")] <- 0
ds3_probit1 <- glm(Damage ~ logPGA, family = binomial(link = "probit"), data = ds3_ex1)
# summary(ds3_probit1)

ds3_ex2 <- dataset_3
ds3_ex2$Damage <- 1
ds3_ex2$Damage[ds3_ex2$damage_grade %in% c("Grade 0", "Grade 1", "Grade 2")] <- 0
ds3_probit2 <- glm(Damage ~ logPGA, family = binomial(link = "probit"), data = ds3_ex2)
# summary(ds3_probit2)

ds3_ex3 <- dataset_3
ds3_ex3$Damage <- 1
ds3_ex3$Damage[ds3_ex3$damage_grade %in% c("Grade 0", "Grade 1", "Grade 2", "Grade 3")] <- 0
ds3_probit3 <- glm(Damage ~ logPGA, family = binomial(link = "probit"), data = ds3_ex3)
# summary(ds3_probit3)

ds3_ex4 <- dataset_3
ds3_ex4$Damage <- 1
ds3_ex4$Damage[ds3_ex4$damage_grade %in% c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4")] <- 0
ds3_probit4 <- glm(Damage ~ logPGA, family = binomial(link = "probit"), data = ds3_ex4)
# summary(ds3_probit4)

ex0_lp_mean <- predict.glm(ds3_probit0, newdata = plot_pga, type = "link")
frag_ex0 <- pnorm(ex0_lp_mean, lower.tail = TRUE)

ex1_lp_mean <- predict.glm(ds3_probit1, newdata = plot_pga, type = "link")
frag_ex1 <- pnorm(ex1_lp_mean, lower.tail = TRUE)

ex2_lp_mean <- predict.glm(ds3_probit2, newdata = plot_pga, type = "link")
frag_ex2 <- pnorm(ex2_lp_mean, lower.tail = TRUE)

ex3_lp_mean <- predict.glm(ds3_probit3, newdata = plot_pga, type = "link")
frag_ex3 <- pnorm(ex3_lp_mean, lower.tail = TRUE)

ex4_lp_mean <- predict.glm(ds3_probit4, newdata = plot_pga, type = "link")
frag_ex4 <- pnorm(ex4_lp_mean, lower.tail = TRUE)

frag_probit3 <- data.frame("PGA" = rep(plot_pga$PGA, 5), "Mean" = c(frag_ex0, frag_ex1, frag_ex2, frag_ex3, frag_ex4),
                           "damage_grade" = rep(c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4"), each = nrow(plot_pga)))
                           
# Fit ordinal regression to Dataset 3.


dataset_3$damage_grade <- ordered(dataset_3$damage_grade, 
                                  levels = c("Grade 0", "Grade 1", 
                                             "Grade 2", "Grade 3",
                                             "Grade 4", "Grade 5"))
ds3_ordinal <- polr(damage_grade ~ logPGA, data = dataset_3, 
                    method = "probit", Hess = TRUE)

summary(ds3_ordinal)

ds3_lp_pred <- expand.grid("logPGA" = log(PGA_list),
                           "damage_grade" = c("Grade 0", "Grade 1", "Grade 2",
                                              "Grade 3","Grade 4"))
# Exclude Grade 5 because fragility curve for probability of exceedance.
ds3_lp_pred$b0_id <- as.numeric(ds3_lp_pred$damage_grade)
ds3_lp_pred$b0 <- ds3_ordinal$zeta[ds3_lp_pred$b0_id]
ds3_ordinal_lp_mean <- ds3_lp_pred$b0 - ds3_ordinal$coefficients*ds3_lp_pred$logPGA
ds3_ordinal_mean <- pnorm(ds3_ordinal_lp_mean, lower.tail = TRUE) 
# Currently, 1- exceedance probability.
ds3_ordinal_mean <- 1 - ds3_ordinal_mean

ds3_ordinal_frag <- data.frame("PGA" = rep(PGA_list, 5),
                               "damage_grade" = rep(c("Grade 0", "Grade 1", "Grade 2",
                                                      "Grade 3","Grade 4"), each = length(PGA_list)),
                               "Mean" = ds3_ordinal_mean)

wide_dataset_1 <- unique(dataset_1[, c("PGA", "logPGA", "damage_grade", "Data")])
wide_dataset_1$Count <- NA
for (k in 1:nrow(wide_dataset_1)){
  wide_dataset_1$Count[k] <- sum(dataset_1$PGA == wide_dataset_1$PGA[k] & dataset_1$damage_grade == wide_dataset_1$damage_grade[k])
}

# Plot stacked barcharts of Dataset 1 with fragility curves
# from ordinal and probit regressions of Dataset 3 (Training data) on top.
# plot5 <- ggplot(data = wide_dataset_1) +
#   geom_bar(aes(fill = damage_grade, y = Count, x = PGA), width = 0.1, alpha = 0.3, position = "fill", stat = "identity") + geom_line(data = ds3_ordinal_frag[ds3_ordinal_frag$damage_grade != "Grade 0", ], aes(x = PGA, y = Mean, color = damage_grade), stat = "identity") + theme_classic() + ggtitle("(a) Ordinal regression") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_viridis_d(begin = 0.2, end = 1) + labs(y = "Probability of exceedance", fill = "Damage grade") + scale_fill_viridis_d(begin = 0, end = 1) + guides(color = FALSE)
# plot6 <- ggplot(data = wide_dataset_1) + geom_bar(aes(fill = damage_grade, y = Count, x = PGA), width = 0.1, alpha = 0.3, position = "fill", stat = "identity") + geom_line(data = frag_probit3[frag_probit3$damage_grade != "Grade 0", ], aes(x = PGA, y = Mean, color = damage_grade), stat = "identity") + theme_classic() + ggtitle("(b) Probit regression") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_viridis_d(begin = 0.2, end = 1) + labs(y = "Probability of exceedance", fill = "Damage grade") + scale_fill_viridis_d(begin = 0, end = 1) + guides(color = FALSE)
# grid_arrange_shared_legend(plot5, plot6, nrow = 1, ncol = 2)



#K-L divergence measure:
KL_probit <- rep(NA, nrow(pga_range))
KL_ordinal <- rep(NA, nrow(pga_range))

no.buildings.test <- no.buildings
test_KL_probit <- matrix(NA, nrow = 100, ncol = nrow(pga_range))
test_KL_ordinal <- matrix(NA, nrow = 100, ncol = nrow(pga_range))

for (l in 1:100){
  
  # Generate a test set:
  set.seed(l+10) 
  test_count <- apply(ds1_prob, MARGIN = 1, FUN = function(x){sample(c("Grade 0", "Grade 1", 
                                                                       "Grade 2", "Grade 3",
                                                                       "Grade 4", "Grade 5"), 
                                                                     size = no.buildings.test, 
                                                                     prob = x, replace = TRUE)})
  
  test_set <- data.frame("PGA" = rep(pga_range$PGA, each = no.buildings.test), 
                         "logPGA" = rep(pga_range$logPGA, each = no.buildings.test), 
                         "damage_grade" = c(test_count[, 1], test_count[, 2], test_count[, 3]))
  test_set$Data <- "Test set"
  
  for (i in 1:nrow(pga_range)){
    
    probit_prob <- matrix(NA, ncol = 6, nrow = nrow(pga_range))
    
    probit_excp <- frag_probit3[frag_probit3$PGA == pga_range$PGA[i], ]
      
    # Use Porter method I for comparing probabilities (use higher exceedance probabilities at crossing):
    for(j in nrow(probit_excp):2){
      if(probit_excp$Mean[j]>probit_excp$Mean[j-1]){
        probit_excp$Mean[j-1] <- probit_excp$Mean[j]
      } 
    }
      
    probit_prob[i, ] <- excp_to_prob(probit_excp)
    
    ordinal_excp <- ds3_ordinal_frag[ds3_ordinal_frag$PGA == pga_range$PGA[i], ]
    probit_prob <- excp_to_prob(probit_excp)
    names(probit_prob) <- c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4", "Grade 5")
    ordinal_prob <- excp_to_prob(ordinal_excp)
    names(ordinal_prob) <- c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4", "Grade 5")
    
    
    test_subset <- test_set[test_set$PGA == pga_range$PGA[i], ]
    test_prob <- c(sum(test_subset$damage_grade == "Grade 0"), 
                   sum(test_subset$damage_grade == "Grade 1"), 
                   sum(test_subset$damage_grade == "Grade 2"), 
                   sum(test_subset$damage_grade == "Grade 3"), 
                   sum(test_subset$damage_grade == "Grade 4"), 
                   sum(test_subset$damage_grade == "Grade 5"))
    test_prob <- test_prob/sum(test_prob)
    
    m_probit_prob <- (test_prob + probit_prob)/2
    m_ordinal_prob <- (test_prob + ordinal_prob)/2
    
    valid_id <- which(test_prob>0 & probit_prob>0)
    valid_id_2 <- which(test_prob>0 & ordinal_prob>0)
    
    test_KL_probit[l, i] <- 
      
      sum(test_prob[valid_id]*log2(test_prob[valid_id]/probit_prob[valid_id]))
      
      # sum(test_prob[valid_id]*log2(test_prob[valid_id]/m_probit_prob[valid_id]))/2 + sum(probit_prob[valid_id]*log2(probit_prob[valid_id]/m_probit_prob[valid_id]))/2
    
    test_KL_ordinal[l, i] <- 
      
      sum(test_prob[valid_id_2 ]*log2(test_prob[valid_id_2 ]/ordinal_prob[valid_id_2]))
      
    }
}


KL_table <- data.frame("Model" = c("Probit", "Ordinal"), 
                       "PGA_0.1" = round(c(mean(test_KL_probit[, 1]), mean(test_KL_ordinal[, 1])), digits = 3), 
                       "PGA_0.4" = round(c(mean(test_KL_probit[, 2]), mean(test_KL_ordinal[, 2])), digits = 3), 
                       "PGA_0.7" = round(c(mean(test_KL_probit[, 3]), mean(test_KL_ordinal[, 3])), digits = 3), 
                       "Average" = round(c(mean(test_KL_probit), mean(test_KL_ordinal)), digits = 3))
KL_table$x <- x

mega_KL_table <- rbind(mega_KL_table, KL_table)

print(paste("Round ", xi, "/", length(x_vec), " done.", sep = ""))

}

mega_KL_table <- mega_KL_table[-1, ]
mega_KL_table$Average <- rowMeans(mega_KL_table[, c("PGA_0.1", "PGA_0.4", "PGA_0.7")])
KL_table_wide <- gather(data = mega_KL_table, key = "PGA", value = "KL", PGA_0.1, PGA_0.4, PGA_0.7, Average)
KL_table_wide$PGA[KL_table_wide$PGA != "Average"] <- substr(KL_table_wide$PGA[KL_table_wide$PGA != "Average"], start = 5, stop = 7)
KL_table_wide$xplot <- KL_table_wide$x*100

# New facet labels for PGA:
PGA.labs <- c("PGA = 0.1", "PGA = 0.4", "PGA = 0.7", "Average")
names(PGA.labs) <- c(0.1, 0.4, 0.7, "Average")

jpeg(file = "D:/Documents/Ordinal_Fragility_Curves/Graphics/fig3.jpeg", height = 1500, width = 2000, res = 300)
ggplot() + geom_line(data = KL_table_wide, aes(x = xplot, y = KL, color = Model)) + facet_wrap(~PGA, labeller = labeller(PGA = PGA.labs)) + theme_bw() + labs(x = "x") + theme(strip.background =element_rect(fill="white"))
dev.off()
