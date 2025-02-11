setwd("G:/Meine Ablage/Master/analysis/vcfs/pheno/")

library(ggplot2)

# data from Filiz et al., 2009

#group6 <- c("BdTR10N",	"BdTR11C",	"BdTR11D",	"BdTR11F",	"BdTR11G",	"BdTR11H",	"BdTR11I")
group1 <- c("BdTR1A",	"BdTR1B",	"BdTR1E",	"BdTR1F",	"BdTR1G",	"BdTR1H",	"BdTR1J",	"BdTR1K",	
            "BdTR1M",	"BdTR1N",	"BdTR2B",	"BdTR2C",	"BdTR2D",	"BdTR2G",	"BdTR2H",	"BdTR2J",
            "BdTR2K",	"BdTR2M",	"BdTR2N",	"BdTR2P",	"BdTR2R",	"BdTR2S")

group1.1 <- c("BdTR1E",	"BdTR1H",	"BdTR1J",	"BdTR1K",	
              "BdTR1M",	"BdTR1N") #only central

group1.2 <- c("BdTR2B",	"BdTR2C",	"BdTR2D",	"BdTR2G",	"BdTR2H",	"BdTR2J",
              "BdTR2K",	"BdTR2M",	"BdTR2N",	"BdTR2P",	"BdTR2R",	"BdTR2S")


Average_height_cm = c("30.0–41.0", "29.5–40.5", "28.7–39.8", "29.0–39.5", "28.2–41.0", "29.5–40.8", "29.0–39.9", "28.9–40.0", "29.1–41.0", "29.9–40.7", "28.9–40.5", "29.0–40.0", "29.1–39.9", "29.6–40.1", "30.0–41.0", "29.7–41.8", "29.2–40.2", "28.5–39.0", "28.0–40.8", "30.9–40.0", "30.0–41.7", "29.4–40.1", "29.0–39.2", "29.9–42.2", "30.0–41.2", "30.0–42.4", "31.1–39.5", "30.0–39.0", "30.7–41.1", "29.9–40.0")

`Germination percentage` = c(94, 94, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
data <- data.frame(
  Genotype = c("BdTR1A", "BdTR1B", "BdTR1C", "BdTR1D", "BdTR1E", "BdTR1F", "BdTR1G", "BdTR1H", "BdTR1I", "BdTR1J", "BdTR1K", "BdTR1M", "BdTR1N", "BdTR2A", "BdTR2B", "BdTR2C", "BdTR2D", "BdTR2E", "BdTR2F", "BdTR2G", "BdTR2H", "BdTR2I", "BdTR2J", "BdTR2K", "BdTR2M", "BdTR2N", "BdTR2O", "BdTR2P", "BdTR2R", "BdTR2S"),
  Longitude = c("27828’36.81@E", "26855’53.29@E", "2882’24.71@E", "28835’6.42@E", "2881’52.75@E", "27828’36.81@E", "26855’53.29@E", "2882’24.71@E", "28834’59.02@E", "28835’6.75@E", "29840’40.96@E", "30815’8.93@E", "30814’44.11@E", "30847’19.07@E", "31819’52.01@E", "31853’29.40@E", "31818’33.71@E", "32824’16.47@E", "31853’5.68@E", "32859’7.32@E", "32825’56.46@E", "33832’16.37@E", "32859’17.24@E", "3484’18.34@E", "3485’40.38@E", "33831’10.58@E", "32825’56.46@E", "33832’16.37@E", "32859’17.24@E", "3484’18.34@E"),
  Elevation = c(124, 141, 363, 612, 986, 124, 141, 363, 841, 513, 1007, 1076, 1034, 932, 667, 864, 1301, 1012, 1288, 1596, 787, 872, 1192, 1142, 1406, 1013, 787, 872, 1192, 1142),
  #  `Leaf color` = as.factor(c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)),
  #  `Leaf hairiness` = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
  `Seed Width` = c(0.123, 0.123, 0.120, 0.118, 0.124, 0.124, 0.121, 0.124, 0.121, 0.114, 0.127, 0.121, 0.119, 0.112, 0.112, 0.103, 0.114, 0.122, 0.108, 0.115, 0.119, 0.119, 0.112, 0.111, 0.113, 0.117, 0.119, 0.110, 0.115, 0.116),
  `Seed Length` = c(0.644, 0.665, 0.661, 0.622, 0.666, 0.629, 0.616, 0.636, 0.606, 0.654, 0.652, 0.673, 0.641, 0.593, 0.585, 0.591, 0.586, 0.611, 0.589, 0.609, 0.611, 0.617, 0.586, 0.597, 0.572, 0.610, 0.580, 0.607, 0.611, 0.561),
  `Seed production (weeks)` = c(8, 8, 8, 8, 8, 8, 9, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11),
  `Average biomass (g)` = c(0.89, 0.88, 0.89, 0.90, 0.90, 0.90, 0.87, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.79, 0.78, 0.78, 79.8, 0.80, 0.78, 0.78, 0.78, 0.78, 0.79, 0.80, 0.80, 0.80, 0.80, 0.90, 0.79, 0.80),
  `Seed yield` = c(64, 63, 63, 65, 63, 64, 64, 64, 64, 63, 62, 63, 63, 183, 183, 183, 186, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184),
  `DNA content (pg/2C)` = c(0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.695, 0.695, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696, 0.696),
  `Germination percentage` = c(94, 94, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
)
data$Average_height_cm <- sapply(strsplit(Average_height_cm, "–"), function(x) mean(as.numeric(x)))
data$`Germination percentage`<- as.numeric(`Germination percentage`)
data$group <- c(rep("Group 1.1",13), rep("Group 1.2",17))


#write.table(data, file = "pheno_data", quote = F, sep = ",")

####

## all of group 1

data3 <- data[data$Genotype %in% c("BdTR1A",	"BdTR1B",	"BdTR1E",	"BdTR1F",	"BdTR1G",	"BdTR1H",	"BdTR1J",	"BdTR1K",	
                                   "BdTR1M",	"BdTR1N",	"BdTR2B",	"BdTR2C",	"BdTR2D",	"BdTR2G",	"BdTR2H",	"BdTR2J",
                                   "BdTR2K",	"BdTR2M",	"BdTR2N",	"BdTR2P",	"BdTR2R",	"BdTR2S"),]

              
par(mfrow=c(3,2), mar=c(2, 4, 4, 2), cex = 0.85)

wilcox.test(Elevation ~ group, data = data3, exact = F) #yes
boxplot(Elevation ~ group, data = data3,  
        main=expression(atop(NA,atop(bold("Elevation of collection site (m)"), "("*italic(W)*" = 18, "*italic(P-val)* " = 0.006153)"))),
        boxwex = 0.4, xlab = "", ylab= "Elevation (m)", par(cex.lab = 1.2), cex.main = 1.5)


#wilcox.test(Seed.Width ~ group, data = data3, exact = F) #yes
#boxplot(Seed.Width ~ group, data = data3, 
#        main=expression(atop(NA,atop(bold("Avg. Seed Width (cm)"), "("*italic(W)*" = 114, "*italic(P-val)* " = 0.0004036)"))), 
#        boxwex = 0.4, xlab = "", ylab= "Avg. seed width (cm)", cex.main = 1.5)

wilcox.test(Seed.yield ~ group, data = data3, exact = F) #yes
boxplot(Seed.yield ~ group, data = data3, 
        main=expression(atop(NA,atop(bold("Seed yield"), "("*italic(W)*" = 0, "*italic(P-val)* " = 4.094e-05)"))), 
        boxwex = 0.4, xlab = "", ylab= "Seed yield", cex.main = 1.5)


wilcox.test(Seed.Length ~ group, data = data3, exact = F) #yes
boxplot(Seed.Length ~ group, data = data3,
        main=expression(atop(NA,atop(bold("Avg. Seed length (cm)"), "("*italic(W)*" = 120, "*italic(P-val)* " = 8.654e-05)"))), 
        boxwex = 0.4, xlab = "", ylab= "Avg. seed length (cm)", cex.main = 1.5)

wilcox.test(Seed.production..weeks. ~ group, data = data3, exact = F) #yes
boxplot(Seed.production..weeks. ~ group, data = data3,
        main=expression(atop(NA,atop(bold("Seed production (weeks)"), "("*italic(W)*" = 0, "*italic(P-val)* " = 1.274e-05)"))), 
        boxwex = 0.4, xlab = "", ylab= "Seed production (weeks)", par(cex.lab = 1.1),cex.main = 1.5)

wilcox.test(Germination.percentage ~ group, data = data3, exact = F) #yes
boxplot(Germination.percentage ~ group, data = data3, 
        main=expression(atop(NA,atop(bold("Germination percentage"), "("*italic(W)*" = 0, "*italic(P-val)* " = 1.021e-05)"))), 
        boxwex = 0.4, xlab = "", ylab= "Germination percentage", par(cex.lab = 1.1), cex.main = 1.5)


#wilcox.test(Average.biomass..g. ~ group, data = data3, exact=F) # no, extreme outlier
#boxplot(Average.biomass..g. ~ group, data = data3) 


#wilcox.test(Average_height_cm ~ group, data = data3, exact=F)# yes
#boxplot(Average_height_cm ~ group, data = data3,
#        main=expression(atop(NA,atop(bold("Average height in cm"), "("*italic(W)*" = 27, "*italic(P-val)* " = 0.03177)"))), 
#        boxwex = 0.4, xlab = "", ylab= "Avg. Height (cm)", cex.main = 1.5)

#######################

## ONLY central of group 1.1

data2 <- data[data$Genotype %in% c("BdTR1E",	"BdTR1H",	"BdTR1J",	"BdTR1K", "BdTR1M",	"BdTR1N",
                                   "BdTR2B",	"BdTR2C",	"BdTR2D",	"BdTR2G",	"BdTR2H",	"BdTR2J", "BdTR2K",	"BdTR2M",	"BdTR2N",	"BdTR2P",	"BdTR2R",	"BdTR2S"),]

par(mfrow=c(3,2), mar=c(2, 4, 4, 2), cex = 0.85)

wilcox.test(Elevation ~ group, data = data2, exact = F) #no
boxplot(Elevation ~ group, data = data2,   
        main=expression(atop(NA,atop(bold("Elevation of collection site (m)"), "("*italic(W)*" = 18, "*italic(P-val)* " = 0.1009)"))),
        boxwex = 0.4, xlab = "", ylab= "Elevation (m)", par(cex.lab = 1.2), cex.main = 1.5)


#wilcox.test(Seed.Width ~ group, data = data2, exact = F) #no
#boxplot(Seed.Width ~ group, data = data2, 
#        main=expression(atop(NA,atop(bold("Avg. Seed Width (cm)"), "("*italic(W)*" = 66, "*italic(P-val)* " = 0.005604)"))), 
#        boxwex = 0.4, xlab = "", ylab= "Avg. seed width (cm)", cex.main = 1.5)

wilcox.test(Seed.yield ~ group, data = data2, exact = F) #yes
boxplot(Seed.yield ~ group, data = data2, 
        main=expression(atop(NA,atop(bold("Seed yield"), "("*italic(W)*" = 0, "*italic(P-val)* " = 0.0003498)"))), 
        boxwex = 0.4, xlab = "", ylab= "Seed yield", cex.main = 1.5)

wilcox.test(Seed.Length ~ group, data = data2, exact = F) #yes
boxplot(Seed.Length ~ group, data = data2,
        main=expression(atop(NA,atop(bold("Avg. Seed length (cm)"), "("*italic(W)*" = 72, "*italic(P-val)* " = 0.0008737)"))), 
        boxwex = 0.4, xlab = "", ylab= "Avg. seed length (cm)", cex.main = 1.5)

wilcox.test(Seed.production..weeks. ~ group, data = data2, exact = F) #yes
boxplot(Seed.production..weeks. ~ group, data = data2,
        main=expression(atop(NA,atop(bold("Seed production (weeks)"), "("*italic(W)*" = 0.5, "*italic(P-val)* " = 3.637e-07)"))), 
        boxwex = 0.4, xlab = "", ylab= "Seed production (weeks)", cex.main = 1.5)

wilcox.test(Germination.percentage ~ group, data = data2, exact = F) #yes
boxplot(Germination.percentage ~ group, data = data2, 
        main=expression(atop(NA,atop(bold("Germination percentage"), "("*italic(W)*" = 0, "*italic(P-val)* " = 4.786e-05)"))), 
        boxwex = 0.4, xlab = "", ylab= "Germination percentage", cex.main = 1.5)


#wilcox.test(Average.biomass..g. ~ group, data = data2, exact=F) # no, extreme outlier
#boxplot(Average.biomass..g. ~ group, data = data2) 


#wilcox.test(Average_height_cm ~ group, data = data2, exact=F)# no
#boxplot(Average_height_cm ~ group, data = data2)
