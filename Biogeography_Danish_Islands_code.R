### Island biogeography revisited: How dispersal characteristics shape the relationship of insular plant diversity with area and isolation on Danish islands

### Walentowitz et al.

### last changes applied on 08. June 2021

### Data preparation ###

# clearn work space and memory
gc()
rm(list = ls())

# libraries
require(tidyverse)
require(XLConnect)
require(rms)
require(gridExtra)
require(vegan)

# functions
star <- function(x){a<-""
                  if(x<0.05){a<-"*"}	
                  if(x<0.01){a<-"**"}
                  if(x<0.001){a<-"***"}	
                  a}

lmp <- function(modelobject){if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
                             f <- summary(modelobject)$fstatistic
                             p <- pf(f[1],f[2],f[3],lower.tail=F)
                             attributes(p) <- NULL
                             return(p)}

# Get data
data    <- loadWorkbook("Dispersal_Danish_Islands_data.xlsx")
env     <- readWorksheet(data, sheet = "env", header = TRUE)
spec    <- readWorksheet(data, sheet = "spec", header = TRUE)
trait   <- readWorksheet(data, sheet = "trait", header = TRUE)

# explanation of columns

# table: env
# nr = island specific ID
# Name = island name
# mainland = name of nearest mainland
# distance = distance of island to the mainland
# area = island area
# saltwater = if island is located in saltwater (1) or sweetwater (0)

# table: species
# names = species names
# other column names: correspond to island id (1 indicates species presence and 0 absence)

# table: trait
# X = project specific species ID
# name = species name
# following columns: dispersal syndrom and seed mass of species


# Check data
str(trait) # check data type
sum(!trait$name==spec$name) # should be zero, is zero
trait$seed.mass <-as.numeric(trait$seed.mass) # make numeric

### Complementary calculations

# spec nr. per island
env$nspec.new   <-colSums(spec[2:ncol(spec)]>0)

# calculate community seed mass per island
env$mean.seed.weight<-NA
for(i in colnames(spec)[2:ncol(spec)]){
                                      env[paste("X",env$nr,sep="")==i,
                                      "mean.seed.weight"] <- mean(trait[spec[,i]>0,
                                      "seed.mass"],na.rm=T)
                                      }

# calculate number of dispersal syndrome per islands
head(trait)
trait$autochor    <-as.numeric(trait$autochor)
trait$meteorochor <-as.numeric(trait$meteorochor)
trait$nautochor   <-as.numeric(trait$nautochor)
trait$zoochor     <-as.numeric(trait$zoochor)
trait[,c("autochor","meteorochor", "nautochor", "zoochor")]<-trait[,c("autochor","meteorochor", "nautochor", "zoochor")]/rowSums(trait[,c("autochor","meteorochor", "nautochor", "zoochor")],na.rm=T)

# sums of species per island and dispersal syndrom
head(trait)
env$n.autochor    <-NA
env$n.meteorochor <-NA
env$n.nautochor   <-NA
env$n.zoochor     <-NA

for(i in colnames(spec)[2:ncol(spec)]){
  env[paste("X",env$nr,sep="")==i,"n.autochor"]   <- sum(trait[spec[,i]>0,"autochor"],na.rm=T)
  env[paste("X",env$nr,sep="")==i,"n.meteorochor"]<- sum(trait[spec[,i]>0,"meteorochor"],na.rm=T)
  env[paste("X",env$nr,sep="")==i,"n.nautochor"]  <- sum(trait[spec[,i]>0,"nautochor"],na.rm=T)
  env[paste("X",env$nr,sep="")==i,"n.zoochor"]    <- sum(trait[spec[,i]>0,"zoochor"],na.rm=T)
}

# percentage dispersal syndrome of floral communities per island
env$perc.autochor     <-NA
env$perc.meteorochor  <-NA
env$perc.nautochor    <-NA
env$perc.zoochor      <-NA
env[,c("perc.autochor","perc.meteorochor", "perc.nautochor", "perc.zoochor")]<-env[,c("n.autochor","n.meteorochor", "n.nautochor", "n.zoochor")]/rowSums(env[,c("n.autochor","n.meteorochor", "n.nautochor", "n.zoochor")])*100


### Correlation area and isolation ###
cor.test(env$area,env$distance) # 0.3, thus correlation negligible

### Plot area and isolation against spec richness with lm ###

# linear model number of species ~ area (log-log space)
lm_area <- lm(round(log10(nspec.new)) ~ round(log10(area)),data=env) # poisson; log
summary(lm_area)

# plot (figure 2a)
plot_area <- ggplot(env, aes(x = log10(area), y = log10(nspec.new))) +
  geom_point() +
  geom_line(data = fortify(lm_area),
            aes(x = round(log10(env$area)), y = .fitted))+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 2, y = 3,
           label= paste0("Adj. R² = ",
                         round(summary(lm_area)$adj.r,
                               digits = 2) , " ",
                         star(summary(lm_area)$coefficients[8])))+
  labs(x = "Area (ha, log10)", y = "Species richness (log10)")

plot_area

# check residuals for normality and heteroscedasticity
hist(lm_area$residuals) # histogram of residuals
cor(x = lm_area$residuals, y = env$area, method = "pearson") # cor. of residuals with area


# linear model number of species ~ isolation (log-log space)
lm_iso <- lm(round(log10(nspec.new)) ~ round(log10(distance)),data=env[env$distance > 0,])
summary(lm_iso)

# plot (figure 2b)
plot_iso <- ggplot(env[env$distance > 0,],
                   aes(x = log10(distance), y = log10(nspec.new))) +
  geom_point() +
  geom_line(data = fortify(lm_iso),
            aes(x = round(log10(env$distance[env$distance > 0])),
                y = .fitted))+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 4, y = 3,
           label= paste0("Adj. R² = ",
                         round(summary(lm_iso)$adj.r, digits = 2), " ",
                         star(summary(lm_iso)$coefficients[8])))+
  labs(x = "Isolation (m, log10)", y = "Species richness (log10)")

plot_iso

hist(lm_iso$residuals)
cor(x = lm_iso$residuals, y = env$distance[env$distance > 0], method = "pearson") # cor. of residuals with area

### Plot spec richness and dispersal syndromes ###

# create subset of saltwater islands with at least 10 m distance from the mainland for these calculations to keep the species pool consistent

env2        <-subset(env,env$saltwater==1 & env$distance>10)

# build linear regressions
m1          <-lm(log10(n.zoochor) ~ log10(area),data= env2)
m2          <-lm(log10(n.meteorochor) ~ log10(area),data= env2)
m3          <-lm(log10(n.nautochor) ~ log10(area),data= env2)
m4          <-lm(log10(n.autochor) ~ log10(area),data= env2)

# plot (figure 3)
plot(log10(n.zoochor) ~ log10(area),
     data= env2,xlab="Area (ha, log)",
     ylab="Species richness (log)"
     ,col="goldenrod", pch = 19)

points(x = log10(env2$area), y = log10(env2$n.meteorochor),
       col = "darkslategrey", pch = 19)
points(x = log10(env2$area), y = log10(env2$n.nautochor),
       col = "darkcyan", pch = 19)
points(x = log10(env2$area), y = log10(env2$n.autochor),
       col = "coral2", pch = 19)

abline(m1, col = "goldenrod")
abline(m2, col = "darkslategrey")
abline(m3, col = "darkcyan")
abline(m4, col = "coral2")

mtext(text=bquote(Zoochory: italic(R)^2 == .(round(summary(m1)$adj.r, 2)) * .(star(lmp(m1)))),
      line=-1.5,cex=0.8,
      col="goldenrod",
      side = 3, at=-1.5, adj = 0)
mtext(text=bquote(Anemochory: italic(R)^2 == .(round(summary(m2)$adj.r, 2)) * .(star(lmp(m2)))),
      line=-2.5,cex=0.8,
      col="darkslategrey",
      side = 3, at=-1.5, adj = 0)
mtext(text=bquote(Hydrochory: italic(R)^2 == .(round(summary(m3)$adj.r, 2)) * .(star(lmp(m3)))),
      line=-3.5,cex=0.8,
      col="darkcyan",
      side = 3, at=-1.5, adj = 0)
mtext(text=bquote(Autochory: italic(R)^2 == .(round(summary(m4)$adj.r, 2)) * .(star(lmp(m4)))),
      line=-4.5,cex=0.8,
      col="coral2",
      side = 3, at=-1.5, adj = 0)

## % dispersal syndrome in relationship with area and isolation ###

# models and plot (figure 4)

# a) area - seed weight (H1)
fit_a <- glm(round(mean.seed.weight)~log10(area),data=env2,family=poisson())
summary(fit_a)
m_a_pseudoR2  <- lrm(round(env2$mean.seed.weight)~log10(env2$area))

plot_a <- ggplot(env2, aes(log10(area), mean.seed.weight)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 1, y = 60,
           label= paste0("Pseud R² = ", round(m_a_pseudoR2$stats[["R2"]][1], digits = 2), " ", star(summary(fit_a)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Mean seed mass (mg)")

# b) area - zoochory
fit_b <- glm(round(perc.zoochor)~log10(area),data=env2,family=poisson())
summary(fit_b)
m_b_pseudoR2  <- lrm(round(env2$perc.zoochor)~log10(env2$area))

plot_b <- ggplot(env2, aes(log10(area), perc.zoochor)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 1, y = 40,
           label= paste0("Pseud R² = ", round(m_b_pseudoR2$stats[["R2"]][1], digits = 2), " ", star(summary(fit_b)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Zoochory (%)")

# c) area - hydrochory
fit_c <- glm(round(perc.nautochor)~log10(area),data=env2,family=poisson())
summary(fit_c)
m_c_pseudoR2  <- lrm(round(env2$perc.nautochor)~log10(env2$area))

plot_c <- ggplot(env2, aes(log10(area), perc.nautochor)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 1, y = 52,
           label= paste0("Pseud R² = ", round(m_c_pseudoR2$stats[["R2"]][1], digits = 2), " ", star(summary(fit_c)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Hydrochory (%)")

# d) area - anemochory
fit_d <- glm(round(perc.meteorochor)~log10(area),data=env2,family=poisson())
summary(fit_d)
m_d_pseudoR2  <- lrm(round(env2$perc.meteorochor)~log10(env2$area))

plot_d <- ggplot(env2, aes(log10(area), perc.meteorochor)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 1, y = 15,
           label= paste0("Pseud R² = ", round(m_d_pseudoR2$stats[["R2"]][1], digits = 2), " ", star(summary(fit_d)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Anemochory (%)")

# e) area - autochory
fit_e <- glm(round(perc.autochor)~log10(area),data=env2,family=poisson())
summary(fit_e)
m_e_pseudoR2  <- lrm(round(env2$perc.autochor)~log10(env2$area))

plot_e <- ggplot(env2, aes(log10(area), perc.autochor)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Area (ha, log)", y = "Autochory (%)")

# f) isolation - seed weight (H1)
fit_f <- glm(round(mean.seed.weight)~log10(distance),data=env2,family=poisson())
summary(fit_f)
m_f_pseudoR2  <- lrm(round(env2$mean.seed.weight)~log10(env2$distance))

plot_f <- ggplot(env2, aes(log10(distance), mean.seed.weight)) +
  geom_point() +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 2.3, y = 50,
           label= paste0("Pseud R² = ", round(m_f_pseudoR2$stats[["R2"]][1], digits = 3), " ", star(summary(fit_f)$coefficients[8])))+
  labs(x = "Distance (m, log)", y = "Mean seed mass (mg)")

# g) isolation - zoochory
fit_g <- glm(round(perc.zoochor)~log10(distance),data=env2,family=poisson())
summary(fit_g)

plot_g <- ggplot(env2, aes(log10(distance), perc.zoochor)) +
  geom_point()+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Distance (m, log))", y = "Zoochory (%)")

# h) isolation - hydrochory
fit_h <- glm(round(perc.nautochor)~log10(distance),data=env2,family=poisson())
summary(fit_h)

plot_h <- ggplot(env2, aes(log10(distance), perc.nautochor)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Distance (m, log)", y = "Hydrochory (%)")

# i) isolation - anemochory
fit_i <- glm(round(perc.meteorochor)~log10(distance),data=env2,family=poisson())
summary(fit_i)

plot_i <- ggplot(env2, aes(log10(distance), perc.meteorochor)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Distance (m, log)", y = "Anemochory (%)")

# j) isolation - autochory
fit_j <- glm(round(perc.autochor)~log10(distance),data=env2,family=poisson())
summary(fit_j)

plot_j <- ggplot(env2, aes(log10(distance), perc.autochor)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Distance (m, log)", y = "Autochory (%)")

# Combine plots
x11()
grid.arrange(plot_a, plot_b, plot_c, plot_d, plot_e,
             plot_f, plot_g, plot_h, plot_i, plot_j,
             ncol=5)
dev.off()

### Appendix ###

## Figure S1
# Species Area Relationship

nspec <- env$nspec.new
area <- env$area

# general plot

plot(nspec ~ area,
     xlab = "Island Area (ha)", ylab = "Number of Species",
     ylim = c(1, max(nspec)),
     col = "black",
     pch = 16)

xtmp <- seq(min(area), max(area), len=151)

## The Arrhenius model
marr <- nls(nspec ~ SSarrhenius(area, k, z))
marr
lines(xtmp, predict(marr, newdata=data.frame(area = xtmp)), lwd=2, col = "dodgerblue4")

## Gleason: log-linear
mgle <- nls(nspec ~ SSgleason(area, k, slope))
lines(xtmp, predict(mgle, newdata=data.frame(area=xtmp)),
      lwd=2, col= "salmon")

## loglog
mloglog <- lm(log(nspec) ~ log(area))
mloglog
lines(xtmp, exp(predict(mloglog, newdata=data.frame(area=xtmp))),
      lwd=2, col = "goldenrod")

## Gitay: quadratic of log-linear
mgit <- nls(nspec ~ SSgitay(area, k, slope))
lines(xtmp, predict(mgit, newdata=data.frame(area=xtmp)), 
      lwd=2, col = "yellowgreen")

## Lomolino: using original names of the parameters (Lomolino 2000):
mlom <- nls(nspec ~ SSlomolino(area, Smax, A50, Hill))
mlom
lines(xtmp, predict(mlom, newdata=data.frame(area=xtmp)), 
      lwd=2, col = "tomato4")

## Michaelis-Menten
mmic <- nls(nspec ~ SSmicmen(area, slope, Asym))
lines(xtmp, predict(mmic, newdata = data.frame(area=xtmp)),
      lwd =2, col = "cadetblue4")

legend("bottomright", c("Arrhenius", "Gleason", "log-log linear", "Gitay", 
                        "Lomolino", "Michaelis-Menten"),
       col=c("dodgerblue4","salmon","goldenrod","yellowgreen","tomato4","cadetblue4"),
       lwd=c(2,2,2,2,2,2), 
       lty=c(1,1,1,1,1,1))

## compare models (AIC)
allmods <- list(Arrhenius = marr,
                Gleason = mgle,
                logloglin = mloglog,
                Gitay = mgit, 
                Lomolino = mlom,
                MicMen= mmic)
sapply(allmods, AIC)

# par(mfrow = c(3,2))
hist(marr$m$resid())
hist(mgle$m$resid())
hist(mloglog$residuals)
hist(mgit$m$resid())
hist(mlom$m$resid())
hist(mmic$m$resid())
