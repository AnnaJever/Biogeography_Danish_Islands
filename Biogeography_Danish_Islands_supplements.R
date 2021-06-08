### Island biogeography revisited: How dispersal characteristics shape the relationship of insular plant diversity with area and isolation on Danish islands

### Walentowitz et al.

### last changes applied on 08. June 2021

### Supplementary materials

# use data, packages and additional calculations from "Biogeography_Danish_Islands_code.R"

### Plot spec richness and dispersal syndromes ###

# create subset of saltwater islands with at least 10 m distance from the mainland for these calculations to keep the species pool consistent

# analysis related to dispersal syndromes (figure3) are repeated including all islands (saltwater and sweet water islands), but with >10 m distance from the mainland

env3 <- env[env$distance > 10,]

## % dispersal syndrome in relationship with area and isolation ###

# a) area - seed weight (H1)
fit_a <- glm(round(mean.seed.weight)~log10(area),data=env3,family=poisson())
summary(fit_a)
m_a_pseudoR2  <- lrm(round(env3$mean.seed.weight)~log10(env3$area))

plot_a <- ggplot(env3, aes(log10(area), mean.seed.weight)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 2, y = 75,
           label= paste0("Pseud R² = ", round(m_a_pseudoR2$stats[["R2"]][1], digits = 2), " ", star(summary(fit_a)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Mean seed mass (mg)")

# b) area - zoochory
fit_b <- glm(round(perc.zoochor)~log10(area),data=env3,family=poisson())
summary(fit_b)
m_b_pseudoR2  <- lrm(round(env3$perc.zoochor)~log10(env3$area))

plot_b <- ggplot(env3, aes(log10(area), perc.zoochor)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 2, y = 42,
           label= paste0("Pseud R² = ", round(m_b_pseudoR2$stats[["R2"]][1], digits = 2), " ", star(summary(fit_b)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Zoochory (%)")

# c) area - hydrochory
fit_c <- glm(round(perc.nautochor)~log10(area),data=env3,family=poisson())
summary(fit_c)
m_c_pseudoR2  <- lrm(round(env3$perc.nautochor)~log10(env3$area))

plot_c <- ggplot(env3, aes(log10(area), perc.nautochor)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 2, y = 52,
           label= paste0("Pseud R² = ", round(m_c_pseudoR2$stats[["R2"]][1], digits = 2), " ", star(summary(fit_c)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Hydrochory (%)")

# d) area - anemochory
fit_d <- glm(round(perc.meteorochor)~log10(area),data=env3,family=poisson())
summary(fit_d)

plot_d <- ggplot(env3, aes(log10(area), perc.meteorochor)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Area (ha, log)", y = "Anemochory (%)")

# e) area - autochory
fit_e <- glm(round(perc.autochor)~log10(area),data=env3,family=poisson())
summary(fit_e)
m_e_pseudoR2  <- lrm(round(env3$perc.autochor)~log10(env3$area))

plot_e <- ggplot(env3, aes(log10(area), perc.autochor)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 2, y = 20,
           label= paste0("Pseud R² = ", round(m_e_pseudoR2$stats[["R2"]][1], digits = 2), " ", star(summary(fit_e)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Autochory (%)")

# f) isolation - seed weight (H1)
fit_f <- glm(mean.seed.weight ~ log10(distance),data=env3,family=poisson())
summary(fit_f)
m_f_pseudoR2  <- lrm(round(env3$mean.seed.weight)~log10(env3$distance))

plot_f <- ggplot(env3, aes(log10(distance), mean.seed.weight)) +
  geom_point() +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 3.5, y = 75,
           label= paste0("Pseud R² = ", round(m_f_pseudoR2$stats[["R2"]][1], digits = 3), " ", star(summary(fit_f)$coefficients[8])))+
  labs(x = "Distance (m, log)", y = "Mean seed mass (mg)")

# g) isolation - zoochory
fit_g <- glm(round(perc.zoochor)~log10(distance),data=env3,family=poisson())
summary(fit_g)

plot_g <- ggplot(env3, aes(log10(distance), perc.zoochor)) +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  geom_point() +
  labs(x = "Distance (m, log)", y = "Zoochory (%)")

# h) isolation - hydrochory
fit_h <- glm(round(perc.nautochor)~log10(distance),data=env3,family=poisson())
summary(fit_h)

plot_h <- ggplot(env3, aes(log10(distance), perc.nautochor)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Distance (m, log)", y = "Hydrochory (%)")

# i) isolation - anemochory
fit_i <- glm(round(perc.meteorochor)~log10(distance),data=env3,family=poisson())
summary(fit_i)
m_i_pseudoR2  <- lrm(round(env3$n.meteorochor)~log10(env3$distance))

plot_i <- ggplot(env3, aes(log10(distance), perc.meteorochor)) +
  geom_point() +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 3.5, y = 25,
           label= paste0("Pseud R² = ", round(m_i_pseudoR2$stats[["R2"]][1], digits = 3), " ", star(summary(fit_i)$coefficients[8])))+
  labs(x = "Distance (m, log)", y = "Anemochory (%)")

# j) isolation - autochory
fit_j <- glm(round(perc.autochor)~log10(distance),data=env3,family=poisson())
summary(fit_j)

plot_j <- ggplot(env3, aes(log10(distance), n.autochor)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Distance (m, log)", y = "Autochory (%)")
