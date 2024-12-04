library(ggplot2)
library(patchwork)
###### AIS data ##########
scr_folder <- "mixPFC_src"
source(paste0(scr_folder, "mixPFC_K.R"))
source(paste0(scr_folder, "utility.R"))
set.seed(2023)
# the real data we considered here is the ais data in package dr
# Chiaromonte, Cook and Li (2002)
ais <- dr::ais
sex <- ais$Sex + 1
y <- ais$LBM
x <- cbind(log(ais$SSF), log(ais$Wt), log(ais$Hg), log(ais$Ht),
           log(ais$WCC), log(ais$RCC), log(ais$Hc), log(ais$Ferr))

# pfc
FUN <- function(y) c(y, y^2)
MU_mat <- MU(x, y, type = "pfc", FUN = FUN, H = 5, sig_fit = TRUE, center_F_marg = TRUE)
beta1_pfc <- svd(solve(MU_mat$M) %*% MU_mat$sigma_fit)$u[, 1, drop = FALSE]
beta2_pfc <- svd(solve(MU_mat$M) %*% MU_mat$sigma_fit)$u[, 2, drop = FALSE]
plot(x %*% beta1_pfc, y, 
     pch = ifelse(sex == 1, 16, 17), col = ifelse(sex == 1, "blue", "red"),
     xlab = expression(beta[1]^T*X), main = "Directions by PFC")


plot(x %*% beta2_pfc, y, 
     pch = ifelse(sex == 1, 16, 17), col = ifelse(sex == 1, "blue", "red"),
     xlab = expression(beta[2]^T*X), main = "Directions by PFC")

### plot PFC and two OLS lines ###
FUN <- function(y) c(y, y^2)
MU_mat <- MU(x, y, type = "pfc", FUN = FUN, H = 5, sig_fit = TRUE, center_F_marg = TRUE)
beta1_pfc <- svd(solve(MU_mat$M) %*% MU_mat$sigma_fit)$u[, 1, drop = FALSE]

### treat this as a mixture problem, where sex is the cluster
# generate initial value 
K <- 2
FUN <- function(y) c(y, y^2)
kmean_x <- kmeans(x, K, nstart = 20)$cluster
gams_init <- gen_gams_init_K(kmean_x, true_comp_ind = FALSE)
er_init <- error_rate_K(kmean_x, sex, 2)$error_rate

# run mixpfc with the generated initial values
mixpfc_fit <- mixPFC_K(x, y, FUN = FUN, 
                       n_components = K, gams_init = gams_init,
                       lam_factor = rep(0, K), weight_init = colMeans(gams_init), 
                       center_F_marg = TRUE, compute_loglk = FALSE,
                       max_iter = 50, tol = 1e-3)

beta1 <- svd(mixpfc_fit$B_mats[[1]])$u[, 1, drop=FALSE]
beta2 <- -svd(mixpfc_fit$B_mats[[2]])$u[, 1, drop=FALSE]

sex_est <- est_component(mixpfc_fit$gams)
er_rate <- error_rate_K(sex_est, sex, 2)
sex_est <- er_rate$pred_comp

x1 <- x[sex_est == 1, ]
x2 <- x[sex_est == 2, ]
y1 <- y[sex_est == 1]
y2 <- y[sex_est == 2]


# ----- use ggplot for better quality ----- #
ais_df <- data.frame(dir1 = x %*% beta1,
                     dir2 = x %*% beta2,
                     Y = y, 
                     Gender = factor(sex_est, levels = c(1, 2), labels = c("male", "female")))

p1 <- ggplot(ais_df, aes(x = dir1, y = Y, color = Gender, shape = Gender)) +
  geom_point() +
  geom_smooth(data = subset(ais_df, Gender == "male"),
              method = "lm", se = FALSE, color = rgb(215/255, 25/255, 28/255)) +
  labs(x = bquote(hat(beta)[1]^T * X), y = "Lean Body Mass") + 
  scale_color_manual(values = c("male" = rgb(215/255, 25/255, 28/255),   # RGB for male
                                "female" = rgb(43/255, 131/255, 186/255))) + # RGB for female
  scale_shape_manual(values = c("male" = 1, "female" = 17)) +
  guides(color = "none", shape = "none")  # This line removes the color legend

p2 <- ggplot(ais_df, aes(x = dir2, y = Y, color = Gender, shape = Gender)) +
  geom_point() +
  geom_smooth(data = subset(ais_df, Gender == "female"),
              method = "lm", se = FALSE, color = rgb(43/255, 131/255, 186/255)) +
  scale_color_manual(values = c("male" = rgb(215/255, 25/255, 28/255),   # RGB for male
                                "female" = rgb(43/255, 131/255, 186/255)),
                     guide = "none") + # RGB for female
  scale_shape_manual(values = c("male" = 1, "female" = 17),
                     guide = "none") +
  labs(x = bquote(hat(beta)[2]^T * X), y = "Lean Body Mass")  +
  #guides(color = "none", shape = "none")  # This line removes the color legend
  guides(color = guide_legend(override.aes = list(shape = c(1, 17))))  # Combine the color and shape legends

ais_df_pfc <- data.frame(dir1 = x %*% beta1_pfc, y = y,
                         Gender = factor(sex, levels = c(1, 2), labels = c("male", "female")))
p3 <- ggplot(ais_df_pfc,
             aes(x = dir1, y = y, color = Gender, shape = Gender)) +
  geom_point() +
  geom_smooth(data = subset(ais_df_pfc, Gender == "male"),
              method = "lm", se = FALSE, color = rgb(215/255, 25/255, 28/255)) +
  geom_smooth(data = subset(ais_df_pfc, Gender == "female"),
              method = "lm", se = FALSE, color = rgb(43/255, 131/255, 186/255)) +
  scale_color_manual(values = c("male" = rgb(215/255, 25/255, 28/255),   # RGB for male
                                "female" = rgb(43/255, 131/255, 186/255)),
                     guide = "none") + # RGB for female
  scale_shape_manual(values = c("male" = 1, "female" = 17),
                     guide = "none") +
  labs(x = bquote(hat(beta)^T * X), y = "Lean Body Mass")  +
  guides(color = "none", shape = "none")  # This line removes the color legend
  #guides(color = guide_legend(override.aes = list(shape = c(1, 17))))  # Combine the color and shape legends

p1 + p2 + p3


