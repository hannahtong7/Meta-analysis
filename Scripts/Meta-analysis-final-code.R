#==============================
#Setup and Load Libraries
#==============================
rm(list = ls())
setwd("~/Documents/temp/meta/Files")

library(dplyr)
library(metafor)
library(ggplot2)
library(kableExtra)
library(gridExtra)

#==============================
#Data Import and Filtering
#==============================
raw_data <- read.csv("Meta-analysis.csv")

#Filter for high emission scenarios (DeltaT between 3 and 5°C)
data_filtered <- raw_data %>% 
  filter(DeltaT >= 3, DeltaT <= 5)

#Save filtered data
write.csv(data_filtered, "Meta-analysis_filtered_3to5.csv", row.names = FALSE)

#=======================================
#Data Quality Checks & Pseudoreplication
#=======================================
#Check for repeated measurements (pseudoreplication) per Strain & DeltaT
strain_repeats <- data_filtered %>%
  group_by(Strain, DeltaT) %>%
  tally() %>%
  filter(n > 1)
print(strain_repeats)

#Count number of data points per taxon
taxa_counts <- data_filtered %>%
  group_by(Taxa) %>%
  tally(name = "num_rows")
print(taxa_counts)

#==============================
#Calculate Effect Sizes (lnRR) and Variance
#==============================
data_effect <- data_filtered %>%
  mutate(
    lnRR = log(Treatment_mean / Control_mean),
    lnRR_variance = (Treatment_SD^2 / (Treatment_n * Treatment_mean^2)) +
      (Control_SD^2 / (Control_n * Control_mean^2))
  )

#Identify and remove problematic rows (e.g., where means <= 0 or NA)
data_clean <- data_effect %>%
  filter(Treatment_mean > 0, Control_mean > 0)

#(Re)calculate lnRR and variance (ensuring consistency)
data_clean <- data_clean %>%
  mutate(
    lnRR = log(Treatment_mean / Control_mean),
    lnRR_variance = (Treatment_SD^2 / (Treatment_n * Treatment_mean^2)) +
      (Control_SD^2 / (Control_n * Control_mean^2))
  )

#Quick summaries and save the cleaned dataset
summary(data_clean$lnRR)
summary(data_clean$lnRR_variance)
write.csv(data_clean, "Meta-analysis_filtered_3to5_with_lnRR.csv", row.names = FALSE)
head(data_clean)

#==============================
#Meta-Analysis Models
#==============================
#Basic models for overall effect
model0 <- rma(yi = lnRR, vi = lnRR_variance, data = data_clean, method = "REML")
model1 <- rma.mv(yi = lnRR, V = lnRR_variance, 
                 random = ~ 1 | DOI, 
                 data = data_clean, method = "REML")
model2 <- rma.mv(yi = lnRR, V = lnRR_variance, 
                 random = ~ 1 | DOI/Strain, 
                 data = data_clean, method = "REML")
model2

#Compare models using AIC (lower is better)
AIC(model0, model1, model2)
#Choose the nested random effects model (model2) as the overall model
model_overall <- model2
summary(model_overall)

par(mfrow = c(1,1))

#----------------------------------
#Order data by lnRR and Create Forest Plot
#----------------------------------
clean_data_sorted <- data_clean %>%
  arrange(Taxa, lnRR) %>%
  mutate(Strain_Unique = make.unique(as.character(Strain)))

#Create a forest plot using the sorted strain names
forest(model_overall, 
       slab = clean_data_sorted$Strain_Unique, 
       order = "obs", 
       shade = "zebra",
       xlab = "Effect Size (lnRR)", 
       header = "Strain",              #Column header for the strain labels
       mlab = "Overall Effect", 
       col = "#CCCCFF")

#----------------------------------
#Plot with Taxon-Based Colors
#----------------------------------
#Define taxon color palette
taxon_colors <- c(
  "Cyanobacteria"    = "#1f77b4",
  "Diatom"           = "#636363",
  "Cocolithophore"   = "#ff7f0e",
  "Dinoflagellate"    = "#2ca02c"
)

#Assign colors to each row in clean_data_sorted
row_colours <- taxon_colors[clean_data_sorted$Taxa]

#Plot with taxon-based colors
par(mfrow = c(1,1))
par(xpd = TRUE, mar = c(5, 4, 4, 9))  #Give room for a legend 

forest(model_overall, 
       slab   = clean_data_sorted$Strain_Unique, 
       order  = "obs", 
       shade  = "zebra",
       xlab   = "Effect Size (lnRR)", 
       header = "Strain",              
       mlab   = "Overall Effect", 
       col    = "#CCCCFF",  #CI line color
       colout = row_colours,#fill color by taxon
       cex    = 0.8)

legend("bottomright", inset = c(-0.09, -0.13), 
       legend = names(taxon_colors), 
       fill = taxon_colors, 
       border = "black", 
       box.col = "white", 
       cex = 0.8)

#Forest Plot for Overall Model (Alternate Plot)
forest(model_overall, slab = data_clean$Strain, 
       xlab = "Effect Size (lnRR)", mlab = "Overall Effect")
#This plot shows individual strain effect sizes and the overall summary effect.

#==============================
#Moderator Model: Taxa and Latitude
#==============================
#Clean the Latitude column
data_clean$Latitude <- gsub("−", "-", data_clean$Latitude)  #Replace Unicode minus signs
data_clean$Latitude[data_clean$Latitude == ""] <- NA           #Replace empty strings with NA
data_clean$Latitude <- as.numeric(data_clean$Latitude)         #Convert to numeric

#Remove rows with missing Latitude values and center the Latitude variable
data_clean <- data_clean %>% filter(!is.na(Latitude))
data_clean$Latitude_centered <- scale(data_clean$Latitude, center = TRUE, scale = FALSE)

#Fit the moderator model including both Taxa and centered Latitude
model_moderators <- rma.mv(yi = lnRR, V = lnRR_variance,
                           mods = ~ Taxa + Latitude_centered,
                           random = ~ 1 | DOI/Strain,
                           data = data_clean, method = "REML")
#Display the summary of the moderator model
summary(model_moderators)

#==============================
#Publication Bias Analysis
#==============================
par(mfrow = c(1,1))
par(mar = c(5, 5, 4, 2))

#Create a funnel plot for the overall model
funnel(model_overall, level = c(90, 95, 99), lty = 1, 
       shade = c("#595959", "#7F7F7F", "#A6A6A6"), 
       back = "white", legend = TRUE,
       xlab = "Effect Size (lnRR)", ylab = "Precision (1/SE)", cex.lab = 1.2)
abline(v = model_overall$b[1], col = "#0072B2", lwd = 2, lty = 2)

#Egger's regression test for funnel plot asymmetry
data_clean$sei <- sqrt(data_clean$lnRR_variance)
regtest(x = data_clean$lnRR, sei = data_clean$sei, model = "rma", predictor = "sei")

#==============================
#Taxon-Specific Analyses
#==============================
unique_taxa <- unique(data_clean$Taxa)
taxa_models <- list()
par(mfrow = c(2, 2))  #2x2 layout for 4 taxa

for (taxon in unique_taxa) {
  taxon_data <- subset(data_clean, Taxa == taxon)
  
#Ensure unique strain names and sort by lnRR
  taxon_data <- taxon_data %>%
    arrange(lnRR) %>%
    mutate(Strain_Unique = make.unique(as.character(Strain)))
  
#Fit model for each taxon
  model_taxon <- rma.mv(yi = lnRR, V = lnRR_variance,
                        random = ~ 1 | DOI/Strain,
                        data = taxon_data, method = "REML")
  taxa_models[[taxon]] <- model_taxon
  
#Forest plot for each taxon with formatting
  forest(model_taxon,
         slab   = taxon_data$Strain_Unique,
         order  = "obs",
         shade  = "zebra",
         xlab   = "Effect Size (lnRR)",
         header = "Strain",
         mlab   = paste("Overall Effect –", taxon),
         col    = "#CCCCFF",
         cex    = 0.8)
}

#Create a summary table of taxon-specific effect sizes
taxa_results <- data.frame(
  Taxon = names(taxa_models),
  Effect_Size = sapply(taxa_models, function(model) coef(model)["intrcpt"]),
  CI_Lower = sapply(taxa_models, function(model) model$ci.lb[1]),
  CI_Upper = sapply(taxa_models, function(model) model$ci.ub[1]),
  Sample_Size = sapply(taxa_models, function(model) model$k)
)

#Display the summary table 
taxa_results %>%
  kable("html", digits = 3, caption = "Taxon-specific Effect Sizes (lnRR)") %>%
  kable_styling("striped", full_width = FALSE)

#Plot taxon-specific estimates with 95% CIs
custom_colors <- c(
  "Cyanobacteria" = "#1f77b4",      
  "Diatom" = "#636363",             
  "Cocolithophore" = "#ff7f0e",      
  "Dinoflagellate" = "#2ca02c"       
)
taxa_results <- taxa_results %>%
  arrange(Effect_Size) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon))

ggplot(taxa_results, aes(x = Effect_Size, y = Taxon, color = Taxon, fill = Taxon)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.15, linewidth = 0.8) +
  geom_point(shape = 21, size = 5, stroke = 0.4, color = "black") +
  scale_color_manual(values = custom_colors, guide = "none") +
  scale_fill_manual(values = custom_colors, guide = "none") +
  scale_x_continuous(breaks = seq(-0.20, 0.9, by = 0.1), limits = c(-0.15, 0.9)) +
  theme_minimal(base_size = 13) +
  labs(x = expression("Effect Size (lnRR)"), y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10)),
        plot.margin = margin(10, 15, 10, 10))

#==============================
#Latitudinal Analysis and Plotting - Phytoplankton Overall
#==============================
#Ensure Abs_Latitude is calculated
data_clean$Abs_Latitude <- abs(data_clean$Latitude)

#Fit meta-regression model for Abs_Latitude
model_latitude_only <- rma.mv(
  yi = lnRR,
  V = lnRR_variance,
  mods = ~ Abs_Latitude,
  random = ~ 1 | DOI/Strain,
  data = data_clean,
  method = "REML"
)

#Predict across latitude range
lat_range <- seq(min(data_clean$Abs_Latitude), max(data_clean$Abs_Latitude), length.out = 100)
predictions <- predict(model_latitude_only, newmods = lat_range)
pred_df <- data.frame(
  Abs_Latitude = lat_range,
  lnRR_pred = predictions$pred
)

#Simple ggplot of overall relationship
ggplot(data_clean, aes(x = Abs_Latitude, y = lnRR)) +
  geom_point(color = "black", size = 2) +
  geom_line(data = pred_df, aes(x = Abs_Latitude, y = lnRR_pred), color = "black", linewidth = 1) +
  labs(x = "Absolute Latitude (°)", y = "Effect Size (lnRR)") +
  theme_classic()

#Generate predictions with confidence intervals
lat_range <- seq(min(data_clean$Abs_Latitude), max(data_clean$Abs_Latitude), length.out = 100)
predictions <- predict(model_latitude_only, newmods = lat_range)
pred_df <- data.frame(
  Abs_Latitude = lat_range,
  lnRR_pred = predictions$pred,
  lower_CI = predictions$ci.lb,
  upper_CI = predictions$ci.ub
)

#Make a refined plot
ggplot(data_clean, aes(x = Abs_Latitude, y = lnRR)) +
  geom_point(alpha = 0.7, size = 2.5, color = "black") +
  geom_line(data = pred_df, aes(x = Abs_Latitude, y = lnRR_pred), color = "#1f77b4", linewidth = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df, aes(x = Abs_Latitude, ymin = lower_CI, ymax = upper_CI), alpha = 0.2, fill = "#1f77b4", inherit.aes = FALSE) +
  labs(x = "Absolute Latitude (°)", y = "Effect Size (lnRR)",
       title = "Growth Response of Phytoplankton to Warming Across Latitude") +
  theme_classic(base_size = 13)

#Show full model summary and extract key values
summary(model_latitude_only)
slope_est <- model_latitude_only$b["Abs_Latitude"]
slope_se  <- model_latitude_only$se["Abs_Latitude"]
slope_p   <- model_latitude_only$pval["Abs_Latitude"]
slope_ci_lb <- model_latitude_only$ci.lb["Abs_Latitude"]
slope_ci_ub <- model_latitude_only$ci.ub["Abs_Latitude"]
cat("Abs_Latitude slope estimate =", slope_est, "±", slope_se,
    ", p-value =", slope_p,
    ", 95% CI [", slope_ci_lb, ",", slope_ci_ub, "]\n")
cat("QM statistic =", model_latitude_only$QM, "  df =", model_latitude_only$QMdf, "  p-value =", model_latitude_only$QMp, "\n")
model_no_mod <- rma.mv(
  yi = lnRR, V = lnRR_variance,
  random = ~ 1 | DOI/Strain,
  data = data_clean, method = "REML"
)
R2 <- (model_no_mod$QE - model_latitude_only$QE) / model_no_mod$QE
cat("Proportion of heterogeneity explained by Abs_Latitude =", R2, "\n")

#==============================
#Latitudinal Analysis and Plotting - by Taxa
#==============================
#Create empty lists to store models and predictions
abs_lat_models <- list()
abs_lat_pred_list <- list()

#Loop over each taxon
for (t in unique(data_clean$Taxa)) {
  taxon_data <- filter(data_clean, Taxa == t)
  
#Fit meta-regression model with Abs_Latitude as moderator
  model_abs_lat <- rma.mv(yi = lnRR, V = lnRR_variance,
                          mods = ~ Abs_Latitude,
                          random = ~ 1 | DOI/Strain,
                          data = taxon_data, 
                          method = "REML")
  abs_lat_models[[t]] <- model_abs_lat
  cat("Taxon:", t, "\n")
  print(summary(model_abs_lat))
  cat("\n----------------\n")
  
#Define range and generate predictions
  lat_range <- seq(min(taxon_data$Abs_Latitude), max(taxon_data$Abs_Latitude), length.out = 100)
  pred <- predict(model_abs_lat, newmods = lat_range)
  
  temp_df <- data.frame(
    Taxa = t,
    Abs_Latitude = lat_range,
    pred = pred$pred,
    ci.lb = pred$ci.lb,
    ci.ub = pred$ci.ub
  )
  abs_lat_pred_list[[t]] <- temp_df
}

#Combine predictions and plot
abs_lat_pred_df <- bind_rows(abs_lat_pred_list)
p_abs_lat <- ggplot(data_clean, aes(x = Abs_Latitude, y = lnRR)) +
  geom_point(aes(color = Taxa), size = 2, alpha = 0.7) +
  geom_line(data = abs_lat_pred_df, aes(x = Abs_Latitude, y = pred, color = Taxa), size = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = abs_lat_pred_df, aes(x = Abs_Latitude, ymin = ci.lb, ymax = ci.ub, fill = Taxa), alpha = 0.2, color = NA, inherit.aes = FALSE) +
  facet_wrap(~ Taxa, scales = "free_x") +
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  labs(title = "", x = "Absolute Latitude (°)", y = "Effect Size (lnRR)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")
print(p_abs_lat)

#Same plot with standardised X axis
data_clean <- data_clean %>% mutate(Abs_Latitude = abs(Latitude))
global_lat_range <- seq(min(data_clean$Abs_Latitude), max(data_clean$Abs_Latitude), length.out = 100)
abs_lat_models <- list()
abs_lat_pred_list <- list()

for (t in unique(data_clean$Taxa)) {
  taxon_data <- filter(data_clean, Taxa == t)
  
  model_abs_lat <- rma.mv(yi = lnRR, V = lnRR_variance,
                          mods = ~ Abs_Latitude,
                          random = ~ 1 | DOI/Strain,
                          data = taxon_data, method = "REML")
  abs_lat_models[[t]] <- model_abs_lat
  cat("Taxon:", t, "\n")
  print(summary(model_abs_lat))
  cat("\n----------------\n")
  
  taxon_min <- min(taxon_data$Abs_Latitude)
  taxon_max <- max(taxon_data$Abs_Latitude)
  lat_range <- global_lat_range[global_lat_range >= taxon_min & global_lat_range <= taxon_max]
  pred <- predict(model_abs_lat, newmods = lat_range)
  
  temp_df <- data.frame(
    Taxa = t,
    Abs_Latitude = lat_range,
    pred = pred$pred,
    ci.lb = pred$ci.lb,
    ci.ub = pred$ci.ub
  )
  abs_lat_pred_list[[t]] <- temp_df
}

abs_lat_pred_df <- bind_rows(abs_lat_pred_list)
p_abs_lat <- ggplot(data_clean, aes(x = Abs_Latitude, y = lnRR)) +
  geom_point(aes(color = Taxa), size = 2, alpha = 0.7) +
  geom_line(data = abs_lat_pred_df, aes(x = Abs_Latitude, y = pred, color = Taxa), size = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = abs_lat_pred_df, aes(x = Abs_Latitude, ymin = ci.lb, ymax = ci.ub, fill = Taxa), alpha = 0.2, color = NA, inherit.aes = FALSE) +
  facet_wrap(~ Taxa, scales = "fixed") +
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  labs(title = "Absolute Latitude Effects on lnRR by Taxa",
       x = "Absolute Latitude (°)", y = "Effect Size (lnRR)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")
print(p_abs_lat)

#==============================
#Consideration of Control Temperature
#==============================
#Overall meta-regression model using Control_temp as moderator
model_control <- rma.mv(yi = lnRR, V = lnRR_variance,
                        mods = ~ Control_temp,
                        random = ~ 1 | DOI/Strain,
                        data = data_clean, method = "REML")
summary(model_control)

model_control_taxa <- rma.mv(yi = lnRR, V = lnRR_variance,
                             mods = ~ Control_temp + Taxa,
                             random = ~ 1 | DOI/Strain,
                             data = data_clean, method = "REML")
summary(model_control_taxa)

#Control temperature by taxa:
control_models <- list()
unique_taxa <- unique(data_clean$Taxa)
for (taxon in unique_taxa) {
  taxon_data <- subset(data_clean, Taxa == taxon)
  model_control_taxon <- rma.mv(yi = lnRR, V = lnRR_variance,
                                mods = ~ Control_temp,
                                random = ~ 1 | DOI/Strain,
                                data = taxon_data, method = "REML")
  control_models[[taxon]] <- model_control_taxon
  cat("Taxon:", taxon, "\n")
  print(summary(model_control_taxon))
  cat("\n----------------\n")
}

ct_pred_list <- list()
for (t in unique(data_clean$Taxa)) {
  taxon_data <- filter(data_clean, Taxa == t)
  ct_range <- seq(min(taxon_data$Control_temp), max(taxon_data$Control_temp), length.out = 100)
  pred <- predict(control_models[[t]], newmods = ct_range)
  temp_df <- data.frame(
    Taxa = t,
    Control_temp = ct_range,
    pred = pred$pred,
    ci.lb = pred$ci.lb,
    ci.ub = pred$ci.ub
  )
  ct_pred_list[[t]] <- temp_df
}
ct_pred_df <- bind_rows(ct_pred_list)

p_ct <- ggplot(data_clean, aes(x = Control_temp, y = lnRR)) +
  geom_point(aes(color = Taxa), size = 2, alpha = 0.7) +
  geom_line(data = ct_pred_df, aes(x = Control_temp, y = pred, color = Taxa), size = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = ct_pred_df, aes(x = Control_temp, ymin = ci.lb, ymax = ci.ub, fill = Taxa), alpha = 0.2, color = NA, inherit.aes = FALSE) +
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  facet_wrap(~ Taxa, scales = "free_x") +
  labs(title = "Control Temperature Effects on lnRR by Taxa",
       x = "Control Temperature (°C)", y = "Effect Size (lnRR)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

par(mfrow = c(1,2))
print(p_ct)
print(p_abs_lat)
library(gridExtra)
grid.arrange(p_ct, p_abs_lat, ncol = 2)

#==============================
#Explanation of Key Choices
#==============================
# 1. Including both 'Taxa' and 'Latitude' as moderators in a single model (model_moderators)
#    helps disentangle their individual effects on lnRR.

# 2. Plotting lnRR on the y-axis against Latitude on the x-axis and using predict() from the moderator
#    model provides a more accurate overlay of the fitted relationship than geom_smooth().

#==============================
#Calculate and Summarise Heterogeneity Explained by Moderators
#==============================
model_no_mod <- rma.mv(yi = lnRR, V = lnRR_variance,
                       random = ~ 1 | DOI/Strain,
                       data = data_clean, method = "REML")
summary(model_no_mod)

model_taxa <- rma.mv(yi = lnRR, V = lnRR_variance,
                     mods = ~ Taxa,
                     random = ~ 1 | DOI/Strain,
                     data = data_clean, method = "REML")
summary(model_taxa)

model_latitude <- rma.mv(yi = lnRR, V = lnRR_variance,
                         mods = ~ Abs_Latitude,
                         random = ~ 1 | DOI/Strain,
                         data = data_clean, method = "REML")
summary(model_latitude)

model_control <- rma.mv(yi = lnRR, V = lnRR_variance,
                        mods = ~ Control_temp,
                        random = ~ 1 | DOI/Strain,
                        data = data_clean, method = "REML")
summary(model_control)

model_taxa_lat <- rma.mv(yi = lnRR, V = lnRR_variance,
                         mods = ~ Taxa + Abs_Latitude,
                         random = ~ 1 | DOI/Strain,
                         data = data_clean, method = "REML")
summary(model_taxa_lat)

model_control_taxa <- rma.mv(yi = lnRR, 
                             V = lnRR_variance,
                             mods = ~ Control_temp + Taxa,
                             random = ~ 1 | DOI/Strain,
                             data = data_clean, 
                             method = "REML")
summary(model_control_taxa)

R2_taxa <- (model_no_mod$QE - model_taxa$QE) / model_no_mod$QE
R2_lat <- (model_no_mod$QE - model_latitude$QE) / model_no_mod$QE
R2_control <- (model_no_mod$QE - model_control$QE) / model_no_mod$QE
R2_taxa_lat <- (model_no_mod$QE - model_taxa_lat$QE) / model_no_mod$QE
R2_control_taxa <- (model_no_mod$QE - model_control_taxa$QE) / model_no_mod$QE

moderator_summary <- data.frame(
  Moderator     = c("Taxa", "Abs_Latitude", "Control_temp", "Taxa + Abs_Latitude", "Control_temp + Taxa"),
  QM            = c(model_taxa$QM, model_latitude$QM, model_control$QM, model_taxa_lat$QM, model_control_taxa$QM),
  df            = c(model_taxa$QMdf, model_latitude$QMdf, model_control$QMdf, model_taxa_lat$QMdf, model_control_taxa$QMdf),
  p_value       = c(model_taxa$QMp, model_latitude$QMp, model_control$QMp, model_taxa_lat$QMp, model_control_taxa$QMp),
  QE_residual   = c(model_taxa$QE, model_latitude$QE, model_control$QE, model_taxa_lat$QE, model_control_taxa$QE),
  Pseudo_R2     = c(R2_taxa, R2_lat, R2_control, R2_taxa_lat, R2_control_taxa)
)

moderator_summary %>%
  kable("html", digits = 3, caption = "Moderator Model Summary Table") %>%
  kable_styling("striped", full_width = FALSE)
print(moderator_summary)

#Model with DOI only
model_DOI <- rma.mv(yi = lnRR, V = lnRR_variance,
                    random = ~ 1 | DOI,
                    data = data_clean, method = "ML")

#Model with DOI nested within Strain
model_DOI_Strain <- rma.mv(yi = lnRR, V = lnRR_variance,
                           random = ~ 1 | DOI/Strain,
                           data = data_clean, method = "ML")

#Compare the models
anova(model_DOI, model_DOI_Strain)

ggplot(data_clean, aes(x = lnRR, y = reorder(Strain, lnRR), color = Taxa)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lnRR - sqrt(lnRR_variance), xmax = lnRR + sqrt(lnRR_variance))) +
  facet_wrap(~ Taxa, scales = "free_y") +
  theme_minimal() +
  labs(x = "Effect Size (lnRR)", y = "Strain", title = "Growth Responses by Strain within Taxa")

###########FINAL PLOTS###################

#-----------------------------------------
# 1. Sort Data by Taxa, Then by lnRR
#-----------------------------------------
par(mfrow = c(1,1))
data_clean$Taxa <- factor(data_clean$Taxa, levels = sort(unique(data_clean$Taxa)))
data_sorted <- data_clean %>%
  arrange(Taxa, lnRR) %>%
  mutate(Strain_Unique = make.unique(as.character(Strain)))

#-----------------------------------------
# 2. Fit the Overall Random-Effects Model using sorted data
#-----------------------------------------
model_overall <- rma.mv(
  yi    = lnRR,
  V     = lnRR_variance,
  random= ~ 1 | DOI/Strain,
  data  = data_sorted,     # use the sorted data!
  method= "REML"
)

#-----------------------------------------
# 3. Define Row Positions for Studies and Diamonds
#-----------------------------------------
n_studies <- length(model_overall$yi)
n_taxa <- length(unique(data_sorted$Taxa))
n_diamonds <- n_taxa + 1

study_rows <- seq(from = n_studies, to = 1, by = -1)
diamond_rows <- seq(from = 0, by = -1, length.out = n_diamonds)
ylim_range <- c(min(diamond_rows) - 1, max(study_rows) + 1)

#-----------------------------------------
# 4. Choose Taxon Colors and Match Row Colors
#-----------------------------------------
taxon_colors <- c(
  "Diatom"          = "#636363",
  "Dinoflagellate"  = "#2ca02c", 
  "Cyanobacteria"   = "#1f77b4",
  "Cocolithophore"  = "#ff7f0e"
)
row_colors <- taxon_colors[as.character(data_sorted$Taxa)]

#-----------------------------------------
# 5. Create the Base Forest Plot
#-----------------------------------------
forest(
  x        = model_overall, 
  slab     = data_sorted$Strain_Unique,
  rows     = study_rows,
  ylim     = c(-n_taxa - 2, n_studies + 3),
  refline  = 0,
  xlab     = "Effect Size (lnRR)",
  header   = "Strain",
  shade    = "zebra",
  col      = "#444444",
  colout   = row_colors,
  cex      = 0.6,
  psize    = 1,
  addfit   = FALSE
)

#-----------------------------------------
# 6. Add Subgroup Diamonds for Each Taxon
#-----------------------------------------
unique_taxa <- unique(data_sorted$Taxa)
n_taxa <- length(unique_taxa)
diamond_row <- 0

for (taxon in unique_taxa) {
  t_data <- subset(data_clean, Taxa == taxon)
  t_model <- rma.mv(
    yi     = lnRR,
    V      = lnRR_variance,
    random = ~ 1 | DOI/Strain,
    data   = t_data,
    method = "REML"
  )
  diamond_row <- diamond_row - 1
  addpoly(
    t_model,
    row  = diamond_row,
    mlab = paste0(taxon, " Effect"),
    col  = taxon_colors[taxon]
  )
}

#-----------------------------------------
# 7. Add the Overall Effect Diamond
#-----------------------------------------
diamond_row <- diamond_row - 1
addpoly(
  model_overall,
  row  = diamond_row,
  mlab = "Overall Effect",
  col  = "#AAAAFF"
)

#-----------------------------------------
# 8. Optional: Add a Legend
#-----------------------------------------
legend(
  "bottoml", 
  inset  = c(+0.82, -0.07),
  legend = names(taxon_colors),
  fill   = taxon_colors,
  border = "black",
  box.col= "white",
  cex    = 0.6
)

#Restore old par settings if desired
par(op)

#-----------------------------------------
# 9. Arrange Final Plots Side by Side using Custom Theme
#-----------------------------------------
lnRR_range <- range(data_clean$lnRR, na.rm = TRUE)
my_theme <- theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.2),
    panel.grid.minor = element_line(color = "grey90", size = 0.2),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(hjust = 0.5)
  )

p_ct <- p_ct +
  facet_wrap(~ Taxa, scales = "free_x") +
  coord_cartesian(ylim = lnRR_range) +
  labs(
    title = "B",
    x = "Control Temperature (°C)",
    y = "Effect Size (lnRR)"
  ) +
  my_theme

p_abs_lat <- p_abs_lat +
  facet_wrap(~ Taxa, scales = "free_x") +
  coord_cartesian(ylim = lnRR_range) +
  labs(
    title = "A",
    x = "Absolute Latitude (°)",
    y = "Effect Size (lnRR)"
  ) +
  my_theme

grid.arrange(p_abs_lat, p_ct, ncol = 2, widths = c(1, 1))

version


