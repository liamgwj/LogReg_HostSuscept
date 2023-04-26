# Liam GW Johnson 2023-01-25
# create figures using output from '01_plr-simulation.r'
# analysis from manuscript:
# "Logistic regression is ineffective for predicting host susceptibility based on phylogenetic distance"
# written for R version 4.2.1


# load packages
library(data.table) # version 1.14.2
library(dplyr) # version 1.0.10
library(ggplot2) # version 3.3.6
library(gridExtra) # version 2.3

# don't use scientific notation for plot axes
options(scipen = 999)

# read in simulation outputs
sim_output <- fread("sim-output-v2.csv")

# set factors as factors
sim_output$alpha <- factor(sim_output$alpha, levels = c(3,1,0.1))
sim_output$prop_modified <- factor(sim_output$prop_modified)
sim_output$trait_full <- factor(sim_output$trait_full)
sim_output$modification <- factor(sim_output$modification, levels = c("Remove zeroes",
    "Remove ones", "Change ones to zeroes", "Change zeroes to ones"))

# loop through levels of alpha and proportion modified
#for(i in 1:length(unique(sim_output$alpha))){
#    for(j in 1:length(unique(sim_output$prop_modified))){

# specify alpha and modification proportion to plot
#plot_alpha <- unique(sim_output$alpha)[i]
#plot_propmod <- unique(sim_output$prop_modified)[j]

plot_alpha <- 1

plot_propmod <- 0.5

# generate individual plots

# full model predictions vs true
fig1 <- sim_output %>%

    filter(modification != "Change zeroes to ones" & prop_modified == plot_propmod) %>%

        ggplot(aes(x = trait_full, y = predictions_full)) +
               geom_violin() +
               geom_point() +
               theme_bw() +
               theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 15)) +
               theme(strip.text.x = element_text(size = 15),
                     strip.background = element_rect(fill = NA, color = "black")) +
               facet_grid(cols = vars(alpha)) +
               labs(y = "Predicted probability of susceptibility",
                    x = "True susceptibility status") +
               theme(legend.position = "none") +
               stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
                            colour = "black", width = 0.7)

for(i in c(3,1,0.1)){
tmp <- sim_output %>%
    filter(modification != "Change zeroes to ones" & prop_modified == plot_propmod &
           alpha == i & trait_full == 0) %>% select(predictions_full)
print(i)
print(mean(tmp$predictions_full))

tmp <- sim_output %>%
    filter(modification != "Change zeroes to ones" & prop_modified == plot_propmod &
           alpha == i & trait_full == 1) %>% select(predictions_full)

print(mean(tmp$predictions_full))

}

for(i in c("Remove zeroes", "Remove ones", "Change ones to zeroes")){

tmp <- sim_output %>%
    filter(modification == i & prop_modified == plot_propmod &
           alpha == 1 & trait_full == 0) %>% select(predictions_modified)
print(i)
print(mean(tmp$predictions_modified))

tmp <- sim_output %>%
    filter(modification == i & prop_modified == plot_propmod &
           alpha == 1 & trait_full == 1) %>% select(predictions_modified)

print(mean(tmp$predictions_modified))

}



# modified data models' predictions vs true
pvt <- sim_output %>%

    filter(modification != "Change zeroes to ones" & prop_modified == plot_propmod &
           alpha == plot_alpha) %>%

        ggplot(aes(x = trait_full, y = predictions_modified)) +
               geom_violin() +
               geom_point() +
               theme_bw() +
               theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 15)) +
               theme(strip.text.x = element_text(size = 15),
                     strip.background = element_rect(fill = NA, color = "black")) +
               facet_grid(cols = vars(modification)) +
               labs(y = "Predicted probability of susceptibility",
                    x = "True susceptibility status") +
               theme(legend.position = "none") +
               stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
                            colour = "black", width = 0.7)


# prediction vs prediction
pvp <- sim_output %>%

    filter(modification != "Change zeroes to ones" & prop_modified == plot_propmod &
           alpha == plot_alpha) %>%

        ggplot(aes(x = predictions_full, y = predictions_modified)) +
               geom_abline(slope = 1, intercept = 0) +
               geom_point(alpha = 0.05) +
               theme_bw() +
               theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 15)) +
               facet_grid(cols = vars(modification)) +
               theme(strip.text = element_blank()) +
               xlim(0, 1) + ylim(0, 1) +
               theme(legend.position = "none") +
               labs(x = "Predicted value from full-data model",
                    y = "Predicted value from modified-data model")


# write to file
png(paste0("a", plot_alpha, "pm", plot_propmod, ".png"),
    type = "cairo", width = 440, height = 300, units = 'mm', res = 300)

grid.arrange(tag_facet(pvt, tag_pool = c("A","B","C"), size = 5, open = "", close = ""),
             tag_facet(pvp, tag_pool = c("D","E","F"), size = 5, open = "", close = ""),
             ncol = 1, top = "Data modification")

dev.off()

#}}


png("fig1.png", type = "cairo", width = 440, height = 160, units = 'mm', res = 300)

grid.arrange(tag_facet(fig1, tag_pool = c("A","B","C"), size = 5, open = "", close = ""), ncol = 1, top = "Phylogenetic signal (alpha)")

dev.off()


# supplement - all combinations

for(i in c("Remove zeroes", "Remove ones", "Change ones to zeroes")){

sup_pvt <- sim_output %>%

    filter(modification == i) %>%

        ggplot(aes(x = trait_full, y = predictions_modified)) +
               geom_violin() +
               geom_point() +
               theme_bw() +
               theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 15)) +
               theme(strip.text.x = element_text(size = 15),
                     strip.background = element_rect(fill = NA, color = "black")) +
               facet_grid(cols = vars(alpha), rows = vars(prop_modified)) +
               labs(y = "Predicted probability of susceptibility",
                    x = "True susceptibility status") +
               theme(legend.position = "none") +
               stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
                            colour = "black", width = 0.7)

sup_pvp <- sim_output %>%

    filter(modification == i) %>%

        ggplot(aes(x = predictions_full, y = predictions_modified)) +
               geom_abline(slope = 1, intercept = 0) +
               geom_point(alpha = 0.05) +
               theme_bw() +
               theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 15)) +
               facet_grid(cols = vars(alpha), rows = vars(prop_modified)) +
               theme(strip.text = element_blank()) +
               xlim(0, 1) + ylim(0, 1) +
               theme(legend.position = "none") +
               labs(x = "Predicted value from full-data model",
                    y = "Predicted value from modified-data model")

png(paste0("sup_fig", i, ".png"), type = "cairo", width = 440, height = 800, units = 'mm', res = 300)

grid.arrange(sup_pvt, sup_pvp, ncol = 1, top = paste0("Data modification: ", i))

dev.off()

}




