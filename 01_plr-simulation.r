# Liam GW Johnson 2022-01-24
# analysis from manuscript:
# "Logistic regression is ineffective for predicting host susceptibility based on phylogenetic distance"
# uses custom functions from '00_plr-functions.r'
# written for R version 4.2.1

# load packages
library(foreach) # version 1.5.2
library(doParallel) # version 1.0.17
library(data.table) # version 1.14.2

# set up for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
doParallel::registerDoParallel(cl)

# set simulation parameters to iterate through
alphas <- c(3, 1, 0.1)
modifications <- c("Remove zeroes", "Remove ones", "Change zeroes to ones",
                   "Change ones to zeroes")
props_modified <- c(0.1, 0.5, 0.75)


# define one simulation run ----------------------------------------------------------

plr.sim <- function(alpha, modification, prop_modified){

    # source subroutine functions
    # located here so that functions will be available to all parallel processes
    source("00_plr-functions.r")

    # simulate a phylogeny
    phy <- sim.phy()

    # simulate susceptibility trait on the phylogeny
    trait_full <- sim.trait.on.phy(phy, alpha)

    # calculate phylogenetic distance metric
    logpd <- calc.pd.metric(phy, trait_full)

    ## make object 'logpd' available in the global environment
    ## needed to allow parallelization due to weird behaviour of 'foreach()'
    .GlobalEnv$logpd <- logpd

    # fit a logistic regression model on the trait data
    lr_full <- fit.lr.model(trait_full, logpd)

    # use lr model to predict susceptibility of tips
    predictions_full <- lr.pred(logpd, lr_full)

    # modify the trait data to introduce some errors
    trait_modified <- modify.data(trait_full, prop_modified, modification)

    # fit a second lr model on the modified data
    lr_modified <- fit.lr.model(trait_modified, logpd)

    # predict using modified-data lr model
    predictions_modified <- lr.pred(logpd, lr_modified)

    ## arrange simulation outputs
    sim_out <- data.frame(cbind(logpd, trait_full, trait_modified,
                                predictions_full, predictions_modified))
    sim_out$alpha <- alpha
    sim_out$modification <- modification
    sim_out$prop_modified <- prop_modified
    sim_out$phy <- i
    sim_out$tip <- rownames(sim_out)

    sim_out <- sim_out[,c(6:10,1:5)]

    sim_out$full_intercept <- lr_full$coefficients["(Intercept)"]
    sim_out$full_int_stdError <- coef(summary(lr_full))[1,2]
    sim_out$full_int_p <- coef(summary(lr_full))[1,4]
    sim_out$full_slope <- lr_full$coefficients["logpd"]
    sim_out$full_slope_stdError <- coef(summary(lr_full))[2,2]
    sim_out$full_slope_p <- coef(summary(lr_full))[2,4]

    sim_out$modified_intercept <- lr_modified$coefficients["(Intercept)"]
    sim_out$modified_int_stdError <- coef(summary(lr_modified))[1,2]
    sim_out$modified_int_p <- coef(summary(lr_modified))[1,4]
    sim_out$modified_slope <- lr_modified$coefficients["logpd"]
    sim_out$modified_slope_stdError <- coef(summary(lr_modified))[2,2]
    sim_out$modified_slope_p <- coef(summary(lr_modified))[2,4]

    return(sim_out)

}


# run simulations --------------------------------------------------------------------

for(j in 1:length(alphas)){

    alpha <- alphas[j]

    for(k in 1:length(modifications)){

        modification <- modifications[k]

        for(l in 1:length(props_modified)){

            prop_modified <- props_modified[l]


            sim_out_1comb <- foreach::foreach(i = 1:1000, .combine = rbind) %dopar%{

                sim_out <- plr.sim(alpha, modification, prop_modified)

                sim_out

            }

            # combine outputs from this run with those from prior runs
            if(j + k + l == 3){ sim_out_all <- sim_out_1comb }else{
                sim_out_all <- data.table::rbindlist(list(sim_out_all, sim_out_1comb))}

}}}


# write outputs to file
data.table::fwrite(sim_out_all, "sim-output.csv", row.names = FALSE)


# stop cluster
stopCluster(cl)


