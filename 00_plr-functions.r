# Liam GW Johnson 2022-01-24
# define functions for use in '01_plr-simulation.r'
# analysis from manuscript:
# "Logistic regression is ineffective for predicting host susceptibility based on phylogenetic distance"
# written for R version 4.2.1

# load packages
library(ape) # version 5.6.2
library(phylolm) # version 2.6.2


# simulate a phylogeny -------------------------------------------------------------

sim.phy <- function(){

    tr <- ape::rphylo(n = 500,
                      birth = 0.3,
                      death = 0.2,
                      T0 = 100, # total depth in MYr
                      fossils = FALSE,
                      eps = 0.001 )

    return(tr)

}


# simulate susceptibility trait on a phylogeny -------------------------------------

sim.trait.on.phy <- function(phy, alpha){

    n_cond <- FALSE; while(!n_cond){

        chr <- phylolm::rbinTrait(phy = phy,
                                  beta = -3,
                                  alpha = alpha )

        # beta controls the mean trait value; beta = 0 gives about 50/50 1/0,
        # beta = -1.4 gives about 20/80 1/0

        # alpha controls the phylogenetic signal - towards 0 the trait is more
        # clustered, towards 3 more dispersed

        # check if number of '1' character states falls within specified bounds
        if(sum(chr) <= 25 & sum(chr) >= 6){ n_cond <- TRUE } }

    return(chr)

}


# calculate phylogenetic distance metric ------------------------------------------

calc.pd.metric <- function(phy, chr){

    # get all pairwise phylogenetic distances
    pd <- ape::cophenetic.phylo(phy)

    # select distances from each tip to each host
    hostdist_pd <- pd[names(chr)[which(chr == 1)],]

    # for hosts, set zero-distance to self to NA
#    hostdist_pd[which(hostdist_pd==0)] <- NA

    # calculate log mean distance from each tip to all hosts
    logpd <- log10(colMeans(hostdist_pd, na.rm = TRUE) + 1) # does +1 need to be here?

    return(logpd)

}


# fit a logistic regression model --------------------------------------------------

fit.lr.model <- function(chr, pd){

    lr_out <- glm(chr ~ logpd,
                  family = binomial(link = "logit"),
                  na.action = "na.omit")

    return(lr_out)

}


# use lr model to predict susceptibility of tips -----------------------------------

lr.pred <- function(pd, mod_out){

    # generate predictions using model coefficients
    val_preds <- mod_out$coefficients[1] + mod_out$coefficients[2] * pd

    # convert to predicted probability
    val_prob <- exp(val_preds) / (1 + exp(val_preds))

    # if values in val_preds are >710, exp() gives Inf and Inf/Inf gives NaN
    # convert these NAs to 1's
    val_prob[which(is.na(val_prob))] <- 1

    return(val_prob)

}


# select and apply modifications to trait data -------------------------------------

modify.data <- function(chr, prop_modified, modification_type){

    # specify modification to be made
    if(modification_type == "Remove zeroes"){ from_val <- 0; to_val <- NA }
    if(modification_type == "Remove ones"){ from_val <- 1; to_val <- NA }
    if(modification_type == "Change zeroes to ones"){ from_val <- 0; to_val <- 1 }
    if(modification_type == "Change ones to zeroes"){ from_val <- 1; to_val <- 0 }

    # identify values not to be changed
    keep_ind <- which(chr == c(1,0)[which(c(1,0) != from_val)] )

    # sample the remaining values and convert 'prop_modified' of them to new value
    chr[sample(which(chr==from_val), (length(chr[-keep_ind])*prop_modified))] <- to_val

    return(chr)

}


