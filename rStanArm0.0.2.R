require(tidyverse)
require(coda)
require(broom)
require(bayesplot)
require(MCMCpack)
require(R2jags)
require(rstan)
require(rstanarm)
require(brms)
require(dplyr)
library(magrittr)
library(dplyr)
library(forcats)
library(tidyr)
library(purrr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(ggridges)
library(cowplot)

## do not exclude these options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## read in data
IaData<-read.csv(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Data/Neurophysiology/Ia/Ia_Final_1.csv", header=TRUE, sep=",", check.names = FALSE) # read in data and filter
DataIaRaw<-IaData[,c(2,5,9:25,29:43,117)]

CanData<-read.csv(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Ia Neurons/CanonicalVariablesAll.csv", header=TRUE, sep=",", check.names = FALSE) # read in data and filter
head(CanData)
CanData<-CanData[,c(-1)]
head(CanData)
##little wrangling
CanData$cancer<-as.factor(CanData$cancer)
CanData$chemo<-as.factor(CanData$chemo)



## check data visually 
boxplot(LD1 ~ cancer * chemo, CanData)




#####Build a model and perform model comparison on different predictors
post1 <- stan_glm(LD1 ~ cancer, 
                  data = CanData,
                  prior_intercept = normal(0, 10),   ## weakly informative Gaussian prior on the intercept
                  prior = normal(0, 10, autoscale = FALSE),  ## weakly informative Gaussian prior on the treatment effect
                  prior_aux = student_t(3, 0, 1, autoscale = FALSE), ## student-t prior on the variance 
                  adapt_delta = .99,
                  iter = 4000, 
                  warmup = 400,
                  chains = 4, 
                  thin = 2)
post2 <- update(post1, formula = . ~ chemo)
post3 <- update(post1, formula = . ~ cancer + chemo)
(post4 <- update(post1, formula = . ~ cancer * chemo))

#visualize models and fit
png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Bayesian Inference/4Way Model Comparison.png",res = 250, width = 1000, height =  1000 ,bg = "transparent")

color_scheme_set("gray")
bayesplot_grid(
  pp_check(post1), pp_check(post2), pp_check(post3), pp_check(post4) ,
  xlim = c(-10,10),  
  ylim = c(0,.3), 
  grid_args = list(ncol = 2)
)

dev.off()

#### Leave-one-out cross validation of models
loo1 <- loo(post1, save_psis = TRUE)
loo2 <- loo(post2, save_psis = TRUE)
loo3 <- loo(post3, save_psis = TRUE)
loo4 <- loo(post4, save_psis = TRUE)
(comp <- compare_models(loo1, loo2, loo3, loo4, detail = TRUE))
compare_models(loo1,loo2)  ### get the difference in SE and ELPD to prove model differences are substnatial 
compare_models(loo1,loo3)  ### get the difference in SE and ELPD to prove model differences are substnatial 
compare_models(loo1,loo4)  ### get the difference in SE and ELPD to prove model differences are substnatial 
compare_models(loo2,loo4)  ### get the difference in SE and ELPD to prove model differences are substnatial 
compare_models(loo3,loo4)  ### get the difference in SE and ELPD to prove model differences are substnatial 

launch_shinystan(post4)

################# once the best model is chosen lets subject it to more stringent validation

##really check the model fit in a more controlled manner, choose number of draws from posterior 
##*** this can be used as an alternative to the pp_check above for multple model comparison
posterior <- as.matrix(post4)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("cancer1", "chemo1", "cancer1:chemo1", "sigma"),
           prob = 0.8) + plot_title

png(file="/Users/nickhousley/Desktop/Projects /Papers/GT_2018_CancerChemoInteraction_Paper_1/GT_2018_CancerChemoInteraction_Paper_1/Figures/Bayesian Inference/ModelComparisonPost1.png",res = 250, width = 1000, height =  500 ,bg = "transparent")

color_scheme_set("gray")
ppc_dens_overlay(y = post1$y,
                 yrep = posterior_predict(post1, draws = 500))+xlim(-10,10)
dev.off()


### check for divergent transitions 
color_scheme_set("gray")
mcmc_scatter(
  as.matrix(post4),
  pars = c("cancer1", "cancer1:chemo1"), 
  np = nuts_params(post4), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.95)
)


## NUTS ENERGY DIAGNOSIS
lp_cp <- log_posterior(post4)
np_cp <- nuts_params(post4)
posterior_cp <- as.array(post4)

color_scheme_set("darkgray")
mcmc_parcoord(posterior_cp, np = np_cp)

color_scheme_set("gray")
np <- nuts_params(post4)
mcmc_nuts_energy(np) + ggtitle("NUTS Energy Diagnostic")

mcmc_nuts_divergence(post4,chain = 4)
mcmc_pairs(posterior_cp, np = np_cp,  pars = c("cancer1", "cancer1:chemo1", "chemo1", ("(Intercept)")))



color_scheme_set("red")
mcmc_nuts_divergence(np_cp, lp_cp,chain = 4)

## PSIS diagnostics Pareto smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO) 
print(loo1)
plot(loo1, label_points = TRUE)

### autocorrelation
stan_ac(post4)

### check to make sure there is adequate mixing
color_scheme_set("gray")
mcmc_trace(as.array(post4))

### check to make sure we have reached stationary distribution rhat (gelman rubin diagnositc)
stan_rhat(post4, bins=100)+scale_x_continuous(limit = c(.99, 1.01))

### which measures of central tendency should we use. Visually inspect posteriors and see if mean, median, mode are appropriate.
mcmc_dens(as.array(post4))

### prior to posterior improvement in precision

posterior_vs_prior(post4, color_by = "vs", group_by = TRUE,
                   facet_args = list(scales = "free_y"))



## posterior distributions
plot_title <- ggtitle("Posterior distributions",
                      "with means and 80% intervals")
mcmc_areas(as.matrix(post4), regex_pars = "Intercept|^cancer|^chemo|sigma", prob = 0.8)+plot_title

### check posterior predictive accuracy 
y_pred = posterior_predict(post4)
newdata = CanData %>% cbind(t(y_pred)) %>% gather(key = "Rep",
                                               value = "LD1", -cancer,-chemo,-LD1)
ggplot(newdata) +
  geom_violin(aes(y = LD1, x = cancer, fill = "Model"),alpha = 0.5) +
  geom_violin(data = CanData, aes(y = LD1, x = cancer,fill = "Obs"), alpha = 0.5) +
  geom_point(data = CanData, aes(y = LD1,x = cancer), position = position_jitter(width = 0.1, height = 0),
             color = "black")

ggplot(newdata) +
  geom_violin(aes(y=LD1, x=chemo, fill='Model', group=chemo, color=cancer), alpha=0.5)+
  geom_violin(data = CanData, aes(y = LD1, x = cancer,fill = "Obs"), alpha = 0.5) +
  geom_point(data=CanData, aes(y=LD1, x=chemo, group=chemo,color=cancer))



### alternative posterior predictive check
grid = CanData %>%
  data_grid(cancer, chemo)

fits = grid %>%
  add_fitted_draws(post4)

preds = grid %>%
  add_predicted_draws(post4)

CanData %>%
  ggplot(aes(y = cancer:chemo, x = LD1)) +
  stat_intervalh(aes(x = .prediction), data = preds) +
  stat_pointintervalh(aes(x = .value), data = fits, .width = c(.66, .95), position = position_nudge(y = -0.2)) +
  geom_point() +
  scale_color_brewer()


##Summary statistics
summary(post4)

tidyMCMC(post4$stanfit, conf.int = TRUE, conf.method = "HPDinterval",
         rhat = TRUE, ess = TRUE)


##R^2
mcmc <- as.matrix(post4)
Xmat = model.matrix(~cancer * chemo, CanData)
wch = c(which(colnames(mcmc) == "(Intercept)"), grep("^cancer1|^chemo1", colnames(mcmc)))
coefs = mcmc[, wch]
fit = coefs %*% t(Xmat)
resid = sweep(fit, 2, post4$y, "-")
var_f = apply(fit, 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")



#2-factor ANOVA frequentist check
summary(lm(LD1 ~ cancer * chemo, CanData))










