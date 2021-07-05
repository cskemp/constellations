# Find recurring constellations taking geographic regions into account
# Running this script generated the warnings pasted at the bottom of the file. I didn't pursue this further because all of the asterisms in question fall well outside the top 61 reported in Table S2

library(ggplot2)
library(tidyverse)
library(here)
library(brms)

regions <- read_csv(here("data", "macrogroups.csv"))

clusterscores <- read_csv(here("output","results", "commongroups_df.csv"),
                          col_types = cols(index = col_integer(),
                                           cluster = col_character(),
                                           .default = col_double())) %>%
  pivot_longer(boorong:western_stellarium, names_to="culture", values_to = "score") %>%
  left_join(regions, by=c("culture")) %>%
  mutate(cat = cut(score, breaks = seq(-0.1, 1.0, 0.1), ordered_result=TRUE, labels=FALSE)) %>%
  group_by(index, cluster) %>%
  nest()


mmodel <-  function(data, index) {
  fit <- brm(cat ~ 1 + (1|macrogroup),
              data = data,
              family = cumulative(),
              control = list(adapt_delta = 0.99)
             )

  ff <- fitted(fit, re_formula = NA)
  ff_single <- ff[1,1,]
  centers <- c(0, seq(0.05, 0.95, 0.1))
  pred <- sum(centers * ff_single)
  print(pred)

  line <- sprintf("%d, %f", index, pred)
  # maintain incremental file in case R crashes
  incfile <- here("output","results", "commongroups_inc.csv")
  write(line,file=incfile,append=TRUE)

  return(pred)

}

mmscores <- clusterscores %>%
  mutate(mmscore = map2_dbl(data, index, mmodel))

mmwrite <- mmscores %>%
  select(-data) %>%
  write_csv(here("output","results", "commongroups_mmscore_df.csv"))



#Warning messages:
#1: Problem with `mutate()` input `mmscore`.
#ℹ Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#ℹ Input `mmscore` is `map2_dbl(data, index, mmodel)`.
#ℹ The error occurred in group 288: index = 287, cluster = "(1569, 1524)".
#2: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#3: Problem with `mutate()` input `mmscore`.
#ℹ Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#ℹ Input `mmscore` is `map2_dbl(data, index, mmodel)`.
#ℹ The error occurred in group 305: index = 304, cluster = "(6500, 6830, 6674, 6803, 6549, 6550, 6870, 6971, 6843)".
#4: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#5: Problem with `mutate()` input `mmscore`.
#ℹ Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#ℹ Input `mmscore` is `map2_dbl(data, index, mmodel)`.
#ℹ The error occurred in group 317: index = 316, cluster = "(7130, 7342, 7350)".
#6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#7: Problem with `mutate()` input `mmscore`.
#ℹ Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#ℹ Input `mmscore` is `map2_dbl(data, index, mmodel)`.
#ℹ The error occurred in group 333: index = 332, cluster = "(7878, 7910, 7818, 7631, 7769, 7514, 7903)".
#8: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#9: Problem with `mutate()` input `mmscore`.
#ℹ Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#ℹ Input `mmscore` is `map2_dbl(data, index, mmodel)`.
#ℹ The error occurred in group 368: index = 367, cluster = "(4737, 4931, 4772, 4647, 4681, 4555, 4718, 4368, 4434, 4950, 4792, 4315, 4285, 4702)".
#10: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#11: Problem with `mutate()` input `mmscore`.
#ℹ Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess
#ℹ Input `mmscore` is `map2_dbl(data, index, mmodel)`.
#ℹ The error occurred in group 447: index = 446, cluster = "(953, 1427)".
#12: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#Running the chains for more iterations may help. See
#http://mc-stan.org/misc/warnings.html#tail-ess

