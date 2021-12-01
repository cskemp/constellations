library(ggplot2)
library(tidyverse)
library(here)
library(brms)

# We originally used a mixed model approach to compute the weighted scores in Table 1. This script includes the function we used (mmodel) just in case it's useful in the future.

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


macrogroupavg <-  function(data, index) {
  avgbymacrogroup <- data %>% 
    group_by(macrogroup) %>% 
    summarize(score=mean(score))
   
  pred <- mean(avgbymacrogroup$score)
  
  line <- sprintf("%d, %f", index, pred)
  incfile <- here("output","results", "commongroups_inc.csv")
  write(line,file=incfile,append=TRUE)
  
  return(pred)
}


mmodel <-  function(data, index) {
  fit <- brm(cat ~ 1 + (1|macrogroup),
              data = data,
              family = cumulative(),
              control = list(adapt_delta = 0.99),
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
  #mutate(mmscore = map2_dbl(data, index, mmodel))
  mutate(mmscore = map2_dbl(data, index, macrogroupavg))

mmwrite <- mmscores %>%
  select(-data) %>%
  write_csv(here("output","results", "commongroups_mmscore_df.csv"))
