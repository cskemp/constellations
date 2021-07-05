# Plot model scores

library(ggplot2)
library(scales)
library(tidyverse)
library(here)

mscores <- read_csv(here("output","results", "allmodelscores.csv"))
                    
mscorepic <- mscores %>%
  mutate(scorefn= recode(scorefn, adj_rand = "adjusted Rand index",
                          F= "F score",
                          F10= "F10 score")) %>% 
  mutate(scorefn= factor(scorefn, levels=c("F10 score","F score",   "adjusted Rand index"))) %>% 
  mutate(model = factor(model, levels=c("model_n","model_n_nodilation",   "model_n_nobright", "model_n_nodistance", "model_code_13", "model_knn", "model_knn_thresh", "model_onebig", "model_singleton"))) %>% 
  mutate(model = recode(model, model_n = "GC", 
                         model_n_nodilation = "GC (no scaling)",
                         model_n_nobright = "GC (no brightness)",
                         model_n_nodistance = "GC (no proximity)",
                         model_code_13 = "CODE",
                         model_knn = "k means",
                         model_knn_thresh = "k means (thresh)",
                         model_onebig = "one cluster",
                         model_singleton = "singleton")) %>% 
  ggplot(aes(model, score)) +
  geom_bar(stat="identity") +
  facet_wrap(~scorefn, scales='free_y') +
  theme_classic(base_size=10) +
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("")


plot(mscorepic)

ggsave(here("output","figures", "modelscores.pdf"), plot = mscorepic, width = 6, height=3)


