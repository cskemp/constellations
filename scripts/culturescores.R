  # Plot histograms showing scores for all asterisms in a sky culture

library(ggplot2)
library(scales)
library(tidyverse)
library(stringr)
library(here)
library(ggrepel)

library("scales")

#https://stackoverflow.com/questions/10558918/ggplot-2-facet-grid-free-y-but-forcing-y-axis-to-be-rounded-to-nearest-whole-n

integer_breaks <- function(n = 3, ...) {
  breaker <- pretty_breaks(n, ...)
  function(x) {
    breaks <- breaker(x)
    breaks[breaks == floor(breaks)]
  }
}

# Make Figure 3

regions <- read_csv(here("data", "macrogroups.csv")) %>% 
  mutate(culture = str_remove(culture, "_stellarium"))
  

cscores_orig <- read_csv(here("output","results", "perculturescores.csv"), col_names = c('culture', 'num', 'score') ) %>% 
  left_join(regions, by = "culture")


cscores <- cscores_orig %>% 
  mutate(culture =  recode(culture, "arabic_moon_stations"="arabic", "vanuatu_netwar"="lenekel", "mulapin"="babylonian", "india"="indian", "indomalay" = "Indo-Malay", "maya" = "mayan")) %>% 
  mutate(macrogroup = recode(macrogroup, "asia"="(AS)", "australia" = "(AU)", "nthamerica" = "(N)", "oceania" = "(O)", "sthamerica" = "(S)", "western" = "(W)")) %>%  
  mutate(culture = paste(culture, macrogroup))


cscores_ord <- cscores %>% 
  group_by(culture) %>% 
  summarize(avscore = mean(score))

perculturepic <- cscores %>%
  mutate(culture = str_to_title(culture)) %>% 
  mutate(culture = fct_reorder(culture, score,mean)) %>% 
  group_by(culture) %>%
  ggplot(aes(score)) +
  geom_histogram(bins=10) +
  scale_y_continuous(breaks=integer_breaks(n=3),limits = c(0, NA)) +
  scale_x_continuous(breaks=pretty_breaks(3)) + 
  facet_wrap(~culture, scales='free_y', ncol = 5) +
  labs(x="model score")+
  theme_classic(base_size=10) +
  theme(strip.background = element_blank())
plot(perculturepic)

ggsave(here("output","figures", "perculturehistogram.pdf"), plot = perculturepic, width = 5.625, height=4.5)


# Make Figure S14 

cscores_vsothers<- read_csv(here("output","results", "culturesvsotherscores.csv"), col_names = c('culture', 'num', 'score') ) %>% 
  mutate(culture =  recode(culture, "arabic_moon_stations"="arabic", "vanuatu_netwar"="lenekel", "mulapin"="babylonian", "india"="indian", "indomalay" = "Indo-Malay", "maya" = "mayan"))

perculture_vsothers <-  cscores_vsothers %>%
  mutate(culture = str_to_title(culture)) %>% 
  mutate(culture = fct_reorder(culture, score,mean)) %>% 
  group_by(culture) %>%
  ggplot(aes(score)) +
  geom_histogram(bins=10) +
  scale_y_continuous(breaks=pretty_breaks(3),limits = c(0, NA)) +
  scale_x_continuous(breaks=pretty_breaks(3)) + 
  facet_wrap(~culture, scales='free_y') +
  labs(x="asterism score")+
  theme_classic(base_size=10) +
  theme(strip.background = element_blank())
  
plot(perculture_vsothers)
ggsave(here("output","figures", "perculturehistogram_vsothers.pdf"), plot = perculture_vsothers, width = 6, height=4)


# Make Figure S2 

conncounts <- read_csv(here("output","results", "conncounts.csv"), col_names = c('culture', 'connected', 'disconnected') ) %>% 
  mutate(culture=str_replace(culture, '_stellarium.*', '')) %>% 
  mutate(culture=str_replace(culture, '_processed.*', '')) %>% 
  mutate(total = connected + disconnected)

cscores_summ <- cscores_orig %>% 
  group_by(culture) %>% 
  summarize(n = n())


# n: all constellations (no brightness)
# connected, disconnected: with brightness threshold

all <- left_join(cscores_summ, conncounts)
all_relabel <-  all %>% 
  mutate(culture =  recode(culture, "arabic_moon_stations"="arabic", "vanuatu_netwar"="lenekel", "mulapin"="babylonian", "india"="indian", "indomalay" = "Indo-Malay", "maya" = "mayan")) %>% 
  mutate(culture = str_to_title(culture)) 
  
connpic <- all_relabel %>%
  ggplot(aes(x=connected, y=disconnected)) +
  geom_smooth(method="lm", level=0, color="lightgray")+
  geom_point() +
  geom_text_repel(data=subset(all_relabel, disconnected > 5), aes(label=culture), size=3, max.iter=30000) +
  theme_classic(base_size=10) +
  coord_fixed(ratio=3)+
  theme(strip.background = element_blank()) +
  labs(x="# connected", y="# disconnected")
plot(connpic)

  ggsave(here("output","figures", "connected_counts.pdf"), plot = connpic, width = 5, height=2)


# Make Figure S1 
  
bs <- read_csv(here("output","data", "allhuman_mags.txt"), col_names = c('culture', 'mag')) %>% 
  mutate(culture=str_replace(culture, '_stellarium.*', '')) %>% 
  mutate(culture=str_replace(culture, '_processed.*', '')) %>% 
  mutate(culture =  recode(culture, "arabic_moon_stations"="arabic", "vanuatu_netwar"="lenekel", "mulapin"="babylonian", "india"="indian", "indomalay" = "Indo-Malay", "maya" = "mayan")) 
allbs <- left_join(bs, cscores_ord)
  
bins = seq(-1.5, 6.75, 0.25)
magpic <- allbs %>%
  mutate(culture = str_to_title(culture)) %>% 
  mutate(culture = fct_reorder(culture, avscore,mean)) %>% 
  ggplot(aes(x=mag)) +
  geom_histogram(breaks=bins) +
  theme_classic(base_size=10) +
  theme(strip.background = element_blank()) +
  geom_vline(xintercept=4.5)+
  scale_y_continuous(breaks=pretty_breaks(3)) + 
  facet_wrap(~culture, scales="free_y")+
  labs(x="magnitude")
plot(magpic)
  
ggsave(here("output","figures", "maghist.pdf"), plot = magpic, width = 6, height=4)

bs_summ <- bs %>% 
  group_by(culture) %>% 
  summarize(dropped = sum(mag>4.5)/n())


# Make Figure S15 
# compare model, human scores

mscore <-  cscores %>% group_by(culture) %>% summarize(mscore = mean(score))
hscore <-  cscores_vsothers %>% group_by(culture) %>% summarize(hscore = mean(score))
allscore <- left_join(mscore, hscore)

modvshumanplot <- allscore %>% 
  mutate(culture = str_to_title(culture)) %>% 
  ggplot(aes(x=mscore, y=hscore)) +
  geom_abline(intercept=0, slope=1, color="lightgray") +
  geom_point()+
  geom_text_repel(aes(label=culture)) +
  labs(x="model score", y="other culture score")+
  coord_fixed(ratio=1) +
  theme_classic(base_size=10) +
  theme(strip.background = element_blank()) 
show(modvshumanplot)
  
ggsave(here("output","figures", "perculturehumanvsmodel.pdf"), plot = modvshumanplot, width = 4.5, height=4.5)

# Compute statistics reported in main text

thresh <- 0.2

ppc <- read_csv(here("output","results", "paperpercent_df.csv"))  %>% 
  select(-index) %>% 
   mutate(humthresh = case_when(humanscore>= thresh ~ "acrossc",
                               TRUE ~ "not_acrossc"),
         modthresh =  case_when(modelscore >= thresh ~ "modelled", 
                                TRUE ~ "not_modelled"))
  
ppcsum <- ppc %>% 
  group_by(humthresh, modthresh) %>% 
  summarize(count = n())

ppcsum_hum <- ppc %>% 
  group_by(humthresh) %>% 
  summarize(count = n())

ppcsum_mod <- ppc %>% 
  group_by(modthresh) %>% 
  summarize(count = n())

print(paste("proportion of asterisms that recur across cultures:",
            as.character(round(pull(ppcsum_hum %>% filter(humthresh == "acrossc"), "count") /sum(ppcsum_hum$count), 2)))) 

print(paste("proportion of asterisms captured by model:",
            as.character(round(pull(ppcsum_mod %>% filter(modthresh== "modelled"), "count") /sum(ppcsum_mod$count), 2)))) 

print(paste("proportion of recurring asterisms captured by model:",
            as.character(round(pull(ppcsum %>% filter(modthresh== "modelled" & humthresh == "acrossc"), "count") / pull(ppcsum_hum %>% filter(humthresh == "acrossc"), "count"), 2)))) 

