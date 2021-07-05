# Plot correlation between model weights and edge weights

library(ggplot2)
library(scales)
library(tidyverse)
library(stringr)
library(here)
library(ggrepel)

models = c("model_n", "model_n_nobright", "model_n_nodistance", "model_n_nodilation") 

mnames= c("GC model", "GC (no brightness)", "GC (no proximity)", "GC (no scaling)") 

greekletters = c("Alp" = "alpha", "Bet" = "beta",  "Gam" = "gamma",  "Del" = "delta",  "Eps"  = "epsilon",   "Zet" = "zeta",  "Eta" = "eta",  "The" = "theta",  "Iot" = "iota",  "Kap" = "kap",  "Lam" = "lambda",  "Mu" = "mu",  "Nu" = "nu",  "Xi" = "xi",  "Omi" = "omicron",  "Pi" = "pi",  "Rho" = "rho",  "Sig" = "sigma",   "Ups" = "upsilon",  "Phi" = "phi",  "Chi" = "chi",  "Psi" = "psi",  "Ome" = "omega")

for (i in seq(4)) {
  mname <- mnames[i]
  model <- models[i]
  eweights <- read_csv(here("output","results", paste0("modelhumanedges_", model, "[3.5, 4.0, 4.5].csv")), col_names = c('h1', 'h2', 'name1', 'name2', 'human', 'model') )  %>% 
    mutate_at(c("name1","name2"), funs(str_replace(., "^\\d+(?=[:upper:][[:lower:][:digit:]]+[:upper:])", ""))) %>% 
    mutate_at(c("name1","name2"), funs(str_replace(., "(?<=[:lower:])\\d+(?=[:upper:])", ""))) %>% 
    mutate_at(c("name1","name2"), funs(str_replace(., "(?<=[:lower:]|\\d)(?=[:upper:])", " * ' ' * '"))) %>% 
    mutate_at(c("name1","name2"), funs(str_replace(., "$", "'"))) %>% 
    mutate(first = if_else(name1 < name2, name1, name2)) %>% 
    mutate(second = if_else(name1 < name2, name2, name1)) %>% 
    mutate(name1= first, name2=second) %>% 
    mutate_at(c("name1","name2"), funs(str_replace_all(., greekletters))) %>% 
    mutate(name = paste(name1, name2, sep=' * ", " * ')) %>% 
    mutate(modelprop= model/max(model), humanprop = human/max(human))  %>%       mutate(name = ifelse(modelprop > 0.594 | humanprop > 0.69, name, "''"))
  
  cc <- eweights %>%
      summarize(r2 = cor(human, model))
  
  corrplot <- eweights %>% 
    ggplot(aes(x=model, y = human)) +
    geom_point() +
    geom_text_repel(data =subset(eweights, modelprop > 0.45 | humanprop > 0.5), 
                    aes(x=model, y=human, label=name),
                    parse=TRUE,
                    segment.alpha = 0,
                    max.iter=50000,
                    #hjust=0,
                    point.padding=0.03,
                    color="lightgrey",
                    size = 3) +
    theme_classic(base_size=10) +
    coord_cartesian(clip = "off") +
    theme(aspect.ratio=1) +
    xlim(0, 1.25*max(eweights$model)) +
    labs(x=mname) +
    theme(strip.background = element_blank()) +
    #geom_text_repel(data =subset(eweights, modelprop > 0.4 | humanprop > 0.5), aes(label=name),
   #ggtitle(paste0("'corr' == ", round(cc$r2,2)), parse = TRUE)
   ggtitle(paste0("r= ", round(cc$r2,2)))+
   theme(plot.title=element_text(size=10, hjust=0.5)) 
  
    
  plot(corrplot)
  
  ggsave(here("output","figures", paste0(model, "_corr.eps")), plot = corrplot, width = 3.2, height=3.2)
    ggsave(here("output","figures", paste0(model, "_corr.pdf")), plot = corrplot, width = 3.2, height=3.2)
}  
    
