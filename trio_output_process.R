#script for processing trio-GTA model results 
#author: Laura Hegemann

library(tidyverse)
library(ggsci)

setwd("Z:/projects/Neurodev_trioGCTA/results")

results <- list.files(pattern = "std.csv")

results

att_adhd <- read.csv(results[1]) %>%
  mutate(outcome = "att_adhd")
hyp_adhd <- read.csv(results[2]) %>%
  mutate(outcome = "hyp_adhd")
lang <- read.csv(results[3]) %>%
  mutate(outcome = "language")
motor <- read.csv(results[4]) %>%
  mutate(outcome = "motor")
rrb <- read.csv(results[5]) %>%
  mutate(outcome = "rrb")
sci <- read.csv(results[6]) %>%
  mutate(outcome = "sci")

result_df <- rbind(att_adhd,hyp_adhd) %>%
  rbind(rrb) %>%
  rbind(sci) %>%
  rbind(lang) %>%
  rbind(motor) %>%
  mutate(best_fit = case_when((outcome == "sci" | outcome == "language" | outcome == "motor") & model == "direct" ~ "yes", 
                              (outcome == "rrb" | outcome == "att_adhd"| outcome == "hyp_adhd") & model == "nocov" ~ "yes", 
                              TRUE ~"no")) %>%
  #filter(best_fit == "yes") %>%
  filter(model != "null") %>%
  filter(variable != "mp") %>%
  filter(variable != "e") %>%
  mutate(value = round(value,3),
         Percent_variance = value*100, 
         model = factor(model, levels = c("full", "nocov", "direct")),
         outcome = case_when(outcome == "att_adhd" ~ "Attention", 
                                                  outcome == "hyp_adhd" ~ "Hyperactivity", 
                                                  outcome == "rrb" ~ "RRBI",
                                                  outcome == "sci" ~ "Social & communication", 
                                                  outcome == "motor" ~ "Motor", 
                                                  outcome == "language" ~ "Language"), 
         variable = case_when(variable == "m" ~ "\n Maternal \n", 
                              variable == "p" ~ "\nPaternal \n", 
                              variable == "o" ~ "\nChild \n",
                              variable == "om" ~ "Maternal/Child\nCovariance\n",
                              variable == "op" ~ "Paternal/Child\nCovariance\n"),
         label = case_when(best_fit == "yes" ~ "*"))
 

ggplot(result_df, aes(fill=variable, y=Percent_variance, x=model, color = "white")) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ outcome, switch = "y") +  
  #scale_alpha_discrete(range = c(0.68, 1)) +
  scale_color_manual(values = c("white")) +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"),
        panel.border = element_rect(linetype = "dashed", color = "grey"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 21, color = "black"),
        axis.text.x = element_text(size =18, color = "black", face = "bold"),
        axis.title = element_text(size = 24),
        axis.title.x = element_text(vjust = -1),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  labs(y = "Percent Variance",
       title = "",
       fill = "Genetic Effects",
       x = "Model") +
  guides(alpha = "none", 
         color = "none",
         fill = "none"
  )
 

fits_files <-  as_tibble(list.files(pattern = "fit.csv")) %>%
  rename("File" = value) %>%
  mutate(outcome = str_remove(File, "_(\\w+)(?:_\\w+)?\\.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
             mutate(model = c("full", "nocov", "direct", "null"), 
                    outcome = .y)
        
         ))

fit_data <- bind_rows(fits_files$data)

lrt_files <- as_tibble(list.files(pattern = "lrtest.csv")) %>%
  rename("File" = value ) %>%
  mutate(outcome = str_remove(File, "_(\\w+)(?:_\\w+)?\\.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
                       mutate(model = c("full", "nocov", "direct", "null"), 
                              outcome = .y)
                     
         ))

lrt_data <- bind_rows(lrt_files$data)

full <- lrt_data %>% filter(model == "full") %>% mutate(p_adj = NA_real_)
nocov_p <- lrt_data %>% filter(model == "nocov") %>% mutate(p_adj = p.adjust(.$p..Chisq., method = "fdr", n = 6))
direct_p <- lrt_data %>% filter(model == "direct") %>% mutate(p_adj = p.adjust(.$p..Chisq., method = "fdr", n = 6))
null_p <- lrt_data %>% filter(model == "null") %>% mutate(p_adj = p.adjust(.$p..Chisq., method = "fdr", n = 6))

lrt_data2 <- bind_rows(full, nocov_p, direct_p, null_p)

cor_files <- as_tibble(list.files(pattern = "cor.csv")) %>%
  rename("File" = value ) %>%
  mutate(outcome = str_remove(File, "_(\\w+)(?:_\\w+)?\\.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
                       mutate(outcome = .y)
                     
         ))

cor_data <- bind_rows(cor_files$data)

library(xlsx)
write.xlsx(fit_data, file = "N:/durable/projects/neurodev_general_genetic/results/fits_trioGCTA.xlsx")
write.xlsx(lrt_data, file = "N:/durable/projects/neurodev_general_genetic/results/lrt_trioGCTA.xlsx")
write.xlsx(cor_data, file = "N:/durable/projects/neurodev_general_genetic/results/cor_trioGCTA.xlsx")

result_df_all <- rbind(att_adhd,hyp_adhd) %>%
  rbind(rrb) %>%
  rbind(sci) %>%
  rbind(lang) %>%
  rbind(motor) %>%
  mutate(value = round(value,3),
         outcome = case_when(outcome == "att_adhd" ~ "Attention", 
                             outcome == "hyp_adhd" ~ "Hyperactivity", 
                             outcome == "rrb" ~ "RRBI",
                             outcome == "sci" ~ "Social & communication", 
                             outcome == "motor" ~ "Motor", 
                             outcome == "language" ~ "Language")) %>%
  rename("Genetic Effect" = "variable", 
         "Estimate" = "value") %>%
  pivot_wider(names_from = "Genetic Effect", values_from = Estimate) %>%
  mutate(best_fit = case_when((outcome == "Social & communication" | outcome == "Language" | outcome == "Motor") & model == "direct" ~ "yes", 
                       (outcome == "RRBI" | outcome == "Attention"| outcome == "Hyperactivity") & model == "nocov" ~ "yes", 
                       TRUE ~"no"))

write.xlsx(result_df_all, file = "N:/durable/projects/neurodev_general_genetic/results/estimates_TrioGCTA.xlsx")

