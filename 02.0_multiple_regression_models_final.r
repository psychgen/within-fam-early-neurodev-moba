#script for trio-PGS models 
#author: Laura Hegemann


library(tidyverse)
library(ggsci)
library(lmtest)
library(sandwich)
library(furrr)
library(foreign)
library(sensemakr)
library(phenotools)

#adding a total score for SCQ to existing data set using phenotools

phenovars <- available_variables(source = "moba")
head(phenovars)

pheno_var <- c("SCQ")
subscale <- c( "full")

selected_var <- phenovars %>%
  filter(.,measure %in% pheno_var) %>%
  filter(.,subscale %in% subscale) %>%
  filter(.,grepl("3",questionnaire))

var_names <- as.vector(selected_var$var_name)

data_phenotools <- curate_dataset(variables_required= var_names,
                                  pheno_data_root_dir="N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/",
                                  PDB="2306",
                                  completion_threshold=0.5,
                                  return_items=FALSE,
                                  consistent_items=FALSE,
                                  transformations=NULL)

data_scq <- data_phenotools$moba %>%
  select(c("preg_id","BARN_NR","m_id", "f_id", SCQ_scale = "scq_full_c_3yr"))


#reading in and combining phenotype and genotype data
data_pheno <- readRDS("N:/durable/projects/neurodev_general_genetic/data/data_scales.rds") %>%
  left_join(data_scq) %>%
  mutate(SCI_scale = scale(sqrt(SCI_scale + 1)), #transform because of positive skew
         att_ADHD_scale = scale(att_ADHD_scale), 
         hyp_ADHD_scale = scale(hyp_ADHD_scale),
         RRB_scale = scale(RRB_scale),
         SCQ_scale = scale(SCQ_scale), 
         lang_scale = scale(sqrt(lang_scale + 1)),
         motor_scale = scale(sqrt(motor_scale + 1)),
         rec_lang_scale = case_when(rec_lang_scale == 0 ~ 0, # collapse top categories due to few endorsements, interpret as reporting any rec lang / gross motor difficulties vs no reported difficulties
                                    rec_lang_scale > 0 ~ 1),
         gross_motor_scale = case_when(gross_motor_scale == 0 ~ 0, 
                                       gross_motor_scale > 0 ~ 1), 
         fine_motor_scale = factor(fine_motor_scale, levels = c("0", "1", "2", "3", "4"), ordered = TRUE), 
         exp_lang_scale = factor(as.character(exp_lang_scale), levels = c("0", "1", "2", "3", "4", "5", "6"), ordered = TRUE),
         sex = case_when(sex == "Male" ~"Male", 
                         sex == "Female" ~ "Female"))

data_pgs <- readRDS("N:/durable/projects/neurodev_general_genetic/data/data_pgs_ldpred.rds")

#adding age covariate 

pheno_data_root_dir <- "N:/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/"
age <- read.spss(paste0(pheno_data_root_dir, "PDB2306_Q6_3yrs_v12.sav"), to.data.frame = TRUE) %>%
mutate(age = round(((ALDERRETUR_S6) / 365),2)) %>% 
  select(age, preg_id = PREG_ID_2306, BARN_NR)


data_pheno <- data_pheno %>%
  merge.data.frame(age, by.x = c("preg_id", "BARN_NR"))

#formatting prs data for analyses - one row per trio 
child_prs <- data_pgs %>%
  filter(!is.na(preg_id)) %>%
  select(preg_id,BARN_NR, ends_with("_pgs_res")) %>% 
  rename_with(., .fn = ~paste0(., ".C"), .cols =  ends_with("_pgs_res"))

m_prs <- data_pgs %>%
  filter(!is.na(m_id)) %>%
  select(m_id, ends_with("_pgs_res")) %>% 
  rename_with(., .fn = ~paste0(., ".M"), .cols = ends_with("_pgs_res"))

f_prs <- data_pgs %>%
  filter(!is.na(f_id)) %>%
  select(f_id, ends_with("_pgs_res")) %>% 
  rename_with(., .fn = ~paste0(., ".F"), .cols = ends_with("_pgs_res"))


dat <- data_pheno %>%
  mutate(BARN_NR = as.character(BARN_NR)) %>%
  left_join(., child_prs, by = c("preg_id", "BARN_NR")) %>%
  left_join(., m_prs, by = "m_id") %>%
  left_join(., f_prs, by = "f_id")

#complete trios 
data <- dat %>% 
  filter_at(vars(contains("_pgs_res")), all_vars(!is.na(.))) %>%
  filter_at(vars(contains("_scale")), any_vars(!is.na(.))) 
  
multi <- data %>% #twins/triplets IDs 
  group_by(preg_id) %>%
  filter(n() > 1) %>%
  ungroup()

data <- data %>% #keep one child from pregnancies with multiples
  distinct(preg_id, .keep_all = TRUE) %>%
  filter(age < 4) #remove age outliers
  

nrow(data) #total

#incomplete trios
 
 dat %>% 
  filter_at(vars(contains("_scale")), any_vars(!is.na(.))) %>%
  nrow()

#creating data set with no siblings  
length(unique(data$m_id)) 
data_nosib <- data %>%
  distinct(m_id, .keep_all = TRUE)

#defining terms for models
outcomes <- data %>% 
  select(ends_with("_scale")) %>% 
  select(-c(ADHD_scale, exp_lang_scale, gross_motor_scale, rec_lang_scale, fine_motor_scale)) %>% #only subscale, remove ones that are not linear regression
  colnames()

outcomes_sens <- data  %>% 
  select(c(SCQ_scale)) %>% #SCQ total score sensitivity analyses
  colnames()

pgs <- data %>%
  select(contains("_pgs_res")) %>%
  colnames(.)

covariates <- c("sex", "age") 




#defining functions for multiple regression analysis for use with map()

makeRegression <- function(outcome, predictors, covars) { #function to make regression model syntax 
  formula <- as.formula(
    paste(outcome, 
          paste(c(predictors, covariates), collapse = " + "), 
          sep = " ~ "))
  return(formula)
} 

runRegression <- function(model, data){ #function that runs regression to be able to easily use with map()
  model_obj <- eval(bquote(lm(.(model), data = data)))  
  return(model_obj)
}

coefTable <- function(results_obj, pgs, data){ #creates dataframe with coefficients from a smmry.lm object
  df <- as.data.frame(summary(results_obj)[["coefficients"]]) %>% #create df
    rownames_to_column(var = "predictor") %>%
    rename(p.value = "Pr(>|t|)", 
           std.error = "Std. Error", 
           t.value = "t value") 
  df <- df %>% 
      mutate(lower_ci = Estimate - 1.96*std.error,
             higher_ci = Estimate + 1.96*std.error,
             partial.R2 = NULL,
             sig_0.05_before.adj = if_else(p.value < 0.05, TRUE, FALSE),
             outcome = as.vector(as.character(summary(results_obj)[["terms"]][[2]]))
             ) %>%
    relocate(c(higher_ci, lower_ci), .after = Estimate) %>%
    select(-c(t.value))
  
  outcome <- results_obj[["terms"]][[2]] #getting correct outcome for partial r2 
    
    return(df)
}

isSig <- function(table, col) { 
  df <- table %>%
    filter(predictor != "(Intercept)") %>%
    filter(get(col) == TRUE)
  names <- df$predictor
  return(names)
}
#running analysis
#first: linear multiple regression for scales that meet assumptions interpret r2 only

#all

#predictors for direct effects (child)
pgs_d = data %>%
  select(contains(".C")) %>%
  colnames(.) 

#predictors for indirect effects (mother and father)
pgs_i = data %>%
  select(contains(".M"),contains(".F")) %>%
  colnames(.) 


results_regression_df <- tribble(~outcome, ~model, #defining df col to hold outcome name and models
                      outcomes,map(outcomes, makeRegression, predictors = pgs, covars = covariates)) %>% #making model for each outcome
  unnest(cols = c(outcome,model)) %>% #unnesting lists so each item is a row
  mutate_at("outcome", str_replace,".outcome", "") %>%
  mutate(result_obj = map(model, runRegression, data = data), #running models
         r2_direct = map_dbl(result_obj, group_partial_r2, covariates = pgs_d),
         r2_indirect = map_dbl(result_obj, group_partial_r2, covariates = pgs_i),
         percent_indirect = r2_indirect*100,
         percent_direct = r2_direct*100,
         ratio_dir_indirect = round((percent_direct / percent_indirect), 2)
         ) %>%
  mutate_if(is.list,setNames, .$outcome) #naming list cols to match outcome 

#writing resutls 
writexl::write_xlsx(results_regression_df, "N:/durable/projects/neurodev_general_genetic/results/partialr2_ldpred_full.xlsx") 

#creating plots to check assumptions 

plots <- function(result_obj, outcome){
  
  label <- paste(outcome)
  plot(result_obj, main = paste0(label))
}

par(mfrow = c(2, 2))

map2(results_regression_df$result_obj, results_regression_df$outcome, plots)


#sensitivity analysis running above models in data with no siblings 
results_regression_df_nosib <- tribble(~outcome, ~model, #defining df col to hold outcome name and models
                                 outcomes,map(outcomes, makeRegression, predictors = pgs, covars = covariates)) %>% #making model for each outcome
  unnest(cols = c(outcome,model)) %>% #unnesting lists so each item is a row
  mutate_at("outcome", str_replace,".outcome", "") %>%
  mutate(result_obj = map(model, runRegression, data = data_nosib), #running models
         r2_direct = map_dbl(result_obj, group_partial_r2, covariates = pgs_d),
         r2_indirect = map_dbl(result_obj, group_partial_r2, covariates = pgs_i),
         percent_indirect = r2_indirect*100,
         percent_direct = r2_direct*100,
         ratio_dir_indirect = round((percent_direct / percent_indirect), 2)
  ) %>%
  mutate_if(is.list,setNames, .$outcome) #naming list cols to match outcome 

#writing data
writexl::write_xlsx(results_regression_df_nosib, "N:/durable/projects/neurodev_general_genetic/results/partialr2_ldpred_nosib.xlsx") 


###### regressions for outcome on each PGS trait separately ########

#setting up dataframe for map() to run regressions each row is a separate regression model (i.e. SCI_scale ~ autism_pgs_res.C + autism_pgs_res.F + autism_pgs_res.M)
n_pgstrait <- 5
n_outcomes <- 11

#listing outcomes
outcomes_2 <- data %>% 
  select(ends_with("_scale")) %>% 
  select(-c(ADHD_scale)) %>% 
  colnames() 

outcomes_2 <- unlist(lapply(outcomes_2 , function(x) rep(x, n_pgstrait)))

#repeating 5 different regression models for each outcome 
PGS <- rep(c("autism_pgs_res.C + autism_pgs_res.F + autism_pgs_res.M", 
         "ADHD_pgs_res.C + ADHD_pgs_res.F + ADHD_pgs_res.M",
         "EA_pgs_res.C + EA_pgs_res.F + EA_pgs_res.M", 
         "cog_pgs_res.C + cog_pgs_res.F + cog_pgs_res.M",
         "dyslexia_pgs_res.C + dyslexia_pgs_res.M + dyslexia_pgs_res.F"), n_outcomes)
 
#setting up df for runinng resgression 
results_regression_indv <-  tribble(~outcome, ~predictor, #defining df col to hold outcome name and predictor
                                    outcomes_2, PGS) %>%
                                      unnest(cols = c(outcome ,predictor)) %>%
  mutate(model = paste0(outcome," ~ ", predictor, " + sex + age"),
         predictor = case_when(grepl("autism", predictor) == TRUE ~ "autism pgs", 
                               grepl("ADHD", predictor) == TRUE ~ "ADHD pgs",
                               grepl("dyslexia", predictor) == TRUE ~ "Dyslexia pgs",
                               grepl("cog", predictor) == TRUE ~ "Intell pgs",
                               grepl("EA", predictor) == TRUE ~ "EA pgs"),
         reg_type = case_when(grepl("fine_motor|exp_lang", outcome)== TRUE ~ "Ordinal", 
                              grepl("gross_motor|rec_lang", outcome)== TRUE ~ "Logistic",
                              TRUE ~ "Linear")) 
 

#functions for map() to run models 

#function to run regression - run correct regression depending on in a linear/logistic/or ordinal is needed 

runRegression_2 <- function(model, type = c("Linear", "Logistic", "Ordinal"), data){ 
  if(type == "Linear") {
    model_obj <- eval(bquote(lm(.(model), data = data))) 
  } 
  
  if(type == "Logistic"){
    model_obj <- eval(bquote(glm(.(model), data = data, family = "binomial")))
    } 
  
  if(type == "Ordinal") {
    model_obj <- eval(bquote(MASS::polr(.(model), data = data, Hess=TRUE)))
  }
  
  return(model_obj)
}

#new coeftable
 coefTable_2 <- function(est_clustered,results_obj){
   df <- as.data.frame(summary(results_obj)[["coefficients"]])%>% #create df
     rownames_to_column(var = "predictor") %>%
     rename(std.error = "Std. Error") %>%
     filter(grepl("pgs|sex|(Intercept)|age", predictor)) %>%
     mutate(Estimate= est_clustered[1:nrow(est_clustered),1], 
            std.error = est_clustered[1:nrow(est_clustered),2], 
            torz.value = est_clustered[1:nrow(est_clustered),3],
            p.value = est_clustered[1:nrow(est_clustered),4])
   
   row_df <- nrow(df)
   
 
   
    df <- df %>% #add cols if values reach a sig threshold 
     mutate(lower_ci = Estimate - 1.96*std.error,
            higher_ci = Estimate + 1.96*std.error,
            sig_0.05_before.adj = if_else(p.value < 0.05, TRUE, FALSE),
            outcome = as.vector(as.character(summary(results_obj)[["terms"]][[2]]))
     ) %>%
     relocate(c(higher_ci, lower_ci), .after = Estimate) %>%
     select(c(predictor, Estimate, higher_ci, lower_ci, std.error,torz.value, p.value,  starts_with("sig"), outcome))
   
   }

#running regression models and estimateding clustered SE by m_id 
results_regression_indv <- results_regression_indv %>%
  mutate(results_obj = future_map2(model, reg_type, runRegression_2, data = data),
         est_clustered = map(results_obj, coeftest, vcov = vcovCL, cluster = data$m_id))%>%
  mutate_if(is.list,setNames, .$outcome)


results_regression_indv <- results_regression_indv %>%
  mutate(coef_table = map2(est_clustered, results_obj, coefTable_2), #generate coef results table
         sig_0.05_before.adj = map(coef_table, isSig, "sig_0.05_before.adj"))


# figures and tables to interpret results

#scq sensitivity 
scq_full_table = results_regression_indv %>%
  filter(outcome == "SCQ_scale") %>%
  .$coef_table %>%
  bind_rows() %>%
  filter(grepl("pgs", predictor))  %>%
  separate(predictor, into = c("predictor", "effect"), sep = "_pgs_res.") %>%
  group_by(effect) 

writexl::write_xlsx(scq_full_table, "N:/durable/projects/neurodev_general_genetic/results/regression_beta_SCQfull_ldpred.xlsx")


#results for main analyses 

ggplot_table <- results_regression_indv %>%
  filter(outcome != "SCQ_scale")  


ggplot_table <- ggplot_table$coef_table %>%
 bind_rows() %>%
  filter(grepl("pgs", predictor))  %>%
  separate(predictor, into = c("predictor", "effect"), sep = "_pgs_res.") %>%
  group_by(effect) %>%
  mutate(p.value.adj_new = p.adjust(p.value, method = "fdr", n = 50)) %>%
  ungroup() %>%
  mutate(sig_new = p.value.adj_new < 0.05) %>%
  mutate(predictor = case_when(predictor == "dyslexia" ~ "Dyslexia",
                               predictor == "EA" ~ "EA", 
                               predictor == "cog" ~ "Cognitive Ability", 
                               predictor == "autism" ~ "Autism", 
                               predictor == "ADHD" ~ "ADHD", 
                            TRUE ~ predictor), 
         outcome = case_when(outcome == "att_ADHD_scale" ~ "Inattention", 
                              outcome == "exp_lang_scale" ~ "Expressive language",
                              outcome == "fine_motor_scale" ~ "Fine motor",
                              outcome == "gross_motor_scale" ~ "Gross motor",
                              outcome == "hyp_ADHD_scale" ~ "Hyperactivity",
                              outcome == "rec_lang_scale" ~ "Receptive language" ,
                              outcome == "RRB_scale" ~ "RRBI",
                              outcome == "SCI_scale" ~ "Social & communication", 
                              outcome == "motor_scale" ~ "Motor", 
                             outcome == "lang_scale" ~ "Language"), 
         outcome = ordered(outcome, levels = c("Inattention", "Hyperactivity","Expressive language",
                                               "Receptive language","Language", "Motor", "Fine motor","Gross motor",
                                               "RRBI", "Social & communication")), 
         scale = case_when(outcome == "Inattention" | outcome == "Hyperactivity" ~ "CBCL-ADHD", 
                           outcome == "Fine motor" | outcome == "Gross motor" | outcome == "Motor" ~ "ASQ-Motor",
                           outcome == "Expressive language" | outcome == "Receptive language"| outcome == "Language" ~ "ASQ-Communication", 
                           outcome == "RRBI" | outcome == "Social & communication" ~ "SCQ scale"),
         effect = case_when(effect == "C" ~ "Direct", 
                            effect == "M" ~ "Indirect - Maternal", 
                            effect == "F" ~ "Indirect - Paternal"),
         label = case_when(
           p.value.adj_new > 0.05 ~ "",
           p.value.adj_new > 0.01 ~ "*",
           p.value.adj_new > 0.001 ~ "**",
           !is.na(p.value.adj_new) ~ "***",
           TRUE ~ NA_character_
         ))

library(xlsx)
writexl::write_xlsx(ggplot_table, "N:/durable/projects/neurodev_general_genetic/results/regression_beta_indv_ldpred.xlsx")

#log estimates
ggplot_table %>%
  filter(!(outcome %in% c("Inattention","Hyperactivity", "RRBI","Language", "Motor","Social & communication"))) %>%
  filter(predictor %in% c("ADHD","Autism", "EA", "Cognitive Ability", "Dyslexia")) %>%
  ggplot() +
  geom_pointrange(aes(
    y = exp(Estimate),
    ymin = exp(lower_ci),
    ymax = exp(higher_ci), 
    x = predictor, 
    group = effect, 
    shape = effect,
    color = predictor, 
   alpha = sig_new
   ), 
  size = 1.4,
  position = position_dodge(width = 1)) + 
  scale_alpha_discrete(range = c(0.3,1)) +
  facet_grid(rows = vars(predictor), cols = vars(outcome), scales = "free",  space = "free") +
  geom_text(aes(label = label, x = predictor, group = effect, y = exp(Estimate)), vjust = -0.25, size = 8, position = position_dodge(width = 1)) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  scale_y_continuous(limits = c(0.8, 1.4),breaks=seq(0.8, 1.4, 0.2), trans = "log2", expand = c(-0.05,0.2)) +
  #scale_x_discrete(expand=c(0.1, 0))
  coord_flip() +
  theme_bw() +
  theme(strip.text.y = element_blank(),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"),
        panel.border = element_rect(linetype = "dashed", color = "grey"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 21, color = "black"),
        axis.text.x = element_text(size =18, color = "black", face = "bold"),
        axis.title = element_text(size = 24),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 1),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  scale_color_aaas() +
  labs(y = "Estimate (Odds Ratio)",
       title = "Indirect and Direct Effects of Neurodevelopmental Conditions and Related Triats PGS on Subdomains of Language and Motor Development",
       color = "Polygenic score",
       x = "Polygenic score", 
       shape = "Genetic Effect") +
  guides(alpha = "none") 
  
#beta estimates related traits
ggplot_table %>%
  filter(outcome %in% c("Inattention","Hyperactivity", "RRBI", "Social & communication", "Language", "Motor")) %>%
  filter(predictor %in% c("EA", "Cognitive Ability")) %>%
  ggplot() +
  geom_pointrange(aes(
    ymin = lower_ci,
    ymax = higher_ci, 
    x = predictor, 
    y = Estimate, 
    group = effect, 
    shape = effect,
    color = predictor), 
    size = 1.4,
    position = position_dodge(width = 1)) + 
  geom_text(aes(label = label, x = predictor, group = effect, y = Estimate), vjust = -0.25, size = 8, position = position_dodge(width = 1)) +
  facet_grid(rows = vars(predictor), cols = vars(outcome), scales = "free",  space = "free") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  scale_x_discrete(expand=c(0.1, 0))+
  coord_flip() +
  theme_bw() +
  theme(strip.text.y = element_blank(),
        strip.text.x = element_text(size = 16, face = "bold", color = "black"),
        panel.border = element_rect(linetype = "dashed", color = "grey"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size =18, color = "black"),
        axis.title = element_text(size = 24),
        #axis.title.x = element_text(vjust = -10),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  scale_color_d3() +
  labs(y = "Estimate (Std. Beta)",
       title = "Direct and Indirect Genetic Effects of Educational Attainment and Cognitive Ability PGS on Early Neurodevelopmental Traits",
       color = "Polygenic score",
       x = "Polygenic score", 
       shape = "Genetic Effect") +
  guides(alpha = "none"#, 
         #shape = "none", 
         #color = "none"
         )

###beta estimates ND cond 

ggplot_table %>%
  filter(outcome %in% c("Inattention","Hyperactivity", "RRBI", "Social & communication", "Language", "Motor")) %>%
  filter(predictor %in% c("ADHD","Autism")) %>%
  ggplot() +
  geom_pointrange(aes(
    ymin = lower_ci,
    ymax = higher_ci, 
    x = predictor, 
    y = Estimate, 
    group = effect, 
    shape = effect,
    color = predictor), 
    size = 1.4,
    position = position_dodge(width = 1)) + 
  geom_text(aes(label = label, x = predictor, group = effect, y = Estimate), vjust = -0.25, size = 8, position = position_dodge(width = 1)) +
  facet_grid(rows = vars(predictor), cols = vars(outcome), scales = "free",  space = "free") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_y_continuous(limits = c(-0.08, 0.08)) +
  scale_x_discrete(expand=c(0.1, 0))+
  coord_flip() +
  theme_bw() +
  theme(strip.text.y = element_blank(),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"),
        panel.border = element_rect(linetype = "dashed", color = "grey"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.y = element_text(size = 21, color = "black"),
        axis.text.x = element_text(size =18, color = "black", face = "bold"),
        axis.title = element_text(size = 24),
        #axis.title.x = element_text(vjust = -10),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  scale_color_aaas() +
  labs(y = "Estimate (Std. Beta)",
       title = "Direct and Indirect Genetic Effects of ADHD and Autism PGS on Early Neurodevelopmental Traits",
       color = "Polygenic score",
       x = "Polygenic score", 
       shape = "Genetic Effect") +
  guides(alpha = "none"#, 
         #shape = "none", 
         #color = "none"
  )





  

