
singletriple <- within_data %>% 
  group_by(study) %>% 
  mutate(measure = case_when(
    str_detect(measure, "isk con 30|isk con 60|isk con 90") ~ "quad_slow",
    TRUE ~ measure
  )) %>%
 filter(all(c("single hop", "quad_slow") %in% measure)) %>% 
  filter(measure %in% c("single hop", "quad_slow")) %>% 
  ungroup



bivariate_regress <- function(data, var1, var2, rho = 0.8){
  
   data_filter <- data %>%
    group_by(study) %>% 
    filter(all(c(var1, var2) %in% measure)) %>% 
    filter(measure %in% c(var1, var2)) %>% 
    ungroup %>%
    group_by(study, measure) %>%
    slice(which.min(abs(timepoint_mean - as.numeric(12)))) %>% # take one effect per study, preferably closest to 12 months if multiple
    ungroup() %>%
    filter(timepoint_mean <24) %>%
    mutate(measure = factor(measure, levels = c(var1, var2)))
   
   V <- vcalc(vi,
              cluster = study,
              type = measure,
              rho = rho,
              data = data_filter,
              checkpd = TRUE)
   
   mv <- rma.mv(yi,
                    V,
                    mods = ~ measure - 1, # remove intercept to generate estimate for each.
                    random = list(~ measure | cohort),
                    struct = c("UN"),
                    data = data_filter,
                    cvvc = "varcov",
                    control=list(nearpd=TRUE)) # need this to use  later with matreg
   
  return(mv)
  
}

singletriple <- within_data %>% 
  group_by(study) %>% 
  filter(all(c("single hop", "triple hop") %in% measure)) %>% 
  filter(measure %in% c("single hop", "triple hop")) %>% 
  ungroup

  
  str_detect(measure, "isk con 30|isk con 60|isk con 90") ~ "Slow Isokinetic"
  
singletriple <- singletriple %>%
  group_by(study, measure) %>%
  slice(which.min(abs(timepoint_mean - as.numeric(12)))) %>% # take one effect per study, preferably closest to 12 months if multiple
  ungroup() %>%
  filter(timepoint_mean <24)

singletripleV <- vcalc(vi,
                     cluster = study,
                     type = measure,
                     rho = 0.9,
                     data = singletriple,
                     checkpd = TRUE)


singletriple_mv <- rma.mv(yi,
                         singletripleV,
                         mods = ~ measure - 1, # remove intercept to generate estimate for each.
                         random = list(~ measure | cohort),
                         struct = c("UN"),
                         data = singletriple,
                         cvvc = "varcov",
                         control=list(nearpd=TRUE)) # need this to use  later with matreg


corplot_function <- function(model){
  data <- model$data
  coef <- exp(coef(model))
  
  reg <- matreg(y = 2, x = 1, R = model$G, cov = TRUE, means = coef, V = model$vvc)
  ci <- as.data.frame(ellipse(model$vb, centre = coef, level = 0.95)) %>% rename(x = 1, y = 2)
  pred <- as.data.frame(ellipse(model$G, centre = coef, level = 0.95)) %>% rename(x = 1, y = 2)
  
  xaxis <- model$g.levels.f[[1]][1]
  yaxis <- model$g.levels.f[[1]][2]
  
  plot <- data %>% 
    select(study, measure, yi, vi) %>% # select only the estimates
    pivot_longer(-c(study, measure), names_to = "var", values_to = "val") %>%
    mutate(measure = as.numeric(factor(measure))) %>%
    pivot_wider(id_cols = c(study), names_from = c("var", measure), values_from = "val") %>%
    ggplot(aes(x = exp(yi_1), y = exp(yi_2))) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_abline(intercept = 0.25, slope = 0.75, linetype = "dotted", colour = "grey") +
    geom_abline(intercept = 0.5, slope = 0.5, linetype = "dotted", colour = "grey", alpha = 0.7) +
    geom_vline(xintercept = 1, colour = "grey") + 
    geom_hline(yintercept = 1, colour = "grey") +
    geom_point() +
    geom_abline(intercept = reg$tab$beta[1], slope = reg$tab$beta[2], colour = "red") +
    geom_point(x = coef[1], y = coef[2], inherit.aes = FALSE, colour = "red", size = 3) +
    geom_polygon(data = ci, mapping = aes(x = x, y = y), alpha = 0.15, inherit.aes = FALSE, fill = "red") +
    geom_polygon(data = pred, mapping = aes(x = x, y = y), alpha = 0.15, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0.6, 1.1), ylim = c(0.6, 1.1), clip = "off") +
    scale_x_continuous(breaks = c(0.7, 0.8, 0.9, 1.0), labels = c("30%", "20%", "10%", "0%")) +
    scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 1.0), labels = c("30%", "20%", "10%", "0%")) +
    labs(x = xaxis, y = yaxis) +
    annotate(geom = "text", x = 0.6, y = 0.6, label = "1x", family = "Karla", colour = "grey", vjust = 0, angle = 45) + 
    annotate(geom = "text", x = 0.6, y = 0.7, label = "1.5x", family = "Karla", colour = "grey", vjust = 0, angle = 37) + 
    annotate(geom = "text", x = 0.6, y = 0.8, label = "2x", family = "Karla", colour = "grey", vjust = 0, angle = 30) +
    theme_mgpub() +
    theme(panel.grid.major.x = element_line(colour = "#F1F1F1"),
          panel.grid.major.y = element_line(colour = "#F1F1F1"))
  
  #return(plot)
  return(list(reg, plot))
  
}






# These studies report the raw for injured mean and non-injured, so to calcualte LSI, need to take non-inj / inj
invert6m <- sixmhop_data_pre %>%
  filter(study %in% c("Ebert 2018", "Norte 2020", "Ebert 2019", "Meierbachtol 2017", "Logerstedt 2013", "Moran 2022a",
                      "Sonesson 2021b", "Sonesson 2021a")) %>%
  escalc(ni = acl_n, m1i = noninj_mean, sd1i = noninj_sd, m2i = inj_mean, sd2i = inj_sd, ri = ri, 
         data = ., measure = "ROMC") %>%
  mutate(yi = case_when(
    is.na(yi) & !is.na(lsi_mean) ~ log(lsi_mean),
    TRUE ~ yi),
    vi = case_when(
      is.na(vi) & !is.na(lsi_mean) ~ lsi_sd^2/(acl_n * lsi_mean^2),
      TRUE ~ vi))

# These studies only report LSI, take the inverse of the LSI
invert6mlsi <- sixmhop_data_pre %>%
  filter(study %in% c("Norte 2019b", "Arhos 2020a", "Barnett 2020", "Reinke 2011")) %>%
  select(-c(yi:ci.ub)) %>%
  mutate(lsi_mean = 1/lsi_mean) %>%
  mutate(yi = log(lsi_mean),
         vi = lsi_sd^2/(acl_n * lsi_mean^2))

# Join back Together
sixmhop_data <- bind_rows(sixmhop_data_pre %>% filter(!study %in% c("Ebert 2018", "Norte 2020", "Ebert 2019", "Meierbachtol 2017", "Logerstedt 2013", "Moran 2022a",
                                                                    "Sonesson 2021b", "Sonesson 2021a", "Norte 2019b", "Arhos 2020a", "Barnett 2020", "Reinke 2011")),
                          invert6m) %>%
  bind_rows(., invert6mlsi) 

