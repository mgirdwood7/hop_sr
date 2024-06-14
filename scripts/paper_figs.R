## Hop figures
## Requires all models to be run/loaded in environment

library(ggpubr)
library(mgfunctions)


custom_trans <- scales::trans_new(
  "custom_trans",
  transform = function(x) ifelse(x <= 75, log(x), x),
  inverse = function(x) ifelse(x <= 75, exp(x), x)
)

singlehop <-  mv_plotdetails(singlehop_mv, include_pi = TRUE) + 
  scale_size(limits = c(1, 5000), 
             trans = custom_trans, 
             range = c(1, 7), 
             breaks = c(50, 100, 500, 1000, 5000)) +
  labs(title = "\nSingle Forward Hop", subtitle = " ") +
  coord_cartesian(xlim = c(0,120), ylim = c(0.6, 1.2)) +
  theme(plot.title = element_text(size = 16, family = "BarlowSemiCondensed-SemiBold"))

triplehop <-  mv_plotdetails(triplehop_mv, include_pi = TRUE) + 
  scale_size(limits = c(1, 5000), 
             trans = custom_trans, 
             range = c(1, 7), 
             breaks = c(50, 100, 500, 1000, 5000)) +
  labs(title = "\nTriple Forward Hop", subtitle = " ") +
  coord_cartesian(xlim = c(0,120), ylim = c(0.6, 1.2)) +
  theme(plot.title = element_text(size = 16, family = "BarlowSemiCondensed-SemiBold"))

triplexhop <-  mv_plotdetails(triplexhop_mv, include_pi = TRUE) + 
  scale_size(limits = c(1, 5000), 
             trans = custom_trans, 
             range = c(1, 7), 
             breaks = c(50, 100, 500, 1000, 5000)) +
  labs(title = "\nTriple Crossover Hop", subtitle = " ") +
  coord_cartesian(xlim = c(0,120), ylim = c(0.6, 1.2)) +
  theme(plot.title = element_text(size = 16, family = "BarlowSemiCondensed-SemiBold"))

sidehop <-  mv_plotdetails(sidehop_mv, include_pi = TRUE) + 
  scale_size(limits = c(1, 5000), 
             trans = custom_trans, 
             range = c(1, 7), 
             breaks = c(50, 100, 500, 1000, 5000)) +
  labs(title = "\nSide Hop", subtitle = " ") +
  coord_cartesian(xlim = c(0,120), ylim = c(0.6, 1.2)) +
  theme(plot.title = element_text(size = 16, family = "BarlowSemiCondensed-SemiBold"))

verticalhop <-  mv_plotdetails(verticalhop_mv, include_pi = TRUE) + 
  scale_size(limits = c(1, 5000), 
             trans = custom_trans, 
             range = c(1, 7), 
             breaks = c(50, 100, 500, 1000, 5000)) +
  labs(title = "\nVertical Hop", subtitle = " ") +
  coord_cartesian(xlim = c(0,120), ylim = c(0.6, 1.2)) +
  theme(plot.title = element_text(size = 16, family = "BarlowSemiCondensed-SemiBold"))

sixmhop <-  mv_plotdetails(sixmhop_mv, include_pi = TRUE) + 
  scale_size(limits = c(1, 5000), 
             trans = custom_trans, 
             range = c(1, 7), 
             breaks = c(50, 100, 500, 1000, 5000)) +
  labs(title = "\nSix-metre Timed Hop", subtitle = " ") +
  coord_cartesian(xlim = c(0,120), ylim = c(0.6, 1.2)) +
  theme(plot.title = element_text(size = 16, family = "BarlowSemiCondensed-SemiBold"))

ggarrange(singlehop, triplehop, triplexhop, sixmhop, sidehop, verticalhop, common.legend = TRUE, legend = "right", nrow = 3, ncol = 2)

ggsave("output/plots/hop_within_all2.png", device = png, width = 9, height = 12)
ggsave("output/plots/hop_within_all2.jpeg", device = jpeg, width = 9, height = 12, dpi = 300)



### Relationships
####


# Use this function for each of the comparisons
bivariate_regress <- function(data, var1, var2, rho = 0.8){
  
  data_filter <- data %>%
    group_by(study) %>% 
    filter(all(c(var1, var2) %in% measure)) %>% # filter studies that use the two measures
    filter(measure %in% c(var1, var2)) %>%  # filter to only those measures
    ungroup %>%
    group_by(study, measure) %>%
    slice(which.min(abs(timepoint_mean - as.numeric(12)))) %>% # take one effect per study, preferably closest to 12 months if multiple
    ungroup() %>%
    filter(timepoint_mean <24,
           !is.na(vi)) %>% # remove the long
    mutate(measure = factor(measure, levels = c(var1, var2)))
  
  # need to create covariance matrix
  V <- vcalc(vi,
             cluster = study,
             type = measure,
             rho = rho, # use specified rho
             data = data_filter, 
             checkpd = TRUE)
  
  mv <- rma.mv(yi,
               V,
               mods = ~ measure - 1, # remove intercept to generate estimate for each type of effect
               random = list(~ measure | cohort),
               struct = c("UN"),
               data = data_filter,
               cvvc = "varcov",
               control=list(nearpd=TRUE)) # need this to use  later with matreg
  
  return(mv)
  
}

list_to_matrix <- function(lst) {
  matrix(unlist(lst), nrow = 2)
}

flipmatrix <- function(matrixin){
  matrixin <- list_to_matrix(matrixin)
  matrixout <- matrix(c(matrixin[4], matrixin[3], matrixin[2], matrixin[1]), nrow = 2, byrow = TRUE)
  return(matrixout)
}


## Function conducts regression of the true effects for within and between person analysis and 
## Plots the result
## Input is a bivariate rma.mv model 

corplot_function <- function(model, flip = FALSE){
  data <- model$data
  coef <- exp(coef(model))
  xaxis <- model$g.levels.f[[1]][2]
  yaxis <- model$g.levels.f[[1]][1]
  
  xname <- paste0("yi_", xaxis)
  yname <- paste0("yi_", yaxis)
  
  
  
  data2 <- model$data %>%
    select(study, yi, vi, measure) %>%
    pivot_wider(id_cols = study, names_from = measure,
                values_from = c(yi, vi), 
                names_glue = "{.value}_{measure}") %>%
    mutate(covar = vector("list", length = nrow(.)))
  
  data2$covar <- blsplit(model$V, cluster = data$study)
  
  if(flip == TRUE){
    data2 <- data2 %>%
      mutate(covar2 = map(covar, flipmatrix)) %>%
      rowwise() %>%
      mutate(centre = list(c(exp(!!rlang::sym(xname)), exp(!!rlang::sym(yname))))) %>%
      ungroup() %>%
      mutate(ell = map2(.x = covar2, .y = centre, ~ellipse(matrix(unlist(.x), nrow = 2), centre = unlist(.y)))) %>%
      ungroup() %>%
      mutate(ellx = map(ell, as.data.frame)) %>%
      unnest(ellx)
  } else {
  data2 <- data2 %>%
    mutate(covar2 = map(covar, list_to_matrix)) %>%
    rowwise() %>%
    mutate(centre = list(c(exp(!!rlang::sym(xname)), exp(!!rlang::sym(yname))))) %>%
    ungroup() %>%
    mutate(ell = map2(.x = covar, .y = centre, ~ellipse(matrix(unlist(.x), nrow = 2), centre = unlist(.y)))) %>%
    ungroup() %>%
    mutate(ellx = map(ell, as.data.frame)) %>%
    unnest(ellx)
  }
  
  reg <- matreg(y = 1, x = 2, R = model$G, cov = TRUE, means = coef, V = model$vvc)
  ci <- as.data.frame(ellipse(model$vb, centre = coef, level = 0.95)) %>% rename(x = 1, y = 2)
  pred <- as.data.frame(ellipse(model$G, centre = coef, level = 0.95)) %>% rename(x = 1, y = 2)
  
  
  
  plot <- data %>% 
    select(study, measure, yi, vi) %>% # select only the estimates
    pivot_longer(-c(study, measure), names_to = "var", values_to = "val") %>%
    mutate(measure = as.numeric(factor(measure))) %>%
    pivot_wider(id_cols = c(study), names_from = c("var", measure), values_from = "val") %>%
    ggplot(aes(x = exp(yi_2), y = exp(yi_1))) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    #geom_abline(intercept = 0.5, slope = 0.5, linetype = "dotted", colour = "grey", alpha = 0.7) +
    geom_vline(xintercept = 1, colour = "grey") + 
    geom_hline(yintercept = 1, colour = "grey") +
    geom_point() +
    geom_polygon(data = data2, aes(x = x, y = y, group = study), inherit.aes = FALSE, alpha = 0.07) +
    geom_abline(intercept = reg$tab$beta[1], slope = reg$tab$beta[2], colour = "red") +
    #geom_point(x = coef[1], y = coef[2], inherit.aes = FALSE, colour = "red", size = 3) +
    #geom_polygon(data = ci, mapping = aes(x = x, y = y), alpha = 0.15, inherit.aes = FALSE, fill = "red") +
    #geom_polygon(data = pred, mapping = aes(x = x, y = y), alpha = 0.15, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0.6, 1.1), ylim = c(0.6, 1.1), clip = "off") +
    scale_x_continuous(breaks = c(0.7, 0.8, 0.9, 1.0), labels = c("30%", "20%", "10%", "0%")) +
    scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 1.0), labels = c("30%", "20%", "10%", "0%")) +
    labs(x = xaxis, y = yaxis) +
    #annotate(geom = "text", x = 0.6, y = 0.6, label = "1x", family = "Karla", colour = "grey", vjust = 0, angle = 45) + 
    #annotate(geom = "text", x = 0.6, y = 0.8, label = "2x", family = "Karla", colour = "grey", vjust = 0, angle = 30) +
    theme_mgpub() +
    theme(panel.grid.major.x = element_line(colour = "#F1F1F1"),
          panel.grid.major.y = element_line(colour = "#F1F1F1"))
  
  #return(plot)
  return(list(reg, plot))
  
}

within_data <- within_data %>% filter(study != "Raoul 2018") %>%  rowwise() %>% 
  mutate(yi = case_when(
    study == "Ebert 2018" & measure %in% c("6m timed hop", "triple hop", "triple crossover hop") ~ log(lsi_mean),
    TRUE ~ yi),
    vi = case_when(
      study == "Ebert 2018" & measure %in% c("6m timed hop", "triple hop", "triple crossover hop") ~ lsi_sd^2/(acl_n * lsi_mean^2),
      TRUE ~ vi)
  ) %>% ungroup

singletriple <- bivariate_regress(within_data, "single hop", "triple hop", rho = 0.7)
singletripleplot <- corplot_function(singletriple)[[2]] +
  labs(x = "Triple Forward\nDeficit (%)", y = "Single Forward\nDeficit (%)", title = "Single Forward  ~ Triple Forward", subtitle = " ") +
  theme(plot.title = element_text(size = 14, family = "BarlowSemiCondensed-SemiBold"))


singletriplex <- bivariate_regress(within_data, "single hop", "triple crossover hop", rho = 0.7)
singletriplexplot <- corplot_function(singletriplex)[[2]] +
  labs(x = "Triple Crossover\nDeficit (%)", y = "Single Forward\nDeficit (%)", title = "Single Forward ~ Triple Crossover", subtitle = " ") +
  theme(plot.title = element_text(size = 14, family = "BarlowSemiCondensed-SemiBold"))

triples <- bivariate_regress(within_data, "triple hop", "triple crossover hop", rho = 0.8)
triplesplot <- corplot_function(triples)[[2]] +
  labs(y = "Triple Forward\nDeficit (%)", x = "Triple Crossover\nDeficit (%)", title = "Triple Forward ~ Triple Crossover", subtitle = " ") +
  theme(plot.title = element_text(size = 14, family = "BarlowSemiCondensed-SemiBold"))

sixmhopdata <- bind_rows(sixmhop_data, within_data %>% filter(study %in% sixmhop_data$study, measure == "single hop"))
singlesix <- bivariate_regress(sixmhopdata, "single hop", "6m timed hop", rho = 0.6)
singlesixplot <- corplot_function(singlesix)[[2]] +
  labs(y = "Single Forward\nDeficit (%)", x = "Six-metre Timed\nDeficit (%)", title = "Single Forward ~ Six-metre Timed", subtitle = " ") +
  theme(plot.title = element_text(size = 14, family = "BarlowSemiCondensed-SemiBold"))

ggarrange(singletripleplot, singletriplexplot, triplesplot, singlesixplot, nrow = 2, ncol = 2)
ggsave("output/plots/forwardhop.png", device = png, width = 8, height = 8)
ggsave("output/plots/forwardhop.jpeg", device = jpeg, width = 8, height = 8, dpi = 300)


triplex6mdata <- bind_rows(sixmhop_data, within_data %>% filter(study %in% sixmhop_data$study, measure == "triple crossover hop"))
triplex6m <- bivariate_regress(triplex6mdata, "triple crossover hop", "6m timed hop", rho = 0.6)
triplex6mplot <- corplot_function(triplex6m)[[2]] +
  labs(y = "Triple Crossover\nDeficit (%)", x = "Six-metre Timed\nDeficit (%)", title = "Triple Crossover ~ Six-metre Timed", subtitle = " ") +
  theme(plot.title = element_text(size = 14, family = "BarlowSemiCondensed-SemiBold"))

triple6mdata <- bind_rows(sixmhop_data, within_data %>% filter(study %in% sixmhop_data$study, measure == "triple hop"))
triple6m <- bivariate_regress(triple6mdata, "triple hop", "6m timed hop", rho = 0.6)
triple6mplot <- corplot_function(triple6m)[[2]] +
  labs(y = "Triple Forward\nDeficit (%)", x = "Six-metre Timed\nDeficit (%)", title = "Triple Forward ~ Six-metre Timed", subtitle = " ") +
  theme(plot.title = element_text(size = 14, family = "BarlowSemiCondensed-SemiBold"))


singlevertical <- bivariate_regress(within_data, "single hop", "vertical hop", rho = 0.6)
singleverticalplot <- corplot_function(singlevertical, flip = TRUE)[[2]] +
  labs(y = "Single Forward\nDeficit (%)", x = "Vertical Hop\nDeficit (%)", title = "Single Forward ~ Vertical Hop", subtitle = " ") +
  theme(plot.title = element_text(size = 14, family = "BarlowSemiCondensed-SemiBold"))

singleside <- bivariate_regress(within_data, "single hop", "side hop", rho = 0.6)
singlesideplot <- corplot_function(singleside)[[2]] +
  labs(y = "Single Forward\nDeficit (%)", x = "Side Hop\nDeficit (%)", title = "Single Forward ~ Side Hop", subtitle = " ") +
  theme(plot.title = element_text(size = 14, family = "BarlowSemiCondensed-SemiBold"))

ggarrange(singleverticalplot, singlesideplot, nrow = 1, ncol = 2)
ggsave("output/plots/otherhop.png", device = png, width = 8, height = 4)
ggsave("output/plots/otherhop.jpeg", device = jpeg, width = 8, height = 4, dpi = 300)



ggarrange(triplex6mplot, triple6mplot, singleverticalplot, singlesideplot, nrow = 2, ncol = 2)
ggsave("output/plots/otherhop2.png", device = png, width = 8, height = 8)
ggsave("output/plots/otherhop2.jpeg", device = jpeg, width = 8, height = 8, dpi = 300)



