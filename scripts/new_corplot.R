
corplot_function <- function(model, flip = FALSE){
  data <- model$data
  coef <- exp(coef(model))
  xaxis <- model$g.levels.f[[1]][2]
  yaxis <- model$g.levels.f[[1]][1]
  
  xname <- paste0("yi_", xaxis)
  yname <- paste0("yi_", yaxis)
  
  
  matrixflip <- function(matrixin){
    newmatrix <- matrix(c(matrixin[4], matrixin[3], matrixin[2], matrixin[1]), nrow = 2, byrow = TRUE)
    return(newmatrix)
  }
  
  data2 <- model$data %>%
    select(study, yi, vi, measure) %>%
    pivot_wider(id_cols = study, names_from = measure,
                values_from = c(yi, vi), 
                names_glue = "{.value}_{measure}") %>%
    mutate(covar = vector("list", length = nrow(.)))
  
  data2$covar <- blsplit(model$V, cluster = data$study)
  
  if(flip == TRUE){
    data2 <- data2 %>%
      mutate(covar2 = map(covar, list_to_matrix)) %>%
      mutate(covar2 = map(covar2, matrixflip)) %>%
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
      mutate(ell = map2(.x = covar2, .y = centre, ~ellipse(matrix(unlist(.x), nrow = 2), centre = unlist(.y)))) %>%
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
