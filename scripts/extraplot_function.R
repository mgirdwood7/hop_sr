
mv_plotdetails <- function(model, logscale = FALSE, xlimit = 120, include_pi = FALSE, showgraft = NULL, linerange = FALSE){
  
  knot_pos <- NULL # set to 0 for logic to work later
  
  moderator <- as.character(model$call$mods[2]) # get the moderator variable
  data <- model$data # get the data
  
  
  total_k  <- model$k
  total_studies <- data %>% filter(!is.na(vi)) %>% summarise(n = length(unique(study))) %>% as.numeric()
  total_n <- data %>% group_by(study) %>% arrange(timepoint) %>% slice(1) %>% ungroup()  %>% filter(!is.na(vi)) %>% summarise(n = sum(acl_n)) %>% as.numeric
  
  last_timepoint <- round(max(data$timepoint_mean),0)
  
  data <- data %>%
    mutate(hsgraft = case_when(
      acl_graft_group == "hs" ~ "hs",
      TRUE ~ "other"
    ),
    quadgraft = case_when(
      acl_graft_group == "quad" ~ "quad",
      TRUE ~ "other"
    ))
  
  
  if (str_detect(moderator, "rcs")) { # if a rcs is used, extract the information
    spline <- as.list(model$formula.mods[[2]][[3]]) 
    knot_pos <- unlist(spline[sapply(spline, is.numeric)]) # this gets the vector of knot positions or no. of knots
  }
  
  # depending on the type od moderator (i.e. log, linear, rcs...) calculate the predicted curve points
  mod_type <- if (str_detect(moderator, "log")) {
    
    #log
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods = log(seq(1,120, length = 120)),
                                 transf = exp
    )) %>% mutate(x = row_number())
    
  } else if (str_detect(moderator, "poly")) {
    
    # poly
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods=unname(poly(seq(1,120, length = 120), degree = 2, raw=TRUE)),
                                 transf = exp
    )) %>% mutate(x = row_number())
    
  } else if (length(knot_pos) == 3) {
    
    # 3 knot rcs
    knots <- attr(rcs(model.matrix(model)[,2], knot_pos), "parms")
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods=rcspline.eval(seq(1,120, length = 120), knots, inclx=TRUE),
                                 transf = exp)) %>% 
      mutate(x = row_number())
    
  } else if (length(knot_pos) == 1) { # i.e. cases where only the number of knots is given
    
    # 3 knot rcs
    knots <- attr(rcs(data$timepoint_mean, as.numeric(knot_pos)), "parms") # use model data to get knot positions
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods=rcspline.eval(seq(1,120, length = 120), knots, inclx=TRUE),
                                 transf = exp)) %>% 
      mutate(x = row_number())
    
  } else if (length(knot_pos) == 4) {
    
    # 4 knot rcs
    knots <- attr(rcs(model.matrix(model)[,2], knot_pos), "parms")
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods=rcspline.eval(seq(1,120, length = 120), knots, inclx=TRUE),
                                 transf = exp)) %>% 
      mutate(x = row_number())
    
  } else {
    # Linear
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods = seq(1,120, length = 120),
                                 transf = exp)) %>% mutate(x = row_number())
  }
  
  # calculate ci points
  ci <-  pointsfunction(points)
  pi <-  pi_pointsfunction(points)
  
  
  ci <- ci %>% filter(x < last_timepoint + 6)
  points <- points %>% filter(x < last_timepoint + 6)
  pi <- pi %>% filter(x < last_timepoint + 6)
  
  if (logscale == FALSE & is.null(showgraft) & linerange == FALSE){
    
    # plot data and predicted model
    plot <- data %>%
      ggplot(aes(x = timepoint_mean, y = exp(yi), group = cohort)) +
      geom_hline(yintercept = 1, colour = "dark grey") +
      geom_line(alpha = 0.8, colour = "grey") + 
      geom_point(aes(size = acl_n), alpha = 0.3) + 
      scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2), labels = c("-40%", "-20%", "0%", "+20%")) +
      #scale_y_continuous(labels = scales::percent, limits = c(-0.75, 0.25)) +
      #scale_x_continuous(trans = 'log10', limits = c(3,150)) +
      scale_x_continuous(breaks = c(0, 12, 24, 60, 120), labels = c(0, 1, 2, 5, 10)) +
      #scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120), labels = c("", 1, 2, "", "", 5, "", "", "", "", 10)) +
      coord_cartesian(xlim = c(0, xlimit), ylim = c(0.5, 1.2)) +
      #scale_size(limits = c(1, 4500), trans = custom_trans, range = c(1, 10)) +
      labs(x = "Time since surgery (Years)", y = "Percentage Deficit", colour = "Measure", size = "Participants (n)") +
      #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
      #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
      geom_line(data = points, aes(x = x, y = pred), colour = "red", linewidth = 1.3, inherit.aes = FALSE) +
      geom_polygon(data = ci, aes(x = x, y = ci.ub), fill = "red",  alpha = 0.2, inherit.aes = FALSE) +
      annotate(geom = "text", x = 120, y = 0.62, label = paste0("k = ", total_k, " (", total_studies, " studies)"), family = "Karla", hjust = 1) +
      annotate(geom = "text", x = 120, y = 0.67, label = paste0("n = ", total_n), family = "Karla", hjust = 1) +
      theme_mgpub() +
      theme(panel.grid.major.x = element_line(linewidth = rel(0.5), linetype = 2))
    
  } else if (linerange == TRUE) {
    
    plot <- data %>%
      ggplot(aes(x = timepoint_mean, y = exp(yi), group = cohort)) +
      geom_point(alpha = 0.3) + 
      geom_linerange(data = data, mapping = aes(x = timepoint_mean, ymin = ci.lb, ymax = ci.ub), alpha = 0.3,  inherit.aes = FALSE) +
      geom_line(alpha = 0.1, linewidth = 3) + 
      scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2), labels = c("-40%", "-20%", "0%", "+20%")) +
      #scale_y_continuous(labels = scales::percent, limits = c(-0.75, 0.25)) +
      #scale_x_continuous(trans = 'log10', limits = c(3,150)) +
      coord_cartesian(xlim = c(0,100), ylim = c(0.5, 1.2)) +
      scale_size(range = c(0, 10)) +
      labs(x = "Time since surgery (Months)", y = "Percentage Deficit", colour = "Measure", size = "Participants (n)") +
      #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
      #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
      geom_line(data = points, aes(x = x, y = pred), colour = "red", linewidth = 1.3, inherit.aes = FALSE) +
      geom_polygon(data = ci, aes(x = x, y = ci.ub), fill = "red",  alpha = 0.1, inherit.aes = FALSE) +
      theme_mgpub()
    
  } else if (showgraft == "quad"){
    
    plot <- data %>%
      ggplot(aes(x = timepoint_mean, y = exp(yi), group = cohort, colour = quadgraft)) +
      geom_point(aes(size = acl_n, colour = quadgraft), alpha = 0.3) + 
      scale_fill_manual(values = c("red", "blue")) + 
      geom_line(alpha = 0.8) + 
      scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2), labels = c("-40%", "-20%", "0%", "+20%")) +
      #scale_y_continuous(labels = scales::percent, limits = c(-0.75, 0.25)) +
      #scale_x_continuous(trans = 'log10', limits = c(3,150)) +
      coord_cartesian(xlim = c(0,xlimit), ylim = c(0.5, 1.2)) +
      scale_size(range = c(0, 10)) +
      labs(x = "Time since surgery (Months)", y = "Percentage Deficit", size = "Participants (n)") +
      #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
      #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
      geom_line(data = points, aes(x = x, y = pred), colour = "black", linewidth = 1.3, inherit.aes = FALSE) +
      geom_polygon(data = ci, aes(x = x, y = ci.ub), fill = "black",  alpha = 0.1, inherit.aes = FALSE) +
      theme_mgpub()
    
  } else if (showgraft == "hs"){
    
    plot <- data %>%
      ggplot(aes(x = timepoint_mean, y = exp(yi), group = cohort, colour = hsgraft)) +
      geom_point(aes(size = acl_n, colour = hsgraft), alpha = 0.3) + 
      geom_line(alpha = 0.8) + 
      scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2), labels = c("-40%", "-20%", "0%", "+20%")) +
      #scale_y_continuous(labels = scales::percent, limits = c(-0.75, 0.25)) +
      #scale_x_continuous(trans = 'log10', limits = c(3,150)) +
      coord_cartesian(xlim = c(0,xlimit), ylim = c(0.5, 1.2)) +
      scale_size(range = c(0, 10)) +
      labs(x = "Time since surgery (Months)", y = "Percentage Deficit", colour = "Measure", size = "Participants (n)") +
      #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
      #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
      geom_line(data = points, aes(x = x, y = pred), colour = "red", linewidth = 1.3, inherit.aes = FALSE) +
      geom_polygon(data = ci, aes(x = x, y = ci.ub), fill = "red",  alpha = 0.1, inherit.aes = FALSE) +
      theme_mgpub()
    
  } else {
    
    plot <- data %>%
      ggplot(aes(x = timepoint_mean, y = exp(yi), group = cohort)) +
      geom_hline(yintercept = 1, colour = "dark grey") +
      geom_line(alpha = 0.8, colour = "grey") + 
      geom_point(aes(size = acl_n), alpha = 0.3) + 
      scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2), labels = c("-40%", "-20%", "0%", "+20%")) +
      #scale_y_continuous(labels = scales::percent, limits = c(-0.75, 0.25)) +
      scale_x_continuous(trans = 'log10', limits = c(2,150), breaks = c(3, 12, 36, 120), labels = c(0.25, 1, 3, 10)) +
      #coord_cartesian(xlim = c(0,100), ylim = c(-0.75, 0.25)) +
      scale_size(range = c(0, 10)) +
      labs(x = "Time since surgery (Years)", y = "Percentage Deficit", colour = "Measure", size = "Participants (n)") +
      #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
      #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
      geom_line(data = points, aes(x = x, y = pred), colour = "red", linewidth = 1.3, inherit.aes = FALSE) +
      geom_polygon(data = ci, aes(x = x, y = ci.ub), fill = "red",  alpha = 0.2, inherit.aes = FALSE) +
      theme_mgpub()
    
  } 
  
  if (include_pi == TRUE) {
    plot <- plot +
      geom_polygon(data = pi, aes(x = x, y = pi.ub), fill = "black",  alpha = 0.06, inherit.aes = FALSE)
  }
  
  return(plot)
  
}


custom_trans <- scales::trans_new(
  "custom_trans",
  transform = function(x) ifelse(x <= 75, log(x), x),
  inverse = function(x) ifelse(x <= 75, exp(x), x)
)
