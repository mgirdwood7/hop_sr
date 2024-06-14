# Functions

# Function for combining groups and adding as new columns
combine_groups <- function(data, x){
  
  # creates the column calls for the formula. Use of sym + {{}} to work with dplyr later
  n1 <- sym("n_1")
  n2 <- sym("n_2")
  m1 <- sym(paste0(x, "_mean_1")) 
  m2 <- sym(paste0(x, "_mean_2"))
  sd1 <- sym(paste0(x, "_sd_1"))
  sd2 <- sym(paste0(x, "_sd_2"))

  # names of new variables (no suffix)
  mean_name = paste0(x, "_mean")
  sd_name = paste0(x, "_sd")
  n_name = "n"
  
  # new data frame
  # need to use walrus := with !!names
  data %>%
    mutate(!!n_name := {{ n1 }} + {{ n2 }},
           !!mean_name := ({{ n1 }} * {{ m1 }} + {{ n2 }} * {{ m2 }})/({{ n1 }} + {{ n2}}),
           !!sd_name := sqrt((({{ n1 }} - 1) * {{ sd1 }}^2 + ({{ n2 }} - 1) * {{ sd2 }}^2 + ((({{ n1 }} * {{ n2 }})/({{ n1 }} + {{ n2 }})) * ({{ m1 }}^2 + {{ m2 }}^2 - 2 * {{ m1 }} *
                                                                                                {{ m2 }})))/({{ n1 }} + {{ n2 }} - 1))) 
}



# Combine means - provide with a vectors of n, means to combine
combine_mean <- function(ncol, meancol){
  n <- ncol
  mean <- meancol
  
out <- sum({{n}}*{{mean}})/sum({{n}})
return(out)
}

# Combine sds - provide with a vectors of n, means, sd to combine
combine_sd <- function(ncol, meancol, sdcol){
  n <- ncol
  mean <- meancol
  sd <- sdcol
  
  mean_combined <- sum({{n}}*{{mean}})/sum({{n}})
  out <- sqrt((sum(({{n}}-1)*{{sd}}^2) + sum({{n}}*({{mean}}-mean_combined)^2))/(sum({{n}})-1))
  return(out)
}
  
  

# function for creating a tidy output for meta, like broom::tidy, as not available for class: meta objects.
tidymeta <- function(x){
  tibble(
    "estimate" = x$TE.random,
    "std.error" = x$seTE.random,
    "ci.lb" = x$lower.random,
    "ci.ub" = x$upper.random,
    "statistic" = x$statistic.random,
    "p.value" = x$pval.random,
    "pred" = x$seTE.predict,
    "pi.lb" = x$lower.predict,
    "pi.ub" = x$upper.predict,
    "i.squared" = x$I2,
    "h" = x$H,
    "tau.squared" = x$tau2,
    "tau.squared.se" = x$se.tau2,
    "cochran.qe" = x$Q,
    "p.value.cochran.qe" = x$pval.Q,
    "df.residual" = x$df.Q,
    "nobs" = x$k)
}


# Function to extract information from notes section of data entry form - (e.g. median, iqr etc)
extract_string <- function(x){
  # locate position of each element
  inj_pos <- str_locate(x, "(?<!-)Inj|(?<!-)inj|(?<!-)Injured") # need to make sure "inj" not preceded by non-inj etc
  noninj_pos <- str_locate(x, "Non-inj|non-inj|Noninj|noninj")
  lsi_pos <- str_locate(x, "LSI|lsi")
  
  # extract string based on end of "inj" string location, and start of next string location
  inj_sub <- if (is.na(lsi_pos[1]) & is.na(noninj_pos[1])) { # if only inj information, then take from end of inj location to end of string
    substring(x, inj_pos[2]+1, nchar(x)) 
  } else if(is.na(noninj_pos[1])) { # if no non-inj, then take from end of inj locaiton to start of lsi
    substring(x, inj_pos[2]+1, lsi_pos[1]-1)
  } else { # otherwise split based on positions as normal
    substring(x, inj_pos[2]+1, noninj_pos[1]-1)
  }
  
  noninj_sub <- if(is.na(lsi_pos[1])) {
    substring(x, noninj_pos[2]+1, nchar(x)) 
  } else {
    substring(x, noninj_pos[2]+1, lsi_pos[1]-1) 
  }
  
  lsi_sub <- substring(x, lsi_pos[2]+1, nchar(x)) 
  
  # remove any punctiation and trim whitespace
  inj_sub <- inj_sub %>% str_replace(",|;", "") %>% str_trim()
  noninj_sub <- noninj_sub %>% str_replace(",|;", "") %>% str_trim()
  lsi_sub <- lsi_sub %>% str_replace(",|;", "") %>% str_trim()
  
  return(data.frame(inj_sub, noninj_sub, lsi_sub))
}


# Similar function to above: extract information from notes section of data entry form - (e.g. median, iqr etc)
extract_string_data <- function(x, name){
  # locate position of each element
  med_pos <- str_locate(x, "MD|Median|median") # 
  iqr_pos <- str_locate(x, "IQR|iqr")
  range_pos <- str_locate(x, "Range|range")
  
  
  med <- if (is.na(iqr_pos[1]) & is.na(range_pos[1])) {
      substring(x, med_pos[2]+1, nchar(x))
    } else if (is.na(iqr_pos[1])) {
      substring(x, med_pos[2]+1, range_pos[1]-1)
    } else {
      substring(x, med_pos[2]+1, iqr_pos[1]-1)
    }
  
  iqr <- if (is.na(range_pos[1])) {
      substring(x, iqr_pos[2]+1, nchar(x))
    } else {
      substring(x, iqr_pos[2]+1, range_pos[1]-1)
    } 
  
  range <- substring(x, range_pos[2]+1, nchar(x)) 
  
  # remove any punctuation and symbols and remove all white space
  med <- med %>% str_replace(",|;|\\(|\\)", "") %>% str_trim()
  iqr <- iqr %>% str_replace(",|;|\\(|\\)", "") %>% str_replace("to", "-") %>% str_replace_all(" ", "")
  range <- range %>% str_replace(",|;|\\(|\\)", "") %>% str_replace("to", "-") %>% str_replace_all(" ", "")
  
  return(data.frame(med, iqr, range) %>% rename_with(., ~paste0(name, .x))) # paste name of data before each name
}


# Plot function

mv_plot <- function(mv, data, type){
  
  if(type == "linear") {points <- data.frame(predict(robust(mv, cluster = cohort, clubSandwich = TRUE), newmods = seq(1,200, length = 200))) %>% mutate(x = row_number())}
  
  if(type == "poly") {
    points <- data.frame(predict(robust(mv, cluster = cohort, clubSandwich = TRUE), newmods = unname(poly((seq(1,200, length = 200)), degree=2, raw=TRUE)))) %>% 
      mutate(x = row_number())}
  
  if(type == "spline"){
    knots <- attr(rcs(model.matrix(robust(mv, cluster = cohort, clubSandwich = TRUE))[,2], 3), "parms")
    points <- data.frame(predict(robust(mv, cluster = cohort, clubSandwich = TRUE), newmods = rcspline.eval((seq(1,200, length = 200)), knots, inclx = TRUE))) %>%
      mutate(x = row_number())}
  
  if(type == "log") {
    points <- data.frame(predict(robust(mv, cluster = cohort, clubSandwich = TRUE), newmods = log(seq(1,200, length = 200)))) %>% mutate(x = row_number())
  }
  
  if(type == "factor") {
    points <- data.frame(predict(robust(mv, cluster = cohort, clubSandwich = TRUE), newmods = diag(1, mv$tau2s, mv$tau2s))) %>% mutate(x = c(3,6,9,12,24,48,96))
  }

#get the ci points
poly <- points %>% 
  select(ci.ub, x) %>% 
  bind_rows(., points %>% 
              select(ci.lb, x) %>% 
              rename(ci.ub = ci.lb) %>% 
              arrange(desc(x))) # reverse them so the path isnt crossed

#get the pi points
#poly2 <- points %>% 
#  select(pi.ub, x) %>% 
#  bind_rows(., points %>% 
#             select(pi.lb, x) %>% 
#              rename(pi.ub = pi.lb) %>% 
#              arrange(desc(x))) # reverse them so the path isnt crossed

#plot

if(type == "factor") {
  data %>%
    ggplot(aes(x = as.numeric(as.character(timepoint_cut)), y = yi, group = interaction(cohort, group))) +
    geom_point(aes(size = acl_n), colour = "black", alpha = 0.3) + 
    geom_line(colour = "dark grey", alpha = 0.8) + 
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Time since surgery (Months)", y = "Percentage Deficit") +
    scale_size(range = c(0, 10)) +
    #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
    #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
    geom_line(data = points, aes(x = x, y = pred), colour = "orange", inherit.aes = FALSE) +
    geom_polygon(data = poly, aes(x = x, y = ci.ub), fill = "orange", alpha = 0.3, inherit.aes = FALSE) 
  #geom_polygon(data = poly2, aes(x = x, y = pi.ub), colour = "grey", alpha = 0.1, inherit.aes = FALSE)
} else {
data %>%
  ggplot(aes(x = timepoint_mean, y = yi, group = interaction(cohort, group))) +
  geom_point(aes(size = acl_n), colour = "black", alpha = 0.3) + 
  geom_line(colour = "dark grey", alpha = 0.8) + 
  scale_y_continuous(labels = scales::percent) +
  scale_size(range = c(0, 10)) +
  labs(x = "Time since surgery (Months)", y = "Percentage Deficit") +
  #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
  #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
  geom_line(data = points, aes(x = x, y = pred), colour = "orange", inherit.aes = FALSE) +
  geom_polygon(data = poly, aes(x = x, y = ci.ub), fill = "orange", alpha = 0.3, inherit.aes = FALSE) 
  #geom_polygon(data = poly2, aes(x = x, y = pi.ub), colour = "grey", alpha = 0.1, inherit.aes = FALSE)
}
}


pointsfunction <- function(data){
  data %>%
    select(ci.ub, x) %>% 
    bind_rows(., data %>% 
                select(ci.lb, x) %>% 
                rename(ci.ub = ci.lb) %>% 
                arrange(desc(x)))
}

pi_pointsfunction <- function(data){
  data %>%
    select(pi.ub, x) %>% 
    bind_rows(., data %>% 
                select(pi.lb, x) %>% 
                rename(pi.ub = pi.lb) %>% 
                arrange(desc(x)))
}



####
## Function to help model selection
mod_selection <- function(model) {
  results <- list() # create empty list
  plots <- list()
  
  data <- model$data # get data from model
  timepoint <- sort(data$timepoint_mean) # get the timepoint information
  
  
  # list of different moderator calls
  calls <- list("~timepoint_mean", 
                "~log(timepoint_mean)", 
                "~poly(timepoint_mean, degree = 2, raw = TRUE)", 
                "~rcs(timepoint_mean, 3)",  
                "~rcs(timepoint_mean, 4)"
  )
  
  for (i in calls) { # for each possible moderator format fit the model
    model_res <- update(model, as.formula(i)) # take the input model and add the relevant moderator
    fit_stats <- fitstats(model_res) # provide  fitstats
    results[[i]] <- fit_stats # add to list
    plots[[i]] <- mv_plotfunction(model_res)
  }
  
  results <- results %>%
    map_df(~ data.frame(t(.), row.names = NULL), .id = "mod") %>% # convert list to a dataframe 
    mutate(mod = c("Linear", "Log", "Poly (2)", "3 knot RCS", "4 knot RCS"), .before = 1)
  
  allplots <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], nrow = 2, ncol = 3, legend = "none")
  
  #return(results)
  return(list(results, allplots))
}


####
## Function for plotting mv object


mv_plotfunction <- function(model, logscale = FALSE){
  
  knot_pos <- NULL # set to 0 for logic to work later
  
  moderator <- as.character(model$call$mods[2]) # get the moderator variable
  data <- model$data # get the data
  
  last_timepoint <- round(max(data$timepoint_mean),0)
  
  total_k  <- model$k
  total_studies <- data %>% filter(!is.na(vi)) %>% summarise(n = length(unique(study))) %>% as.numeric()
  total_n <- data %>% distinct(cohort, .keep_all = TRUE)  %>% filter(!is.na(vi)) %>% summarise(n = sum(acl_n)) %>% as.numeric
    
    sum(data$acl_n[!is.na(data$vi)])
  
  if (str_detect(moderator, "rcs")) { # if a rcs is used, extract the information
    spline <- as.list(model$formula.mods[[2]][[3]]) 
    knot_pos <- unlist(spline[sapply(spline, is.numeric)]) # this gets the vector of knot positions or no. of knots
  }
  
  # depending on the type od moderator (i.e. log, linear, rcs...) calculate the predicted curve points
  mod_type <- if (str_detect(moderator, "log")) {
    
    #log
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods = log(seq(1,100, length = 100)),
                                 transf = exp
                                 )) %>% mutate(x = row_number())
    
  } else if (str_detect(moderator, "poly")) {
    
    # poly
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods=unname(poly(seq(1,100, length = 100), degree = 2, raw=TRUE)),
                                 transf = exp
                                 )) %>% mutate(x = row_number())
    
  } else if (length(knot_pos) == 3) {
    
    # 3 knot rcs
    knots <- attr(rcs(model.matrix(model)[,2], knot_pos), "parms")
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods=rcspline.eval(seq(1,100, length = 100), knots, inclx=TRUE),
                                 transf = exp)) %>% 
      mutate(x = row_number())
    
  } else if (length(knot_pos) == 1) { # i.e. cases where only the number of knots is given
    
    # 3 knot rcs
    knots <- attr(rcs(data$timepoint_mean, as.numeric(knot_pos)), "parms") # use model data to get knot positions
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods=rcspline.eval(seq(1,100, length = 100), knots, inclx=TRUE),
                                 transf = exp)) %>% 
      mutate(x = row_number())
    
  } else if (length(knot_pos) == 4) {
    
    # 4 knot rcs
    knots <- attr(rcs(model.matrix(model)[,2], knot_pos), "parms")
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods=rcspline.eval(seq(1,100, length = 100), knots, inclx=TRUE),
                                 transf = exp)) %>% 
      mutate(x = row_number())
    
  } else {
    # Linear
    points <- data.frame(predict(robust(model, cluster = cohort, clubSandwich = TRUE), 
                                 newmods = seq(1,100, length = 100),
                                 transf = exp)) %>% mutate(x = row_number())
  }
  
  # calculate ci points
  ci <-  pointsfunction(points)
  
  ci <- ci %>% filter(x < last_timepoint + 6)
  points <- points %>% filter(x < last_timepoint + 6)
  
  
  if (logscale == FALSE){
    
    # plot data and predicted model
    plot <- data %>%
      ggplot(aes(x = timepoint_mean, y = exp(yi), group = cohort)) +
      geom_point(aes(size = acl_n), alpha = 0.3) + 
      geom_line(alpha = 0.8) + 
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
      annotate(geom = "text", x = 100, y = 0.5, label = paste0("k = ", total_k, " (", total_studies, " studies)"), family = "Karla", hjust = 1) +
      annotate(geom = "text", x = 100, y = 0.55, label = paste0("n = ", total_n), family = "Karla", hjust = 1) +
      theme_mgpub()
    
  } else {
    
    plot <- data %>%
      ggplot(aes(x = timepoint_mean, y = exp(yi), group = cohort)) +
      geom_point(aes(size = acl_n), alpha = 0.3) + 
      geom_line(alpha = 0.8) + 
      #scale_y_continuous(labels = scales::percent, limits = c(-0.75, 0.25)) +
      scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2), labels = c("-40%", "-20%", "0%", "+20%")) +
      scale_x_continuous(trans = 'log10', limits = c(2,150)) +
      #coord_cartesian(xlim = c(0,100), ylim = c(-0.75, 0.25)) +
      scale_size(range = c(0, 10)) +
      labs(x = "Time since surgery (Months)", y = "Percentage Deficit", colour = "Measure", size = "Participants (n)") +
      #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
      #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
      geom_line(data = points, aes(x = x, y = pred), colour = "red", linewidth = 1.3, inherit.aes = FALSE) +
      geom_polygon(data = ci, aes(x = x, y = ci.ub), fill = "red",  alpha = 0.1, inherit.aes = FALSE) +
      theme_mgpub()
    
  }
  
  return(plot)
  
}



# This function returns some information around the predicted fit of the model 
# Estimates at 1, 2 and 5 years
# Last Data point
# zero_crossing point

predict_details <- function(model) {
  
  knot_pos <- NULL # set to 0 for logic to work later
  
  moderator <- as.character(model$call$mods[2]) # get the moderator variable
  data <- model$data # get the data
  
  total_k <- length(!is.na(data$vi))
  total_studies <- as.numeric(data %>% filter(!is.na(vi)) %>% summarise(length(unique(cohort))))
  total_n <- sum(data$acl_n[!is.na(data$vi)])
  
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
  
  
  cutpoints <- points %>%
    filter(x %in% c(12, 24, 60)) %>% # get the specific timepoint data
    select(pred, ci.lb, ci.ub) %>%
    mutate(across(where(is.numeric), ~round(.x*100, 1))) %>% # transform to percentage
    mutate(timepoint = c("1 year", "2 years", "5 years"), .before = 1) %>%
    mutate(new = paste0(pred, " (", ci.lb, " to ", ci.ub, ")")) %>% 
    select(timepoint, new) %>% 
    pivot_wider(names_from = timepoint, values_from = new) 
  
  logpoints <- points %>% # need to transform back to logscale to be able to use zero_crossings function
    mutate(across(where(is.numeric), ~log(.)))
  
  lastdata <- round(max(data$timepoint_mean)/12,1) # get the highest timepoint in fitted data
  plotlastdata <- lastdata + 6 # add 5 months to this for predicted fit purposes
  
  lastpoints <- logpoints %>% filter(x < plotlastdata + 1) # filter predicted data based on this
  
  zerocrossing <- round(modelbased::zero_crossings(lastpoints$ci.ub)[1]/12,1) # get the zerocrossing point
  
  if (!is.na(zerocrossing)) { 
    zerocrossing <- paste0(zerocrossing, " years") # Add months string if not NA
  }
  
  table <- cutpoints %>% 
    bind_cols('Zero crossing' = zerocrossing) %>% 
    bind_cols('Last Data Point' = paste0(lastdata, " years"))
  
  return(table)
  
}


## Function conducts regression of the true effects for within and between person analysis and 
## Plots the result
## Input is a bivariate rma.mv model 
corplot_function <- function(model){
  data <- model$data # get data
  coef <- exp(coef(model)) # get co-efficients for matreg as well as plot
  
  reg <- matreg(y = 2, x = 1, R = model$G, cov = TRUE, means = coef, V = model$vvc) # regress true effects
  ci <- as.data.frame(ellipse(model$vb, centre = coef, level = 0.95)) %>% rename(x = 1, y = 2) # calculate CI ellipse
  pred <- as.data.frame(ellipse(model$G, centre = coef, level = 0.95)) %>% rename(x = 1, y = 2) # calculate pi (95% coverage) ellipse
  
  plot <- data %>% 
    select(study, type, measure_2, yi, vi) %>% # select only the estimates
    pivot_longer(-c(study, type, measure_2), names_to = "var", values_to = "val") %>%
    pivot_wider(id_cols = c(study, measure_2), names_from = c(type, "var"), values_from = "val") %>%
    ggplot(aes(x = exp(casecontrol_yi), y = exp(within_yi))) + 
    geom_abline(intercept = 0, slope = 1, colour = "grey") + # plot line of perfect agreement
    geom_abline(intercept = 0.25, slope = 0.75, linetype = "dotted", colour = "grey") + # line
    geom_abline(intercept = 0.5, slope = 0.5, linetype = "dotted", colour = "grey", alpha = 0.7) + # line showing 50% difference
    geom_vline(xintercept = 1, colour = "grey") + # plot 0 line
    geom_hline(yintercept = 1, colour = "grey") +
    geom_point() +
    geom_abline(intercept = reg$tab$beta[1], slope = reg$tab$beta[2], colour = "red") + # regression coefficient from matreg
    geom_point(x = coef[1], y = coef[2], inherit.aes = FALSE, colour = "red", size = 3) +
    geom_polygon(data = ci, mapping = aes(x = x, y = y), alpha = 0.15, inherit.aes = FALSE, fill = "red") + # confidence ellipse
    geom_polygon(data = pred, mapping = aes(x = x, y = y), alpha = 0.15, inherit.aes = FALSE) + # pi ellipse
    coord_cartesian(xlim = c(0.4, 1.15), ylim = c(0.4, 1.15), clip = "off") + # allow clipping
    scale_x_continuous(breaks = c(0.5, 0.75, 1.0), labels = c("50%", "25%", "0%")) +
    scale_y_continuous(breaks = c(0.5, 0.75, 1.0), labels = c("50%", "25%", "0%")) +
    labs(x = "Between Person Deficit", y = "Within Person Deficit") +
    theme_mgpub()
  
  #return(plot)
  return(list(reg, plot))
  
}


## Function for graft plots

graft_plot <- function(model, graft = "quad"){
  
  name <- if (graft == "quad") {
    "Extensor"
  } else if (graft == "hs") {
    "Flexor"
  } else {
    "name"
  }
  
  knot_pos <- NULL # set to 0 for logic to work later
  
  moderator <- as.character(model$call$mods[2]) # get the moderator variable
  data <- model$data # get the data
  
  last_timepoint <- data %>% filter(!is.na(vi)) %>% summarise(round(max(timepoint_mean),0)) %>% as.numeric
  
  if (str_detect(moderator, "rcs")) { # if a rcs is used, extract the information
    spline <- as.list(model$formula.mods[[2]][[3]][[3]]) 
    knot_pos <- unlist(spline[sapply(spline, is.numeric)]) # this gets the vector of knot positions or no. of knots
    
  }
  
  
  
  # depending on the type od moderator (i.e. log, linear, rcs...) calculate the predicted curve points
  mod_type <- if (str_detect(moderator, "log")) {
    
    #log
    line1 <- data.frame(predict(model, 
                                newmods = cbind(rep(0, 100), log(seq(1,100, length = 100)), rep(0, 100)), 
                                transf = exp)) %>% mutate(x = row_number())
    
    line2 <- data.frame(predict(model, 
                                newmods = cbind(rep(1, 100), log(seq(1,100, length = 100)), rep(1, 100)*log(seq(1,100, length = 100))), 
                                transf = exp)) %>% mutate(x = row_number())
    
    
  }  else if (length(knot_pos) == 3) {
    
    # 3 knot rcs
    line1 <- data.frame(predict(model, newmods = 
                                  matrix(c(rep(0, 100), 
                                           rcspline.eval(seq(1, 100, length = 100), knot_pos, inclx = TRUE), 
                                           rcspline.eval(seq(1, 100, length = 100), knot_pos, inclx = TRUE)*0),
                                         nrow = 100, 
                                         ncol = 5), transf = exp)) %>% mutate(x = row_number())
    line2 <- data.frame(predict(model, newmods = 
                                  matrix(c(rep(1, 100), 
                                           rcspline.eval(seq(1, 100, length = 100), knot_pos, inclx = TRUE), 
                                           rcspline.eval(seq(1, 100, length = 100), knot_pos, inclx = TRUE)*1),
                                         nrow = 100, 
                                         ncol = 5), transf = exp)) %>% mutate(x = row_number())
    
    
  } else if (length(knot_pos) == 4) {
    
    # 4 knot rcs
    line1 <- data.frame(predict(model, newmods = 
                                  matrix(c(rep(0, 100), 
                                           rcspline.eval(seq(1, 100, length = 100), knot_pos, inclx = TRUE), 
                                           rcspline.eval(seq(1, 100, length = 100), knot_pos, inclx = TRUE)*0),
                                         nrow = 100, 
                                         ncol = 7), transf = exp)) %>% mutate(x = row_number())
    line2 <- data.frame(predict(model, newmods = 
                                  matrix(c(rep(1, 100), 
                                           rcspline.eval(seq(1, 100, length = 100), knot_pos, inclx = TRUE), 
                                           rcspline.eval(seq(1, 100, length = 100), knot_pos, inclx = TRUE)*1),
                                         nrow = 100, 
                                         ncol = 7), transf = exp)) %>% mutate(x = row_number())
    
  }
  
  
  # calculate ci points
  ci1 <-  pointsfunction(line1)
  ci2 <-  pointsfunction(line2)
  
  ci1 <- ci1 %>% filter(x < last_timepoint + 6)
  ci2 <- ci2 %>% filter(x < last_timepoint + 6)
  line1 <- line1 %>% filter(x < last_timepoint + 6)
  line2 <- line2 %>% filter(x < last_timepoint + 6)
  
  plot <- data %>%
    ggplot(aes(x = timepoint_mean, y = exp(yi), group = interaction(cohort, graft_group2), colour = graft_group2)) +
    geom_point(aes(size = acl_n), alpha = 0.1) + 
    geom_line(alpha = 0.1) + 
    scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2), labels = c("-40%", "-20%", "0%", "+20%")) +
    scale_colour_manual(values = c("red", "blue"), labels = c(name, "Other"), name = "Graft") +
    #scale_y_continuous(labels = scales::percent, limits = c(-0.75, 0.25)) +
    #scale_x_continuous(trans = 'log10', limits = c(3,150)) +
    coord_cartesian(xlim = c(0,100), ylim = c(0.5, 1.2)) +
    scale_size(range = c(0, 10)) +
    labs(x = "Time since surgery (Months)", y = "Percentage Deficit", colour = "Measure", size = "Participants (n)") +
    #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
    #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
    geom_line(data = line1, aes(x = x, y = pred), colour = "red", linewidth = 1.3, inherit.aes = FALSE) +
    geom_line(data = line2, aes(x = x, y = pred), colour = "blue", linewidth = 1.3, inherit.aes = FALSE) +
    geom_polygon(data = ci1, aes(x = x, y = ci.ub), fill = "red",  alpha = 0.1, inherit.aes = FALSE) +
    geom_polygon(data = ci2, aes(x = x, y = ci.ub), fill = "blue",  alpha = 0.1, inherit.aes = FALSE) +
    theme_mgpub()
  
  return(plot)
  
}



corplot_function_2 <- function(model){
  data <- model$data
  coef <- exp(coef(model))
  
  data2 <- model$data %>%
    select(study, yi, vi, type) %>%
    pivot_wider(id_cols = study, names_from = type,
                values_from = c(yi, vi), 
                names_glue = "{.value}_{type}") %>%
    mutate(covar = vector("list", length = nrow(.)))
  
  data2$covar <- blsplit(model$V, cluster = data$cohort)
  
  
  data2 <- data2 %>%
    mutate(covar2 = map(covar, list_to_matrix)) %>%
    rowwise() %>%
    mutate(centre = list(c(exp(yi_casecontrol), exp(yi_within)))) %>%
    ungroup() %>%
    mutate(ell = map2(.x = covar, .y = centre, ~ellipse(matrix(unlist(.x), nrow = 2), centre = unlist(.y)))) %>%
    ungroup() %>%
    mutate(ellx = map(ell, as.data.frame)) %>%
    unnest(ellx)
  
  
  reg <- matreg(y = 2, x = 1, R = model$G, cov = TRUE, means = coef, V = model$vvc)
  ci <- as.data.frame(ellipse(model$vb, centre = coef, level = 0.95)) %>% rename(x = 1, y = 2)
  pred <- as.data.frame(ellipse(model$G, centre = coef, level = 0.95)) %>% rename(x = 1, y = 2)
  
  plot <- data %>% 
    select(study, type, measure_2, yi, vi) %>% # select only the estimates
    pivot_longer(-c(study, type, measure_2), names_to = "var", values_to = "val") %>%
    pivot_wider(id_cols = c(study, measure_2), names_from = c(type, "var"), values_from = "val") %>%
    ggplot(aes(x = exp(casecontrol_yi), y = exp(within_yi))) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_abline(intercept = 0.25, slope = 0.75, linetype = "dotted", colour = "grey") +
    geom_abline(intercept = 0.5, slope = 0.5, linetype = "dotted", colour = "grey", alpha = 0.7) +
    geom_vline(xintercept = 1, colour = "grey") + 
    geom_hline(yintercept = 1, colour = "grey") +
    geom_point() +
    geom_polygon(data = data2, aes(x = x, y = y, group = study), inherit.aes = FALSE, alpha = 0.07) +
    geom_abline(intercept = reg$tab$beta[1], slope = reg$tab$beta[2], colour = "red", linewidth = 1) +
    #geom_point(x = coef[1], y = coef[2], inherit.aes = FALSE, colour = "red", size = 3) +
    #geom_polygon(data = ci, mapping = aes(x = x, y = y), alpha = 0.15, inherit.aes = FALSE, fill = "red") +
    #geom_polygon(data = pred, mapping = aes(x = x, y = y), alpha = 0.15, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0.4, 1.15), ylim = c(0.4, 1.15), clip = "off") +
    scale_x_continuous(breaks = c(0.5, 0.75, 1.0), labels = c("50%", "25%", "0%")) +
    scale_y_continuous(breaks = c(0.5, 0.75, 1.0), labels = c("50%", "25%", "0%")) +
    annotate(geom = "text", x = 0.4, y = 0.4, label = "1x", family = "Karla", colour = "grey", vjust = 0, angle = 45) + 
    annotate(geom = "text", x = 0.4, y = 0.55, label = "1.5x", family = "Karla", colour = "grey", vjust = 0, angle = 30) + 
    annotate(geom = "text", x = 0.4, y = 0.7, label = "2x", family = "Karla", colour = "grey", vjust = 0, angle = 20) +
    labs(x = "Between Person Deficit", y = "Within Person Deficit") +
    theme_mgpub()
  
  #return(plot)
  return(list(reg, plot))
  
}

list_to_matrix <- function(lst) {
  matrix(unlist(lst), nrow = 2)
}
