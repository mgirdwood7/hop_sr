

#### Quad Within ####
# Data for Quads
hop <- within_data %>%
  filter(str_detect(measure, "triple hop")) %>%
  filter(timepoint_mean >0.2) %>% # no pre-operative data
  mutate(es_id = row_number())

# Reducing timepoint down to categories (3, 6, 12, 24, 48, 96 months post)
hop_cat <- hop %>%
  mutate(timepoint_cut = cut(timepoint_mean, 
                             breaks = c(0, 4.5, 9, 18, 36, 72, Inf), 
                             labels = c(3, 6, 12, 24, 48, 96))) %>%
  group_by(cohort, timepoint_cut) %>%
  # if multiple timepoints allocated to same category, take the closest to the categorical timepoint
  slice(which.min(abs(timepoint_mean - as.numeric(as.character(timepoint_cut))))) %>% 
  ungroup()


# not including measure as almost all studies measure only one fast/iso/slow so measure
# is indistinguishable from study level effects.

hopV <- vcalc(vi = vi,
               cluster = cohort,
               time1 = timepoint_mean,
               phi = 0.85,
               checkpd = TRUE,
               data = hop_cat)

hop_mv <- rma.mv(yi, hopV, 
                  mods = ~log(timepoint_mean), 
                  data = hop_cat, 
                  random = list(~ timepoint_mean|cohort),
                  struct = "CAR")

hop_mv_3 <- rma.mv(yi, hopV, 
                 mods = ~rcs(timepoint_mean, 3), 
                 data = hop_cat, 
                 random = list(~ timepoint_mean|cohort),
                 struct = "CAR")


hop_mv_4 <- rma.mv(yi, hopV, 
                   mods = ~rcs(timepoint_mean, 4), 
                   data = hop_cat, 
                   random = list(~ timepoint_mean|cohort),
                   struct = "CAR")

hop_mv_5 <- rma.mv(yi, hopV, 
                   mods = ~rcs(timepoint_mean, 5), 
                   data = hop_cat, 
                   random = list(~ timepoint_mean|cohort),
                   struct = "CAR")

hopknots <- attr(rcs(model.matrix(hop_mv_4)[,2], 5), "parms")

hop_pred <- data.frame(predict(hop_mv_4, newmods=rcspline.eval(seq(1,100, length = 100), hopknots, inclx=TRUE))) %>% mutate(x = row_number(), measure_2 = "Isometric")


hop_pred <- data.frame(predict(hop_mv, newmods = log(seq(1,100, length = 100)))) %>% mutate(x = row_number(), measure_2 = "Isometric")
hop_ci <- pointsfunction(hop_pred) %>% mutate(measure_2 = "Slow Isokinetic")

hop_cat %>%
  ggplot(aes(x = timepoint_mean, y = yi, group = interaction(cohort, group))) +
  geom_point(aes(size = acl_n), alpha = 0.3) + 
  geom_line(alpha = 0.8) + 
  scale_y_continuous(labels = scales::percent, limits = c(-0.75, 0.25)) +
  scale_x_continuous(limits = c(0, 200)) +
  #scale_x_continuous(trans = 'log10', limits = c(3,150)) +
  scale_size(range = c(0, 10)) +
  labs(title = "Hop - Between people", x = "Time since surgery (Months)", y = "Percentage Deficit",
        size = "ACL group n") +
  #geom_smooth(aes(x = timepoint_mean, y = yi), method = "lm", colour = "green", inherit.aes = FALSE) +
  #geom_abline(intercept = -0.2198, slope = 0.0014, colour = "red") +
  geom_line(data = hop_pred, aes(x = x, y = pred), colour = "red", inherit.aes = FALSE) +
  geom_polygon(data = hop_ci, aes(x = x, y = ci.ub), fill = "red", alpha = 0.1, inherit.aes = FALSE) +
  theme_mgpub()
