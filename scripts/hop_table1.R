

demo2 <- demo %>%
  filter((study == "Bodkin 2020" & name == "Repeat measure group") | 
           (study == "Dempsey 2019" & name == "All included") | # use combined group for table only
           (study == "Kuenze 2019b" & name %in% c("Male", "Female")) | # this one is same as in analysis
           (study == "Shibata 2019" & name == "All") | # use the all group
           (study == "Blakeney 2018" & name %in% c("G6M", "Healthy Control")) | # same as analysis
           (study == "Siney 2010" & name == "All") | # use the combined group for table only
           (study == "Ebert 2021a" & name == "Entire cohort") |
           (study == "Schwery 2022" & name %in% c("Male", "Female")) | # use different combination for this
           (study == "Hogberg 2023" & name == "12 months") | 
           (group_indep != 0 | is.na(group_indep))) # all the other indep group studies

# combined data so one line per study
demo3 <- demo2 %>%
  filter(str_detect(name, "Healthy Control", negate = TRUE), # not healthy control groups
         str_detect(graft, "Contralateral|contralateral", negate = TRUE),
         str_detect(graft, "ACLD", negate = TRUE)) %>%  # remove contralateral graft groups
  group_by(study) %>%
  #filter(n() > 1) %>%
  mutate(n_new = sum(n),
         female_new = sum(female),
         age_mean_new = combine_mean(n, age_mean),
         age_sd_new = combine_sd(n, age_mean, age_sd),
         bmi_mean_new = combine_mean(n, bmi_mean),
         bmi_sd_new = combine_sd(n, bmi_mean, bmi_sd),
         group_new = paste(name, collapse = " + "),
         graft_new = case_when( # combine graft names - if all the same then preserve the original, if different then mixed
           length(unique(graft)) == 1 ~ graft,
           TRUE ~ paste(unique(graft), collapse = " + ")
         )) %>%
  distinct(study, .keep_all = TRUE) %>%
  ungroup() %>%
  select(-c(n, female, age_mean, age_sd, bmi_mean, bmi_sd, name, graft, group, tegner, group_indep, group_indep_notes, timepoint_combo, combination_notes)) %>%
  rename_with(~str_replace(., "_new", ""), ends_with("_new")) # rename the new columns, removing "_new"


# function to paste things togehter so that NAs become blank ""
paste_na <- function(...,sep=" ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}


demo4 <- demo3 %>%
  mutate(age_strings = map(age_other, ~extract_string_data(.x, name = "age_"))) %>% # create new string variable using extract_string function
  unnest(age_strings) %>% # unnest 3 columns created (inj, non-inj and lsi)
  mutate(bmi_strings = map(bmi_note, ~extract_string_data(.x, name = "bmi_"))) %>% # create new string variable using extract_string function
  unnest(bmi_strings) %>% # unnest 3 columns created (inj, non-inj and lsi)
  mutate(age_med = ifelse(is.na(age_med), NA, paste("Median", age_med)), # add words to describe other statistics for table
         age_iqr = ifelse(is.na(age_iqr), NA, paste("IQR", age_iqr)),
         age_range = ifelse(is.na(age_range), NA, paste("Range", age_range)),
         bmi_med = ifelse(is.na(bmi_med), NA, paste("Median", bmi_med)),
         bmi_iqr = ifelse(is.na(bmi_iqr), NA, paste("IQR", bmi_iqr))) %>%
  select(-c(age_other, bmi_note)) %>% # remove old other columns 
  mutate(across(where(is.numeric), ~round(.x, 2))) %>% # round numeric data to 2 digis
  mutate(age_other = paste_na(age_med, age_iqr, age_range), # paste together all other statistics
         bmi_other = paste_na(bmi_med, bmi_iqr, bmi_range),
         age_sd = ifelse(is.na(age_sd), NA, paste0("(", age_sd, ")")),
         bmi_sd = ifelse(is.na(bmi_sd), NA, paste0("(", bmi_sd, ")")),
         ) %>% 
  select(-c(age_med, age_iqr, age_range, bmi_med, bmi_iqr)) %>%
  mutate(age = paste_na(age_mean, age_sd, age_other),
         bmi = paste_na(bmi_mean, bmi_sd, bmi_other)) %>%
  select(study, n, female, age, bmi, group, graft) %>%
  arrange(study)


details <- data %>%
  filter(str_detect(measure, "single hop|triple hop|triple crossover hop|side hop|vertical hop|6m timed hop|medial hop|lateral hop|square hop"),
         study %in% demo$study) %>%
  select(study, measure, timepoint_mean, units) %>%
  mutate(measure = str_replace(measure, "single hop", "Single forward hop"),
         measure = str_replace(measure, "triple hop", "Triple hop"),
         measure = str_replace(measure, "triple crossover hop", "Triple crossover hop"),
         measure = str_replace(measure, "side hop", "Side hop"),
         measure = str_replace(measure, "vertical hop", "Vertical hop"),
         measure = str_replace(measure, "6m timed hop", "6m timed forward hop"),
    measure = str_replace(measure, "medial hop", "Medial hop"),
    measure = str_replace(measure, "lateral hop", "Lateral hop"),
    measure = str_replace(measure, "square hop", "Square hop"),
    units = str_replace(units, "lsi", "LSI")) %>%
  mutate(timepoint_mean = round(timepoint_mean/12, 1)) %>%
  group_by(study) %>%
  distinct(measure, .keep_all = TRUE) %>%
  summarise(hop = paste(measure, collapse = "\n"),
            units = paste(units, collapse = "\n"))

# Timepoint

timepoint <- data %>%
  filter(str_detect(measure, "single hop|triple hop|triple crossover hop|side hop|vertical hop|6m timed hop|medial hop|lateral hop|square hop"),
         study %in% demo$study) %>%
  select(study, timepoint, timepoint_mean, units, n) %>%
  group_by(study) %>%
  distinct(timepoint_mean, .keep_all = TRUE) %>%
  ungroup() %>%
  group_by(study, timepoint) %>%
  mutate(time_new = round(combine_mean(n, timepoint_mean)/12, 1)) %>%
  ungroup() %>%
  group_by(study) %>%
  summarise(time = paste(sort(unique(time_new[time_new>0.16])), collapse = ", "))

details2 <- left_join(demo4, left_join(details, timepoint, by = "study"),
                      by = "study") %>%
  mutate(nf = paste0(n, " (", female, ")")) %>%
  select(-c(n, female)) %>%
  select(study, nf, everything())

write_csv(details2, "output/tables/table1.csv")


## Healthy Controls
healthycontrol <- demo %>%
  filter(str_detect(name, "Healthy Control")) %>%
  group_by(study) %>%
  #filter(n() > 1) %>%
  mutate(n_new = sum(n),
         female_new = sum(female),
         age_mean_new = combine_mean(n, age_mean),
         age_sd_new = combine_sd(n, age_mean, age_sd),
         bmi_mean_new = combine_mean(n, bmi_mean),
         bmi_sd_new = combine_sd(n, bmi_mean, bmi_sd),
         group_new = paste(name, collapse = " + "),
         graft_new = case_when( # combine graft names - if all the same then preserve the original, if different then mixed
           length(unique(graft)) == 1 ~ graft,
           TRUE ~ paste(unique(graft), collapse = " + ")
         )) %>%
  distinct(study, .keep_all = TRUE) %>%
  ungroup() %>%
  select(-c(n, female, age_mean, age_sd, bmi_mean, bmi_sd, name, graft, group, tegner, group_indep, group_indep_notes, timepoint_combo, combination_notes)) %>%
  rename_with(~str_replace(., "_new", ""), ends_with("_new")) # rename the new columns, removing "_new"

healthycontrol <- healthycontrol %>%
  mutate(age_strings = map(age_other, ~extract_string_data(.x, name = "age_"))) %>% # create new string variable using extract_string function
  unnest(age_strings) %>% # unnest 3 columns created (inj, non-inj and lsi)
  mutate(bmi_strings = map(bmi_note, ~extract_string_data(.x, name = "bmi_"))) %>% # create new string variable using extract_string function
  unnest(bmi_strings) %>% # unnest 3 columns created (inj, non-inj and lsi)
  mutate(age_med = ifelse(is.na(age_med), NA, paste("Median", age_med)), # add words to describe other statistics for table
         age_iqr = ifelse(is.na(age_iqr), NA, paste("IQR", age_iqr)),
         age_range = ifelse(is.na(age_range), NA, paste("Range", age_range)),
         bmi_med = ifelse(is.na(bmi_med), NA, paste("Median", bmi_med)),
         bmi_iqr = ifelse(is.na(bmi_iqr), NA, paste("IQR", bmi_iqr))) %>%
  select(-c(age_other, bmi_note)) %>% # remove old other columns 
  mutate(across(where(is.numeric), ~round(.x, 2))) %>% # round numeric data to 2 digis
  mutate(age_other = paste_na(age_med, age_iqr, age_range), # paste together all other statistics
         bmi_other = paste_na(bmi_med, bmi_iqr, bmi_range),
         age_sd = ifelse(is.na(age_sd), NA, paste0("(", age_sd, ")")),
         bmi_sd = ifelse(is.na(bmi_sd), NA, paste0("(", bmi_sd, ")")),
  ) %>% 
  select(-c(age_med, age_iqr, age_range, bmi_med, bmi_iqr)) %>%
  mutate(age = paste_na(age_mean, age_sd, age_other),
         bmi = paste_na(bmi_mean, bmi_sd, bmi_other)) %>%
  select(study, n, female, age, bmi) %>%
  arrange(study) %>%
  mutate(nf = paste0(n, " (", female, ")")) %>%
  select(-c(n, female)) %>%
  select(study, nf, everything())

details3 <- left_join(details2, healthycontrol %>% select(study, nf) %>% rename(hc_nf = nf), by = "study")

write_csv(details3, "output/tables/table1_hc.csv")


  