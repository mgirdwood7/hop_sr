
library(DescTools)
data <- read.csv("data/processed/data_20231128.csv")
demographics <- read.csv("data/processed/demographics_20231128.csv")
methods <- read.csv("data/processed/methods_20231128.csv")


cc_studies <- casecontrol %>% filter(str_detect(measure, "single hop|triple hop|triple crossover hop|side hop|vertical hop|6m timed hop|medial hop|lateral hop|square hop")) %>%
  distinct(study, .keep_all = TRUE) %>%
  filter(!study %in% c("Geoghegan 2007", "Karanikas 2004", "Karanikas 2005", "Nicholas 2001", 
                       "Urabe 2002", "Tate 2017", "Dalton 2011", "Hall 2015", "Noehren 2014", "Rahova 2020",
                       "Thomas 2013", "Balki 2019", "Rhatomy", "Sullivan 2022"))

within_studies <- within_data %>% filter(str_detect(measure, "single hop|triple hop|triple crossover hop|side hop|vertical hop|6m timed hop|medial hop|lateral hop|square hop")) %>%
  distinct(study, .keep_all = TRUE) %>%
  filter(!study %in% c("Geoghegan 2007", "Karanikas 2004", "Karanikas 2005", "Nicholas 2001", 
                       "Urabe 2002", "Tate 2017", "Dalton 2011", "Hall 2015", "Noehren 2014", "Rahova 2020",
                       "Thomas 2013", "Balki 2019", "Rhatomy", "Sullivan 2022"))


hop_studies <- bind_rows(cc_studies %>% select(study, study_id, cohort), within_studies %>% select(study, study_id, cohort)) %>%
  distinct(study, .keep_all = TRUE)


demo <- demographics %>% filter(study %in% c(hop_studies$study, "Arundale 2017", "Capin 2019")) %>%
  left_join(., methods %>% select(study, country, cohort), by = "study") 

# Studies
demo %>%
  summarise(n_studies = length(unique(study)),
            n_countries = length(unique(country)))

# Total n - unique participants (removing n where same cohort followed up)
uniquen <- demo %>% 
  mutate(year = as.numeric(str_extract(study, "\\d+"))) %>% 
  mutate(year = ifelse(study == "Capin 2019", 2017, year)) %>%
  arrange(cohort, year) %>%
  group_by(cohort) %>%
  filter(year == min(year)) %>%
  ungroup() 

uniquen %>%
  summarise(n = sum(n),
            n_female = sum(female, na.rm = TRUE))


bind_rows(casecontrol, within_data) %>%
  filter(study %in% c(hop_studies$study)) %>%
  filter(timepoint_mean != 0.1) %>%
  mutate(timepoint_mean = round(timepoint_mean, 0)) %>%
  group_by(study, timepoint_mean) %>%
  distinct(timepoint_mean, .keep_all = TRUE) %>%
  ungroup() %>%
  group_by(timepoint_mean) %>%
  summarise(freq = n()) %>% arrange(desc(freq))

demo %>%
  summarise(
    n_countries = unique(country))

timepointdata <- within_data %>%
  group_by(study) %>%
  mutate(n_timepoints = case_when(
    any(timepoint_mean< 0.2) ~ n_timepoints - 1,
    TRUE ~ n_timepoints
  )) %>%
  ungroup() %>%
  distinct(study, .keep_all = TRUE) %>% 
  group_by(cohort) %>%
  mutate(cohort_timepoints = sum(n_timepoints, na.rm = TRUE)) %>%
  ungroup() 

timepointdata %>%
  filter(study %in% hop_studies$study) %>%
  summarise(mode = Mode(n_timepoints, na.rm = TRUE),
            median = median(n_timepoints, na.rm = TRUE),
            morethantwo = sum(n_timepoints > 1, na.rm = TRUE))

timepointdata %>%
  filter(study %in% hop_studies$study) %>%
  filter(cohort_timepoints > 1) %>%
  summarise(mode = Mode(n_timepoints, na.rm = TRUE),
            median = median(n_timepoints, na.rm = TRUE),
            max = max(n_timepoints, na.rm = TRUE),
            q1 = quantile(n_timepoints, 0.25, na.rm = TRUE),
            q3 = quantile(n_timepoints, 0.57, na.rm = TRUE))

demo %>%
  group_by(study) %>%
  summarise(n = sum(n),
            n_f = sum(female)) %>%
  ungroup %>%
  mutate(n_m = n - n_f) %>%
  summarise(all_m = sum(n_f == 0, na.rm = TRUE),
            all_f = sum(n_m == 0, na.rm = TRUE))

demo %>% summarise(age = paste(min(age_mean, na.rm = TRUE), " - ", max(age_mean, na.rm = TRUE)))

demo %>% summarise(age = paste(min(bmi_mean, na.rm = TRUE), "-", max(bmi_mean, na.rm = TRUE)))


demographics %>%
  filter(study %in% c(hop_studies$study)) %>% 
  filter(!is.na(graft)) %>%
  group_by(study) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  summarise(median = median(n, na.rm = TRUE),
            max = max(n, na.rm = TRUE),
            q1 = quantile(n, 0.25, na.rm = TRUE),
            q3 = quantile(n, 0.75, na.rm = TRUE))

demographics %>%
  filter(study %in% c(uniquen$study)) %>% 
  mutate(graft_group = case_when(
    str_detect(graft, "Contralateral") ~ "other",
    str_detect(graft, "LARS") ~ "other",
    str_detect(graft, "HS") ~ "hs",
    str_detect(graft, "BPTB|PT|QT") ~ "quad",
    is.na(graft) ~ NA_character_,
    str_detect(graft, "mixed|Mixed") ~ "mixed",
    TRUE ~ "other"
  )) %>%
  group_by(graft_group) %>%
  summarise(sum(n))

demographics %>%
  filter(study %in% c(hop_studies$study)) %>% 
  filter(!is.na(graft)) %>%
  summarise(median = median(bmi_mean, na.rm = TRUE),
            max = max(bmi_mean, na.rm = TRUE),
            min = min(bmi_mean, na.rm = TRUE),
            q1 = quantile(bmi_mean, 0.25, na.rm = TRUE),
            q3 = quantile(bmi_mean, 0.75, na.rm = TRUE))

demographics %>%
  filter(study %in% c(hop_studies$study)) %>% 
  filter(!is.na(graft)) %>%
  summarise(median = median(age_mean, na.rm = TRUE),
            max = max(age_mean, na.rm = TRUE),
            min = min(age_mean, na.rm = TRUE),
            q1 = quantile(age_mean, 0.25, na.rm = TRUE),
            q3 = quantile(age_mean, 0.75, na.rm = TRUE))
