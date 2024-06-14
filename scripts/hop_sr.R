# Hop test analysis

library(metafor)
library(rms)
library(clubSandwich)

within_data <- read_csv("data/processed/within_data.csv") %>%
  filter(!study %in% c("Geoghegan 2007", "Karanikas 2004", "Karanikas 2005", "Nicholas 2001", 
                       "Urabe 2002", "Tate 2017", "Dalton 2011", "Hall 2015", "Noehren 2014", "Rahova 2020",
                       "Thomas 2013", "Balki 2019", "Rhatomy", "Sullivan 2022"))

# filtering studies that were included in hip review (where no limit on publication year), but need to remove forthis review
casecontrol <- read_csv("data/processed/casecontrol.csv") %>%
  filter(!study %in% c("Geoghegan 2007", "Karanikas 2004", "Karanikas 2005", "Nicholas 2001", 
                       "Urabe 2002", "Tate 2017", "Dalton 2011", "Hall 2015", "Noehren 2014", "Rahova 2020",
                       "Thomas 2013", "Balki 2019", "Rhatomy", "Sullivan 2022"))



## Single Leg Hop 

singlehop_data <- within_data %>%
  rename(acl_graft_group = graft_group) %>%
  filter(str_detect(measure, "single hop")) %>%
  filter(timepoint_mean >0.2, # no pre-operative data
         timepoint_mean < 120,
         str_detect(graft, "contralateral|Contralateral", negate = TRUE)) %>% 
  mutate(es_id = row_number()) 

hopV <- impute_covariance_matrix(vi = singlehop_data$vi, 
                                  cluster = singlehop_data$cohort,
                                  ti = singlehop_data$timepoint_mean,
                                  ar1 = 0.85,
                                  check_PD = TRUE,
                                  smooth_vi = TRUE,
                                  return_list = FALSE)

# Empty model for determining best fit
singlehop_empty <- rma.mv(yi, hopV, 
                        data = singlehop_data, 
                        random = list(~ timepoint_mean|cohort),
                        struct = c("CAR"))


singlehop_mv <- rma.mv(yi, hopV, 
                       mods = ~rcs(timepoint_mean, 4),
                          data = singlehop_data, 
                          random = list(~ timepoint_mean|cohort),
                          struct = c("CAR"))


## Triple  Hop

triplehop_data <- within_data %>%
  rename(acl_graft_group = graft_group) %>%
  filter(str_detect(measure, "triple hop")) %>%
  filter(timepoint_mean >0.2, # no pre-operative data
         timepoint_mean < 120,
         str_detect(graft, "contralateral|Contralateral", negate = TRUE)) %>% 
  mutate(es_id = row_number())


triplehopV <- impute_covariance_matrix(vi = triplehop_data$vi, 
                                 cluster = triplehop_data$cohort,
                                 ti = triplehop_data$timepoint_mean,
                                 ar1 = 0.85,
                                 check_PD = TRUE,
                                 smooth_vi = TRUE,
                                 return_list = FALSE)

# Empty model for determining best fit
triplehop_empty <- rma.mv(yi, triplehopV, 
                          data = triplehop_data, 
                          random = list(~ timepoint_mean|cohort),
                          struct = c("CAR"))


triplehop_mv <- rma.mv(yi, triplehopV, 
                       mods = ~log(timepoint_mean),
                       data = triplehop_data, 
                       random = list(~ timepoint_mean|cohort),
                       struct = c("CAR"))





## Triple Crossover Hop

triplexhop_data <- within_data %>%
  rename(acl_graft_group = graft_group) %>%
  filter(str_detect(measure, "triple crossover hop")) %>%
  filter(timepoint_mean >0.2, # no pre-operative data
         timepoint_mean < 120,
         str_detect(graft, "contralateral|Contralateral", negate = TRUE)) %>% 
  mutate(es_id = row_number()) 


triplexhopV <- impute_covariance_matrix(vi = triplexhop_data$vi, 
                                       cluster = triplexhop_data$cohort,
                                       ti = triplexhop_data$timepoint_mean,
                                       ar1 = 0.85,
                                       check_PD = TRUE,
                                       smooth_vi = TRUE,
                                       return_list = FALSE)

# Empty model for determining best fit
triplexhop_empty <- rma.mv(yi, triplexhopV, 
                          data = triplexhop_data, 
                          random = list(~ timepoint_mean|cohort),
                          struct = c("CAR"))


triplexhop_mv <- rma.mv(yi, triplexhopV, 
                       mods = ~log(timepoint_mean),
                       data = triplexhop_data, 
                       random = list(~ timepoint_mean|cohort),
                       struct = c("CAR"))



## Side Hop

sidehop_data <- within_data %>%
  rename(acl_graft_group = graft_group) %>%
  filter(str_detect(measure, "side hop")) %>%
  filter(timepoint_mean >0.2, # no pre-operative data
         timepoint_mean < 120,
         str_detect(graft, "contralateral|Contralateral", negate = TRUE)) %>% 
  mutate(es_id = row_number()) 

sidehopV <- impute_covariance_matrix(vi = sidehop_data$vi, 
                                        cluster = sidehop_data$cohort,
                                        ti = sidehop_data$timepoint_mean,
                                        ar1 = 0.85,
                                        check_PD = TRUE,
                                        smooth_vi = TRUE,
                                        return_list = FALSE)

# Empty model for determining best fit
sidehop_empty <- rma.mv(yi, sidehopV, 
                           data = sidehop_data, 
                           random = list(~ timepoint_mean|cohort),
                           struct = c("CAR"))


sidehop_mv <- rma.mv(yi, sidehopV, 
                        mods = ~rcs(timepoint_mean, 3),
                        data = sidehop_data, 
                        random = list(~ timepoint_mean|cohort),
                        struct = c("CAR"))



## Vertical Hop

verticalhop_data <- within_data %>%
  rename(acl_graft_group = graft_group) %>%
  filter(str_detect(measure, "vertical hop")) %>%
  filter(timepoint_mean >0.2, # no pre-operative data
         timepoint_mean < 120,
         str_detect(graft, "contralateral|Contralateral", negate = TRUE)) %>% 
  mutate(es_id = row_number()) 

verticalhopV <- impute_covariance_matrix(vi = verticalhop_data$vi, 
                                     cluster = verticalhop_data$cohort,
                                     ti = verticalhop_data$timepoint_mean,
                                     ar1 = 0.85,
                                     check_PD = TRUE,
                                     smooth_vi = TRUE,
                                     return_list = FALSE)

# Empty model for determining best fit
verticalhop_empty <- rma.mv(yi, verticalhopV, 
                        data = verticalhop_data, 
                        random = list(~ timepoint_mean|cohort),
                        struct = c("CAR"))

# Log fit best
verticalhop_mv <- rma.mv(yi, verticalhopV, 
                     mods = ~log(timepoint_mean),
                     data = verticalhop_data, 
                     random = list(~ timepoint_mean|cohort),
                     struct = c("CAR"))




## 6m Timed Hop

sixmhop_data_pre <- within_data %>%
  rename(acl_graft_group = graft_group) %>%
  filter(str_detect(measure, "6m timed hop")) %>%
  filter(timepoint_mean >0.2, # no pre-operative data
         timepoint_mean < 120,
         str_detect(graft, "contralateral|Contralateral", negate = TRUE)) %>% 
  mutate(es_id = row_number()) 
# Some studies invert 6m hop data results (as 6m hop score = better = lower score)
# Have individually checked all studies, some need to be switched

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

sixmhopV <- impute_covariance_matrix(vi = sixmhop_data$vi, 
                                         cluster = sixmhop_data$cohort,
                                         ti = sixmhop_data$timepoint_mean,
                                         ar1 = 0.85,
                                         check_PD = TRUE,
                                         smooth_vi = TRUE,
                                         return_list = FALSE)

# Empty model for determining best fit
sixmhop_empty <- rma.mv(yi, sixmhopV, 
                            data = sixmhop_data, 
                            random = list(~ timepoint_mean|cohort),
                            struct = c("CAR"))

# Log fit best.
sixmhop_mv <- rma.mv(yi, sixmhopV, 
                         mods = ~log(timepoint_mean),
                         data = sixmhop_data, 
                         random = list(~ timepoint_mean|cohort),
                         struct = c("CAR"))



##### Case Control ####



singlehopcont_data <- casecontrol %>%
  rename(timepoint_mean = acl_timepoint_mean) %>%
  filter(str_detect(measure, "single hop")) %>%
  filter(timepoint_mean >0.2, # no pre-operative data
         timepoint_mean < 120) %>% 
  mutate(es_id = row_number()) 

hopcontV <- impute_covariance_matrix(vi = singlehopcont_data$vi, 
                                 cluster = singlehopcont_data$cohort,
                                 ti = singlehopcont_data$timepoint_mean,
                                 ar1 = 0.85,
                                 check_PD = TRUE,
                                 smooth_vi = TRUE,
                                 return_list = FALSE)

# Empty model for determining best fit
singlehopcont_empty <- rma.mv(yi, hopcontV, 
                          data = singlehopcont_data, 
                          random = list(~ timepoint_mean|cohort),
                          struct = c("CAR"))


singlehopcont_mv <- rma.mv(yi, hopcontV, 
                       mods = ~log(timepoint_mean),
                       data = singlehopcont_data, 
                       random = list(~ timepoint_mean|cohort),
                       struct = c("CAR"))


forest_cc_function <- function(model){
  forest.meta(model,
              sortvar = acl_timepoint_mean,
              common = FALSE,
              prediction = TRUE,
              #xlab = label,
              #xlab.pos = xlab,
              #xlim = xlim,
              smlab = "",
              leftcols = c("study", "acl_timepoint_mean", "acl_n"),
              leftlabs = c("Study", "Months\npost ACLR", "n\nACLR"),
              rightcols = c("effect.ci"),
              rightlabs = c("RoM [95% CI]"),
              just.addcols = "left",
              addrows.below.overall = 1,
              digits.addcols = 1,
              print.pval.Q = FALSE,
              fontfamily = "Karla",
              ff.predict = 1,
              ref = 1,
              col.diamond = "black",
              #print.subgroup.labels = TRUE, 
              #subgroup.name = "Timepoint",
              #test.subgroup = FALSE, 
              #subgroup.hetstat = FALSE, 
              col.random = "grey")
}

## Triple Hop

# All except Casp have triple crossover hop
# Use triple crossover hop + Casp data
# Remove 2nd Patterson timepoint
# Order by timepoint?

triplehopcont_data <- casecontrol %>%
  filter(str_detect(measure, "triple hop|triple crossover hop")) %>%
  filter(study == "Casp 2021" | measure == "triple crossover hop",
         acl_timepoint_mean < 30) %>%
  metagen(TE = yi, seTE = sei, studlab = study, data = ., sm = "ROM")



# Only 4 studies measure triple or triple crossover hop case control
# Only 1 with multiple timepoint (Patterson)
# All at 12 months or less (Patterson has second timepoint at 5 year)

## Side Hop
# Only 4 studies, 1 with multiple timepoint (Patterson)
# 3 <24m, 2 > 50m.

sidehopcont_data <- casecontrol %>%
  filter(str_detect(measure, "side hop")) %>%
  filter(!(study == "Patterson 2020a" & acl_timepoint_mean == 12.0)) %>%
  #filter(acl_timepoint_mean < 30) %>%
  metagen(TE = yi, seTE = sei, studlab = study, data = ., sm = "ROM")


## Vertical Hop
# 5 studies
# Varied timepoints...

verticalhopcont_data <- casecontrol %>%
  filter(str_detect(measure, "vertical")) %>%
  metagen(TE = yi, seTE = sei, studlab = study, data = ., sm = "ROM")


## 6m timed hop
# 2 studies only

sixmhopcont_data <- casecontrol %>%
  filter(str_detect(measure, "6m timed hop")) %>%
  escalc(n1i = con_n, n2i = acl_n, m1i = con_inj_mean, sd1i = con_inj_sd, m2i = acl_inj_mean, sd2i = acl_inj_sd, 
         data = ., measure = "ROM") %>%
  metagen(TE = yi, seTE = sei, studlab = study, data = ., sm = "ROM")

otherhops <- within_data %>% 
  filter(str_detect(measure, "square hop|stair hop|lateral hop|medial hop|4-jump test|speedy hop test"))

