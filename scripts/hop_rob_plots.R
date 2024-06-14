# Rob tables

library(tidyverse)

rob <- read_csv("data/processed/rob_final.csv")

# Need to filter data for each paper.

rob <- rob %>% filter(study %in% demo$study)

rob2 <- rob %>% mutate(rob_analysis = case_when(
  outcome_c != "Low" & blinding_assessor_c != "Low" & selection_c != "Low" ~ "High",
  TRUE ~ "Low"))

rob <- rob %>%
  arrange(study) %>%
  mutate(row = row_number(),
         plotgroup = cut(row, breaks = 4, labels = FALSE)) %>%
  select(-row) %>%
  arrange(plotgroup, desc(study))

rob <- rob %>%
  select(study, ends_with("_c"), plotgroup) %>%
  rename_with(~str_replace(.x, "_c", ""), ends_with("_c")) %>%
  pivot_longer(., -c(study, plotgroup), names_to = "Domain", values_to = "Value") %>%
  mutate(val2=if_else(Value == "Low", "+", if_else(Value == "High", "-", if_else(is.na(Value), "", "?")))) %>%
  mutate(Domain=factor(Domain,levels=unique(Domain))) %>%
  mutate(Value = factor(Value, levels = c("Low", "Unclear", "High", NA))) %>%
  mutate(study = factor(study, levels = unique(study)),
         plotgroup = factor(plotgroup))

#Plot with traffic light table
ggplot(data = rob, aes(y = study, x = Domain)) +
  geom_tile(color="black", fill="white", size = 0.7) +
  geom_point(aes(color=as.factor(val2)), size=8) +
  geom_text(aes(label = val2), size = 8) +
  scale_x_discrete(position = "top", labels = c("Sequence Generation", "Allocation", "Patient Blinding", "Therapist Blinding", "Assessor Blinding",
                                                "Outcome Measurement", "Selection", "Attrition", "Analysis")) +
  #scale_y_discrete(limits=rev(levels(as.factor(rob$study)))) +
  scale_color_manual(values = c("-" = "#BF0000",
                                "+" = "#02C100",
                                "?" = "#E2DF07"), na.value = "white") +
  theme_minimal() +
  #coord_equal() +
  theme_mgpub() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black", angle = 60, hjust=0),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 4,
        strip.text = element_blank(),
        strip.background = element_blank()) +
  facet_wrap(~ plotgroup, scale = 'free', nrow = 1)

ggsave("output/plots/hop_rob.png", width = 20, height = 25, limitsize = FALSE)



# Summary ROB Plot with %
ggplot(data = rob) +
  geom_bar(mapping = aes(x = Domain, fill = Value), width = 0.7, position = "fill", color = "black") +
  coord_flip(ylim = c(0, 1)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual("Risk of Bias",
                    labels = c(" Low risk of bias ",
                               " Unclear risk of bias ", 
                               " High risk of bias"),
                    values = c(High = "#BF0000",
                               Low = "#02C100", 
                               Unclear = "#E2DF07"),
                    na.value = "white") +
  scale_x_discrete(limits=rev(levels(as.factor(rob$Domain))), 
                   labels = rev(c("Sequence Generation", "Allocation", "Patient Blinding", "Therapist Blinding", "Assessor Blinding",
                              "Outcome Measurement", "Selection", "Attrition", "Analysis"))) + 
  scale_y_continuous(labels = scales::percent) +
  theme_mgpub() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(size = 14))

ggsave("output/plots/rob_summary.png", width = 10, height = 5, )


## Funnel Plot
## Hip abduction k=12

png("output/plots/funnelplot.png", height = 800, pointsize = 25, width = 1000)
par(family = "Karla")
hip_within_meta %>%
  filter(analysis_group == "hip abd") %>%
  select(rma) %>%
  pluck(1) %>%
  pluck(1) %>%
  metafor::funnel(.,
                  studlab = TRUE)
dev.off()

hip_within_es %>%
  filter(analysis_group == "hip abd") %>%
  mutate(y = yi/sei, x = 1/sei) %>%
  lm(y ~ x, data = .) %>%
  summary()

fastdata %>%
  mutate(y = yi/sei, x = 1/sei) %>%
  lm(y ~ x, data = .) %>%
  summary()
