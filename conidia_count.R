library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbreak)
library(scales)

conidia_count <- read.csv("./conidia_count.csv")

conidia_count_long <- conidia_count %>%
  pivot_longer(
    cols = -replicate,
    names_to = "condition",
    values_to = "conidia_count"
  ) %>%
  mutate(
    condition = recode(condition,
                       "full.strength.PDA" = "PDA",
                       "half.strength.PDA" = "Half-PDA",
                       "one.tenth.strength.PDA" = "OneTenth-PDA"
    )
  ) %>%
  mutate(condition = factor(condition,
                            levels = c("PDA", "Half-PDA", "OneTenth-PDA")))
conidia_summary <- conidia_count_long %>%
  group_by(condition) %>%
  summarize(
    mean_n = mean(conidia_count),
    sd_n = sd(conidia_count),
    .groups = "drop")

t.test(conidia_count ~ condition,
       data = filter(conidia_count_long, condition %in% c("PDA", "Half-PDA")))

t.test(conidia_count ~ condition,
       data = filter(conidia_count_long, condition %in% c("PDA", "OneTenth-PDA")))

t.test(conidia_count ~ condition,
       data = filter(conidia_count_long, condition %in% c("Half-PDA", "OneTenth-PDA")))

ggplot(conidia_summary, aes(x = condition, y = mean_n)) +
  geom_bar(stat = "identity", width = 0.65, fill = "#4682B4") +
  geom_errorbar(aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n), width = 0.1) +
  scale_y_continuous(
    limits = c(0, 100000000),
    breaks = seq(0, 100000000, 20000000)
  ) +
  labs(
    x = "Condition",
    y = "Mean conidia count",
    title = "Conidia count across PDA conditions"
  ) +
  theme_minimal(base_size = 16)

ggplot(conidia_summary, aes(x = condition, y = mean_n)) +
  geom_bar(stat = "identity", width = 0.65, fill = "#4682B4") +
  geom_errorbar(
    aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n),
    width = 0.1
  ) +
  scale_y_break(c(5000000, 40000000), scales = 0.5) +
  scale_y_continuous(
    labels = label_scientific(),
    breaks = c(
      0, 2000000, 5000000,
      40000000, 60000000, 80000000, 100000000, 120000000)
  ) +
  labs(
    x = "Condition",
    y = "Conidia count"
  ) +
  theme_minimal(base_size = 20)

#ggsave("conidia_counts.png", width = 7, height = 5, bg = "white")
