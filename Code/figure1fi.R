library(ggplot2)
library(dplyr)
library(readxl)
df <- read_excel("I:/CarbonF/CarbonSink_China/Data/rate_AGBeq.xlsx", col_types = c("text", "numeric", "numeric"))
## group carbon growth rate by ID
summary_rate <- df %>%
  group_by(ID) %>%
  summarise(
    mean = mean(C_rate),
    sd = sd(C_rate),
    n = n()
  )
ggplot(summary_rate, aes(x = factor(ID), y = mean)) +
  geom_col(width = 0.6, fill = "grey80", color = "black") +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.15,
                linewidth = 0.6) +
  geom_jitter(data = df,
              aes(x = factor(ID), y = C_rate),
              inherit.aes = FALSE,
              width = 0.08,
              height = 0,
              size = 2.5,
              alpha = 0.05) +
  labs(x = " ",
       y = "Carbon growth rate (MgC/ha/yr)") +
  theme_classic(base_size = 30)
ggsave("I:/CarbonF/CarbonSink_China/Output/Figure1i.jpg",
       plot = last_plot(),
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)
## group AGBeq by ID
summary_AGBeq <- df %>%
  group_by(ID) %>%
  summarise(
    mean = mean(AGBeq),
    sd = sd(AGBeq),
    n = n()
  )
ggplot(summary_AGBeq, aes(x = factor(ID), y = mean)) +
  geom_col(width = 0.6, fill = "grey80", color = "black") +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.15,
                linewidth = 0.6) +
  geom_jitter(data = df,
              aes(x = factor(ID), y = AGBeq),
              inherit.aes = FALSE,
              width = 0.08,
              height = 0,
              size = 2.5,
              alpha = 0.05) +
  labs(x = " ",
       y = "AGBeq (Mg/ha)") +
  theme_classic(base_size = 30)
ggsave("I:/CarbonF/CarbonSink_China/Output/Figure1f.jpg",
       plot = last_plot(),
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)

