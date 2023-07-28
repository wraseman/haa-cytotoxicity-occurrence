# Figure 1: CTI Values ####

index.tox %>%
  mutate(Color.Group = factor(Color.Group, levels = c("HAA5", "HAA9", "TTHM", "Other"))) %>%
  mutate(Group = factor(Group, levels = c("HAA", "THM", "HAN", "HAM", "HAL", "HNM", "HC"))) %>%
  ggplot(aes(x = Group, y = CTI)) +
  geom_hline(yintercept = c(0.1, 1, 10, 100, 1000), alpha = 0.3) +
  geom_point(aes(color = Color.Group), position = position_dodge(width = 0.3), size = 1.5, alpha = 0.7) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
  annotation_logticks(sides = "l", alpha = 0.3) +
  scale_color_manual(values = c(brewer.pal(4, "Set1")[1], brewer.pal(4, "Set1")[2], brewer.pal(4, "Set1")[3], "Black")) +
  labs(x = "", y = expression(Cytotoxicity~Index*','~CTI~(mu*M^-1~x10^3)), color = "") +
  theme(legend.position = "bottom", axis.title.x = element_blank())


ggsave("Output/Figures/Manuscript Revision/Figure1.png", width = 4.5, height = 3, units = "in", dpi = 1000)



# Figure 2: Concentrations and CTIs ####

HAA9.rank <- df.sum %>%
  filter(HAA.Group == "HAA9") %>%
  mutate(p.rank = percent_rank(c.mass_HAA9)) %>%
  select(PWSID, p.rank, c.mass_HAA9)

HAA6Br.rank <- df.sum %>%
  filter(HAA.Group == "HAA6Br") %>%
  mutate(p.rank = percent_rank(c.mass_HAA6Br)) %>%
  select(PWSID, p.rank)

CAT.rank <- df.sum %>%
  filter(HAA.Group == "HAA9") %>%
  mutate(p.rank = percent_rank(CAT_HAA9)) %>%
  select(PWSID, p.rank)

CAT.rank.HAA6Br <- df.sum %>%
  filter(HAA.Group == "HAA6Br") %>%
  mutate(p.rank = percent_rank(CAT_HAA9)) %>%
  select(PWSID, p.rank)

species.quantiles.conc <-  df %>%
  filter(HAA.Group == "HAA9") %>%
  group_by(group, DBP, subclass) %>%
  summarize(`95th` = quantile(c.mass, 0.95, na.rm = T),
            `98th` = quantile(c.mass, 0.98, na.rm = T),
            `99th` = quantile(c.mass, 0.99, na.rm = T)
  ) %>%
  pivot_longer(-c(DBP, group, subclass), names_to = "Quantile", values_to = "c.mass")

species.quantiles.CAT <- df %>%
  filter(HAA.Group == "HAA9") %>%
  left_join(., CAT.rank, by = "PWSID") %>%
  filter(p.rank > 0.5) %>%
  group_by(PWSID) %>%
  mutate(CAT.frac = CAT / sum(CAT, na.rm = T)) %>%
  group_by(group, DBP, subclass) %>%
  summarize(`95th` = quantile(CAT.frac, 0.95, na.rm = T),
            `98th` = quantile(CAT.frac, 0.98, na.rm = T),
            `99th` = quantile(CAT.frac, 0.99, na.rm = T)
  ) %>%
  pivot_longer(-c(DBP, group, subclass), names_to = "Quantile", values_to = "CAT.frac")

pct.detected <- df %>% 
  filter(HAA.Group == "HAA9") %>%
  drop_na(c.mass) %>%
  group_by(DBP, group, subclass) %>%
  summarize(n = n()) %>%
  mutate(Pct.Detected = round(n / 4747*100, 0)) %>%
  mutate(Pct.Detected = paste(Pct.Detected, "%", sep = ""))


p2a <- df %>%
  filter(HAA.Group == "HAA9") %>%
  mutate(DBP = fct_reorder(DBP, order)) %>%
  mutate(BSF = n.br / halogens) %>%
  ggplot(aes(x = DBP, y = c.mass)) +
  geom_boxplot(aes(color = group), outlier.alpha = 0, outlier.size = 0.7) +
  geom_point(data = species.quantiles.conc, aes(shape = as.factor(Quantile), color = group), size = 2) +
  geom_text(data = pct.detected, aes(y = 40, label = Pct.Detected), size = rel(2.5)) +
  scale_y_continuous(limits = c(0,40), breaks = seq(0,40,5)) +
  facet_grid(cols = vars(subclass), scales = "free", space = "free", labeller = as_labeller(c("mHAA" = "Mono-HAAs", "dHAA" = "Di-HAAs", "tHAA" = "Tri-HAAs"))) +
  labs(shape = "Percentile", color = "Regulatory\nGroup", x = "", y = expression(atop(Concentration~at~Max, paste(LAA~Location~(mu*g/L))))) +
  scale_color_brewer(palette = "Set2")


p2b <-  index %>%
  filter(class == "HAA") %>%
  mutate(DBP = fct_reorder(DBP, order)) %>%
  mutate(BSF = n.br / halogens) %>%
  ggplot(aes(x = DBP, y = CTI)) +
  geom_col(aes(fill = group)) +
  facet_grid(cols = vars(subclass), scales = "free", space = "free", labeller = as_labeller(c("mHAA" = "Mono-HAAs", "dHAA" = "Di-HAAs", "tHAA" = "Tri-HAAs"))) +
  labs(shape = "Percentile", fill = "Regulatory\nGroup", x = "", y = expression(atop(Cytotoxicity~Index, paste((mu*M^-1~x~10^3))))) +
  scale_y_continuous(limits = c(0,120), breaks = seq(0,120,20)) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")

p2c <- df %>%
  filter(HAA.Group == "HAA9") %>%
  left_join(., CAT.rank, by = "PWSID") %>%
  ungroup() %>%
  filter(p.rank >= 0.5) %>%
  mutate(DBP = fct_reorder(DBP, order)) %>%
  mutate(BSF = n.br / halogens) %>%
  group_by(PWSID) %>%
  mutate(CAT.frac = CAT / sum(CAT, na.rm = T)) %>%
  mutate(DBP = fct_reorder(DBP, order)) %>%
  mutate(BSF = n.br / halogens) %>% #select(PWSID, STG2_MCL, DBP, CAT, CAT.frac) %>% write_csv("Output/Data/Figure2c data_CAT fractions.csv")
  ggplot(aes(x = DBP, y = CAT.frac)) +
  geom_boxplot(aes(color = group), outlier.alpha = 0, outlier.size = .7) +
  geom_point(data = species.quantiles.CAT, aes(shape = as.factor(Quantile), color = group), size = 2) +
  facet_grid(cols = vars(subclass), scales = "free", space = "free", labeller = as_labeller(c("mHAA" = "Mono-HAAs", "dHAA" = "Di-HAAs", "tHAA" = "Tri-HAAs"))) +
  labs(shape = "Percentile", color = "Regulatory\nGroup", x = "", y = "Contributing Fraction to CAT")+
  scale_color_brewer(palette = "Set2")


p2a / p2b / p2c  +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Output/Figures/Manuscript Revision/Figure2.png", width = 6, height = 7, units = "in", dpi = 1000)


# Figure 3: Regulatory Groups vs CAT ####

## Plot A (HAA9)

p3a1 <- df.sum %>%
  filter(HAA.Group == "HAA9") %>%
  mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
  mutate(bin = case_when(
    CAT_HAA9_rank < 0.5 ~ "0-50",
    CAT_HAA9_rank < 0.9 ~ "50-90",
    CAT_HAA9_rank < 0.99 ~ "90-99",
    T ~ "99+"
    )
  ) %>% 
  ggplot(aes(x = c.mass_HAA9, y = CAT_HAA9, color = bin)) +
  geom_point(alpha = 0.3) +
  geom_vline(data = . %>% summarize(c.mass_HAA9 = quantile(c.mass_HAA9, c(0.95, 0.98))), aes(xintercept = c.mass_HAA9, linetype = as.factor(round(c.mass_HAA9, 1)))) +
  labs(x = expression(HAA9~Max~LAA~Concentration~(mu*g/L)), y = "Calculated Additive Toxicity", color = "HAA9 CAT\nPercentile", linetype = "HAA Concentration\nPercentile", subtitle = "All HAA9 Species") +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100, 10)) +
  scale_y_continuous(limits = c(0,6), breaks = seq(0,6,1)) +
  scale_linetype_discrete(labels = c("53.2" = "95th", "60.7" = "98th"))

  p3a2 <- df.sum_no.mHAA %>%
    filter(HAA.Group == "HAA9") %>%
    mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
    mutate(bin = case_when(
      CAT_HAA9_rank < 0.5 ~ "0-50",
      CAT_HAA9_rank < 0.9 ~ "50-90",
      CAT_HAA9_rank < 0.99 ~ "90-99",
      T ~ "99+"
      )
    ) %>% 
  ggplot(aes(x = c.mass_HAA9, y = CAT_HAA9, color = bin)) +
    geom_point(alpha = 0.3) +
    geom_vline(data = . %>% summarize(c.mass_HAA9 = quantile(c.mass_HAA9, c(0.95, 0.98))), aes(xintercept = c.mass_HAA9, linetype = as.factor(round(c.mass_HAA9, 1)))) +
    labs(x = expression(HAA9~Max~LAA~Concentration~(mu*g/L)), y = "Calculated Additive Toxicity", color = "HAA9 CAT\nPercentile", linetype = "HAA Concentration\nPercentile", subtitle = "Excluding mono-HAAs") +
    scale_x_continuous(limits = c(0,100), breaks = seq(0,100, 10)) +
    scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1.2,0.2)) +
    scale_linetype_discrete(labels = c("51.6" = "95th", "58.6" = "98th"))
  

## Plot B (HAA6Br)

p3b1 <- df.sum %>%
  filter(HAA.Group == "HAA6Br") %>%
  mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
  mutate(bin = case_when(
    CAT_HAA9_rank < 0.5 ~ "0-50",
    CAT_HAA9_rank < 0.9 ~ "50-90",
    CAT_HAA9_rank < 0.99 ~ "90-99",
    T ~ "99+"
    )
  ) %>%
  ggplot(aes(x = c.mass_HAA6Br, y = CAT_HAA9, color = bin)) +
  geom_point(alpha = 0.3) +
  geom_vline(data = . %>% summarize(c.mass_HAA6Br = quantile(c.mass_HAA6Br, c(0.95, 0.98))), aes(xintercept = c.mass_HAA6Br, linetype = as.factor(round(c.mass_HAA6Br, 1)))) +
  labs(x = expression(HAA6Br~Max~LAA~Concentration~(mu*g/L)), y = "Calculated Additive Toxicity", color = "HAA9 CAT\nPercentile", linetype = "HAA Concentration\nPercentile", subtitle = "All HAA6Br Species") +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100, 10)) +
  scale_y_continuous(limits = c(0,6), breaks = seq(0,6,1)) +
  scale_linetype_discrete(labels = c("20.6" = "95th", "27" = "98th"))

p3b2 <- df.sum_no.mHAA %>%
  filter(HAA.Group == "HAA6Br") %>%
  mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
  mutate(bin = case_when(
    CAT_HAA9_rank < 0.5 ~ "0-50",
    CAT_HAA9_rank < 0.9 ~ "50-90",
    CAT_HAA9_rank < 0.99 ~ "90-99",
    T ~ "99+"
    )
  ) %>%
  ggplot(aes(x = c.mass_HAA6Br, y = CAT_HAA9, color = bin)) +
  geom_point(alpha = 0.3) +
  geom_vline(data = . %>% summarize(c.mass_HAA6Br = quantile(c.mass_HAA6Br, c(0.95, 0.98))), aes(xintercept = c.mass_HAA6Br, linetype = as.factor(round(c.mass_HAA6Br, 1)))) +
  labs(x = expression(HAA6Br~Max~LAA~Concentration~(mu*g/L)), y = "Calculated Additive Toxicity", color = "HAA9 CAT\nPercentile", linetype = "HAA Concentration\nPercentile", subtitle = "Excluding mono-HAAs") +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100, 10)) +
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1.2,0.2)) +
  scale_linetype_discrete(labels = c("19.7" = "95th", "26" = "98th"))

(p3a1 | p3b1) / (p3a2 | p3b2) +
   plot_layout(guides = "collect") +
   plot_annotation(tag_levels = "a") &
   theme(plot.tag = element_text(face = 'bold'))

 ggsave("Output/Figures/Manuscript Revision/Figure3.png", width = 7, height = 5, units = "in", dpi = 1000)
 
 
# Figure 4: HAA species vs HAA9 ####
 df %>%
   filter(HAA.Group == "HAA9") %>%
   left_join(., HAA9.rank, by = "PWSID") %>%
   group_by(PWSID) %>%
   mutate(c.mass.frac = c.mass / sum(c.mass, na.rm = T)) %>%
   mutate(CAT.frac = CAT / sum(CAT, na.rm = T)) %>%
   select(PWSID, DBP, c.mass, c.mass.frac, CAT, CAT.frac, p.rank, c.mass_HAA9) %>% #filter(p.rank > 0.995) %>% group_by(DBP) %>% summarize(across(c.mass:CAT.frac, ~ median(., na.rm = T)))   
   mutate(p.rank = signif(p.rank, digits = 2)) %>%
   group_by(p.rank) %>% mutate(c.mass_HAA9 = mean(c.mass_HAA9, na.rm = T)) %>%
   mutate(group = case_when(
     DBP %in% c("DCAA", "TCAA") ~ "A",
     DBP %in% c("BCAA", "BDCAA", "DBCAA") ~ "B",
     T ~ "C"
   )) %>%
   group_by(PWSID, p.rank, group, c.mass_HAA9) %>%
   summarize(across(c.mass:CAT.frac, ~sum(., na.rm = T))) %>%
   group_by(group, p.rank, c.mass_HAA9) %>%
   summarize(c.mass.median = median(c.mass, na.rm = T),
             c.mass.frac = median(c.mass.frac, na.rm = T),
             CAT = median(CAT, na.rm = T),
             CAT.frac = median(CAT.frac, na.rm = T),
             c.mass.upper = quantile(c.mass, 0.9, na.rm  = T),
             c.mass.lower = quantile(c.mass, 0.1, na.rm = T),
             c.mass.max = max(c.mass, na.rm = T)
    ) %>%
   ggplot(aes(x = p.rank, y = c.mass.median)) +
   geom_text(data = tibble(group = c("A", "B", "C"), yval = c(65, 40, 60), p.rank = c(0,0,0), label = c("a","b","c")),
               aes(y = yval, label = label), fontface = "bold") +
   geom_line(size = 1) +
   geom_line(aes(y = c.mass.upper), alpha = 0.5) +
   geom_line(aes(y = c.mass.lower), alpha = 0.5) +
   geom_line(aes(y = c.mass.max), alpha = 0.5, linetype = "dotdash") +
   facet_wrap(vars(group), ncol = 1, scales = "free_y", labeller = as_labeller(c("A" = "DCAA + TCAA", "B" = "BCAA + BDCAA + DBCAA", "C" = "Other HAAs"))) +
   labs(x = "PWS HAA9 Percent Rank", y = expression(DBP~Concentration~(mu*g/L))) +
   scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
   theme(strip.text = element_text(size = rel(1)))
 
 
 ggsave("Output/Figures/Manuscript/Figure4.png", width = 4, height = 6, units = "in", dpi = 1000)
 
# Figure 5: HAA species vs HAA6Br ####
 
df %>%
   filter(HAA.Group == "HAA6Br") %>%
   left_join(., HAA6Br.rank, by = "PWSID") %>%
   filter(n.br > 0) %>%
   group_by(PWSID) %>%
   mutate(c.mass.frac = c.mass / sum(c.mass, na.rm = T)) %>%
   mutate(CAT.frac = CAT / sum(CAT, na.rm = T)) %>%
   select(PWSID, DBP, c.mass, c.mass.frac, CAT, CAT.frac, p.rank) %>% #filter(p.rank > 0.995) %>% group_by(DBP) %>% summarize(across(c.mass:CAT.frac, ~ median(., na.rm = T)))   
   mutate(p.rank = signif(p.rank, digits = 2)) %>%
   ungroup() %>%
   mutate(group = case_when(
     DBP %in% c("BCAA", "BDCAA") ~ "A",
     DBP %in% c("DBAA", "DBCAA", "TBAA") ~ "B",
     T ~ "C"
    )
   ) %>%
   group_by(PWSID, p.rank, group) %>%
   summarize(across(c.mass:CAT.frac, ~sum(., na.rm = T))) %>%
   group_by(group, p.rank) %>%
   summarize(c.mass.median = median(c.mass, na.rm = T),
             c.mass.frac = median(c.mass.frac, na.rm = T),
             CAT = median(CAT, na.rm = T),
             CAT.frac = median(CAT.frac, na.rm = T),
             c.mass.upper = quantile(c.mass, 0.9, na.rm  = T),
             c.mass.lower = quantile(c.mass, 0.1, na.rm = T),
             c.mass.max = max(c.mass, na.rm = T)
   ) %>%
   ggplot(aes(x = p.rank, y = c.mass.median)) +
   geom_text(data = tibble(group = c("A", "B", "C"), yval = c(30, 60, 6.5), p.rank = c(0,0,0), label = c("a","b","c")),
             aes(y = yval, label = label), fontface = "bold") +   geom_line(size = 1) +
   geom_line(aes(y = c.mass.upper), alpha = 0.5) +
   geom_line(aes(y = c.mass.lower), alpha = 0.5) +
   geom_line(aes(y = c.mass.max), alpha = 0.5, linetype = "dotdash") +
   facet_wrap(vars(group), ncol = 1, scales = "free_y", labeller = as_labeller(c("A" = "BCAA + BDCAA", "B" = "DBAA + DBCAA + TBAA", "C" = "MBAA"))) +
   labs(x = "PWS HAA6Br Percent Rank", y = expression(DBP~Concentration~(mu*g/L))) +
   scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
   theme(strip.text = element_text(size = rel(1)))
 
 
 ggsave("Output/Figures/Manuscript/Figure5.png", width = 4, height = 6, units = "in", dpi = 1000)
 
 
# Figure 6: Source water quality vs treatment used ####
 
 # Figure 6a Chloramines
p6a <- df.sum %>% ungroup() %>%
   filter(HAA.Group == "HAA9") %>%
   drop_na(Chloramines) %>%
   group_by(TOC.cat, Br.cat, Chloramines) %>%
   summarize(n = n()) %>% 
   group_by(TOC.cat, Br.cat) %>%
   mutate(rate = n/sum(n)) %>% 
   filter(Chloramines == T) %>% #filter(TOC.cat == "> 4")
   ggplot(aes(x = TOC.cat, y = Br.cat)) +
   geom_tile(aes(fill = rate), color = "white") +
   geom_label(aes(label = paste(round(rate*100,0), "%", sep = "")), fill = alpha(c("white"),0.7)) +
   scale_fill_viridis_b(limits = c(0,1), n.breaks = 10) +
   labs(x = expression(Source~TOC~(mg/L)), y = expression(Bromide~(mu*g/L)), title = "% of PWS Using Chloramines") +
   theme(legend.position = "none")
 
 # Figure 6b Advanced Treatment
 # Treatment decisions by TOC and Br  
p6b <- df.sum %>% ungroup() %>%
   filter(HAA.Group == "HAA9") %>%
   drop_na(Adv.Treat) %>% #filter(Bromide_ppb > 50 & TOC_ppm > 2) %>% group_by(Adv.Treat) %>% summarize(n = n())
   group_by(TOC.cat, Br.cat, Adv.Treat) %>%
   summarize(n = n()) %>% 
   group_by(TOC.cat, Br.cat) %>%
   mutate(rate = n/sum(n)) %>% 
   filter(Adv.Treat == T) %>% #filter(TOC.cat == "> 4")
   ggplot(aes(x = TOC.cat, y = Br.cat)) +
   geom_tile(aes(fill = rate), color = "white") +
   geom_label(aes(label = paste(round(rate*100,0), "%", sep = "")), fill = alpha(c("white"),0.7)) +
   scale_fill_viridis_b(limits = c(0,1), n.breaks = 10) +
   labs(x = expression(Source~TOC~(mg/L)), y = expression(Bromide~(mu*g/L)), title = "% of PWS Using GAC, IX, or Biofiltration") +
   theme(legend.position = "none") 

p6a + p6b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave("Output/Figures/Manuscript/Figure6.png", width = 7, height = 3, units = "in", dpi = 1000)


# Figure 7: Relationship between source WQ/treatment and Brominated HAA formation (and CAT) ####

p7a <- df.sum %>% ungroup() %>%
  filter(HAA.Group == "HAA6Br") %>%
  drop_na(Chloramines) %>% 
  filter(HAA6Br_high == T) %>% 
  group_by(TOC.cat, Br.cat, Chloramines) %>%
  summarize(value = n()) %>%
  mutate(alpha.cat = ifelse(value == 1, "Dim", "Full")) %>%
  ggplot(aes(x = TOC.cat, y = Br.cat)) +
  geom_tile(aes(fill = value, alpha = alpha.cat), color = "white") +
  geom_label(aes(label = round(value,1)), fill = alpha(c("white"),0.7)) +
  scale_fill_viridis_b(limits = c(0,30), n.breaks = 20) +
  labs(x = expression(Source~TOC~(mg/L)), y = expression(Bromide~(mu*g/L)), title = "# PWS Exceeding 95th Percentile HAA6Br") +
  theme(legend.position = "none", strip.text = element_text(size = rel(1.2))) +
  facet_wrap(vars(Chloramines), labeller = as_labeller(c("FALSE" = "No Chloramines", "TRUE" = "Chloramines"))) +
  scale_alpha_manual(values = c(0, 1))

p7b <- df.sum %>% ungroup() %>%
  filter(HAA.Group == "HAA6Br") %>%
  drop_na(Chloramines) %>%
  filter(CAT_high == T) %>% 
  group_by(TOC.cat, Br.cat, Chloramines) %>%
  summarize(value = n()) %>%
  mutate(alpha.cat = ifelse(value == 1, "Dim", "Full")) %>%
  ggplot(aes(x = TOC.cat, y = Br.cat)) +
  geom_tile(aes(fill = value, alpha = alpha.cat), color = "white") +
  geom_label(aes(label = round(value,1)), fill = alpha(c("white"),0.7)) +
  scale_fill_viridis_b(limits = c(0,30), n.breaks = 20) +
  labs(x = expression(Source~TOC~(mg/L)), y = expression(Bromide~(mu*g/L)), title = "# PWS Exceeding 95th Percentile CAT") +
  theme(legend.position = "none", strip.text = element_text(size = rel(1.2))) +
  facet_wrap(vars(Chloramines), labeller = as_labeller(c("FALSE" = "No Chloramines", "TRUE" = "Chloramines"))) +
  scale_alpha_manual(values = c(0, 1))

p7a / p7b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave("Output/Figures/Manuscript Revision/Figure7.png", width = 7, height = 7, units = "in", dpi = 1000)

