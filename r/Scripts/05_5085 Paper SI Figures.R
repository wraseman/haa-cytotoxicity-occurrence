# Figure S1: CTI and GTI Values ####

# S1: Linear Scale
ps1a <- index %>%
  filter(class == "HAA") %>%
  select(DBP, subclass, group, CTI, GTI) %>%
  pivot_longer(CTI:GTI, names_to = "name", values_to = "value") %>%
  ggplot(aes(x = DBP, y = value, fill = group)) +
  geom_col() +
  geom_vline(xintercept = c(2.5, 5.5), alpha = 0.5) +
  facet_grid(cols = vars(name), labeller = as_labeller(c("CTI" = "Cytotoxicity Index (CTI)", "GTI" = "Genotoxicity Index (GTI)"))) +
  labs(x = "", y = "In Vitro Potency Index Value", fill = "HAA Regulatory\nGroup") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0,120), breaks = seq(0,120,20))

# S2: Log Scale
ps1b <- index %>%
  filter(class == "HAA") %>%
  select(DBP, subclass, group, CTI, GTI) %>%
  pivot_longer(CTI:GTI, names_to = "name", values_to = "value") %>%
  ggplot(aes(x = DBP, y = value, fill = group)) +
  geom_col() +
  geom_vline(xintercept = c(2.5, 5.5), alpha = 0.5) +
  facet_grid(cols = vars(name), labeller = as_labeller(c("CTI" = "Cytotoxicity Index (CTI)", "GTI" = "Genotoxicity Index (GTI)"))) +
  labs(x = "", y = "In Vitro Potency Index Value", fill = "HAA Regulatory\nGroup") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), breaks = c(0.1, 1, 10, 100))

ps1a / ps1b +
  plot_layout(guides = "collect") 

ggsave("Output/Figures/SI/FigureS1.png", width = 6, height = 4, units = "in", dpi = 1000)

# Figure S2: Sensitivity of HAA regulatory approaches to implicating PWS with high, moderate, and low CAT ####


ps2a1 <- tibble(pctile = seq(0.9, 0.99, 0.01)) %>%
  expand_grid(., df.sum %>% filter(HAA.Group == "HAA9")) %>%
  group_by(pctile) %>% 
  mutate(c.mass_HAA9_rank = percent_rank(c.mass_HAA9)) %>%
  mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
  mutate(bin = case_when(
    CAT_HAA9_rank < 0.5 ~ "0-50",
    CAT_HAA9_rank < 0.9 ~ "50-90",
    CAT_HAA9_rank < 0.99 ~ "90-99",
    T ~ "99+"
  )
  ) %>%
  select(pctile, PWSID, c.mass_HAA9_rank, CAT_HAA9_rank, bin) %>%
  group_by(pctile, bin) %>%
  mutate(Exceed.HAA9.conc = ifelse(c.mass_HAA9_rank > pctile, T, F)) %>%
  group_by(pctile, bin, Exceed.HAA9.conc) %>% 
  summarize(n = n()) %>%
  group_by(pctile, bin) %>%
  mutate(n = n/sum(n)) %>%
  filter(Exceed.HAA9.conc == T) %>% #filter(bin == "90-99")
  ggplot(aes(x = pctile, y = n, color = bin)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0.9,0.99), breaks = seq(0.9,0.99,0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1), labels = scales::percent_format()) +
  labs(x = "HAA9 Concentration Percentile", y = "Percent of UCMR4 PWS Implicated", color = "HAA9 CAT Percentile") +
  theme(legend.position = "none")


ps2a2 <- tibble(pctile = seq(0.9, 0.99, 0.01)) %>%
  expand_grid(., df.sum %>% filter(HAA.Group == "HAA9")) %>%
  group_by(pctile) %>% 
  mutate(c.mass_HAA9_rank = percent_rank(c.mass_HAA9)) %>%
  mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
  mutate(bin = case_when(
    CAT_HAA9_rank < 0.5 ~ "0-50",
    CAT_HAA9_rank < 0.9 ~ "50-90",
    CAT_HAA9_rank < 0.99 ~ "90-99",
    T ~ "99+"
  )
  ) %>%
  select(pctile, PWSID, c.mass_HAA9_rank, CAT_HAA9_rank, bin) %>%
  group_by(pctile, bin) %>%
  mutate(Exceed.HAA9.conc = ifelse(c.mass_HAA9_rank > pctile, T, F)) %>%
  group_by(pctile, bin, Exceed.HAA9.conc) %>% 
  summarize(n = n()) %>%
  group_by(pctile, bin) %>%
  filter(Exceed.HAA9.conc == T) %>% #filter(bin == "90-99")
  ggplot(aes(x = pctile, y = n, color = bin)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0.9,0.99), breaks = seq(0.9,0.99,0.01)) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0,400,50)) +
  labs(x = "HAA9 Concentration Percentile", y = "Number of UCMR4 PWS Implicated", color = "HAA9 CAT Percentile") +
  theme(legend.position = "none")

ps2b1 <- tibble(pctile = seq(0.9, 0.99, 0.01)) %>%
  expand_grid(., df.sum %>% filter(HAA.Group == "HAA6Br")) %>%
  group_by(pctile) %>% 
  mutate(c.mass_HAA6Br_rank = percent_rank(c.mass_HAA6Br)) %>%
  mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
  mutate(bin = case_when(
    CAT_HAA9_rank < 0.5 ~ "0-50",
    CAT_HAA9_rank < 0.9 ~ "50-90",
    CAT_HAA9_rank < 0.99 ~ "90-99",
    T ~ "99+"
  )
  ) %>%
  select(pctile, PWSID, c.mass_HAA6Br_rank, CAT_HAA9_rank, bin) %>%
  group_by(pctile, bin) %>%
  mutate(Exceed.HAA6Br.conc = ifelse(c.mass_HAA6Br_rank > pctile, T, F)) %>%
  mutate(pctile = as_factor(pctile), bin = as_factor(bin), Exceed.HAA6Br.conc = as_factor(Exceed.HAA6Br.conc)) %>%
  group_by(pctile, bin, Exceed.HAA6Br.conc, .drop = F) %>% 
  summarize(n = n()) %>%
  group_by(pctile, bin) %>%
  mutate(n = n/sum(n)) %>%
  filter(Exceed.HAA6Br.conc == T) %>% #filter(bin == "90-99")
  mutate(pctile = as.numeric(as.character(pctile))) %>%
  mutate(n = ifelse(n == 0, -1, n)) %>%
  ggplot(aes(x = pctile, y = n, color = bin)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0.9,0.99), breaks = seq(0.9,0.99,0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1), labels = scales::percent_format()) +
  labs(x = "HAA6Br Concentration Percentile", y = "Percent of UCMR4 PWS Implicated", color = "HAA9 CAT Percentile")

ps2b2 <- tibble(pctile = seq(0.9, 0.99, 0.01)) %>%
  expand_grid(., df.sum %>% filter(HAA.Group == "HAA6Br")) %>%
  group_by(pctile) %>% 
  mutate(c.mass_HAA6Br_rank = percent_rank(c.mass_HAA6Br)) %>%
  mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
  mutate(bin = case_when(
    CAT_HAA9_rank < 0.5 ~ "0-50",
    CAT_HAA9_rank < 0.9 ~ "50-90",
    CAT_HAA9_rank < 0.99 ~ "90-99",
    T ~ "99+"
  )
  ) %>%
  select(pctile, PWSID, c.mass_HAA6Br_rank, CAT_HAA9_rank, bin) %>%
  group_by(pctile, bin) %>%
  mutate(Exceed.HAA6Br.conc = ifelse(c.mass_HAA6Br_rank > pctile, T, F)) %>%
  mutate(pctile = as_factor(pctile), bin = as_factor(bin), Exceed.HAA6Br.conc = as_factor(Exceed.HAA6Br.conc)) %>%
  group_by(pctile, bin, Exceed.HAA6Br.conc, .drop = F) %>% 
  summarize(n = n()) %>%
  group_by(pctile, bin) %>%
  #  mutate(n = n/sum(n)) %>%
  filter(Exceed.HAA6Br.conc == T) %>% #filter(bin == "90-99")
  mutate(pctile = as.numeric(as.character(pctile))) %>%
  mutate(n = ifelse(n == 0, -1, n)) %>%
  ggplot(aes(x = pctile, y = n, color = bin)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0.9,0.99), breaks = seq(0.9,0.99,0.01)) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0,400,50)) +
  labs(x = "HAA6Br Concentration Percentile", y = "Number of UCMR4 PWS Implicated", color = "HAA9 CAT Percentile")


ps2a1 / ps2a2 | ps2b1 / ps2b2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Output/Figures/SI/FigureS2.png", width = 7, height = 6, units = "in", dpi = 1000)

# Figure S3: Relationship between HAA9 CAT and MBAA concentration ####
df.sum %>%
  filter(HAA.Group == "HAA9") %>%
  mutate(CAT_HAA9_rank = percent_rank(CAT_HAA9)) %>%
  mutate(bin = case_when(
    CAT_HAA9_rank < 0.5 ~ "0-50",
    CAT_HAA9_rank < 0.9 ~ "50-90",
    CAT_HAA9_rank < 0.99 ~ "90-99",
    T ~ "99+"
    )
  )  %>%
  ggplot(aes(x = c.mass_MBAA, y = CAT_HAA9)) +
  geom_point(aes(color = bin), alpha = 0.3) +
  labs(x = expression(MBAA~Concentration~(mu*g/L)), y = "HAA9 CAT", color = "HAA9 CAT %ile") +
  scale_x_continuous(limits = c(0,8), breaks = seq(0,8, 1)) +
  scale_y_continuous(limits = c(0,6), breaks = seq(0,6,1))

ggsave("Output/Figures/SI/FigureS3.png", width = 4, height = 3, units = "in", dpi = 1000)


# Figure S4: Species contribution to HAA9 as a function of HAA9 ####

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
    DBP %in% c("DCAA", "TCAA") ~ "DCAA & TCAA",
    DBP %in% c("BCAA", "BDCAA", "DBCAA") ~ "BCAA, BDCAA, & DBCAA",
    T ~ "Other HAAs"
  )) %>%
  mutate(group = factor(group, levels = c("DCAA & TCAA", "BCAA, BDCAA, & DBCAA", "Other HAAs"))) %>%
  # group_by(PWSID, p.rank, group, c.mass_HAA9) %>%
  # summarize(across(c.mass:CAT.frac, ~sum(., na.rm = T))) %>%
  group_by(DBP, group, p.rank, c.mass_HAA9) %>%
  summarize(c.mass.median = median(c.mass, na.rm = T),
            c.mass.frac = median(c.mass.frac, na.rm = T),
            CAT = median(CAT, na.rm = T),
            CAT.frac = median(CAT.frac, na.rm = T),
            c.mass.upper = quantile(c.mass, 0.9, na.rm  = T),
            c.mass.lower = quantile(c.mass, 0.1, na.rm = T),
            c.mass.max = max(c.mass, na.rm = T)
  ) %>%
  ggplot(aes(x = p.rank, y = c.mass.median, color = group)) +
  #geom_vline(xintercept = 0.975, linetype = "dashed") +
  geom_line(size = 1) +
  geom_line(aes(y = c.mass.upper), alpha = 0.8, size = 0.3) +
  geom_line(aes(y = c.mass.max), alpha = 0.8, size = 0.3, linetype = "dotdash") +
  facet_wrap(vars(DBP), scales = "free_y") +
 # geom_hline(data = tibble(group = c("BAA", "BCAA", "BDCAA", "CAA", "DBAA", "DBCAA", "DCAA", "TBAA", "TCAA"), yint = c(8, 20, 20, 30, 40, 14, 50, 25, 50)), aes(yintercept = yint), alpha = 0) +
  labs(x = "PWS HAA9 Percent Rank", y = expression(Species~Concentration~(mu*g/L)), color = "") +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2")

ggsave("Output/Figures/SI/FigureS4.png", width = 7, height = 6, units = "in", dpi = 1000)


# Figure S5: Species contribution to HAA6Br as a function of HAA6Br ####
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
    DBP %in% c("BCAA", "BDCAA") ~ "BCAA & BDCAA",
    DBP %in% c("DBAA", "DBCAA", "TBAA") ~ "DBAA, DBCAA, & TBAA",
    T ~ "BAA"
  )
  ) %>%
  mutate(group = factor(group, levels = c("BCAA & BDCAA", "DBAA, DBCAA, & TBAA", "BAA"))) %>%
  # group_by(PWSID, p.rank, group) %>%
  # summarize(across(c.mass:CAT.frac, ~sum(., na.rm = T))) %>%
  group_by(DBP, group, p.rank) %>%
  summarize(c.mass.median = median(c.mass, na.rm = T),
            c.mass.frac = median(c.mass.frac, na.rm = T),
            CAT = median(CAT, na.rm = T),
            CAT.frac = median(CAT.frac, na.rm = T),
            c.mass.upper = quantile(c.mass, 0.9, na.rm  = T),
            c.mass.lower = quantile(c.mass, 0.1, na.rm = T),
            c.mass.max = max(c.mass, na.rm = T)
  ) %>%
  ggplot(aes(x = p.rank, y = c.mass.median, color = group)) +
  #geom_vline(xintercept = 0.975, linetype = "dashed") +
  geom_line(size = 1) +
  geom_line(aes(y = c.mass.upper), alpha = 0.8, size = 0.3) +
  geom_line(aes(y = c.mass.max), alpha = 0.8, size = 0.3, linetype = "dotdash") +
  facet_wrap(vars(DBP), scales = "free_y") +
  # geom_hline(data = tibble(group = c("BAA", "BCAA", "BDCAA", "CAA", "DBAA", "DBCAA", "DCAA", "TBAA", "TCAA"), yint = c(8, 20, 20, 30, 40, 14, 50, 25, 50)), aes(yintercept = yint), alpha = 0) +
  labs(x = "PWS HAA6Br Percent Rank", y = expression(Species~Concentration~(mu*g/L)), color = "") +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2")

ggsave("Output/Figures/SI/FigureS5.png", width = 7, height = 4, units = "in", dpi = 1000)

#Figure S6: Source water and treatment characteristics vs HAA9 ####

df.sum %>% ungroup() %>%
  filter(HAA.Group == "HAA9") %>%
  drop_na(Chloramines) %>% 
  filter(HAA9_high == T) %>% 
  group_by(TOC.cat, Br.cat, Chloramines) %>%
  summarize(value = n()) %>%
  mutate(alpha.cat = ifelse(value == 1, "Dim", "Full")) %>%
  ggplot(aes(x = TOC.cat, y = Br.cat)) +
  geom_tile(aes(fill = value, alpha = alpha.cat), color = "white") +
  geom_label(aes(label = round(value,1)), fill = alpha(c("white"),0.7)) +
  #theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_viridis_b(limits = c(0,30), n.breaks = 20) +
  labs(x = expression(Source~TOC~(mg/L)), y = expression(Bromide~(mu*g/L)), title = "# PWS Exceeding 95th Percentile HAA9") +
  theme(legend.position = "none") +
  facet_wrap(vars(Chloramines), labeller = as_labeller(c("FALSE" = "No Chloramines", "TRUE" = "Chloramines"))) +
  scale_alpha_manual(values = c(0, 1))

ggsave("Output/Figures/SI/FigureS6.png", width = 6, height = 3, units = "in", dpi = 1000)

