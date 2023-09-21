# Rapsodyn -> fichier qui contient la fonction pour la recherche du nombre de jour avec t max par periode

f_tmaxcount <- function(sel.year, sel.trial, sel.nitrogen, sel.genotype, sel.rep, df.plant, df.climat, tmax1, tmax2){

  # date sowing
  date_sowing <- df.plant %>%
    filter(local == sel.trial) %>%
    filter(genotype == sel.genotype) %>%
    filter(rep == sel.rep) %>%
    select(date_sowing)

  # date flowering
  date_flowering <- df.plant %>%
    filter(local == sel.trial) %>%
    filter(genotype == sel.genotype) %>%
    filter(rep == sel.rep) %>%
    select(date_flowering)

  # date harvest
  date_harvest <- df.plant %>%
    filter(local == sel.trial) %>%
    filter(genotype == sel.genotype) %>%
    filter(rep == sel.rep) %>%
    select(date_harvest)

  # Calcul du TT_cum entre sowing et harvest
  df.climat.out <- df.climat %>%
    filter(trial == sel.trial) %>%
    filter(between(date, date_sowing, date_harvest)) %>%
    mutate(TT_cum = cumsum(TT_day))

  # Recherche TT_cum pour sowing, flowering et harvest
  TT_sowing <- df.climat.out %>%
    filter(date == date_sowing[1, 1]) %>%
    select(TT_cum)
  names(TT_sowing) <- "TT_sowing"

  TT_flowering <- df.climat.out %>%
    filter(date == date_flowering[1, 1]) %>%
    select(TT_cum)
  names(TT_flowering) <- "TT_flowering"

  TT_harvest <- df.climat.out %>%
    filter(date == date_harvest[1, 1]) %>%
    select(TT_cum)
  names(TT_harvest)<- "TT_harvest"

  # seed_yield, seed_number, PMG, lipid_content, protein_content
  analysis <- df.plant %>%
    filter(local == sel.trial) %>%
    filter(genotype == sel.genotype) %>%
    filter(rep == sel.rep) %>%
    select(seed_yield, seed_number, PMG, lipid_content, protein_content)


  # Calcul de TT_P300, TT_P600 et TT_P1000, puis recherche de la date correspondante
  
  TT_P300 <- TT_flowering + 300
  date_P300 <- df.climat.out %>%
    filter(TT_cum == nth(TT_cum, which.min(abs(TT_cum - TT_P300[1, 1])))) %>% # recherche la valeur la plus proche
    select(date)
  names(TT_P300) <- "TT_P300"
  names(date_P300) <- "date_P300"

  TT_P600 <- TT_flowering + 600
  date_P600 <- df.climat.out %>%
    filter(TT_cum == nth(TT_cum, which.min(abs(TT_cum - TT_P600[1, 1])))) %>% # recherche la valeur la plus proche
    select(date)
  names(TT_P600) <- "TT_P600"
  names(date_P600) <- "date_P600"

  TT_P1000 <- TT_flowering + 1000
  date_P1000 <- df.climat.out %>%
    filter(TT_cum == nth(TT_cum, which.min(abs(TT_cum - TT_P1000[1, 1])))) %>% # recherche la valeur la plus proche
    select(date)
  names(TT_P1000) <- "TT_P1000"
  names(date_P1000) <- "date_P1000"


  # Calcul du nombre de jours dans chaque intervalle superieur a tmax

  Ptot_tmax1 <- df.climat.out %>%
    filter(between(date, date_flowering[1,1], date_harvest[1,1])) %>%
    filter(temp_max>=tmax1) %>%
    nrow()
  
  Ptot_tmax2 <- df.climat.out %>%
    filter(between(date, date_flowering[1,1], date_harvest[1,1])) %>%
    filter(temp_max>=tmax2) %>%
    nrow()
  
  P300_tmax1 <- df.climat.out %>%
    filter(between(date, date_flowering[1,1], date_P300[1,1])) %>%
    filter(temp_max>=tmax1) %>%
    nrow()

  P300_tmax2 <- df.climat.out %>%
    filter(between(date, date_flowering[1,1], date_P300[1,1])) %>%
    filter(temp_max>=tmax2) %>%
    nrow()

  P600_tmax1 <- df.climat.out %>%
    filter(between(date, date_P300[1,1], date_P600[1,1])) %>%
    filter(temp_max>=tmax1) %>%
    nrow()

  P600_tmax2 <- df.climat.out %>%
    filter(between(date, date_P300[1,1], date_P600[1,1])) %>%
    filter(temp_max>=tmax2) %>%
    nrow()

  P1000_tmax1 <- df.climat.out %>%
    filter(between(date, date_P600[1,1], date_P1000[1,1])) %>%
    filter(temp_max>=tmax1) %>%
    nrow()

  P1000_tmax2 <- df.climat.out %>%
    filter(between(date, date_P600[1,1], date_P1000[1,1])) %>%
    filter(temp_max>=tmax2) %>%
    nrow()

  harvest_tmax1 <- df.climat.out %>%
    filter(between(date, date_P1000[1,1], date_harvest[1,1])) %>%
    filter(temp_max>=tmax1) %>%
    nrow()

  harvest_tmax2 <- df.climat.out %>%
    filter(between(date, date_P1000[1,1], date_harvest[1,1])) %>%
    filter(temp_max>=tmax2) %>%
    nrow()


  df.out <- data.frame(sel.year, sel.trial, sel.nitrogen, sel.genotype, sel.rep, TT_sowing, date_sowing, TT_flowering, date_flowering, TT_P300, date_P300, TT_P600, date_P600, TT_P1000, date_P1000, TT_harvest, date_harvest, tmax1, tmax2, Ptot_tmax1, Ptot_tmax2, P300_tmax1, P300_tmax2, P600_tmax1, P600_tmax2, P1000_tmax1, P1000_tmax2, harvest_tmax1, harvest_tmax2, analysis)

  return(df.out)
}
