# Organisation des donnees Rapsodyn pour Methode AIC (Akmouche 2019)
# Objectif: determiner le nombre de jours avec pic thermique apres la floraison du colza

library(dplyr)

# Parametres
T_base = 0 # 0 pour Colza
tmax1 = 25 # t max 25
tmax2 = 30 # t max 30

# Environnement de travail 
setwd("C:/Users/Magno Lethicia/Documents/7. COSMOS/6. AIC RAPSODYN/RAPSODYN_Data script R")

# Compile la fonction f_tmaxcount
source("rapsodyn_fonctions.R")

# Importation des donnees
df.plant <- read.csv(file = "data/plant.csv", sep = ";", header = T)
df.climat <- read.csv(file = "data/climat.csv", sep = ";", header = T)

# Conversion des dates dans un format manipulable
df.plant <- mutate(df.plant, date_sowing = as.Date(date_sowing, "%d/%m/%Y"))
df.plant <- mutate(df.plant, date_flowering = as.Date(date_flowering, "%d/%m/%Y"))
df.plant <- mutate(df.plant, date_harvest = as.Date(date_harvest, "%d/%m/%Y"))

df.climat$date <-as.Date(df.climat$date, "%d/%m/%Y")

# Calcul le TT_day
df.climat <- df.climat %>%
  mutate(TT_day = ifelse(temp_moy>T_base, temp_moy - T_base, 0))

# Debut de la boucle principale 
if (exists("out.all")){rm(out.all)}
for(i in 1:nrow(df.plant)){
  sel.year <- df.plant$year[i]
  sel.trial <- df.plant$local[i]
  sel.nitrogen <- df.plant$nitrogen[i]
  sel.genotype <- df.plant$genotype[i]
  sel.rep <- df.plant$rep[i]
  out <- f_tmaxcount(sel.year, sel.trial, sel.nitrogen, sel.genotype, sel.rep, df.plant, df.climat, tmax1, tmax2)
  if (exists("out.all")) {out.all <- rbind(out.all, out) } else {out.all <-  out}
  }

View(out.all)
write.csv(out.all,"out.csv", row.names = FALSE, quote = FALSE)
