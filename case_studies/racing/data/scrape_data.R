library(rvest)
library(lubridate)

# Scrape races and filter for Map Randos
N_pages <- 200
suffixes <- c()
race_types <- c()

for (p in 1:N_pages) {
  url <- paste0('https://racetime.gg/smr?page=', p, '#past')
  races <- read_html(url) %>% html_elements("[class='race']")

  page_suffixes <- races %>% html_attr('href')
  suffixes <- c(suffixes, page_suffixes)
  
  page_race_types <- races %>% html_elements('.goal') %>% html_text2()
  race_types <- c(race_types, page_race_types)
}

race_info <- data.frame(suffixes, race_types)
valid_races <- race_info[race_info$race_types == "Map Rando",]

# Scrape individual race info
N_races <- nrow(valid_races)

race_urls <- c()
race_info <- c()
race_datetimes <- c()
race_N_entrants_f <- c()
race_N_entrants_dnf <- c()

race_f_start_idxs <- c()
race_f_end_idxs <- c()
f_idx <- 1

race_entrant_f_names <- c()
race_entrant_f_times <- c()

race_dnf_start_idxs <- c()
race_dnf_end_idxs <- c()
dnf_idx <- 1

race_entrant_dnf_names <- c()

seed_urls <- c()
seeds <- c()
versions <- c()
skill_assumptions <- c()
item_progressions <- c()

# Loop over races in chronological order
for (r in N_races:1) {
  print(paste("Processing race", r))
  
  race_url <- paste0('https://racetime.gg', valid_races$suffixes[r])
  race_html <- read_html(race_url) 
  seed_url <- race_html %>% html_elements('.info') %>% 
              html_elements('a') %>% html_attr('href')

  if (length(seed_url) == 0) next
  if (!grepl('https://maprando.com/seed/', seed_url, fixed=TRUE)) next
    
  seed_html <- read_html(seed_url) 
  seed_names <- seed_html %>% html_element('.card-body') %>% 
                html_elements('.row') %>% html_elements('.col-5') %>%
                html_text()
  seed_config <- seed_html %>% html_element('.card-body') %>% 
                 html_elements('.row') %>% html_elements('.col-7')
  seed <- seed_config[which("Seed name:" == seed_names)] %>% html_text()
  version <- seed_config[which("Version:" == seed_names)] %>% html_text2()
  skill_assumption <- seed_config[which("Skill assumption:" == seed_names)] %>% 
                      html_text2()
  item_progression <- seed_config[which("Item progression:" == seed_names)] %>% 
                      html_text2()
  
  if (skill_assumption != "Hard" | item_progression != "Tricky") next
  
  seed_urls <- c(seed_urls, seed_url)
  seeds <- c(seeds, seed)
  versions <- c(versions, version)
  skill_assumptions <- c(skill_assumptions, skill_assumption)
  item_progressions <- c(item_progressions, item_progression)
  
   
  race_urls <- c(race_urls, race_url)

  info <- race_html %>% html_elements('.info') %>% html_text()
  race_info <- c(race_info, info)
  
  datetime <- (race_html %>% html_elements('.datetime'))[1] %>% 
              html_attr('datetime')
  race_datetimes <- c(race_datetimes, datetime)
  
  entrants <- html_elements(race_html, '.race-entrants') 
  
  statuses <- entrants %>% html_elements('.status') %>% html_text2()
  f_filter <- statuses == "Finished"
  N_entrants_f <- sum(f_filter)
  race_N_entrants_f <- c(race_N_entrants_f, N_entrants_f)
  
  dnf_filter <- statuses == "DNF"
  N_entrants_dnf <- sum(dnf_filter)
  race_N_entrants_dnf <- c(race_N_entrants_dnf, N_entrants_dnf)
  
  names <- entrants %>% html_elements('.name') %>% html_text()
  race_entrant_f_names <- c(race_entrant_f_names, names[f_filter])
  race_entrant_dnf_names <- c(race_entrant_dnf_names, names[dnf_filter])

  times <- entrants %>% html_elements('.finish-time') %>% html_text()
  times <- times[c(rep(c(TRUE, FALSE), N_entrants_f), 
                   rep(FALSE, N_entrants_dnf))]
  times <- period_to_seconds(hms(times))
  race_entrant_f_times <- c(race_entrant_f_times, times)
  
  race_f_start_idxs <- c(race_f_start_idxs, f_idx)
  race_f_end_idxs <- c(race_f_end_idxs, f_idx + N_entrants_f - 1)
  f_idx <- f_idx + N_entrants_f
  
  if (N_entrants_dnf > 0) {
    race_dnf_start_idxs <- c(race_dnf_start_idxs, dnf_idx)
    race_dnf_end_idxs <- c(race_dnf_end_idxs, dnf_idx + N_entrants_dnf - 1)
    dnf_idx <- dnf_idx + N_entrants_dnf
  } else {
    race_dnf_start_idxs <- c(race_dnf_start_idxs, 0)
    race_dnf_end_idxs <- c(race_dnf_end_idxs, 0)
  }
}

# Convert names to factors
all_entrants <- c(race_entrant_f_names, race_entrant_dnf_names)
uniq_entrants <- sort(unique(all_entrants))
N_entrants <- length(uniq_entrants)

race_entrant_f_idxs <- factor(race_entrant_f_names, 
                              levels=uniq_entrants, 
                              labels=1:N_entrants)

race_entrant_dnf_idxs <- factor(race_entrant_dnf_names, 
                                levels=uniq_entrants, 
                                labels=seq_along(uniq_entrants))


# Write data to file
df <- data.frame("level" = 1:N_entrants, "name" = uniq_entrants)
write.csv(df, "entrant_level_defs.csv")

df <- data.frame(race_urls, race_datetimes,
                 race_N_entrants_f, race_N_entrants_dnf,
                 race_f_start_idxs, race_f_end_idxs,
                 race_dnf_start_idxs, race_dnf_end_idxs,
                 seed_urls, seeds, versions,
                 skill_assumptions, item_progressions)
write.csv(df, "race_info.csv")

df <- data.frame(race_entrant_f_idxs,
                 race_entrant_f_times)
write.csv(df, "race_entrant_f_info.csv")


df <- data.frame(race_entrant_dnf_idxs)
write.csv(df, "race_entrant_dnf_info.csv")
