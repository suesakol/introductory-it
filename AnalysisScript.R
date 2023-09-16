##########################################################################
##### Usage events and L2 knowledge of Introductory it construction  #####
#####               Supplementary materials: R script                #####
##########################################################################



# 1. Preamble ----------------------------------------------------------------

# (1) Load necessary packages

# If packages not yet installed, run: 
# install.packages(c("tidyverse", "lme4", "lmerTest", "MuMin", "magrittr", "gridExtra"))

library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(magrittr)
library(gridExtra)




# (2) Create a custom function: Normalized entropy: Hnorm (after Gries 2021)

hnorm <- function(x){
  percentage <- x/sum(x)
  find_hnorm <- -sum(percentage * log(percentage))/log(length(percentage))
  find_hnorm
  }








# 2. Corpus data: Analysis ---------------------------------------------------

# (1) File import

adj_freq <- read_delim(file = "./Data/Adjective_Frequency_Raw.txt", 
                       delim = "\t") |> 
  rename(Frequencies = n)


# (2) Descriptive information 

# - instances of each variant

adj_freq |> 
  group_by(Variant) |> 
  summarize(
    sum = sum(Frequencies)
    )

# - calculate type and token frequencies

adj_sm <- adj_freq |> 
  filter(Variant != "Other")

adj_sm |> 
  distinct(Adjectives) |> 
  count()

adj_sm |> 
  group_by(Variant) |> 
  summarize(
    type = n_distinct(Adjectives),
    n    = sum(Frequencies), 
    ttr  = (type/n) *100,
    ) |> 
  ungroup()


# - calculate normalized entropy

adj_sm |> 
  group_by(Variant) |> 
  summarize(
    hnorm = hnorm(Frequencies)
    )


rm(adj_sm)


# (3) Data preparation and ΔP calculation

# - pivot the table to a wide format (each column = raw frequency)

adj_freq <- adj_freq |>  
  mutate(
    Frequencies = as.double(Frequencies)
    ) |>  
  pivot_wider(
    names_from = Variant, 
    values_from = Frequencies
    ) 


# - Identify adjectives that appeared in *other* patterns
#   NOTE: These *other* patterns have introductory-it as subject of 1st clause
#         and adjective as predicate lemma

adj_not_target <- adj_freq |>  
  filter(if_all(.cols = c(Adj_to, Adj_that),
                .fns  = ~is.na(.)
                )
         ) |>  
  pull(Adjectives)


# - replace NA with zero in adj_freq tibble

adj_freq <- adj_freq |> 
  mutate(across(.cols = where(is.numeric),
                .fns = ~if_else(is.na(.), 0, .)
                )
         )


# - sort frequencies

adj_freq |>
  select(!Other) |> 
  arrange(desc(Adj_to)) |> 
  head(15)

adj_freq |>
  select(!Other) |> 
  arrange(desc(Adj_that)) |> 
  head(15)


# - count shared types

adj_freq |>  
  filter(Adj_that > 0 & Adj_to > 0) |> 
  summarize(
    types = n_distinct(Adjectives)
  )


# - count unique types and hapaxes

adj_freq |> 
  filter(Adj_that > 0 & Adj_to == 0) |> 
  summarize(
    types = n_distinct(Adjectives)
    )

adj_freq |>  
  filter(Adj_that == 1 & Adj_to == 0) |>  
  summarize(
    total = n_distinct(Adjectives)
    )


adj_freq |> 
  filter(Adj_to > 0 & Adj_that == 0) |> 
  summarize(
    types = n_distinct(Adjectives)
    )

adj_freq |>  
  filter(Adj_to == 1 & Adj_that == 0) |>  
  summarize(
    total = n_distinct(Adjectives)
    )


# - ΔP calculation
# --- calculate raw co-occurrence frequencies in a two-by-two table

adj_freq <- adj_freq |>  
  rename_with(.cols = starts_with("Adj_"),
              .fn = ~str_replace(., "Adj_", "freq_adj")
              ) |> 
  rename_with(.fn = tolower) |> 
  rename(freq_other = other) |>  
  mutate(
    total = sum(freq_adjto) + sum(freq_adjthat) + sum(freq_other),
    
    #variant: adj_to
    freq_not_adjto = freq_adjthat + freq_other,
    other_words_adjto = sum(freq_adjto) - freq_adjto,
    all_but_adjto = total - freq_adjto - freq_not_adjto - other_words_adjto,
    
    sum_adjto = sum(freq_adjto),
    sum_not_adjto = sum(freq_adjthat + freq_other),
    sum_adj_total_adjto = freq_adjto + freq_not_adjto,
    sum_notadj_total_adjto = other_words_adjto + all_but_adjto,
    not_adj_in_adjto = sum_adjto - sum_adj_total_adjto,
    not_adj_not_adjto = sum_not_adjto - sum_adj_total_adjto,
    
    #variant: adj_that
    freq_not_adjthat = freq_adjto + freq_other,
    other_words_adjthat = sum(freq_adjthat) - freq_adjthat,
    all_but_adjthat = total - freq_adjthat - freq_not_adjthat - other_words_adjthat,
    
    sum_adjthat = sum(freq_adjthat),
    sum_not_adjthat = sum(freq_adjto + freq_other),
    sum_adj_total_adjthat = freq_adjthat + freq_not_adjthat,
    sum_notadj_total_adjthat = other_words_adjthat + all_but_adjthat,
    not_adj_in_adjthat = sum_adjthat - sum_adj_total_adjthat,
    not_adj_not_adjthat = sum_not_adjthat - sum_adj_total_adjthat,
    
    zero = 0
    )


# --- calculate Delta P values for each variant

adj_freq <- adj_freq |>  
  mutate(
    
    #Adj-to
    adjto_DPwc = (freq_adjto / (freq_adjto + freq_not_adjto)) - 
      (other_words_adjto / (other_words_adjto + all_but_adjto)),
    
    adjto_DPcw = (freq_adjto / (freq_adjto + other_words_adjto)) - 
      (freq_not_adjto / (freq_not_adjto + all_but_adjto)),
    
    adjto_DPcw_upp = (sum_adj_total_adjto / (sum_adj_total_adjto + not_adj_in_adjto)) - 
      (zero / (zero + sum_not_adjto)),
    
    adjto_DPcw_low = (zero / (zero + sum_adjto)) - 
      (sum_adj_total_adjto / (sum_adj_total_adjto + not_adj_not_adjto)),
    
    adjto_DPcw_adjusted = (adjto_DPcw - adjto_DPcw_low) / (adjto_DPcw_upp - adjto_DPcw_low),
    
    
    #Adj-that
    adjthat_DPwc = (freq_adjthat / (freq_adjthat + freq_not_adjthat)) - 
      (other_words_adjthat / (other_words_adjthat + all_but_adjthat)),
    
    adjthat_DPcw = (freq_adjthat / (freq_adjthat + other_words_adjthat)) - 
      (freq_not_adjthat / (freq_not_adjthat + all_but_adjthat)),

    adjthat_DPcw_upp = (sum_adj_total_adjthat / (sum_adj_total_adjthat + not_adj_in_adjthat)) - 
      (zero / (zero + sum_not_adjthat)),
    
    adjthat_DPcw_low = (zero / (zero + sum_adjthat)) - 
      (sum_adj_total_adjthat / (sum_adj_total_adjthat + not_adj_not_adjthat)),
    
    adjthat_DPcw_adjusted = (adjthat_DPcw - adjthat_DPcw_low) / (adjthat_DPcw_upp - adjthat_DPcw_low)
    ) |>  
  mutate(across(.cols = adjto_DPwc:adjthat_DPcw_adjusted,
                .fns  = ~round(., digits = 3)
                )
         )


# --- mutate new column that sums each adjective's total raw frequency
# --- keep relevant columns

adj_freq <- adj_freq |>  
  mutate(
    freq_total = freq_adjto + freq_adjthat + freq_other
    ) |>  
  relocate(freq_total, .after = freq_other) |>  
  select(
    adjectives, freq_adjto, freq_adjthat, freq_other, freq_total, 
    ends_with("_DPwc"),
    ends_with("_DPcw"), 
    ends_with("adjusted")
    )


# --- remove adjectives not in the two target variants

adj_freq <- adj_freq |>  
  filter(!(adjectives %in% adj_not_target)
         )

rm(adj_not_target)


# --- sort adjectives by their Delta P scores

adj_freq |> 
  arrange(desc(adjthat_DPcw), desc(adjthat_DPcw_adjusted)) |>  
  select(
    adjectives, 
    starts_with("freq"), 
    starts_with("adjthat_")
    ) |>  
  head(15)

adj_freq |>  
  arrange(desc(adjto_DPcw), desc(adjto_DPcw_adjusted)) |>  
  select(
    adjectives, 
    starts_with("freq"), 
    starts_with("adjto_")
    ) |> 
  head(15)








# 3. Production data: Analysis --------------------------------------------

# (1) Data preparation

# --- lowercase answers & remove spaces in the column Response
# --- make intro-it variants in the column Variant consistent with corpus file

produc_dat <- read_csv("./Data/Production_introIT.csv") |> 
  mutate(
    Response = str_to_lower(Response),
    Response = str_replace(Response, " ", ""),
    Variant  = if_else(Variant == "Adj-that", "Adj_that", "Adj_to")
    )


# --- convert NA in gender to Not specified

produc_dat <- produc_dat |>  
  mutate(
    Sex = if_else(is.na(Sex) == TRUE, "NotSpec", Sex)
    )


# (2) Summarize demographic data: age, gender, education, and L2 proficiency

produc_dat |>  
  distinct(Participant, .keep_all = TRUE) |>  
  group_by(L1) |>  
  summarize(
    m = mean(Age, na.rm = TRUE),
    sd = sd(Age, na.rm = TRUE)
    ) |>  
  ungroup()

produc_dat |>  
  distinct(Participant, .keep_all = TRUE) |>  
  group_by(L1, Sex) |> 
  count() |> 
  ungroup()

produc_dat |>  
  distinct(Participant, .keep_all = TRUE) |>  
  group_by(L1, Degree) |>  
  count() |> 
  ungroup()

produc_dat |>  
  filter(L1 == "Thai") |>  
  distinct(Participant, .keep_all = TRUE) |>  
  summarize(
    m = mean(TestScore, na.rm = TRUE),
    sd = sd(TestScore, na.rm = TRUE),
    min = min(TestScore, na.rm = TRUE),
    max = max(TestScore, na.rm = TRUE) 
    )


# (3) Summarize experiment information

# --- % of incomplete trials (misspelled, incomplete, non-adj responses)

produc_dat |>  
  filter(Incomplete == 1) |>  
  summarize(
    n = n(),
    pt = ( n/nrow(produc_dat) ) * 100
    )

produc_dat |>  
  filter(Incomplete == 1) |> 
  group_by(L1) |>  
  summarize(n = n()) |> 
  ungroup()


# --- % of past-participle responses

produc_dat |>  
  filter(Incomplete != 1) |>  
  count(L1, Past_participle) |>  
  ungroup() |>  
  select(L1, n) |>  
  group_by(L1) |>  
  mutate(
    percent = (n/sum(n)) * 100
    )


# --- count of remaining total

produc_dat |> 
  filter(Past_participle == 0 & Incomplete == 0) |> 
  group_by(L1) |>  
  summarize(total = n() ) |>  
  ungroup()


# --- average responses per variant and L1

produc_dat |>  
  filter(Past_participle == 0 & Incomplete == 0) |>  
  count(L1, Participant, Variant) |> 
  ungroup() |> 
  group_by(L1, Variant) |> 
  summarize(
    m  = mean(n),
    sd = sd(n)
    ) |> 
  ungroup()


# --- number of types (per variant or in total)

produc_dat |>  
  filter(Past_participle == 0 & Incomplete == 0) |> 
  select(Variant, Response) |> 
  group_by(Variant) |>  
  summarize(
    type = n_distinct(Response)
    ) |>  
  ungroup()

produc_dat |>  
  filter(Past_participle == 0 & Incomplete == 0) |> 
  select(Variant, Response) |>  
  summarize(type_total = n_distinct(Response)) |>  
  ungroup()


# (4) Remove unwanted items and count 

# --- remove participles and incomplete items

dat <- produc_dat |>  
  filter(Past_participle == 0 & Incomplete == 0) |>  
  select(!(Past_participle:Incomplete) )


# --- count Elicited adjectives
# --- in total, 281 types with five most frequent items being
# --- important (53), necessary (29), obvious (26), likely (25), possible (24)

adj_list <- dat |>  
  count(Response, Variant) |>  
  pivot_wider(names_from = Variant, 
              values_from = n
              ) |>  
  mutate(across(.cols = where(is.numeric),
                .fns  = ~as.double(.)
                )
         ) |>  
  mutate(across(.cols = where(is.numeric),
                .fns  = ~if_else(is.na(.), 0, .)
                )
         )

adj_list |>  
  arrange(-Adj_to)

adj_list |>  
  arrange(-Adj_that)


rm(adj_list)


# --- count number of items per participants

dat |>  
  count(Participant) |>  
  mutate(perc = (n/sum(n)) * 100) |>  
  arrange(desc(n)) |>  
  head(3)

dat |>  
  count(Participant) |>  
  mutate(perc = (n/sum(n)) * 100) |>  
  arrange(desc(n)) |>  
  tail(3)








# 4. Mixed-effects modeling -----------------------------------------------

# (1) Join two tibbles: elicitation and corpus data

# NOTE: a column "cue" created to indicate if each adjective is cued by Adj_to 
#       or Adj_that. DPwc scores, not DPcw, used as a test condition below as
#       some items with low frequencies (e.g., abstract) had the same DPcw for 
#       two variants, after rounding

item_dat <- left_join(dat, 
                      adj_freq |>  
                        mutate(cue = if_else(adjto_DPwc > adjthat_DPwc, 
                                             "Adj_to", "Adj_that"
                                             )
                               ), 
                      by = c("Response" = "adjectives") 
                      )


# --- mutate new columns:
# (1) DPcw (delta P w/ construction as cue). If a word is cued by Adj-to, 
#     assign DPcw_adjusted of that word in Adj-to, else Adj-that
# (2) frequencies in "preferred" variant. If a word is cued by Adj-to, assign
#     the frequency of that word in that variant

# This is irrespective of the frames in which participants supplied adjectives

item_dat <- item_dat |> 
  mutate(
    DPcw         = if_else(cue == "Adj_that", 
                           adjthat_DPcw_adjusted, 
                           adjto_DPcw_adjusted
                           ),
    freq_variant = if_else(cue == "Adj_that", 
                           freq_adjthat, 
                           freq_adjto
                           )
    ) |>  
  select(!c(freq_adjto:freq_other, 
            adjto_DPwc:adjthat_DPcw_adjusted
            )
         ) |> 
  relocate(c(freq_variant, DPcw), .after = freq_total)


# --- count answers not attested in corpus

item_dat |> 
  filter(if_all(.cols = freq_total:cue,
                .fns  = ~is.na(.))
         ) |>  
  summarize(n = n())


item_dat |>  
  filter(if_all(.cols = freq_total:cue, 
                .fns  = ~is.na(.)
                )
         ) |>  
  select(Response) |>  
  count(Response) |>  
  arrange(desc(n))


# --- remove responses not in corpus data
# --- use NAs in columns freq_total to cue to identify those responses 

item_dat_sm <- item_dat |>  
  filter(if_all(.cols = freq_total:cue,
                .fns  = ~!is.na(.)
                )
         )


# --- visualize: histograms of DP

item_dat_sm |>  
  filter(cue == "Adj_that") |> 
  ggplot(aes(x = DPcw)) +
  geom_histogram()

item_dat_sm |>  
  filter(cue == "Adj_to") |> 
  ggplot(aes(x = DPcw)) +
  geom_histogram()


# --- mutate outcome variable: "Match"
# --- Did each elicited response "fit" the frame?

item_dat_sm <- item_dat_sm |>  
  mutate(
    match = if_else(Variant == cue, 1, 0)
    )


# --- count number of types & tokens and calculate TTR
# --- 183 types remain; 63 types = hapax legomena

item_dat_sm |>  
  summarize(
    type = n_distinct(Response) 
    )

item_dat_sm |>  
  count(Response) |>  
  ungroup() |>  
  arrange(desc(n)) |>  
  filter(n == 1) |>  
  summarize(
    hapax = n_distinct(Response)
    )

item_dat_sm |> 
  group_by(L1, Variant) |>  
  summarize(
    type = n_distinct(Response),
    n    = n(),
    TTR  = (type/n)
    ) |>  
  ungroup()


# --- obtain crosstabs
# --- all combinations present; no counts of zero
# --- evidently clear that Adj_that had more zero
# --- and this mostly down to answers from Thai L1 participants

item_dat_sm |>  
  count(L1, match)

item_dat_sm |> 
  count(Variant, match)

item_dat_sm |>  
  count(L1, Variant, match) 

item_dat_sm |> 
  count(Degree, match)

item_dat_sm |>  
  count(L1, Degree, match)


# --- visualize data

item_dat_sm |> 
  distinct(Participant, .keep_all = TRUE) |>  
  ggplot(aes(x = AWEQ)) +
  geom_histogram()

# Left skewed: score of 5 most frequent


ggplot(data = item_dat_sm, aes(x = Variant, y = match)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) 

ggplot(data = item_dat_sm, aes(x = Variant, y = match, color = L1)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1)

ggplot(data = item_dat_sm, aes(x = L1, y = match, color = Variant)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  facet_wrap(~ Degree)


# --- mutate new columns
# --- 1) Add order in which words appear
# --- 2) Add type of prompt (is vs. seem)

item_dat_sm <- item_dat_sm |> 
  group_by(Participant, Variant) |>  
  mutate(
    Final_order = row_number()
    ) |>  
  ungroup() |>  
  relocate(Final_order, .before = Response) |>  
  mutate(
    Link_verb  = str_extract(Frame, pattern = ".+?(?=\\.)")
    ) |>  
  relocate(Link_verb, .after = Variant)


# --- visualize: trial number normally distributed looking
# --- (though the shape of responses for Adj-that didn't look normal-ish)
# --- we will not transform this variable

item_dat_sm |>  
  count(Participant, Variant) |>  
  ggplot(aes(x = n)) +
  geom_histogram()


# --- check frequencies of adjectives in production data
# --- combine those hapax legomena into one category "other" in Response_rnd

item_dat_sm <- item_dat_sm |> 
  group_by(Response) |>  
  mutate(n = n()) |>  
  mutate(
    Response_rnd = if_else(n == 1, "other", Response)
    ) |>  
  ungroup() |> 
  relocate(Response_rnd, .after = Response) |>  
  select(!n)


# --- standardize AWEQ
# --- here we used the entire data set to calculate mean() and sd()
# --- another option would be to use distinct participants for calculation
# --- the values obtained would differ only in 2nd or 3rd decimal place

# NOTE: we didn't log(AWEQ, base = 2) --> doesn't change shape of distribution

item_dat_sm <- item_dat_sm |>  
  mutate(
    AWEQ_s = (AWEQ - mean(AWEQ, na.rm = TRUE)) / sd(AWEQ, na.rm = TRUE) 
    ) |> 
  relocate(AWEQ_s, .after = AWEQ)


# --- standardize trial order

item_dat_sm <- item_dat_sm |>  
  mutate(
    Final_s = (Final_order - mean(Final_order)) / sd(Final_order)
    ) |>  
  relocate(Final_s, .after = Final_order)


item_dat_sm |>  
  summarize(across(.cols = c(AWEQ, Final_order),
                   .fns  = mean)
            )

# AWE-Q: 3.52 out of 5
# Average number of answers: 3.55 responses


item_dat_sm <- item_dat_sm |> 
  mutate(across(.cols = c(AWEQ, AWEQ_s, Final_s),
                .fns  = ~round(., digits = 3)
                )
         )


# --- log & standardize logged frequencies
# --- logging DPcw did not change the shape of the distribution

item_dat_sm <- item_dat_sm |>  
  mutate(across(.cols  = c(freq_total, freq_variant),
                .fns   = ~log(., base = 2),
                .names = "{.col}_log"
                )
         ) |>  
  mutate(across(.cols = c(ends_with("_log"), DPcw),
                .fns  = ~ (. - mean(.)) / sd(.),
                .names = "{.col}_s")
         ) |> 
  mutate(across(.cols = c(contains("_log"), DPcw_s),
                .fns  = ~round(., digits = 3)
                )
         ) |>  
  relocate(c(freq_total_log, 
             freq_total_log_s, 
             freq_variant_log, 
             freq_variant_log_s
             ), 
           .after = freq_variant) |>  
  relocate(DPcw_s, .after = DPcw)


# --- visualize data: Spine plots
# --- spine plots to visualize categorical DV and some target continuous IV

spineplot(factor(match) ~ freq_variant_log_s, data = item_dat_sm,
          ylab = "Matched responses", xlab = "Standardized logged freq")

spineplot(factor(match) ~ freq_total_log_s, data = item_dat_sm,
          ylab = "Matched responses", xlab = "Standardized logged total freq")

spineplot(factor(match) ~ DPcw_s, data = item_dat_sm,
          ylab = "Matched responses", xlab = "Standardized DP Cx as cue")

spineplot(factor(match) ~ AWEQ_s, data = item_dat_sm,
          ylab = "Matched responses", xlab = "Standardized writing exp.")

spineplot(factor(match) ~ Final_s, data = item_dat_sm,
          ylab = "Matched responses", xlab = "Standardized answer order")





# (2) Check: Linear mixed-effects model
#     Did English L1 and Thai L1 participants produced roughly same #items?

# 2.1 create a data frame (no. of items as outcome)

num_dat <- item_dat_sm |>  
  count(L1, Participant, Variant) |>  
  ungroup()

num_dat |>  
  group_by(L1, Variant) |>  
  summarize(
    m = mean(n),
    sd = sd(n)
    )


# 2.2 run a linear model

model1 <- lmer(n ~ L1 * Variant + (1 | Participant), 
               contrasts = list(Variant = contr.sum, L1 = contr.sum),
               data = num_dat
               )

summary(model1)

rm(num_dat, model1)

# No. of items generated did not differ by neither L1, variant, nor their
# interaction (p > 0.05)





# (3) Build models


# ------------------------------ RQ 1 ------------------------------ #

# ----- Intercept-only model ----- #

m0.glmer <- glmer(match ~ 1 + (1 | Participant) + (1 | Response_rnd),
                  data = item_dat_sm,
                  family = binomial(link = "logit")
                  )

summary(m0.glmer)

# In intercept-only model deviance = model misfit (Hox et al 2017)
# Deviance 1044.7; SDs of participants and items = 0.34 and 0.24
# There is variation that can be accounted for by predictors

rm(m0.glmer)



# ----- Model fitting ----- #

# We begin by adding level-1 & level-2 predictors and intra-level interactions
# Create sum contrasts

item_dat_sm <- item_dat_sm |>  
  mutate(across(.cols = c(Degree, L1, Variant, Link_verb),
                .fns  = ~factor(.)
                )
         )

contrasts(item_dat_sm$Variant)   <- contr.sum(2)
contrasts(item_dat_sm$Link_verb) <- contr.sum(2)
contrasts(item_dat_sm$Degree)    <- contr.sum(2)
contrasts(item_dat_sm$L1)        <- contr.sum(2)



# ---
# Run a model: m1s with level-1 predictors

# We added level-1 predictors of interest; we had only one association measure,
# DPcw, since our task set-up was such that construction cued adjectives. Below,
# we added frequency/association measures one at a time

m1.1 <- glmer(match ~ 1 + Final_s + Link_verb + Variant +
                (1 | Participant) + (1 | Response_rnd),
              data = item_dat_sm,
              family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 1e5)
                                     )
              )


car::vif(m1.1) |>  
  as_tibble(rownames = "predictors")

# VIF scores: almost equal 1



m1.2 <- update(m1.1, . ~ . + DPcw_s)
summary(m1.2)

car::vif(m1.2) |> 
  as_tibble(rownames = "predictors")

# DPcw_s independently significant, as well as trial order and variant
# VIF scores: still close to 1



m1.3 <- update(m1.2, . ~ . + freq_variant_log_s)
summary(m1.3)

car::vif(m1.3) |> 
  as_tibble(rownames = "predictors")

# DPcw_s remained significant, freq_variant wasn't
# VIF scores of freq_variant and DPcw were 1.12 and 1.05; still close to 1



m1.4 <- update(m1.3, . ~ . + freq_total_log_s)
summary(m1.4)

car::vif(m1.4) |> 
  as_tibble(rownames = "predictors")

# VIF scores of the two frequency measures went above 100
# Since our focus was primarily on freq_variant, we went back to m1.3



# ---
# Run a model: m2s with level-2 predictors

m2.1 <- update(m1.3, . ~ . + Degree + L1 + AWEQ_s)
summary(m2.1)

car::vif(m2.1) |> 
  as_tibble(rownames = "predictors")

# VIFs: close to 1 for all predictors



# ---
# Run a model: intra-level interactions

m3.1 <- update(m2.1, . ~ . 
               + freq_variant_log_s:DPcw_s)
summary(m3.1)

car::vif(m3.1) |> 
  as_tibble(rownames = "predictors")

# VIF of freq_variant_log_s and DPcw_s still around 1 while that of the
# interaction was 1.33



m3.2 <- update(m3.1, . ~ . 
               + Variant:freq_variant_log_s 
               + Variant:DPcw_s)

summary(m3.2)

car::vif(m3.2) |> 
  as_tibble(rownames = "predictors")

# All VIFs still 1 (range: 1.01-1.37)



# ---
# Run a model: Add theoretically relevant varying slopes for predictors of int.
# Participants: Variant, freq_variant_log_s, and DPcw_s 
# Response_rnd (items): L1 

m4.1 <- update(m3.2, . ~ . 
               - (1 | Participant) 
               - (1 | Response_rnd) 
               + (1 + Variant + freq_variant_log_s + DPcw_s | Participant)
               + (1 + L1 | Response_rnd)
               )

summary(m4.1)

# Impression:
# Singular fit reported

# In the two sources of random-effects variation, Participant and Response_rnd, 
# high (and perfect) correlation observed. Plus, a few components of the 
# random-effect structure had small SDs

# Run a principal components analysis (PCA) on the random-effects structure

summary(rePCA(m4.1))

# The 1st row of each matrix (= SD) indicates that (1) we need only 1
# element for Response_rnd and at most 2 for Participant
# We begin with the elements that had the smallest SD; L1 for Response_rnd



m4.2 <- update(m4.1, . ~ .
               - (1 + L1 | Response_rnd)
               + (1 | Response_rnd)
               )

summary(m4.2)


# LRT supported dropping the random slope for L1
# X^2(2) = 0.566, p = 0.754
# AIC and BIC of m4.2 are also smaller (and deviance was only slightly bigger)

anova(m4.2, m4.1, test = "Chisq")



# We simplify random-effects terms for Participant, beginning with smallest SD
# which is DPcw

m4.3 <- update(m4.2, . ~ . 
               - (1 + Variant + freq_variant_log_s + DPcw_s | Participant)
               + (1 + Variant + freq_variant_log_s | Participant)
               )

summary(m4.3)

# Deletion supported by LRT, X^2(4) = 0.593, p = 0.964

anova(m4.3, m4.2, test = "Chisq")



# In the next step, we dropped random slope for freq_variant, which LRT supports, 
# X^2(3) = 1.211, p = 750. Warning about singular fit also disappeared

m4.4 <- update(m4.3, . ~ . 
               - (1 + Variant + freq_variant_log_s | Participant)
               + (1 + Variant | Participant)
               )

anova(m4.4, m4.3, test = "Chisq")



# ---
# Run a model: Add inter-level interaction terms
# SD for Variant is substantial (0.96), so we include inter-level interaction

m5.1 <- update(m4.4, . ~ . 
               + L1:freq_variant_log_s 
               + L1:DPcw_s
               + L1:Variant
               )

summary(m5.1)

car::vif(m5.1) |> 
  as_tibble(rownames = "parameters")

# VIFs are right around 1
# SD of Variant (random slope) came down (from ~1.04 to ~0.97) but SD of
# random intercept for participants went down to zero



# ---
# Test: Which predictors can be dropped? 

(m5.1_dropped <- drop1(m5.1, test = "Chisq"))


# (1) Two interaction terms (freq_variant_log_s:L1 and DPcw_s:L1) had largest 
# p-value. They are dropped. 

m6.1 <- update(m5.1, . ~ . 
               - freq_variant_log_s:L1
               - DPcw_s:L1
               )


# LRT: m5.1 wasn't better than m6.1, X^2(2) = 0.038, p = 0.981
# Thus, we used m6.1

anova(m6.1, m5.1, test = "Chisq")


# Next, though AWEQ_s has the largest p-value among the predictors left, we will 
# not drop it as it is part of our focus. Instead, we drop: 

m6.2 <- update(m6.1, . ~ . 
               - DPcw_s:freq_variant_log_s 
               - Degree
               - Link_verb
               )

anova(m6.2, m6.1, test = "Chisq")

# LRT supports m6.2, X^2(3) = 3.178, p = 0.365. We then probe which predictors 
# can additionally be dropped.


(m6.2_dropped <- drop1(m6.2, test = "Chisq"))



# ---
# Test: Whether the effect of AWEQ is curved, per reviewer's comment

m6.3 <- update(m6.2, . ~ . 
               - AWEQ_s
               + poly(AWEQ_s, 2)
               )

# We didn't find evidence that AWEQ_s had a curved effect
# m6.3 wasn't better, X^2(1) = 0.382, p = 0.536
# Thus, we go ahead and use m6.2 as our "final" model

anova(m6.3, m6.2, test = "Chisq")


car::vif(m6.2) |> 
  as_tibble(rownames = "parameters")

# All VIF scores are right around 1



m.final <- m6.2

rm(m1.1, m1.2, m1.3, m1.4, m2.1, m3.1, m3.2, 
   m4.1, m4.2, m4.3, m4.4,
   m5.1, m5.1_dropped,
   m6.1, m6.2, m6.2_dropped, m6.3)




# ----- Model interpretation ----- #

# Test: whether final model is better than the null model
# The null model --> intercept + random-effect components of final model

m.null <- glmer(match ~ 1 + (1 | Response_rnd) + (1 + Variant | Participant),
                data = item_dat_sm,
                family = binomial(link = "logit")
                )


anova(m.final, m.null, test = "Chisq")

# m.final significantly different from m.null, X^2(9) = 179.81, p < 0.001

AIC(m.final)
AIC(m.null)

rm(m.null)



# ---
# Obtain: marginal and conditional R^2
# Marginal R^2 captures variance explained by fixed effects only; 
# Conditional R^2 --> variance explained by fixed and random effects

MuMIn::r.squaredGLMM(m.final)

# marginal R^2 = 0.36
# conditional R^2 = 0.52



# ---
# Obtain C score

probs <-  1 / (1 + exp( -fitted(m.final) ) )
Hmisc::somers2(probs, item_dat_sm$match)

rm(probs)





# ----- Additional analysis ----- #

# Removing two participants
# Two L1 Thai-L2 English participants (no. 53 & 59) may have been early 
# learners of English. To ensure that results didn't change with or without 
# these two participants, we dropped them and re-ran the analysis

item_dat_sm2 <- item_dat_sm |>  
  filter(!Participant %in% c(53, 59)
         )

m.two_engl2_out <- glmer(match ~ Final_s + Variant + DPcw_s + 
                           freq_variant_log_s + L1 + AWEQ_s + 
                           Variant:freq_variant_log_s + Variant:DPcw_s + 
                           Variant:L1 +
                           (1 + Variant | Participant) + (1 | Response_rnd),
                         data = item_dat_sm2,
                         family = binomial(link = "logit"),
                         control = glmerControl(optimizer = "bobyqa", 
                                                optCtrl = list(maxfun = 2e5)
                                                )
                         )

summary(m.two_engl2_out)

# dropping the two participants did not change the results
# we therefore included these two participants in our analysis

rm(item_dat_sm2, m.two_engl2_out)





# ------------------------------ RQ 2 ------------------------------ #

# ----- Data preparation ----- #

# --- create a data frame with only Thai L1 participants
# --- remove those that did not provide their TOEFL scores
# --- grand mean-center and standardize TOEFL scores

dat_thai <- item_dat_sm |>  
  filter(L1 == "Thai") |>  
  filter(!Participant %in% c(51, 66, 73, 79)
         ) |>  
  mutate(
    TestScore_s = (TestScore - mean(TestScore)) / sd(TestScore)
    ) |>  
  relocate(TestScore_s, .after = TestScore)


# --- standardize continuous variables

dat_thai <- dat_thai |>  
  mutate(
    AWEQ_s = (AWEQ - mean(AWEQ, na.rm = TRUE)) / sd(AWEQ, na.rm = TRUE),
    Final_s = (Final_order - mean(Final_order)) / sd(Final_order)
    ) |>  
  relocate(AWEQ_s, .after = AWEQ) |>
  relocate(Final_s, .after = Final_order)


dat_thai <- dat_thai |>  
  mutate(across(.cols = c(ends_with("_log"), DPcw),
                .fns  = ~ (. - mean(.)) / sd(.),
                .names = "{.col}_s")
         ) |> 
  mutate(across(.cols = c(contains("_log"), 
                          TestScore_s, AWEQ_s, Final_s, DPcw_s),
                .fns  = ~round(., digits = 3)
                )
         ) |>  
  relocate(freq_total_log_s, .after = freq_total_log) |> 
  relocate(freq_variant_log_s, .after = freq_variant_log) |>  
  relocate(DPcw_s, .after = DPcw)


# --- re-code answers with frequency = 1 as "other"

dat_thai <- dat_thai |>  
  group_by(Response) |>  
  mutate(n = n()) |>  
  mutate(
    Response_rnd = if_else(n == 1, "other", Response)
    ) |>  
  ungroup() |>  
  relocate(Response_rnd, .after = Response) |> 
  select(!n)



# ----- Model fitting ----- #

# Begin with the set of predictors from the m.final model
# Here, L1 and related interactions dropped and TestScore_s added

contrasts(dat_thai$Variant)   <- contr.sum(2)
contrasts(dat_thai$Link_verb) <- contr.sum(2)


m.thai1 = glmer(match ~ Final_s + Variant + DPcw_s + freq_variant_log_s + 
                  AWEQ_s + TestScore_s + 
                  Variant:freq_variant_log_s + Variant:DPcw_s + 
                  (1 + Variant | Participant) + (1 | Response_rnd),
                data = dat_thai,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 1e5) 
                                       )
                )


# Singular fit reported

summary(rePCA(m.thai1))

# Principal component analysis on the random-effects structure revealed that
# we only needed 1 element for Participant, and there's virtually no variation
# in response_rnd

m.thai2 <- update(m.thai1, . ~ .
                  - (1 | Response_rnd)
                  )

anova(m.thai2, m.thai1, test = "Chisq")

# LRT supports dropping the term, X^2(1) = 0, p = 1


m.thai3 <- update(m.thai2, . ~ . 
                  - (1 + Variant | Participant)
                  + (1 | Participant)
                  )

anova(m.thai3, m.thai2, test = "Chisq")

# LRT does not support dropping the term
# all VIF are around 1

car::vif(m.thai2) |>
  as_tibble(rownames = "parameters")



# ---
# Test: Which predictors can be dropped? 

(m.thai2_dropped <- drop1(m.thai2, test = "Chisq"))


# trial order, AWE_Q, Variant:freq_variant_log_s can all be dropped
# We will not drop TOEFL scores since it is our focus

m.thai4 <- update(m.thai2, . ~ .
                  - Final_s
                  - AWEQ_s
                  - Variant:freq_variant_log_s
                  )

anova(m.thai4, m.thai2, test = "Chisq")

# LRT supports dropping, X^2(3) = 3.89, p = 0.274


# We use m.thai4 as our "final" model
# All VIF right around 1

car::vif(m.thai4) |>
  as_tibble(rownames = "parameters")

m.th_final <- m.thai4


rm(m.thai1, m.thai2, m.thai2_dropped, m.thai3, m.thai4)


summary(m.th_final)





# ----- Model interpretation ----- #

# Test: whether final model is better than the null model
# The null model --> intercept + random-effect components of final model

m.th_null <- glmer(match ~ 1 + (1 + Variant | Participant),
                   data = dat_thai,
                   family = binomial(link = "logit")
                   )


anova(m.th_final, m.th_null, test = "Chisq")

# m.th_final significantly different from m.null, X^2(5) = 86.87, p < 0.001

AIC(m.th_final)
AIC(m.th_null)

rm(m.th_null)



# ---
# Obtain: marginal and conditional R^2
# Marginal R^2 captures variance explained by fixed effects only; 
# Conditional R^2 --> variance explained by fixed and random effects

MuMIn::r.squaredGLMM(m.th_final)

# marginal R^2 = 0.42
# conditional R^2 = 0.60



# ---
# Obtain: C score

probs <-  1 / (1 + exp( -fitted(m.th_final) ) )

Hmisc::somers2(probs, dat_thai$match)

rm(probs)





# ----- Additional analysis ----- #

# Four participants with missing TOEFL
# NOTE: 1. Use %$% from magrittr to explode variable
#       2. Use pluck() from purr to "pull" the mean value (= 102) from the 
#          calculation and replace NA with that value

dat_thai2 <- item_dat_sm |> 
  filter(L1 == "Thai") |>  
  mutate(TestScore = if_else(is.na(TestScore) == TRUE, 
                             dat |>  
                               distinct(Participant, .keep_all = TRUE) %$%
                               mean(TestScore, na.rm = TRUE) |>  
                               pluck(), 
                             TestScore
                             ) ,
         TestScore_s = (TestScore - mean(TestScore)) / sd(TestScore)
         ) |>  
  relocate(TestScore_s, .after = TestScore)


m.thai_sm <- glmer(match ~ 1 + Variant + DPcw_s + freq_variant_log_s + 
                     TestScore_s + Variant:DPcw_s +
                     (1 + Variant | Participant),
                   data = dat_thai2,
                   family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "bobyqa", 
                                          optCtrl = list(maxfun = 2e5)
                                          )
                   )

summary(m.thai_sm)
rm(dat_thai2, m.thai_sm)





# (4) Summarize results


# ------------------------------ RQ 1 ------------------------------ #

# --- obtain crosstab

item_dat_sm |>  
  count(L1, Variant, match) |>  
  group_by(L1, Variant) |>  
  mutate(
    Prop = (n/sum(n)) * 100 
    )


# --- count cued responses (match = 1)
# --- all answers in L1 English group had relatively high Delta P

item_dat_sm |> 
  filter(L1 == "English") |> 
  group_by(Variant, match, DPcw, freq_variant) |>  
  count(Response) |> 
  ungroup() |>  
  filter(match != 0) |>  
  select(!match) |>  
  group_by(Variant) |>  
  arrange(desc(n)) |>  
  ungroup()

item_dat_sm |> 
  filter(L1 == "Thai") |> 
  group_by(Variant, match, DPcw, freq_variant) |>  
  count(Response) |> 
  ungroup() |>  
  filter(match != 0) |>  
  select(!match) |>  
  group_by(Variant) |>  
  arrange(desc(n)) |>  
  ungroup()



# --- count non-cued answers

item_dat_sm |>
  filter(L1 == "English") |> 
  group_by(Variant, match, DPcw, freq_variant) |>  
  count(Response) |>  
  ungroup() |>  
  filter(match == 0) |>  
  select(!match) |>  
  group_by(Variant) |> 
  arrange(desc(n)) |> 
  ungroup()

item_dat_sm |>  
  filter(L1 == "Thai") |> 
  group_by(Variant, match, DPcw, freq_variant) |>  
  count(Response) |>  
  ungroup() |>  
  filter(match == 0) |>  
  select(!match) |>  
  group_by(Variant) |> 
  arrange(desc(n)) |> 
  ungroup()



# --- extract fixed-effects estimates: fixef(mod_XXX) and 95% CI: confint(mod_XXX)
# --- and calculate odd ratio with an exponent: exp()

confint(m.final, parm = "beta_", method = "Wald") |>  
  as_tibble(rownames = "Parameters") |>  
  mutate(
    Coeff = fixef(m.final)
    ) |>  
  relocate(Coeff, .after = Parameters) # |>  
  mutate(across(where(is.numeric), exp)
         )


 
# ---
# Obtain: Extract prediction from model
# predict(type = "response") gives us probability

item_dat_sm_est <- item_dat_sm |> 
  bind_cols(predict(m.final, type = "response") |> as_tibble()
            ) |>  
  rename(estimate = value)


# --- obtain summary of the estimate
  
item_dat_sm_est |>  
  group_by(Variant) |>  
  summarize(
    m = mean(estimate),
    sd = sd(estimate)
    ) |>  
  ungroup()

item_dat_sm_est |>  
  group_by(L1) |>  
  summarize(
    m = mean(estimate),
    sd = sd(estimate)
    ) |> 
  ungroup()



# ---
# Simulate prediction from the model: (type = "responses") gives us probability
# Use bootMer to simulate prediction in order to calculate prediction interval

set.seed(4566)

bstrap <- lme4::bootMer(m.final, 
                        function(x) predict(x, type = "response"), 
                        nsim = 100, 
                        re.form = NA,
                        parallel = 'multicore'
                        )


# bstrap is a "bootMer" object; we'll get predicted values from this object
# this can be done with: $t --> bstrap$t will pull these values out
# after converting this into a tibble, we'll get 785 columns (= rows in our data)

# NOTE: We assign this object to df_strp

df_strp <- left_join(item_dat_sm |> 
                       #create a column with row numbers (as characters)
                       #facilitate joining the two tibble
                       mutate(Rows = as.character(row_number())
                               ),
                     bstrap$t |> 
                       as_tibble() |>  
                       #convert this tibble into a long format for joining
                       pivot_longer(cols      = everything(),
                                    names_to  = "Sample",
                                    values_to = "Estimate"), 
                     by = c("Rows" = "Sample")
                     )



# ---
# Understand & visualize marginal effects (mean proportion of cued responses) 

# (1) -- L1s x variants interaction
# Calculate model estimates & boostrapped CI

sum_t <- df_strp |>  
  group_by(L1, Variant, Rows) |>  
  summarize(
    est = mean(Estimate)
    ) |>  
  ungroup() |>  
  group_by(L1, Variant) |>  
  summarize(
    m  = mean(est),
    sd = sd(est),
    se = sd/sqrt(n()),
    low = m - qt(1 - ((1 - 0.95) / 2), n() - 1) * se,
    up  = m + qt(1 - ((1 - 0.95) / 2), n() - 1) * se
    ) |> 
  ungroup()


# Plot the interaction

sum_t |>  
  ggplot(aes(x = Variant, y = m, color = L1)) +
  geom_point(position = position_dodge(width = 0.3),
             size  = 2, 
             shape = 15
             ) +
  geom_errorbar(aes(ymin = low, ymax = up), 
                position = position_dodge(width = 0.3),
                linewidth = 0.75,
                width = 0.3
                ) +
  scale_color_manual(values = c("#000000", "#777777"),
                     name = "Group:",
                     labels = c("L1 English", "L2 English"),
                     guide  = guide_legend(direction = "horizontal") 
                     ) +
  scale_x_discrete(labels = c("Adj-that", "Adj-to") ) +
  labs(
    x = "Introductory-it variants", 
    y = "Proportion of cued responses"
    ) +
  theme_bw() +
  theme(legend.position = c(0.25, 0.93), 
        panel.grid      = element_blank(),
        text = element_text(size = 12),
        axis.text       = element_text(size = 12),
        axis.title      = element_text(size = 12, face = "bold"),
        legend.box.background = element_rect(colour = "black")
        )

ggsave("fig1.png", dpi = 300)

rm(sum_t)



# (2) -- Variants x Delta P & variants x frequency interactions

inter_t1 <- df_strp |>  
  mutate(
    DPcw_cat = round(DPcw, digits = 1)
    ) |>  
  filter(DPcw_cat != 0.2) |> 
  group_by(Variant, DPcw_cat, Rows) |>  
  summarize(
    est = mean(Estimate)
    ) |> 
  ungroup() |>  
  group_by(Variant, DPcw_cat) |>  
  summarize(
    m  = mean(est),
    total = n(),
    sd = sd(est),
    se = sd/sqrt(n()),
    low = m - qt(1 - ((1 - 0.95) / 2), n() - 1) * se,
    up  = m + qt(1 - ((1 - 0.95) / 2), n() - 1) * se
    ) |> 
  ungroup() 
  

fig1 <- inter_t1 |>  
  ggplot(aes(x = DPcw_cat, y = m, color = Variant)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up, group = Variant), alpha = 0.3, color = NA) +
  labs(x = "DeltaP scores of adjectives given their attracting variants", y = "") +
  scale_color_manual(values = c("#eb4034", "#365c45"),
                     name = "Variants:",
                     labels = c("Adj-that", "Adj-to"),
                     guide  = guide_legend(direction = "horizontal") 
  ) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.15), 
        panel.grid      = element_blank(),
        text = element_text(size = 11),
        axis.text       = element_text(size = 11),
        axis.title      = element_text(size = 12, face = "bold"),
        axis.title.y   = element_blank(),
        legend.box.background = element_rect(colour = "black")
        )


inter_t2 <- df_strp |> 
  mutate(
    freq_cat = round(freq_variant_log, digits = 0)
    ) |>  
  group_by(Variant, freq_cat, Rows) |>  
  summarize(est = mean(Estimate)) |>  
  ungroup() |>  
  group_by(Variant, freq_cat) |>  
  summarize(
    m  = mean(est),
    total = n(),
    sd = sd(est),
    se = sd/sqrt(n()),
    low = m - qt(1 - ((1 - 0.95) / 2), n() - 1) * se,
    up  = m + qt(1 - ((1 - 0.95) / 2), n() - 1) * se
    ) |> 
  ungroup()


fig2 <- inter_t2 |>  
  ggplot(aes(x = freq_cat, y = m, color = Variant)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up, group = Variant), alpha = 0.3, color = NA) +
  labs(x = "Logged adjective frequencies in preferred variants", y = "Proportion of cued responses") +
  scale_color_manual(values = c("#eb4034", "#365c45"),
                     name = "Variants:",
                     labels = c("Adj-that", "Adj-to"),
                     guide  = guide_legend(direction = "horizontal") 
                     ) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.15), 
        panel.grid      = element_blank(),
        text = element_text(size = 11),
        axis.text       = element_text(size = 11),
        axis.title      = element_text(size = 12, face = "bold"),
        legend.box.background = element_rect(colour = "black")
  )


figure <- gridExtra::grid.arrange(fig2, fig1, ncol = 2) 
ggsave("fig2.png", figure, dpi = 300, width = 15, height = 8)



# ------------------------------ RQ 2 ------------------------------ #

# Extract fixed-effects estimates: fixef(mod_XXX) and 95% CI: confint(mod_XXX)
# And finally calculate odd ratio with an exponent: exp()

confint(m.th_final, parm = "beta_", method = "Wald") |>  
  as_tibble(rownames = "Parameters") |> 
  mutate(
    Coeff = fixef(m.th_final)
    ) |>  
  relocate(Coeff, .after = Parameters) |>  
  mutate(
    across(where(is.numeric), exp)
    )



# ---
# Obtain: Extract prediction from model
# predict(type = "response") gives us probability

dat_thai_est <- dat_thai |> 
  bind_cols(predict(m.th_final, type = "response") |> as_tibble()
            ) |>  
  rename(estimate = value)
