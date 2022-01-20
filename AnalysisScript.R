#############################################################
### Analysis -- Introductory it construction      ###########
### Sakol Suethanapornkul and Sarut Supasiraprapa ###########
#############################################################

library(tidyverse)
library(rsample)
library(ggrepel)

library(magrittr)

library(lme4)
library(lmerTest)


library(effects)
library(ggeffects)










# Corpus analysis: Establishing adjective-variant associations ------------



##### Pre-processing 

# Read file as tibble

corpusIntro <- read_csv("Data/Corpus_introIT.csv")


# 1. Prepare data by extracting adjectives from POS tag column
# (1) str_extract: in each line extract first instance of adjective (POS = JJ) 
# since we only have one _JJ per line, str_extract is sufficient
# matched patterns are vectorized (e.g., stored in a vector)
# regex in "()" captures compound adj eg. far-fetched, well-known, and even
# words like O.K. with period between them by a group here --> [\w|.]
# (2) str_replace: substitute _JJ, _JJR, _JJS [regex = _JJ.*] with nothing

# NOTE: %>% is embedded inside mutate to outline steps taken inside this column

corpusDat <- corpusIntro %>% 
  mutate(Adjectives = 
           str_extract(TaggedQueryItem, "([\\w|.]*)-?([\\w|.]*)_JJ.") %>% 
           str_replace(., "_JJ.*", "") %>% 
           str_to_lower(.)
         )


# OPTIONAL: Extract verbs and adverbs from POS tag column

# (1) since there are 1+ verbs (e.g., has become), we use str_extract_all
# with this, we get a list inside a column "Verbs". we unnest this list with
# unnest_wider() providing column name to the function. The function unnests
# each element to its specific column (e.g., "has" to ...1 and "become" to ...2)
# NOTE: this is a crude way to extract verbs (elements such as modals not included)

# (2) then, we strip POS tags from each columns and unite the two columns
# forming a "Verb" column that contains verbs 

corpusDat <- corpusIntro %>% 
  mutate(Adjectives = 
           str_extract(TaggedQueryItem, "([\\w|.]*)-?([\\w|.]*)_JJ.") %>% 
           str_replace(., "_JJ.*", "") %>% 
           str_to_lower(.),
         Verbs   = str_extract_all(TaggedQueryItem, "([\\w|']*)_V.."),
         Adverbs = str_extract(TaggedQueryItem, "(\\w)*_RB[S|R]?")
         ) %>% 
  unnest_wider(col = Verbs) %>% 
  rename(V1 = ...1,
         V2 = ...2) %>% 
  mutate(across(.cols = V1:V2, ~str_replace(.x, "_V.*", ""))
         ) %>% 
  unite(col = "Verbs", V1:V2, sep = " ", na.rm = TRUE) %>% 
  mutate(Adverbs = str_replace(Adverbs, "_R.*", "")
         )


# Convert superlative adjectives to stem & make spelling consistent
# Use case_when() for multiple if_else patterns and end with TRUE ~ existing vector

corpusDat <- corpusDat %>% 
  mutate(Adjectives = case_when(Adjectives %in% c("best", "better") ~ "good",
                                Adjectives == "cheaper" ~ "cheap",
                                Adjectives == "clearer" ~ "clear",
                                Adjectives == "costlier" ~ "costly",
                                Adjectives %in% c("easier", "easiest") ~ "easy",
                                Adjectives == "fairer" ~ "fair",
                                Adjectives == "farfetched" ~ "far-fetched",
                                Adjectives == "faster" ~ "fast",
                                Adjectives %in% c("harder", "hardest") ~ "hard",
                                Adjectives == "likelier" ~ "likely",
                                Adjectives %in% c("o.k.", "ok") ~ "okay",
                                Adjectives == "odder" ~ "odd",
                                Adjectives == "pleasanter" ~ "pleasant",
                                Adjectives %in% c("safer", "safest") ~ "safe",
                                Adjectives %in% c("simpler", "simplest") ~ "simple",
                                Adjectives == "sweeter" ~ "sweet",
                                Adjectives == "tougher" ~ "tough",
                                Adjectives == "truer" ~ "true",
                                Adjectives %in% c("wiser", "wisest") ~ "wise",
                                Adjectives == "worse" ~ "bad",
                                TRUE ~ Adjectives
                                )
         )





##### Checking frequencies

# Inspect total number of items per variant

corpusDat %>% 
  count(Variant)


# Inspect number of types (total)

corpusDat %>% 
  distinct(Adjectives) %>% 
  count()


# Get type-token ratio of the complete data set

corpusDat %>% 
  group_by(Variant) %>% 
  summarize(type = n_distinct(Adjectives), 
            n = n(), 
            TTR = (type/n)*100
            ) %>% 
  ungroup()


# If words were distributed uniformly, how many instances per word in Adj-that?
# ANS = 60.5

corpusDat %>% 
  filter(Variant == "Adj_that") %>% 
  count(Adjectives) %>%
  summarize(uniform = sum(n)/length(Adjectives)
            )


# Plot the distribution 

TopLable <- corpusDat %>% 
  filter(Variant == "Adj_that") %>%
  count(Adjectives, sort = TRUE) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 5)
  

corpusDat %>% 
  filter(Variant == "Adj_that") %>% 
  count(Adjectives, sort = TRUE) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 100) %>% 
  ggplot(aes(x = ID, y = n)) +
  geom_point(color = "black") +
  geom_line(stat = "identity") +
  ggrepel::geom_text_repel(data = TopLable, aes(label = Adjectives), 
                           show.legend = FALSE, hjust = -0.5) +
  labs(x = "Rank order", y = "Frequencies") +
  geom_hline(yintercept = 60.5, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
        ) +
  ggsave("Adj-that-fulldata.png", dpi = 300)

  
# Repeat the same step with the Adj-to variant

corpusDat %>% 
  filter(Variant == "Adj_to") %>% 
  count(Adjectives, sort = TRUE) %>% 
  summarize(uniform = sum(n)/length(Adjectives)
            )

# ANS = 47.7


# Plot the distribution 

TopLable <- corpusDat %>% 
  filter(Variant == "Adj_to") %>% 
  count(Adjectives, sort = TRUE) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 5)


corpusDat %>% 
  filter(Variant == "Adj_to") %>% 
  count(Adjectives, sort = TRUE) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 100) %>% 
  ggplot(aes(x = ID, y = n)) +
  geom_point(color = "black") +
  geom_line(stat = "identity") +
  ggrepel::geom_text_repel(data = TopLable, aes(label = Adjectives), 
                           show.legend = FALSE, hjust = -0.55) +
  labs(x = "Rank order", y = "Frequencies") +
  geom_hline(yintercept = 47.7, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
        ) +
  ggsave("Adj-to-fulldata.png", dpi = 300)


rm(TopLable)


# Inspect frequency counts in Adj_that and Adj_to

corpusDat %>% 
  group_by(Adjectives) %>% 
  summarize(Adj_that = sum(Variant == "Adj_that"), 
            Adj_to = sum(Variant == "Adj_to")
            ) %>% 
  arrange(-Adj_that) %>% 
  ungroup()


# Create a sample table for DCA with 'necessary' as an example

corpusDat %>% 
  mutate(Necessary = if_else(Adjectives == "necessary", 1, 0) ) %>% 
  group_by(Variant) %>% 
  count(Necessary) %>% 
  ungroup()


# Obtain hapax legomena of each variant

corpusDat %>% 
  filter(Variant == "Adj_that") %>% 
  # filter(Variant == "Adj_to") %>% 
  count(Adjectives) %>% 
  filter(n == 1) %>% 
  summarize(total = n_distinct(Adjectives))


# Relative entropy
# dispersion measure for categorical data (e.g., adjectives in each variant)
# approximates 1 as distributions become more even and 0 for uneven distributions

# Create a function to calculate relative entropy

hrel <- function(x){
  percentage <- x/sum(x)
  find_hrel <- -sum(percentage * log(percentage))/log(length(percentage))
  find_hrel
  }


# Calculate entropy (adj-that & adj-to done separately)

corpusDat %>% 
  filter(Variant == "Adj_that") %>% 
  # filter(Variant == "Adj_to") %>% 
  count(Adjectives, sort = TRUE) %>% 
  summarize(hrel = hrel(n))





##### DCA analysis with script from Gries (2015)

# Write a data set to .txt

write_delim(corpusDat %>% 
              select(Variant, Adjectives) %>% 
              as.data.frame(), 
            file = "./Data/AdjectiveList.txt", 
            delim = "\t")


# Import DCA results back into the current session and create two new variables:
# (1) grand-mean centered coll_strength & (2) standardized coll_strength

dca <- read_delim(file = "./Data/DCA_Corpus.txt", 
                  delim = "\t") %>% 
  mutate(Coll_Str_c = Coll_strength - mean(Coll_strength),
         Coll_Str_s = (Coll_strength - mean(Coll_strength)) / sd(Coll_strength)
         )








# Production results: Analysis -----------------------------------------------


##### Preprocessing

# Read in file; remove white spaces in Response; re-code categories inside Variant

produc_dat <- read_csv("./Data/Production_introIT.csv") %>% 
  mutate(Response = str_to_lower(Response) %>% str_replace(., " ", ""),
         Variant  = if_else(Variant == "Adj-that", "Adj_that", "Adj_to")
         )


# Recode gender, NA to Not specified

produc_dat <- produc_dat %>% 
  mutate(Sex = if_else(is.na(Sex) == TRUE, "NotSpec", Sex)
         )





##### Basic information

# Age and gender

produc_dat %>% 
  distinct(Participant, .keep_all = TRUE) %>% 
  group_by(L1) %>% 
  summarize(m = mean(Age, na.rm = TRUE), 
            sd = sd(Age, na.rm = TRUE)
            ) %>% 
  ungroup()


produc_dat %>% 
  distinct(Participant, .keep_all = TRUE) %>% 
  group_by(L1, Sex) %>% 
  count() %>% 
  ungroup()


# Counts

# 1) % of incomplete trials (misspelled, incomplete, non-adj responses)

produc_dat %>% 
  filter(Incomplete == 1) %>% 
  summarize(n = n(), 
            pt = (n/1057) * 100)


produc_dat %>% 
  filter(Incomplete == 1) %>%
  group_by(L1) %>% 
  summarize(n = n()) %>% 
  ungroup()


# 2) % of past-participle responses

produc_dat %>% 
  filter(Incomplete != 1) %>% 
  group_by(L1, Past_participle) %>% 
  count() %>% 
  ungroup() %>% 
  select(L1, n) %>% 
  group_by(L1) %>% 
  mutate(percent = (n/sum(n)) * 100
         )


# 3) count of remaining total

produc_dat %>% 
  filter(Past_participle == 0 & Incomplete == 0) %>%
  group_by(L1) %>% 
  summarize(total = n() ) %>% 
  ungroup()


# 4) average of responses per variant per L1?

produc_dat %>% 
  filter(Past_participle == 0 & Incomplete == 0) %>% 
  group_by(L1, Participant, Variant) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(L1, Variant) %>% 
  summarize(m = mean(n),
            sd = sd(n)
            ) %>% 
  ungroup()


# 5) number of types (per variant or in total)

produc_dat %>% 
  filter(Past_participle == 0 & Incomplete == 0) %>%
  select(Variant, Response) %>% 
  group_by(Variant) %>% 
  summarize(type = n_distinct(Response)) %>% 
  ungroup()


produc_dat %>% 
  filter(Past_participle == 0 & Incomplete == 0) %>%
  select(Variant, Response) %>% 
  summarize(type_total = n_distinct(Response)) %>% 
  ungroup()





##### Main statistical analysis

# Join tibbles (production data and dca data)

dat <- left_join(produc_dat %>% 
                   filter(Past_participle == 0 & Incomplete == 0) %>% 
                   select(!c(Past_participle, Incomplete)
                          ), 
                 dca %>% 
                   select(Adjectives, Preference, starts_with("Coll")
                          ), 
                 by = c("Response" = "Adjectives")
                 )


# Code whether each response was a distinctive collexeme of a given variant
# Remove 116 items that were not in DCA from analysis; from 901 down to 785 rows

dat %>% 
  filter(is.na(Preference)) %>% 
  summarize(n = n())
  
dat <- dat %>% 
  filter(!is.na(Preference) 
         ) %>% 
  mutate(Match = if_else(Variant == Preference, 1, 0)
         ) %>% 
  relocate(Match, .after = Preference)


# Check no. of types

dat %>% 
  summarize(type = n_distinct(Response) )


dat %>% 
  filter(Match == 1) %>% 
  group_by(L1, Variant) %>% 
  summarize(type = n_distinct(Response), 
            n = n(), 
            TTR = (type/n)*100) %>% 
  ungroup()


# After cases are dropped, how many are left, by L1 and variant

dat %>% 
  group_by(L1, Variant) %>% 
  count() %>% 
  ungroup()


# Mutate new columns
# 1) Add order in which words appear, after responses not in DCA dropped
# 2) Add type of prompt (is vs. seem)

dat <- dat %>% 
  group_by(Participant, Variant) %>% 
  mutate(Final_order = row_number()) %>% 
  ungroup() %>% 
  relocate(Final_order, .after = Frame) %>% 
  mutate(Link_verb  = str_extract(Frame, pattern = ".+?(?=\\.)")
         ) %>% 
  relocate(Link_verb, .after = Variant)


# Check if number of items differed by L1 and variant
# Step 1: Create a dataframe

num_dat <- dat %>% 
  group_by(L1, Participant, Variant) %>% 
  count()


num_dat %>% 
  group_by(L1, Variant) %>% 
  summarize(m = mean(n),
            sd = sd(n)
            )


# Step 2: 

model1 <- lmer(n ~ L1 * Variant + (1 | Participant), 
               contrasts = list(Variant = contr.sum, L1 = contr.sum),
               data = num_dat)

summary(model1)

rm(num_dat, model1)


# Check frequencies of adjectives in production data
# Combine those hapax legomena into one category "other"
# Response_rnd used for random intercept for item

dat <- dat %>% 
  group_by(Response) %>% 
  mutate(n = n()) %>% 
  mutate(Response_rnd = if_else(n == 1, "other", Response)
         ) %>% 
  ungroup() %>% 
  relocate(Response_rnd, .after = Response) %>% 
  select(!n)


# Grand-mean center & scale AWEQ

dat <- dat %>% 
  mutate(AWEQ_c = AWEQ - mean(AWEQ),
         AWEQ_s = (AWEQ - mean(AWEQ)) / sd(AWEQ) 
         ) %>% 
  relocate(c(AWEQ_c, AWEQ_s), .after = AWEQ)
  

# Grand-mean center trial order [Final_order]

dat <- dat %>% 
  mutate(Final_c = Final_order - mean(Final_order),
         Final_s = (Final_order - mean(Final_order)) / sd(Final_order)
         ) %>% 
  relocate(c(Final_c, Final_s), .after = Final_order)




# Intercept-only model -----

# Begin with a plot: correct proportion by participants, not exactly aligned with bernouli trials
# This however shows there is substantial variability across subjects

dat %>% 
  group_by(Participant) %>% 
  summarize(n = n(),
            total = mean(Match),
            se = sd(Match)/sqrt(n)
            ) %>% 
  mutate(se = if_else(is.na(se), 0, se),
         ID = str_c("S",Participant)
         ) %>%
  ungroup() %>% 
  ggplot(aes(x = factor(ID, levels = unique(ID)), y = total)) +
  geom_pointrange(aes(ymin = total - se, ymax = total + se)) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = rel(0.6), angle = 60, vjust = 0.5),
        axis.title = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
        ) #+
  ggsave("Proportion.png", width = 11, height = 11, dpi = 300)


# Run a model

m0 <- glmer(Match ~ 1 + (1 | Participant) + (1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit")
            )


# Calculate Intraclass correlation (ICC) for logistic model

icc <- function(x) {
  cal <- (x^2) / ( (x^2) + ( (pi^2) / 3) )
  cal
  }


VarCorr(m0) %>% 
  as_tibble() %>% 
  mutate(icc = icc(sdcor)
         )


# Since there is no predictor, deviance can be taken to suggest model misfit (Hox et al 2017)
# Deviance 1046.1

  


# Model building ----- 

# We begin by adding level-1 predictors, then level-2 predictors, then intra-level interactions


# Create sum contrasts

dat <- dat %>% 
  mutate(across(.cols = c(Degree, L1, Variant, Link_verb),
                .fns  = ~as.factor(.)
                )
         )

contrasts(dat$Variant)   <- contr.sum(2)
contrasts(dat$Link_verb) <- contr.sum(2)
contrasts(dat$Degree)    <- contr.sum(2)
contrasts(dat$L1)        <- contr.sum(2)


# Run a model: m1 with level-1 predictors

m1 <- glmer(Match ~ 1 + Final_s + Link_verb + Variant + Coll_Str_s + 
                (1 | Participant) + (1 | Response_rnd),
              data = dat,
              family = binomial(link = "logit")
              )

summary(m1)


# Run a model: m2 with level-1 and level-2 predictors 

m2 <- glmer(Match ~ 1 + Final_s + Link_verb + Variant + Coll_Str_s +
              L1 + AWEQ_s + (1 | Participant) + (1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit")
            )

# We did include degree but a model with this predictor wasn't significantly different from m2
# X^2(1) = 0.259, p = 0.611


# Run a model: m3 with intra-level interaction
# An interaction between L1 and AWEQ not significant, nor was each of 2 variables

# a model with an interaction between L1 and AWEQ did not improve fit (when compared to m3)
# X^2(1) = 1.711, p = 0.191

m3 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s +
              L1 + AWEQ_s + (1 | Participant) + (1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit")
            )


# Run a model: Add theoretically relevant varying slopes
# Participants: Variant and Coll_str_s 
# Response_rnd (items): L1

m4 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s + L1 + 
              AWEQ_s + (1 + Variant * Coll_Str_s | Participant) +
              (1 + L1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa", 
                                   optCtrl = list(maxfun = 2e5)
                                   )
            )


m5 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s + L1 + 
              AWEQ_s + (1 + Variant * Coll_Str_s | Participant) +
              (1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa", 
                                   optCtrl = list(maxfun = 2e5)
                                   )
            )


m6 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s + L1 + 
              AWEQ_s + (1 + Variant + Coll_Str_s | Participant) +
              (1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa", 
                                   optCtrl = list(maxfun = 2e5)
                                   )
            )


m7 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s + L1 + 
              AWEQ_s + (1 + Variant | Participant) +
              (1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa", 
                                   optCtrl = list(maxfun = 2e5)
                                   )
            )


m8 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s + L1 + 
              AWEQ_s + (1 + Coll_Str_s | Participant) +
              (1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa", 
                                   optCtrl = list(maxfun = 2e5)
                                   )
            )


# Adding random slopes resulted in overfitting (a singular fit was reported)
# We thus adopted m3

rm(m0, m1, m2, m4, m5, m6, m7, m8)


# Run a model: Add cross-level interaction 
# We had only one cross-level interaction term (L1:Variant).
# Additional interaction terms did not improve model fit.

m9 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s +
              L1 + L1:Variant + AWEQ_s + (1 | Participant) +
              (1 | Response_rnd),
            data = dat,
            family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa", 
                                   optCtrl = list(maxfun = 2e5)
                                   )
            )


m10 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s +
               L1 + L1:Variant + L1:Coll_Str_s + AWEQ_s + 
               (1 | Participant) + (1 | Response_rnd),
             data = dat,
             family = binomial(link = "logit"),
             control = glmerControl(optimizer = "bobyqa", 
                                    optCtrl = list(maxfun = 2e5)
                                    )
             )

# model with the 3-way interaction (L1:Variant:Coll_Str_s) not shown

anova(m3, m9)
anova(m9, m10)      #likelihood ratio test: not significant


rm(m3, m9, m10)


# Run a model: Rerun m9 and draw inference

m_final <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s +
                   L1 + L1:Variant + AWEQ_s + (1 | Participant) +
                   (1 | Response_rnd),
                 data = dat,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa", 
                                        optCtrl = list(maxfun = 2e5)
                                        )
                 )


# Additional Analysis: Removing two participants
# Two Thai L1 participants (no. 53 & 59) may have been early learners of English
# To ensure that results didn't change with or without these two participants
# we dropped the two participants and re-ran the analysis

dat2 <- dat %>% 
  filter(!Participant %in% c(53, 59))

m_subset <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s +
                    L1 + L1:Variant + AWEQ_s + (1 | Participant) +
                    (1 | Response_rnd),
                  data = dat2,
                  family = binomial(link = "logit"),
                  control = glmerControl(optimizer = "bobyqa", 
                                         optCtrl = list(maxfun = 2e5)
                                         )
                  )

# dropping the two participants did not change the results substantially
# L1 became significant though (from being marginally significant)
# The estimate of L1 however rarely changed. We decided to retain participants

rm(dat2, m_subset)




# Model building: L2 proficiency -----

# 1) Create a dataframe with only Thai L1 participants
#    Add average TOEFL score to 4 participants (who didn't provide info)
#    Grand mean-center and standardize TOEFL scores

# NOTE: use %$% from magrittr to explode variable
# NOTE: use pluck() from purr to pluck a single element from a vector
# Since there's only one value in mean(), we can simply call pluck()

dat_thai <- dat %>% 
  filter(L1 == "Thai") %>% 
  mutate(TestScore = if_else(is.na(TestScore) == TRUE, 
                             dat %>% 
                               distinct(Participant, .keep_all = TRUE) %$%
                               mean(TestScore, na.rm = TRUE) %>% pluck(), 
                             TestScore
                             ) 
         )


dat_thai <- dat_thai %>% 
  mutate(TestScore_s = (TestScore - mean(TestScore)) / sd(TestScore)
         ) %>% 
  relocate(TestScore_s, .after = TestScore)


# 2) Begin with the set of predictors from the m_final model
#    Here, L1 is dropped (since there's no longer the other level to compare)
#    Add proficiency to model and cross-level interaction terms

m_th1 <- glmer(Match ~ 1 + Final_s + Link_verb * Variant * Coll_Str_s + AWEQ_s +
                 TestScore_s + TestScore_s:Coll_Str_s + TestScore_s:Coll_Str_s:Variant +
                 (1 | Participant) + (1 | Response_rnd),
               data = dat_thai,
               family = binomial(link = "logit"),
               control = glmerControl(optimizer = "bobyqa", 
                                      optCtrl = list(maxfun = 2e5)
                                      )
               )


# None of the interaction terms were significant (including variant*strength)
# Thus, we simplified the model. Plus, we dropped control variables
# zooming in on main effects we're interested in

m_th2 <- glmer(Match ~ 1 + Variant + Coll_Str_s + AWEQ_s + TestScore_s + 
                 (1 | Participant) + (1 | Response_rnd),
               data = dat_thai,
               family = binomial(link = "logit"),
               control = glmerControl(optimizer = "bobyqa", 
                                      optCtrl = list(maxfun = 2e5)
                                       )
               )


# Additional Analysis: Drop four participants without TOEFL scores

dat_thai2 <- dat %>% 
  filter(L1 == "Thai") %>% 
  filter(!Participant %in% c(51, 66, 73, 79)) %>% 
  mutate(TestScore_s = (TestScore - mean(TestScore, na.rm = T)) / sd(TestScore, na.rm = T)
         ) %>% 
  relocate(TestScore_s, .after = TestScore)


m_th2_s <- glmer(Match ~ 1 + Variant + Coll_Str_s + AWEQ_s + TestScore_s + 
                   (1 | Participant) + (1 | Response_rnd),
                 data = dat_thai2,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa", 
                                        optCtrl = list(maxfun = 2e5)
                                        )
                 )


# No change in terms of significant predictors

rm(dat_thai2, m_th2_s)





##### Reporting

# crosstab

dat %>% 
  count(L1, Variant, Match) %>% 
  group_by(L1, Variant) %>% 
  mutate(Prop = (n/sum(n)) * 100 
         )


# counts of matched responses (match = 1)
# most common answers had high collostructional strength

dat %>%
  # filter(L1 == "Thai") %>%
  group_by(Variant, Match, Coll_strength) %>% 
  count(Response) %>% 
  ungroup() %>% 
  filter(Match != 0) %>% 
  select(!Match) %>% 
  group_by(Variant) %>% 
  arrange(desc(n)) %>% 
  ungroup() # %>% 
  head(20)


# counts of matched responses by AWE-Q scores
# we split participants into groups based on median AWE-Q scores
# (we may also group by L1 but the patterns don't change)

dat %>% 
  #group_by(L1) %>% 
  mutate(AWEQ_group = if_else(AWEQ > median(AWEQ), "High", "Low"
                              )
         ) %>% 
  group_by(AWEQ_group) %>% 
  count(Match) %>% 
  ungroup()


chisq.test(matrix(c(157, 235, 148, 245), nrow = 2, byrow = T))


# Extract fixed-effects estimates --- with fixef(mod_XXX) --- and 95% CI ---with confint(mod_XXX)
# And finally calculate odd ratio with an exponent -- exp --

confint(m_final, parm = "beta_", method = "Wald") %>% 
  as_tibble(rownames = "Parameters") %>% 
  mutate(Coeff = fixef(m_final)
         ) %>% 
  relocate(Coeff, .after = Parameters) #%>% 
  mutate(across(where(is.numeric), exp)
         )


# Extract prediction from model
# predict(type = "response") gives us probability

dat_est <- dat %>% 
  bind_cols(predict(m_final, type = "response") %>% 
              as_tibble()
            ) %>% 
  rename(Estimate = value)


# Obtain summary of the estimate

dat_est %>% 
  group_by(Variant) %>% 
  summarize(m = mean(Estimate),
            sd = sd(Estimate)
            ) %>% 
  ungroup()


dat_est %>% 
  group_by(L1) %>% 
  summarize(m = mean(Estimate),
            sd = sd(Estimate)
            ) %>% 
  ungroup()


# Simulate prediction from the model: (type = "responses") gives us probability
# Use bootMer to simulate prediction in order to calculate prediction interval

bstrap <- lme4::bootMer(m_final, 
                        function(x) predict(x, type = "response"), 
                        nsim = 100, 
                        re.form = NA,
                        parallel = 'multicore'
                        )


# bstrap is a "bootMer" object; we'll need to get predicted values from this object
# this can be done with: $t --> bstrap$t will pull these values out
# once we convert this into a tibble, we'll get 785 columns (which are rows in our data)

# NOTE: We assign this object to dat_strp

dat_strp <- left_join(dat %>%
                        #create a column with row numbers (as characters)
                        #facilitate joining the two tibble
                        mutate(Rows = as.character(row_number())
                               ),
                      bstrap$t %>% 
                        as_tibble() %>% 
                        #convert this tibble into a long format for joining
                        pivot_longer(cols      = everything(),
                                     names_to  = "Sample",
                                     values_to = "Estimate"), 
                      by = c("Rows" = "Sample")
                      )


# Visualize marginal effects (i.e., mean proportion of collexeme responses) 

# L1s x variants interaction
# Calculate model estimates & boostrapped CI

sum_t <- dat_strp %>% 
  group_by(L1, Variant, Rows) %>% 
  summarize(est = mean(Estimate)
            ) %>% 
  ungroup() %>% 
  group_by(L1, Variant) %>% 
  summarize(m  = mean(est),
            sd = sd(est),
            se = sd/sqrt(n()),
            low = m - qt(1 - ((1 - 0.95) / 2), n() - 1) * se,
            up  = m + qt(1 - ((1 - 0.95) / 2), n() - 1) * se
            ) %>% 
  ungroup()


sum_t %>% 
  ggplot(aes(x = Variant, y = m, color = L1)) +
  geom_point(position = position_dodge(width = 0.3),
             size  = 3, 
             shape = 15
             ) +
  geom_errorbar(aes(ymin = low, ymax = up), 
                position = position_dodge(width = 0.3),
                size = 1,
                width = 0.3
                ) +
  scale_color_manual(values = c("#000000", "#777777"),
                     name = "Group:",
                     labels = c("English L1", "Thai L1-English L2"),
                     guide  = guide_legend(direction = "horizontal") 
                     ) +
  scale_x_discrete(labels = c("Adj-that", "Adj-to") ) +
  labs(x = "Introductory-it variants", y = "Proportions of collexeme response") +
  theme_bw() +
  theme(legend.position = c(0.34, 0.95), 
        panel.grid      = element_blank(),
        text = element_text(size = 12),
        axis.text       = element_text(size = 12),
        axis.title      = element_text(size = 12, face = "bold"),
        legend.box.background = element_rect(colour = "black")
        )

ggsave("L1_var.png", dpi = 300)

rm(sum_t)


# Visualize marginal effects (i.e., mean proportion of collexeme responses) 

# variants * strength interaction
# To begin, we create a new tibble that bins Coll_strength by quantile

quantile(dat$Coll_strength, probs = c(0.25, 0.50, 0.75))

dat_strp <- dat_strp %>% 
  mutate(Coll_Str_b = case_when(Coll_strength <  1.771 ~ 1,
                                Coll_strength >= -1.771  & Coll_strength < 14.271 ~ 2,
                                Coll_strength >=  14.271  & Coll_strength < 201.751 ~ 3,
                                Coll_strength >=  201.751 ~ 4)
         )


# Calculate summary statistics by variant and strength

sum_t <- dat_strp %>% 
  group_by(Variant, Coll_Str_b, Rows) %>% 
  summarize(est = mean(Estimate)
            ) %>% 
  ungroup() %>% 
  group_by(Variant, Coll_Str_b) %>% 
  summarize(m  = mean(est),
            sd = sd(est),
            se = sd/sqrt(n()),
            low = m - qt(1 - ((1 - 0.95) / 2), n() - 1) * se,
            up  = m + qt(1 - ((1 - 0.95) / 2), n() - 1) * se
            ) %>% 
  ungroup()


sum_t %>% 
  ggplot(aes(x = Coll_Str_b, y = m) ) +
  geom_line(aes(group = Variant, color = Variant)) +
  geom_ribbon(aes(ymin = low, ymax = up, group = Variant), alpha = 0.3) +
  annotate(geom = "text", x = 1.2, y = 0.45, label = "Adj-that") +
  annotate(geom = "text", x = 1.2, y = 0.7, label = "Adj-to") +
  labs(x = "Collostructional strengths", y = "Mean proportion of collexeme responses") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),
                     labels = c("<2", "2-14", "14-202", ">202") 
                     ) +
  scale_color_manual(values = c("#000000", "#000000")
                     ) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 12),
        axis.text   = element_text(size = 12),
        axis.text.x = element_text(size = 12, hjust = 0.8),
        axis.title  = element_text(size = 12, face = "bold"),
        panel.grid  = element_blank() 
        )

ggsave("Var_strengh.png", dpi = 300)


# Visualize marginal effects (i.e., mean proportion of collexeme responses) 

# strength x variants x link verb
# Calculate model estimates & boostrapped CI

sum_t <- dat_strp %>% 
  group_by(Variant, Link_verb, Coll_Str_b, Rows) %>% 
  summarize(est = mean(Estimate)
            ) %>% 
  ungroup() %>% 
  group_by(Variant, Link_verb, Coll_Str_b) %>% 
  summarize(m  = mean(est),
            sd = sd(est),
            se = sd/sqrt(n()),
            low = m - qt(1 - ((1 - 0.95) / 2), n() - 1) * se,
            up  = m + qt(1 - ((1 - 0.95) / 2), n() - 1) * se
            ) %>% 
  ungroup()
  
  
# Plot 

quantile(dat$Coll_strength, probs = c(0.25, 0.50, 0.75))
axis <- c("<2", "2-14", "14-202", ">202")


sum_t %>% 
  mutate(Variant = if_else(Variant == "Adj_that", "Adj-that", "Adj-to") 
         ) %>%  
  ggplot(aes(x = Coll_Str_b, y = m) ) +
  geom_line(aes(group = 1)) +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.3) +
  labs(x = "Collostructional strengths", y = "Mean proportion of collexeme responses") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),
                     labels = c("<2", "2-14", "14-202", ">202")) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text   = element_text(size = 12),
        axis.text.x = element_text(size = 12, hjust = 0.8),
        axis.title  = element_text(size = 12, face = "bold"),
        panel.grid  = element_blank() 
        ) +
  facet_grid(Link_verb ~ Variant)
  

ggsave("Verb_variant_strength.png", dpi = 300)


rm(sum_t)

#---------------------------------------------
# We could have used ggeffects:ggemmeans() to visualize the three-way interaction

# To do so, we can specify the terms in terms = c() using variable names from the model
# For continuous variables e.g., Coll_strength, we can get marginal effects at specified values

predict_df <- ggeffects::ggemmeans(m_final, 
                                   terms = c("Link_verb", "Variant", "Coll_Str_s [-1.5, 0, 1.5, 3, 4.5]")
                                   )

# We obtain a tibble, with predicted values calculated for 1st term provided (x [Link_verb])
# Two other columns are grouping variables (group [Variants] and facet [Strengths])
# Rename for clarity

predict_df <-  predict_df %>% 
  rename(Verbs     = x,
         Variants  = group,
         Strengths = facet)


# Plot (but first we create a label for x axis -- convert coll_str back to original)

xaxis <- (seq(from = -1.5, to = 4.5, by = 1.5) * sd(dca$Coll_strength)) + mean(dca$Coll_strength)
xaxis <- round(xaxis, digits = 0)  
  

ggplot(data = predict_df, 
       aes(x = Strengths, y = predicted, group = Variants) ) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3 ) +
  # scale_x_discrete(breaks = c(0, 2, 4, 6, 8),
  #                  labels = xaxis ) +
  labs(x = "Collostructional strengths", y = "Mean proportion of collexeme responses") +
  theme_bw() +  
  theme(text = element_text(size = 12),
        axis.text  = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid = element_blank() 
        ) +
  facet_grid(Verbs ~ Variants) 


rm(predict_df)

#---------------------------------------------


# The item: It seems [blank] to --> As strengths increased, mean proportion went down.
# To get a sense of what happened, we filter out responses for Adj-to from dat:

probe <- left_join(dat_est %>% 
                     filter(Variant == "Adj_to") %>% 
                     count(Response, Link_verb, Match),
                   dca %>% 
                     select(Adjectives, starts_with("Coll")
                            ), 
                   by = c("Response" = "Adjectives")
                   )


dat_est %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(Link_verb) %>% 
  count(Match) %>% 
  ungroup()
  
# Filter out highly distinctive collexemes

probe %>% 
  filter(Coll_strength > 10, Match == 0)


# NOTE: A few observations
# (1) high-strength adjectives supplied fewer times (e.g., difficult, hard, important, necessary)
# (2) responses that are collexemes of *Adj-that* were supplied more often
# more cases with match = 0 in seems (30) than is (21)
# (2) for four adjectives--likely, obvious, unlikely, clear-- (which prefer Adj-that)
# they are supplied in Adj-to items a greater number of times with seem than is




# L2 proficiency (with English L2 data only)

# Extract fixed-effects estimates --- with fixef(mod_XXX) --- and 95% CI ---with confint(mod_XXX)
# And finally calculate odd ratio with an exponent -- exp --

confint(m_th2, parm = "beta_", method = "Wald") %>% 
  as_tibble(rownames = "Parameters") %>% 
  mutate(Coeff = fixef(m_th2)
         ) %>% 
  relocate(Coeff, .after = Parameters) 


# L2 proficiency was not significant
# We checked if higher vs. lower proficiency answered "more accurately"

dat_thai %>% 
  mutate(Test_bin = if_else(TestScore > median(TestScore), "High", "Low" )
         ) %>% 
  group_by(Test_bin) %>% 
  count(Match) 


chisq.test(matrix(c(55, 85, 82, 105), nrow = 2, byrow = TRUE))








##### Additional information: Collexeme analysis in production data #####

# Obtain a full list of responses from production data
# NOTE: we use the dataframe "product_dat" rather than the dataframe "dat" 
#       dat consists of adjectives that were attested in the dca
#       we'd like to conduct DCA on all adjectives supplied by participants

fulldat <- produc_dat %>% 
  filter(Past_participle == 0 & Incomplete == 0) %>% 
  select(Variant, Response) %>% 
  arrange(Variant)

write_delim(fulldat, file = "ProductionFullList.txt", delim = "\t")



# Obtain a list of responses from English L2 subjects

thaidat <- produc_dat %>% 
  filter(Past_participle == 0 & Incomplete == 0) %>% 
  filter(L1 == "Thai") %>% 
  select(Variant, Response) %>% 
  arrange(Variant)

write_delim(thaidat, file = "ProductionThaiList.txt", delim = "\t")


rm(fulldat, thaidat)








##### Additional information: Frequency or strength #####

# Read DCA results into the session

dca_thai <- read_delim(file = "./Data/DCA_ProductionThai.txt", delim = "\t") %>% 
  select(!c(starts_with("Freq_exp"), starts_with("DeltaP")) 
         )


# Join two tibbles and filter out rows that do not exist in DCA of corpus data

dca_thai <- left_join(dca_thai, 
                     dca %>% 
                       select(!c(starts_with("Freq_exp"), starts_with("DeltaP"),
                                 Coll_Str_c, Coll_Str_s
                                 )
                              ),
                     by = "Adjectives",
                     suffix = c(".pro", ".cor")
                     ) %>% 
  filter(!is.na(Freq_raw_that.cor))


# Correlation test 

# Begin by +1 to raw frequency to prevent log(0) = inf
# To every single word or to those words that preferences match

dca_thai_cor <- dca_thai %>% 
  # mutate(Match = if_else(Preference.pro == Preference.cor, 1, 0)
  #        ) %>%
  # filter(Match == 1) %>%
  select(Adjectives, Preference.pro, starts_with("Freq_raw"), starts_with("Coll_strength")) %>% 
  mutate(across(.cols = starts_with("Freq_raw"),
                .fns  = ~ . + 1)
         ) %>% 
  mutate(across(.cols = where(is.numeric),
                .fns  = ~log(., base = 2)
                )
         )


# Log frequencies

dca_thai_cor %>% 
  filter(Preference.pro == "Adj_that") %$% 
  cor.test(Freq_raw_that.pro, Freq_raw_that.cor)

# Full set   --> Cor: r(65) = 0.594, p < 0.001
# Match only --> Cor: r(24) = 0.578, p < 0.01


dca_thai_cor %>% 
  filter(Preference.pro == "Adj_to") %$% 
  cor.test(Freq_raw_to.pro, Freq_raw_to.cor)

# Full set   --> Cor: r(47) = 0.748, p < 0.001
# Match only --> Cor: r(38)  = 0.761, p < 0.001 



# Collostructional strength

dca_thai_cor %>% 
  filter(Preference.pro == "Adj_that") %$%
  cor.test(Coll_strength.pro, Coll_strength.cor)

# Full set   --> Cor: r(65) = 0.05, p = 0.663
# Match only --> Cor: r(24) = 0.23, p = 0.241


dca_thai_cor %>% 
  filter(Preference.pro == "Adj_to") %$%
  cor.test(Coll_strength.pro, Coll_strength.cor)

# Full set   --> Cor: r(47) = 0.42, p = 0.002
# Match only --> Cor: r(38) = 0.458, p = 0.002


rm(dca_thai, dca_thai_cor)
