# Usage Events and L2 Constructional Knowledge: Two Variants of Introductory-it construction

Sakol Suethanapornkul (Thammasat University) [[Personal website](https://sakol.netlify.app)]

This repository contains data and an R script for the paper "Usage Events and L2 Constructional Knowledge: A Study of Two Variants of the Introductory-*it* Construction" which is under review.


[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png


## Overview

The following files are inside the **Data** folder:

* **AdjectiveList.txt**: a long-format table of adjectives and the variants with which they are attested. Minus the column-header row, the file contains 39,489 lines. This is divided into 15,862 lines of Adj-that and 23,627 of Adj-to. 
* **Corpus_introIT.csv**: corpus data from COCA. Note that this file contains only instances of Adj-that and Adj-to variants. Instances that were not part of the construction were removed, per steps discussed in the article.
* **DCA_Corpus.txt**: A distinctive collexeme analysis report of the corpus data.
* **DCA_ProductionFull.txt**: A distinctive collexeme analysis report of the elicitation data, English L1 and Thai L1 data combined.
* **DCA_ProductionThai.txt**: A distinctive collexeme analysis report of the elicitation data, Thai L1 data only.
* **Demographics.csv**: Participants' demographic information.
* **Production_introIT.csv**: Elicitation data from English L1 and Thai L1 participants.



## Elicitation data

Elicitation data can be found in the **Production_introIT.csv** file, which consists of the following columns:

* `Participant`: Participant ID
* `Time_min`: Time spent completing the experiment (in minutes)
* `Age`: Participants' age (in years)
* `Sex`: Participants' gender (three levels: Male, Female, or Non-specified)
* `Degree`: Current graduate degree (two levels: MS or Doctoral)
* `Program`: Field of study
* `L1`: Participants' first language (two levels: English or Thai)
* `TestScore`: TOEFL iBT scores
* `AWEQ`: Average scores on the Academic Writing Experience Questionnaire
* `Frame`: Target stimuli (e.g., Is.that = It is _____ that; Seem.that = It seems _____ that, etc.) Numbers refer to order of responses
* `Variant`: Introductory-*it* variant the frame instantiates (two levels: Adj-that or Adj-to)
* `Response`: Answers supplied by the participants
* `Past_participle`: Whether responses are past participle (e.g., *understood*, *argued*)
* `Incomplete`: Whether responses are typos, mispelled words, or non-adjective responses



## Demographic dnformation
Demographic information of study participants can be found in the **Demographics.csv**, which consists of the following columns:


