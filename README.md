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
* `AWEQ`: Average scores on the Academic Writing Experience Questionnaire (AWEQ)
* `Frame`: Target stimuli (e.g., Is.that = It is _____ that; Seem.that = It seems _____ that, etc.) Numbers refer to order of responses
* `Variant`: Introductory-*it* variant the frame instantiates (two levels: Adj-that or Adj-to)
* `Response`: Answers supplied by the participants
* `Past_participle`: Whether responses are past participle (e.g., *understood*, *argued*)
* `Incomplete`: Whether responses are typos, mispelled words, or non-adjective responses



## Demographic information
Demographic information of study participants can be found in the **Demographics.csv**, which consists of the following columns:

* `Participant`: Participant ID
* `AWEQ`: Average scores on the Academic Writing Experience Questionnaire (AWEQ)
* `AWE.1`: Rating on the 1st item of the AWEQ (*Write an email to professors, colleagues, or peers on academic-related topics*)
* `AWE.2`: Rating on the 2nd item of the AWEQ (*Provide written feedback on students' papers (as a course TA or instructor)*)
* `AWE.3`: Rating on the 3rd item of the AWEQ (*Submit an abstract for a conference presentation*)
* `AWE.4`: Rating on the 4th item of the AWEQ (*Write a 20+ page paper*)
* `AWE.5`: Rating on the 5th item of the AWEQ (*Write a short summary of a research paper or project*)
* `AWE.6`: Rating on the 6th item of the AWEQ (*Write a critique of an article or study*)
* `AWE.7`: Rating on the 7th item of the AWEQ (*Write a peer review for a project proposal*)
* `AWE.8`: Rating on the 8th item of the AWEQ (*Review manuscripts for a journal*)
* `AWE.9`: Rating on the 9th item of the AWEQ (*Write a part of the thesis or dissertation*)
* `AWE.10`: Rating on the 10th item of the AWEQ (*Prepare an academic presentation*)
* `AWE.11`: Rating on the 11th item of the AWEQ (*Create a PowerPoint for a formal research presentation*)
* `AWE.12`: Rating on the 12th item of the AWEQ (*Submit a grant proposal for funding*)
* `AWE.13`: Rating on the 13th item of the AWEQ (*Submit a research proposal*)
* `AWE.14`: Rating on the 14th item of the AWEQ (*Write a memo or report*)
* `AWE.15`: Rating on the 15th item of the AWEQ (*Take written minutes for a meeting*)
* `AWE.16`: Rating on the 16th item of the AWEQ (*Write a lab report*)
* `AWE.17`: Rating on the 17th item of the AWEQ (*Write an annotated bibliography*)
* `AWE.18`: Rating on the 18th item of the AWEQ (*Complete a report for broken or missing equipment*)

The following columns pertain specifically to the Thai L1 participants (`NA` for the English L1 subjects):

* `AgeLearn`: The age at which participants started learning English
* `AgeFluent`: The age at which participants became fluent in English
* `AgeRead`: The age at which participants started reading in English
* `AgeFluentRead`: The age at which participants became fluent in reading in English
* `PrimInstruct`: A primary means through which participants learned English
  ++ Class = mainly through class instruction
  ++ Interact = mainly through NS interaction
  ++ Both = a mixture of both

