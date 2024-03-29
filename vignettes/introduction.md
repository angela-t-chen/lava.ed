
# Introduction to the lava.ed R Package

This vignette provides an overview of using the `lava.ed` R package to predict low acuity visits in health data. The package offers four trained models, detailed in Chen, Kuzma, and Friedman (2024), and includes a function called `rule_based_algos`. This function applies the logic of seven previously published ICD-10 code-based algorithms for visit classification. Please refer to the documentation for details on using the `rule_based_algos` function.

The four trained models in the `lava.ed` package are:

1. `lava.l`: Low Acuity Validation Algorithm using logistic regression, trained on all variables from Chen et al. (2024). 
2. `lava.l.sub`: Logistic regression model trained on the "influential subset" of nine variables, including age, any imaging, CBC, CT scan, flu test, IV fluids, pregnancy test, urinalysis, and urine culture.

3. `lava.x`: Low Acuity Validation Algorithm using XGBoost, trained on all variables from Chen et al. (2024). 
4. `lava.x.sub`: XGBoost model trained on the "influential subset" of nine variables, including age, any imaging, CBC, CT scan, flu test, IV fluids, pregnancy test, urinalysis, and urine culture.


# Variables

## All training variables (LAVA-X and LAVA-L)
The variables listed below inform predictions made by the `lava.x` and `lava.l` models. These variable names correspond to those in the NHAMCS dataset, which was used to train the models. When applying these models to other datasets, ensure that variable names are aligned with those listed below.

| NHAMCS variable | Description                     | NHAMCS variable | Description                          |
| --------------- | ------------------------------- | --------------- | ------------------------------------ |
| _age_           | Patient age in years            | _endoint_       | Endotracheal intubation              |
| _sex_           | Sex (Female/Male)               | _flutest_       | Influenza test                       |
| _diag1_         | Diagnosis #1 in ICD-10          | _glucose_       | Glucose, serum                       |
| _anyimage_      | Any imaging                     | _hivtest_       | HIV test                             |
| _bac_           | Blood alcohol concentration     | _incdrain_      | Incision & drainage (I&D)            |
| _bladcath_      | Bladder catheter                | _ivfluids_      | IV fluids                            |
| _bloodcx_       | Blood culture                   | _lactate_       | Lactate                              |
| _bnp_           | Brain natriuretic peptide       | _lft_           | Liver enzymes/Hepatic function panel |
| _bpap_          | BPAP/CPAP                       | _lumbar_        | Lumbar puncture (LP)                 |
| _buncreat_      | Creatinine/renal function panel | _mri_           | MRI                                  |
| _cardenz_       | Cardiac enzymes                 | _nebuther_      | Nebulizer therapy                    |
| _catscan_       | CT scan (any)                   | _othimage_      | Other imaging                        |
| _cbc_           | Complete blood count            | _pregtest_      | Pregnancy/HCG test                   |
| _centline_      | Central line                    | _pttinr_        | Prothrombin time (PT/PTT/INR)        |
| _cpr_           | CPR                             | _skinadh_       | Skin adhesives                       |
| _ctab_          | CT scan – abdomen/pelvis        | _suture_        | Suturing/staples                     |
| _ctchest_       | CT scan – chest                 | _toxscren_      | Toxicology screen                    |
| _cthead_        | CT scan – head                  | _ultrasnd_      | Ultrasound                           |
| _ctother_       | CT scan - other                 | _urine_         | Urinalysis (UA) or urine dipstick    |
| _ddimer_        | D-dimer                         | _urinecx_       | Culture, urine                       |
| _edhiv_         | HIV infection/AIDS              | _woundcx_       | Culture, wound                       |
| _electrol_      | Electrolytes                    | _xray_          | X-ray                                |

## Influential subset of variables (LAVA-L-sub and LAVA-X-sub)
When limited clinical information is available or for a simpler approach, the following variables can be used for prediction with `lava.l.sub` and `lava.x.sub`:

* age
* any imaging
* CBC
* CT scan
* flu test
* IV fluids
* pregnancy test
* urinalysis
* urine culture 


# Input Data

## Data structure
To use the trained models, the R package requires all categorical variables, including the primary diagnosis, to be one-hot-encoded into dummy variables as input. Primary diagnosis variables are included at the level of the first three ICD-10 characters. Dummy columns can be generated using the `dummy_cols` function from the `fastDummies` package.

Each input column should consist of `No/Yes` values indicating whether the categorical variable applies to a given observation. Below is an example of such input format compatible with the `lava.x.sub` and `lava.l.sub` models:

| age|sex_Male |anyimage_Yes |cbc_Yes |catscan_Yes |flutest_Yes |ivfluids_Yes |pregtest_Yes |urine_Yes |urinecx_Yes |
|---:|:--------|:------------|:-------|:-----------|:-----------|:------------|:------------|:---------|:-----------|
|  30|Yes      |Yes          |Yes     |Yes         |No          |Yes          |No           |No        |No          |
|  64|Yes      |No           |Yes     |No          |No          |No           |No           |Yes       |No          |
|  73|No       |No           |No      |No          |No          |No           |No           |No        |No          |
|  69|Yes      |No           |No      |No          |No          |No           |No           |No        |No          |
|  51|Yes      |No           |Yes     |No          |No          |Yes          |No           |No        |No          |
|  56|No       |Yes          |No      |No          |No          |No           |No           |No        |No          |
|  23|Yes      |No           |No      |No          |No          |No           |No           |No        |No          |
|  49|Yes      |Yes          |No      |No          |No          |No           |No           |No        |No          |
|  50|Yes      |No           |No      |No          |No          |No           |No           |No        |No          |
|  53|No       |No           |No      |No          |No          |No           |No           |No        |No          |

## Primary diagnoses
The `lava.x` and `lava.l` models were trained on all ICD-10 codes found in any of the seven rule-based algorithms, totaling 306 ICD codes. Please refer to Appendix Table 2 from Chen et al. (2024) for the list of these codes. To use `lava.x` and `lava.l` models, ensure that all 306 codes are included as dummy variables in the input dataset. 


# Using the models
All four model objects are preloaded with the R package and can be called simply by typing the model name. Once the input data is processed, applying the models for prediction follows the standard procedures. For instance, if the processed dataframe is named `processed_input`:

`pred <- predict(lava.x, newdata = processed_input)`

The default output for the `lava.l` and `lava.l.sub` models will be the logit-transformed probabilities. Turn them back into probabilities by using the inverse logit function:

`prob.predictions <- 1 / (1 + exp(-pred))`.

The default output for the `lava.x` and `lava.x.sub` models will be a vector of predictions, indicating whether each observation is predicted to be low acuity.


## Changing the probability threshold
As described in Chen et al. (2024), while the default probability for visit classification is set at 0.5, researchers may adjust this threshold based on their specific research objectives. The adjustment process is straightforward and involves first predicting the probability of a value being low acuity. Classifications low acuity will then be based on whether the probability exceeds the chosen threshold.

*Example using a 0.3 threshold for classification using `lava.x` or `lava.x.sub`:*

`pred.prob <- predict(lava.l, processed_input, type = "prob")[, "Yes"]`

`pred.prob.class <- as.factor(ifelse(pred.prob >= 0.3, "Yes", "No"))`
