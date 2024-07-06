This repository contains the code for conducting the simulation study and data analysis for the paper: `Statistical Inference on Relative Risk under Covariate-Adaptive Randomization`, by Fengyu Zhao and Feifang Hu.

**Simulation Study**

As described in Section 5, we have three simulation scenarios. There are mainly four files:

-`Scenario_1.R` contains the simulation of Case 1 and Case 2 for the first scenario, which assumes the underlying model is logistic regression

-`Scenario_2.R` contains the simulation of for the second scenario, which assumes the underlying model does not follow logistic regression

-`Re_design.R` contains the code for the re-design of COVID-19 vaccine trial

-`CAR_cpp.cpp` is the file to perform three Covariate adaptive randomiation - Pocock and Simon's marginal procedure (PS), the stratified permuted block design (STRPB), and Hu \& Hu's  procedures (HH)
