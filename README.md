# Regression Standardization for Causal Inference

Goals: create a unified interface for regression standardization to obtain estimates of causal effects such as the average treatment effect, or relative treatment effect. 

1. Should be easy to use for applied practitioners, i.e., as easy as running glm or coxph. 
2. We want to implement modern, theoretically grounded, doubly-robust estimators, and their associated variance estimators. ✓
3. We want it to be extensible for statistical researchers, i.e., possible to implement new estimators and get other models used within the interface. 
4. Robust and clear documentation with lots of examples and explanation of the necessary assumptions. 
5. Test stability of package with simulations. 

## Key features

1. At least the same models supported as stdReg, but with better formatted output, summary functions, and tidy functions with broom. 
2. Options for double robust estimation by specifying a IPTW model. ✓
3. Better support for other causal parameters of interest, e.g., using restricted mean survival and survival probabilities.

## TODO

1. Check the https://cran.r-project.org/web/views/CausalInference.html task view for other packages that do similar things. 

