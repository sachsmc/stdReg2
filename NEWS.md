# stdReg2 1.0.6
* Remove dependence on AF for example dataset


# stdReg2 1.0.5
* Switch to sandwich package for calculation of estimating equations for glms. Was giving incorrect results for noncanonical link functions. Thanks @snhansen

# stdReg2 1.0.4
* minor tweak to documentation 

# stdReg2 1.0.3
* bugfix -- standard error calculation in standardized_coxph for survival probabilities for factors or dummy variables. Thanks @dominicmagirr for the report
* bugfix -- add tidy method for std_custom
* Improved documentation of return objects

# stdReg2 1.0.2
* bugfixes -- warning from model.matrix
* added support for factors in glm methods

# stdReg2 1.0.1

* Documentation improvements
* add progressbar argument to standardize and standardize_level

# stdReg2 1.0.0

* Initial CRAN submission.
