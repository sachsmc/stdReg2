# Birthweight data clustered on the mother.

This dataset is borrowed from "An introduction to Stata for health
reserachers" (Juul and Frydenberg, 2010). The dataset contains data on
189 mothers who have given birth to one or several children. In total,
the dataset contains data on 487 births.

## Usage

``` r
data(clslowbwt)
```

## Format

The dataset is structured so that each row corresponds to one
birth/child. It contains the following variables:

- id:

  the identification number of the mother.

- birth:

  the number of the birth, i.e. "1" for the mother's first birth, "2"
  for the mother's second birth etc.

- smoke:

  a categorical variable indicating if the mother is a smoker or not
  with levels "`0. No`" and "`1. Yes`".

- race:

  the race of the mother with levels "`1. White`", "`2. Black`" or
  "`3. Other`".

- age:

  the age of the mother at childbirth.

- lwt:

  weight of the mother at last menstruational period (in pounds).

- bwt:

  birthweight of the newborn.

- low:

  a categorical variable indicating if the newborn is categorized as a
  low birthweight baby (\<2500 grams) or not with levels "`0. No`" and
  "`1. Yes`".

- smoker:

  a numeric indicator if the mother is a smoker or not. Recoded version
  of the variable "`smoke`" where "`0.No`" is recoded as "0" and
  "`1.Yes`" is recoded as "1".

- lbw:

  a numeric indicator of whether the newborn is categorized as a low
  birthweight baby (\<2500 grams) or not. Recoded version of the
  variable "`low`" where "`0.No`" is recoded as "0" and "`1.Yes`" is
  recoded as "1".

**The following changes have been made to the original data in Juul &
Frydenberg (2010):**

\- The variable "`low`" is recoded into the numeric indicator variable
"`lbw`":

`clslowbwt$lbw <- as.numeric(clslowbwt$low == "1. Yes")`

\- The variable "`smoke`" is recoded into the numeric indicator variable
"`smoker`":

`clslowbwt$smoker <- as.numeric(clslowbwt$smoke == "1. Yes")`

## References

Juul, Svend & Frydenberg, Morten (2010). *An introduction to Stata for
health researchers*, Texas, Stata press, 2010 (Third edition).
