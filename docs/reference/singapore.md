# Case-control study on oesophageal cancer in Chinese Singapore men.

This dataset is borrowed from "Aetiological factors in oesophageal
cancer in Singapore Chinese" by De Jong UW, Breslow N, Hong JG,
Sridharan M, Shanmugaratnam K (1974).

## Usage

``` r
data(singapore)
```

## Format

The dataset contains the following variables:

- Age:

  age of the patient.

- Dial:

  dialect group where 1 represent "`Hokhien/Teochew`" and 0 represent
  "`Cantonese/Other`".

- Samsu:

  a numeric indicator of whether the patient consumes Samsu wine or not.

- Cigs:

  number of cigarettes smoked per day.

- Bev:

  number of beverage at "burning hot" temperatures ranging between 0 to
  3 different drinks per day.

- Everhotbev:

  a numeric indicator of whether the patients ever drinks "burning hot
  beverage" or not. Recoded from the variable "`Bev`".

- Set:

  matched set identification number.

- CC:

  a numeric variable where 1 represent if the patient is a case, 2
  represent if the patient is a control from the same ward as the case
  and 3 represent if the patient is control from orthopedic hospital.

- Oesophagealcancer:

  a numeric indicator variable of whether the patient is a case of
  oesophageal cancer or not.

**The following changes have been made to the data from the original
data in De Jong UW (1974):**

\- The variable "`Bev`" is recoded into the numeric indicator variable
"`Everhotbev`":

`singapore$Everhotbev <- ifelse(singapore$Bev >= 1, 1, 0)`

## References

De Jong UW, Breslow N, Hong JG, Sridharan M, Shanmugaratnam K. (1974).
Aetiological factors in oesophageal cancer in Singapore Chinese. *Int J
Cancer* Mar 15;13(3), 291-303.

[http://faculty.washington.edu/heagerty/Courses/b513/WEB2002/datasets.html](http://faculty.washington.edu/heagerty/Courses/b513/WEB2002/datasets.md)
