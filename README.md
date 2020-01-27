# measerror-msm

Version 0.1.0

Code for generating table 1 from the preprint 'Quantitative bias analysis for a misclassified confounding variable in point-treatment studies: a comparison between marginal structural models and conditional models'.

The R code in /src can be used to generate the table that comprises the results of the simulation study. The table that is produced by the script can be used in LaTeX files. 

The script will use the .rds files available in ./data/processed.

Download the repository by using:
```console
git clone https://github.com/LindaNab/measerror-msm
```

This script is tested on:
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14.6

With attached packages:
xtable_1.8-4
rsimsum_0.6.1

## Project organization

```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── requirements.txt
├── bin                <- Compiled and external code, ignored by git (PG)
│   └── external       <- Any external source code, ignored by git (RO)
├── config             <- Configuration files (HW)
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                <- Source code for this project (HW)

```


## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)

## Citation

Please [cite this project as described here](/CITATION.md).
