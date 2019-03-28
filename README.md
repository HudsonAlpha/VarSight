# VarSight

VarSight is a collection of scripts used to test the application of classification algorithms for identifying variants that will be returned on a clinical report.  VarSight relies on [Codicem](http://envisiongenomics.com/codicem-analysis-platform/) to perform common clinical pipeline pre-processing, namely variant annotation and filtering.  In addition to data available in Codicem, it also uses gene rankings that are based on Human Phenotype Ontology (HPO) terms using both a built in calculation (cosine score) and the external tool [PyxisMap](https://github.com/HudsonAlpha/LayeredGraph).  The classifiers in VarSight ingest a total of 95 features that are run through a feature selection process to get the top 20 features.  Those 20 features are used to test both classification of reported variants and prioritization of reported variants in a larger filtered list.  Results from the latest version of VarSight tests are available in the folder labeled "paper".

## Running VarSight (Beta)
To run the trained models, you need a Codicem JSON file containing pre-annotated and pre-filtered variants.  The filter file used can be found in the `CODI_metadata` subdirectory.  Additionally, a file with one HPO term per line describing the patient is also required.  Finally, the trained classifiers (pickled in the `models` subdirectory) are required.  Given that information, the following command will print the VarSight predicted ranked ordering of variants:
```
python3 scripts/VarSight.py codicemJson hpoFN modelsFN
```

For more information on the VarSight command-line interface, run the help command:
```
python3 scripts/VarSight.py -h
```

## Sub-folders
1. CODI_metadata - Contains metadata files that describe how Codicem data was filtered and how annotations were used by the classifiers.
2. ExomiserTemplates - Contains the templated parameter files used to run Exomiser as a comparison.
3. models - Contains the pickled models after hyperparameter tuning and fitting to the training set.  The models can be used by the VarSight.py script to rank variants from a Codicem filtered variant file.
4. paper - Contains the main LaTeX file for the paper, multiple Jinja2 template files to use for rendering our results tables automatically, and the results JSON file that contains all output from the classifier results.
5. scripts - Contains all Python3 scripts used to test the classifiers and report results.

## Publications

Pre-print on bioRxiv: [Holt, James M. et al. "VarSight: Prioritizing Clinically Reported Variants with Binary Classification Algorithms."](https://www.biorxiv.org/content/10.1101/532440v2)