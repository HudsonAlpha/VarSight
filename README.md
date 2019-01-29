# VarSight

VarSight is a collection of scripts used to test the application of classification algorithms for identifying variants that will be returned on a clinical report.  VarSight relies on [Codicem](http://envisiongenomics.com/codicem-analysis-platform/) to perform common clinical pipeline pre-processing, namely variant annotation and filtering.  In addition to data available in Codicem, it also uses gene rankings that are based on Human Phenotype Ontology (HPO) terms using both a built in calculation (cosine score) and the external tool [PyxisMap](https://github.com/HudsonAlpha/LayeredGraph).  The classifiers in VarSight ingest a total of 50 features that are used to test both classification of reported variants and prioritization of reported variants in a larger filtered list.  Results from v1 of VarSight tests are available in the folder labeled "paper".

## Sub-folders
1. CODI_metadata - Contains metadata files that describe how Codicem data was filtered and how annotations were used by the classifiers.
2. paper - Contains the main LaTeX file for the paper, multiple Jinja2 template files to use for rendering our results tables automatically, and the results JSON file that contains all output from the classifier results.
3. scripts - Contains all Python3 scripts used to test the classifiers and report results.

## Publications

[Pre-print on bioRxiv: Holt, James M. et al. "VarSight: Prioritizing Clinically Reported Variants with Binary Classification Algorithms."](https://www.biorxiv.org/content/10.1101/532440v1)