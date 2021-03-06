# VarSight Scripts
This folder contains a wide variety of scripts that were used to gather data, train the classifiers, test the classifiers, and report the results.

## VarSight main
`VarSight.py` contains the main command-line interface for running VarSight after models have been trained.  Please refer to the main README for more instructions.

## Data Parsing Scripts
The following files are primarily used to pre-process data into easy-to-manage formats in Python.

1. CodiDumpUtil.py - This script contains helper functions for parsing a filtered Codicem JSON. It's main purpose is to load a Codicem file and full out the fields requested for use by the classifiers.
2. ExomiserUtil.py - This script contains helper functions for parsing prioritized variants from Exomiser.
3. HPOUtil.py - This script contains helper functions for calculating the gene rankings based on the cosine score from the Human Phenotype Ontology (HPO) terms.
4. OntologyUtil.py - This script contains helper functions for ontologies.  It is primarily support code for HPOUtil.py.
5. PVPUtil.py - This script contains helper functions for parsing prioritized variants from DeepPVP.
6. PhenGenUtil.py - This script contains helper function for parsing prioritized variants from Phen-Gen.
7. PyxisMapUtil.py - This script contains helper functions for retrieving gene rankings based on the PyxisMap ranks using the HPO terms.
8. SummaryDBUtil.py - This script contains helper functions for parsing the database dump files containing metadata for the Undiagnosed Diseases Network including samples IDs, HPO terms, and reported primary variants.  Note: these files are not available on GitHub due to Personal Health Information (PHI).

## Training/Testing Scripts
The following files perform the core workhorse training, testing, and reporting of results from VarSight:

1. TestLearners.py - This command line tool contains the primary test functions for gathering data, loading data, cleaning/reformatting data, training the models, testings the models, and reporting results.
    1. `gather` - This sub-routine will pre-gather ranks from PyxisMap and HPOUtil for our use during analysis.
    2. `analyze` - This sub-routine will actually perform the analysis and write the results to rendered .tex files that are automatically pulled into the LaTeX paper.  Here's a breakdown of options available in `analyze` mode (original paper results used `-g`):
    ```
    usage: TestLearners.py analyze [-h] [-R] [-P] [-e | -g | -r]

    optional arguments:
    -h, --help         show this help message and exit
    -R, --regenerate   regenerate all test results even if already available
    -P, --path-only    only uses pathogenic and likely pathogenic variants as
                        true positives
    -e, --exact-mode   use the dev-specified hyperparameters (single-execution)
    -g, --grid-mode    perform a grid search for the best hyperparameters (long
                        multi-execution)
    -r, --random-mode  perform a random search for the best hyperparameters
                        (short multi-execution)
    ```

## Paper Scripts
The following files were used to generate data and/or figure for the paper and supplementary documents:
1. SupplementGen.py - This script parses the Codicem filter JSON file and creates some .dot files that `graphviz` converts into figures for the supplementary document.  These figures visualize the filtering process that was by Codicem before returning results to analysts.  We use the variants that pass this filter for training and testing the classifiers that are a part of the core paper.
2. TestLearners.py - Refer to "Training/Testing Scripts" section