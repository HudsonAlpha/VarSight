# Paper Generation Files
This folder contains the content necessary to create the paper from LaTeX.  Note that common class files for LaTeX are not included here.

## Process for Future Revisions
1. Run `python3 scripts/TestLearners.py -R -g` to regenerate the results.  This will update `results.json` and render the templates in this sub-directory.
2. Manually update descriptive text in `main.tex` as appropriate for any core results changes.
3. Run LaTeX to PDF as normal on `main.tex` to update `main.pdf`.

## Results Documents
1. classifier_template.tex - contains the Jinja2 template for creating the performance table for the tested classifiers
2. data_template.tex - contains the Jinja2 template for creating the ranking results from single-value measures, external tools, and the tested classifiers
3. features_template.tex - contains the Jinja2 template for creating the RandomForest feature importance table
4. results.json - contains the raw results data from running `scripts/TestLearners.py`; this JSON is used by the Jinja2 templates in this subdirectory for rendering

## LaTeX Documents
1. main.pdf - the LaTeX generated document; tagged v1 was submitted to pre-print
2. main.tex - the LaTeX file used to generate main.pdf; uses rendered versions of the templates in this sub-directory to automatically populate core results (i.e. tables/figures)