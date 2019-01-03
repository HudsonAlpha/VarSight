'''
Primary analysis script for testing various learning methods
'''

import datetime
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import re

import CodiDumpUtil
import HPOUtil
import OntologyUtil
import PyxisMapUtil
import SummaryDBUtil

from imblearn.ensemble import BalancedRandomForestClassifier, RUSBoostClassifier, EasyEnsembleClassifier

from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier

PYXIS_DATE='12192018'   

def getPyxisMapResults():
    '''
    Script for downloading the PyxisMap results, update the date constant when needed
    '''
    ontFN = '/Users/matt/githubProjects/LayeredGraph/HPO_graph_data/hp.obo'
    annotFN = '/Users/matt/githubProjects/LayeredGraph/HPO_graph_data/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt'
    rootDir = '/Users/matt/githubProjects/VarSight/pyxis_ranks_'+PYXIS_DATE
    if not os.path.exists(rootDir):
        os.makedirs(rootDir)
    
    #we need the altIDMap
    node, edges, parents, altIDMap = OntologyUtil.loadGraphStructure(ontFN)
    hpoOnt = HPOUtil.HPOUtil(ontFN, annotFN)
    
    caseData = SummaryDBUtil.loadSummaryDatabase(altIDMap)
    for sl in sorted(caseData.keys()):
        #TODO: figure out how to properly handle this
        if '/' in sl:
            print('Skipping '+sl+' due to poor formatting of name.')
            continue
        
        #make sure we have HPO terms for the case
        if len(caseData[sl]['hpoTerms']) > 0:
            #pyxismap ranks
            jfn = rootDir+'/'+sl+'_pyxis.json'
            if os.path.exists(jfn):
                print(jfn+' already exists, skipping.')
            else:
                #get the data, make sure no errors were thrown and save it
                pyxisRanks, rankLen, response = PyxisMapUtil.getPyxisMapResults(caseData[sl]['hpoTerms'])
                
                if rankLen > 0:
                    j = {
                        'ranks' : pyxisRanks,
                        'rankLen' : rankLen,
                        'access' : str(datetime.datetime.utcnow())
                    }
                    
                    print('Saving results to '+jfn)
                    fpo = open(jfn, 'w+')
                    json.dump(j, fpo)
                    fpo.close()
                else:
                    print('Failed to retrieve results for '+sl)
                    print(response)
            
            #hpoutil ranks
            jfn = rootDir+'/'+sl+'_hpoutil.json'
            if os.path.exists(jfn):
                print(jfn+' already exists, skipping.')
            else:
                #get the data, make sure no errors were thrown and save it
                hpoRanks, hpoLen = hpoOnt.rank(caseData[sl]['hpoTerms'])
                
                if hpoLen > 0:
                    j = {
                        'ranks' : hpoRanks,
                        'rankLen' : hpoLen,
                        'access' : str(datetime.datetime.utcnow())
                    }
                    
                    print('Saving results to '+jfn)
                    fpo = open(jfn, 'w+')
                    json.dump(j, fpo)
                    fpo.close()
                else:
                    print('Failed to retrieve results for '+sl+' from HPOUtil')

def runAnalysis():
    '''
    Core test function
    '''
    pyxisRootDir = '/Users/matt/githubProjects/VarSight/pyxis_ranks_'+PYXIS_DATE
    
    #load the categorical metadata
    metaFN = '/Users/matt/githubProjects/VarSight/CODI_metadata/fields_metadata.json'
    fp = open(metaFN, 'rt')
    catMeta = json.load(fp)
    fp.close()

    #build some lookup dicts
    indexDict = {}
    indexLen = {}
    ignoreDict = {}
    catLabels = []
    
    #treat all categorical things the same other than where it is sourced from
    labelBlocks = ['allelicInformation', 'genes', 'variantTranscripts']

    for lg in labelBlocks:
        for d in catMeta[lg]:
            #make sure we never get the same category twice (shouldn't happen)
            assert(d['key'] not in indexDict)

            if d['interpret'] == 'single':
                #singles are easy
                catLabels.append(d['key'])
            elif d['interpret'] == 'multiple':
                #multiples, store the number of possiblities in indexLen and the index in indexDict
                indexDict[d['key']] = {}
                indexLen[d['key']] = len(d['accepted'])
                for i, v in enumerate(sorted(d['accepted'])):
                    #we are .lower()-ing everything, make sure that doesn't introduce a conflict
                    assert(v not in indexDict[d['key']])
                    indexDict[d['key']][v.lower()] = i
                ignoreDict[d['key']] = set([v.lower() for v in d.get('ignored', [])])
                catLabels += [d['key']+'_'+v for v in sorted(d['accepted'])]

    #These are the fields used from CODICEM dumps
    floatValues = [
        ('CADD Scaled', -1.0),
        ('phylop conservation', -30.0),
        ('phylop100 conservation', -30.0),
        ('phastcon100 conservation', -1.0),
        ('Mappability', -1.0),
        ('GERP rsScore', -13.0),
        ('Gnomad Exome AF', -1.0),
        ('Gnomad Exome Hom alt allele count', -1.0),
        ('Gnomad Exome Hemi alt allele count', -1.0),
        ('Gnomad Exome total allele count', -1.0),
        ('Gnomad Genome AF', -1.0),
        ('Gnomad Genome Hom alt allele count', -1.0),
        ('Gnomad Genome Hemi alt allele count', -1.0),
        ('Gnomad Genome total allele count', -1.0),
    ]
    fieldsOnly = set([l for l, d in floatValues]) | set([d['key'] for d in catMeta['allelicInformation']])
    
    #(label, reduce function, default)
    geneValues = [
        ('RVIS Score', min, 10.0),
        ('GHIS Score', max, -10.0),
        ('HIS Score', max, -10.0)
    ]
    geneFieldsOnly = set([l for l, r, d in geneValues]) | set([d['key'] for d in catMeta['genes']])
    
    #(label, default)
    seqValues = [
        ("percent reads", -1.0),
        ("total depth", -1.0)
    ]
    seqFieldsOnly = set([l for l, d in seqValues])
    
    #numerical variant-transcript information
    #"ADA Boost Splice Prediction" - looks like 0-1 with "NA" in there; TODO: what is good/bad?
    #"Random Forest Splice Prediction" - looks like 0-1 with "NA"; TODO: what is good/bad?

    #(label, reduce function, default)
    transValues = [
        ("ADA Boost Splice Prediction", max, -1.0),
        ("Random Forest Splice Prediction", max, -1.0)
    ]
    transLabelsOnly = set([l for l, r, d in transValues]) | set([d['key'] for d in catMeta['variantTranscripts']])

    #TODO: handle all categorical inputs we care about
    #possible good way(s) to do it: https://blog.myyellowroad.com/using-categorical-data-in-machine-learning-with-python-from-dummy-variables-to-deep-category-66041f734512
    #OHE - one-hot encoding (aka, each term is a bit)
    #feature hashing method, looks built in to sklearn

    #categorical allelicInformation
    #"Ensembl Regulatory Feature" - could be broken in ~6 boolean variables since each is a combination of categories
    #"Type" - 3 categories: SNV, INSERTION, DELETION
    #"HGMD association confidence" - list of comma-separated "High" or "Low"; also "NA" if nothing is there
    #"HGMD assessment type" - list of comma-separated "DFP", "DM", "FP", "DP" ("R" and "DM?" might show up FYI);
    #    after reading paper it looks like this: DM > DFP > DP > FP
    #"ClinVar Classification" - has many categories (https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/), but I think we only care about B, LB, V, LP, and P
    #    conflicting is another category that we should figure out how to interpret
    #"variant_attribute" - "Low Complexity Region", "NA", or "Simple Repeat"; might be worth including
    #CHECKED: this seems to not be in my files; "protein_alt" - looks boolean, "Protein Altering" or "NA"
    
    #categorical gene information
    #"Essentiality" - "Essential", "NA", "Neutral", and "Non-essential"; can likely be a bucket since there are possibly multiple genes
    #TODO: OMIM and HGMD disease IDs - could be distilled into a count, do we care about these?

    #categorical variant-transcript information
    #"Meta Svm Prediction" - "D", "T", and "NA"; D = damaging, T = tolerated
    #"PolyPhen HV Prediction" - "B", "D", "P", and "NA"; D - probably damaging, P - possibly damaging, B - benign
    #"PolyPhen HD Prediction" - same categories as above
    #"Provean Prediction" - "D", "N", and "NA"; D - damaging, N - neutral
    #"SIFT Prediction" - "D", "T", and "NA"; D - damaging, T - tolerated
    #"Effects" - multiple categories (example: "AA Deletion", "Intergenic", etc.); seems to break LogReg when we added this
    #"Affected Regions" - multiple categories (example: "3' Coding Exon", "3' UTR intron", etc.), do we want this or not?
    
    #we need the altIDMap - not sure we actually need it here since the HPO query was done beforehand, but I suppose it doesn't hurt
    #node, edges, altIDMap = PyxisMapUtil.loadGraphStructure()
    ontFN = '/Users/matt/githubProjects/LayeredGraph/HPO_graph_data/hp.obo'
    nodes, edges, parents, altIDMap = OntologyUtil.loadGraphStructure(ontFN)
    
    #1 - load the summary database
    caseData = SummaryDBUtil.loadSummaryDatabase(altIDMap)
    
    #2 - load each CODI dump that HAS a result from the summary database and reformat the data
    allValues = []
    allCatValues = []
    allCatBreak = {}
    allClassifications = []
    for sl in sorted(caseData.keys()):
        #we only need to load those which actually have a return
        if len(caseData[sl]['primaries']) == 0:
            print('Skipping '+sl+', no primaries')
            continue
        primarySet = set([p['variant'] for p in caseData[sl]['primaries']])
        #hpoTerms = caseData[sl]['hpoTerms']
        
        #if there is no CODICEM dump, we obviously can't do anything either
        codiDump = CodiDumpUtil.loadCodiDump(sl, fieldsOnly, geneFieldsOnly, seqFieldsOnly, transLabelsOnly)
        if len(codiDump) == 0:
            print('Skipping '+sl+', empty or absent CODICEM dump')
            continue
        
        #finally, make sure we have a PyxisMap dump we can use
        pyxisFN = pyxisRootDir+'/'+sl+'_pyxis.json'
        if not os.path.exists(pyxisFN):
            print('Skipping '+sl+', no PyxisMap dump available')
            continue
        fp = open(pyxisFN, 'r')
        pyxisJson = json.load(fp)
        fp.close()
        
        print("Loading "+sl+"...")

        #reformat for easy lookups
        pyxisRanks = pyxisJson['ranks']
        pyxisLen = pyxisJson['rankLen']
        
        #now do hpoUtil also
        hpoFN = pyxisRootDir+'/'+sl+'_hpoutil.json'
        if not os.path.exists(hpoFN):
            print('Skipping '+sl+', no HPOUtil dump available')
            continue
        fp = open(hpoFN, 'r')
        hpoJson = json.load(fp)
        fp.close()
        
        hpoRanks = hpoJson['ranks']
        hpoLen = hpoJson['rankLen']
        
        #now, we have all the pieces, put them together into tests
        values = []
        catValues = []
        catBreak = {}
        classifications = []
        
        print('Primary set: ', primarySet)
        foundPrimaries = set([])
        
        for variant in codiDump:
            vData = []
            catData = []
            
            #determine whether this variant is primary or not
            vName = variant['variant'].replace('Chr', 'chr')
            isPrimary = (1.0 if (vName in primarySet) else 0.0)
            if isPrimary:
                foundPrimaries.add(vName)
            
            #first, add in the best PyxisMap rank
            geneList = variant['genes']
            pyxisBest = min([pyxisRanks.get(g, (pyxisLen+1, 0.0))[0] / pyxisLen for g in geneList])
            #vData.append(1.0 if pyxisBest <= 1.0 else 0.0)
            vData.append(pyxisBest)
            
            #now, add the best HPOUtil rank
            hpoBest = min([hpoRanks.get(g, (hpoLen+1, 0.0))[0] / hpoLen for g in geneList])
            #TODO: if missing data marked
            vData.append(hpoBest)
            
            #now add all float values from CODICEM
            for fvl, invalidValue in floatValues:
                if variant[fvl] == 'NA':
                    #vData.append(0.0)
                    vData.append(invalidValue)
                else:
                    #vData.append(1.0)
                    vData.append(float(variant[fvl]))

            #add things from the gene values
            for gvl, reduceFunc, invalidValue in geneValues:
                bestVal = reduceFunc([float(v) if v != 'NA' else invalidValue for v in variant[gvl]])
                vData.append(bestVal)
            
            #add things from the seq values
            for svl, invalidValue in seqValues:
                if variant[svl] == 'NA':
                    vData.append(invalidValue)
                else:
                    vData.append(float(variant[svl]))
            
            #add things from the transcript values
            for tvl, reduceFunc, invalidValue in transValues:
                bestVal = reduceFunc([float(v) if v != 'NA' else invalidValue for v in variant[tvl]])
                vData.append(bestVal)

            #add categorical information
            #for cDict in catMeta['allelicInformation']:
            for lg in labelBlocks:
                for cDict in catMeta[lg]:
                    fieldKey = cDict['key']
                    fieldType = cDict['interpret']
                    if fieldKey not in catBreak:
                        catBreak[fieldKey] = []

                    if fieldType == 'single':
                        #singles - add a single float value corresponding to the dictionary lookup in the JSON
                        fieldDict = cDict['values']
                        catData.append(float(fieldDict[variant[fieldKey]]))
                        catBreak[fieldKey].append(float(fieldDict[variant[fieldKey]]))
                    elif fieldType == 'multiple':
                        #multiple - create a bin for each category, the array values are bincounts
                        arr = [0.0]*indexLen[fieldKey]
                        if fieldKey in ["ClinVar Classification"]:
                            #split on comma or '/'
                            csvVals = re.split('[,/]', variant[fieldKey].lower())
                        elif fieldKey in ["Ensembl Regulatory Feature"]:
                            #allow missing values here
                            csvVals = variant.get(fieldKey, 'NA').lower().split(';')
                        elif type(variant[fieldKey]) == list:
                            if fieldKey in ["Effects", "Affected Regions"]:
                                #list of CSV values
                                csvVals = []
                                for vfk in variant[fieldKey]:
                                    csvVals += vfk.lower().split(',')
                            else:
                                #we just need to lower it down
                                csvVals = [vfk.lower() for vfk in variant[fieldKey]]
                        else:
                            #split on commas only
                            csvVals = variant[fieldKey].lower().split(',')
                        
                        for v in csvVals:
                            #remove white space around the values
                            if (v.strip() not in ignoreDict[fieldKey]):
                                try:
                                    arr[indexDict[fieldKey][v.strip()]] += 1
                                except:
                                    print(variant[fieldKey], csvVals)
                                    raise Exception("Unexpected key, see above")
                        catData += arr
                        catBreak[fieldKey].append(arr)
                    else:
                        raise Exception("Unexpected interpretation type: "+fieldType)

            #uncomment to see each variant's info
            #print(vName, isPrimary, vData, sep='\t')
            
            #values.append(vData+catData)
            values.append(vData)
            catValues.append(catData)
            classifications.append(isPrimary)
        
        if sum(classifications) != len(primarySet):
            print('Some primaries were not found')
            print('All: ', primarySet)
            print('Found: ', foundPrimaries)
            #exit()
        
        if sum(classifications) > 0:
            allValues.append(values)
            allCatValues.append(catValues)
            for k in catBreak:
                if k not in allCatBreak:
                    allCatBreak[k] = []
                allCatBreak[k] += catBreak[k]
            allClassifications.append(classifications)
        else:
            print('Ignoring all variants from '+sl+', no primaries found in CODICEM dump')
        
        print()
    
    #determine the category mode
    OHE_MODE = 0
    PCA_MODE = 1
    PCA_BREAK_MODE = 2
    currentCatMode = PCA_BREAK_MODE

    featureLabels = (['PyxisMap', 'HPOUtil']+
        [l for l, d in floatValues]+
        #[d['key'] for d in catMeta['allelicInformation']]+
        [l for l, r, d in geneValues]+
        [l for l, d in seqValues]+
        [l for l, r, d in transValues]
    )

    if currentCatMode == OHE_MODE:
        #basic *-hot encoding mode, this was our default originally
        vVals = np.vstack(allValues)
        vCat = np.vstack(allCatValues)
        xFinal = np.hstack((vVals, vCat))
        featureLabels += catLabels

    elif currentCatMode == PCA_MODE:
        vVals = np.vstack(allValues)
        vCat = np.vstack(allCatValues)

        numComp = 20 #set to None to do all
        #in our initial tests, 20 gets 90% of the variance
        
        pca = PCA(numComp)
        pca.fit(vCat)
        print('pca-explained', pca.explained_variance_ratio_)
        print('pca-cumsum', np.cumsum(pca.explained_variance_ratio_))
        pcaCat = pca.transform(vCat)
        xFinal = np.hstack((vVals, pcaCat))
        featureLabels += ['PCA-'+str(x) for x in range(1, pcaCat.shape[1]+1)]

    elif currentCatMode == PCA_BREAK_MODE:
        vVals = np.vstack(allValues)
        numComp = 2
        pcaBlocks = []

        for lg in labelBlocks:
            for cDict in catMeta[lg]:
                fieldKey = cDict['key']
                fieldType = cDict['interpret']

                if fieldType == 'single':
                    #single value, no PCA required
                    pcaBlocks.append(np.array(allCatBreak[fieldKey]).reshape((len(allCatBreak[fieldKey]), 1)))
                    featureLabels.append(fieldKey)
                elif fieldType == 'multiple':
                    rawStackVals = np.vstack(allCatBreak[fieldKey])
                    pca = PCA(numComp)
                    pca.fit(rawStackVals)
                    print(fieldKey, 'pca-explained', pca.explained_variance_ratio_)
                    pcaCat = np.array(pca.transform(rawStackVals))
                    pcaBlocks.append(pcaCat)
                    featureLabels.append(fieldKey+'-PCA1')
                    featureLabels.append(fieldKey+'-PCA2')
                else:
                    raise Exception('Unexpected fieldType')
        xFinal = np.hstack([vVals]+pcaBlocks)
    else:
        raise Exception('Unexpected currentCatMode')

    #3 - train and test the learners
    #xFinal = np.vstack(allValues)
    yFinal = np.array(np.hstack(allClassifications), dtype='int64')

    runClassifiers(xFinal, yFinal, featureLabels)

def runClassifiers(values, classifications, featureLabels):
    '''
    @param values - a matrix with R rows and C columns, where there are "C" features
    @param classifications - an array of length R corresponding to the above values
    '''
    print('Values:', values.shape)
    print('Classifications (and bincount):', classifications.shape, np.bincount(classifications))
    
    #split the data
    train_x, test_x, train_y, test_y = train_test_split(values, classifications, test_size=0.25, stratify=classifications)
    
    #class_weight="balanced" is when there is a major imbalance between the number of true negatives and true positives
    EXACT_MODE=0
    GRID_MODE=1
    RANDOM_MODE=2
    currentMode=EXACT_MODE

    classifiers = [
        ('rf_bal', RandomForestClassifier(random_state=0, class_weight='balanced', max_depth=4, n_estimators=300, min_samples_split=2, max_features='sqrt'), 
        {
            'random_state' : [0],
            'class_weight' : ['balanced'],
            'n_estimators' : [100, 200, 300],
            'max_depth' : [2, 3, 4],
            'min_samples_split' : [2, 3],
            'max_features' : ["sqrt", "log2"]
        }),
        #('svc_bal', SVC(probability=True, class_weight="balanced")),#this one doesn't scale well past 10k samples
        #('mlp', MLPClassifier()),#this one doesn't work, presumably because there is no balanced option
        ('LogReg_bal', LogisticRegression(random_state=0, class_weight='balanced', penalty='l2', C=0.01, solver='newton-cg', max_iter=200), 
        {
            'random_state' : [0],
            'class_weight' : ['balanced'],
            'penalty' : ['l2'],
            'C' : [0.01, 0.1, 1.0, 10.0, 100.0],
            'solver' : ['newton-cg', 'liblinear'],
            'max_iter' : [200]
        }),
        ('imb_rf', BalancedRandomForestClassifier(random_state=0, n_estimators=300, max_depth=4, min_samples_split=2, max_features='sqrt'), 
        {
            'random_state' : [0],
            'n_estimators' : [100, 200, 300],
            'max_depth' : [2, 3, 4],
            'min_samples_split' : [2, 3],
            'max_features' : ["sqrt", "log2"]
        }),
        #('imb_rus', RUSBoostClassifier(random_state=0)), #this one never seems to have a good result on the test set, I think it's overfitting due to boosting
        ('imb_ee', EasyEnsembleClassifier(random_state=0, n_estimators=50), 
        {
            'random_state' : [0],
            'n_estimators' : [10, 20, 30, 40, 50]
        })
    ]
    
    #things to calculate
    aucs = []
    rocs = []
    prs = []
    pr_aucs = []
    
    for clf_label, raw_clf, hyperparam in classifiers:
        #test the classifier?
        print(clf_label)
        print('\ttraining...')

        if currentMode == EXACT_MODE:
            clf = raw_clf
        elif currentMode == GRID_MODE:
            #grid search over hyperparameters
            clf = GridSearchCV(raw_clf, hyperparam, cv=10, scoring='balanced_accuracy', n_jobs=-1)
        elif currentMode == RANDOM_MODE:
            clf = RandomizedSearchCV(raw_clf, hyperparam, cv=10, scoring='balanced_accuracy', n_jobs=-1)
        else:
            raise Exception("Unexpected classifier mode")
        
        #regardless of above approach, we now fit it
        clf.fit(train_x, train_y)
        if currentMode in [GRID_MODE, RANDOM_MODE]:
            print('\tBest params:', clf.best_params_)

        #I was using this wrong, compute_sample_weight is if we have more trustworthy samples in our data; we might use this eventually but not yet
        #clf.fit(train_x, train_y, compute_sample_weight('balanced', train_y))
        
        try:
            #print('\tfeature_important', clf.feature_importances_)
            print('\tfeature_important:')#, clf.best_estimator_.feature_importances_)
            for j, l in enumerate(featureLabels):
                #print('', '', l, np.sum(clf.feature_importances_[2*j:2*(j+1)]), clf.feature_importances_[2*j:2*(j+1)], sep='\t')
                if currentMode == EXACT_MODE:
                    print('', '', l, clf.feature_importances_[j], sep='\t')
                elif currentMode in [GRID_MODE, RANDOM_MODE]:
                    print('', '', l, clf.best_estimator_.feature_importances_[j], sep='\t')
        except:
            pass
        
        #balanced accuracy - an accuracy score that is an average across the classes
        print('\tbalanced_train_acc', balanced_accuracy_score(train_y, clf.predict(train_x)))
        print('\tbalanced_test_acc', balanced_accuracy_score(test_y, clf.predict(test_x)))
        
        #confusion matrix - exactly what it sounds like, 2x2 grid in our case
        print('\tconf_matrix_train', confusion_matrix(train_y, clf.predict(train_x)))
        print('\tconf_matrix_test', confusion_matrix(test_y, clf.predict(test_x)))
        
        #roc_curve stuff - could be misleading due to imbalance
        y_pred_rf = clf.predict_proba(test_x)[:, 1]
        false_positive_rate, true_positive_rate, thresholds = roc_curve(test_y, y_pred_rf)
        roc_auc = auc(false_positive_rate, true_positive_rate)
        aucs.append(roc_auc)
        rocs.append((false_positive_rate, true_positive_rate))
        print('\troc_auc', roc_auc)
        
        #precision-recall curve stuff - should be less misleading
        precision, recall, pr_thresholds = precision_recall_curve(test_y, y_pred_rf)
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)
        prs.append((recall, precision))
        print('\tpr_auc', pr_auc)
        
        if currentMode == EXACT_MODE:
            #cross validation scores
            cv = StratifiedKFold(n_splits=10)
            
            #balanced accuracy - "The balanced accuracy in binary and multiclass classification problems to deal with imbalanced datasets. It is defined as the average of recall obtained on each class."
            bal_cvs_scores = cross_val_score(clf, train_x, train_y, cv=cv, scoring='balanced_accuracy')
            print('\tbal_cvs_scores', bal_cvs_scores)
            print("\tbal_Accuracy: %0.4f (+/- %0.4f)" % (bal_cvs_scores.mean(), bal_cvs_scores.std() * 2))
            
            #f1_weighted - see https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html#sklearn.metrics.f1_score
            f1w_cvs_scores = cross_val_score(clf, train_x, train_y, cv=cv, scoring='f1_weighted')
            print('\tf1w_cvs_scores', f1w_cvs_scores)
            print("\tf1w_Accuracy: %0.4f (+/- %0.4f)" % (f1w_cvs_scores.mean(), f1w_cvs_scores.std() * 2))
        elif currentMode in [GRID_MODE, RANDOM_MODE]:
            #CV is already calculated here so we can just dump it out
            print('\tscoring systems:')
            means = clf.cv_results_['mean_test_score']
            stds = clf.cv_results_['std_test_score']
            for mean, std, params in zip(means, stds, clf.cv_results_['params']):
                print("\t\t%0.3f (+/-%0.03f) for %r"
                    % (mean, std * 2, params))
        print()
    
    #ROC curve
    plotFN = '/Users/matt/githubProjects/VarSight/images/codi_rf_roc.png'
    plt.figure()
    plt.plot([0, 1], [0, 1], 'k--')
    
    for i, (label2, raw_clf, raw_params) in enumerate(classifiers):
        plt.plot(rocs[i][0], rocs[i][1], label=('%s (%0.4f)' % (label2, aucs[i])))
    
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve - all data')
    plt.legend(loc='best')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid()
    plt.savefig(plotFN)
    plt.close()
    
    #precision recall curve
    plotFN = '/Users/matt/githubProjects/VarSight/images/codi_rf_pr.png'
    plt.figure()
    
    for i, (label2, raw_clf, raw_params) in enumerate(classifiers):
        plt.plot(prs[i][0], prs[i][1], label=('%s (%0.4f)' % (label2, pr_aucs[i])))
    
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PR curve - all data')
    plt.legend(loc='best')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid()
    plt.savefig(plotFN)
    plt.close()

if __name__ == '__main__':
    #this will re-run the PyxisMap tests, should update dates though
    #getPyxisMapResults()
    
    #core tests
    runAnalysis()
    