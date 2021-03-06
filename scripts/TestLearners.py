'''
Primary analysis script for testing various learning methods
'''

import argparse as ap
import datetime
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import re

import CodiDumpUtil
import ExomiserUtil
import HPOUtil
import OntologyUtil
import PhenGenUtil
import PVPUtil
import PyxisMapUtil
import SummaryDBUtil

from imblearn.ensemble import BalancedRandomForestClassifier, RUSBoostClassifier, EasyEnsembleClassifier

from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectKBest, SelectFdr, SelectFpr, SelectFwe
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import make_scorer
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

#CONSTANTS
PYXIS_DATE='01042019'   

#mode for analyzing the classifiers
EXACT_MODE=0
GRID_MODE=1
RANDOM_MODE=2

def getPyxisMapResults():
    '''
    Script for pre-computing HPO-based results to the "pyxis_ranks_<PYXIS_DATE>" subfolder specified in the 
    constant above.  This will do two things: 1) hit the PyxisMap endpoint specified in PyxisMapUtil and
    2) perform the cosine score computation from HPOUtil.  Outputs are saved to files so we don't run this 
    every time we want to train a new model.
    '''
    ontFN = '/Users/matt/githubProjects/LayeredGraph/HPO_graph_data/hp.obo'
    annotFN = '/Users/matt/githubProjects/LayeredGraph/HPO_graph_data/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt'
    rootDir = '/Users/matt/githubProjects/VarSight/pyxis_ranks_'+PYXIS_DATE
    if not os.path.exists(rootDir):
        os.makedirs(rootDir)
    
    #we need the altIDMap
    node, edges, parents, altIDMap = OntologyUtil.loadGraphStructure(ontFN)
    hpoOnt = HPOUtil.HPOUtil(ontFN, annotFN)
    
    #load case data and go through them searching for HPO terms
    caseData = SummaryDBUtil.loadSummaryDatabase(altIDMap, False)
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

def loadFormattedData(args):
    '''
    This function is primarily about loading and formatting data prior to training/testing any classifiers.
    @param args - arguments from the command line
    @return - tuple (xFinal, yFinal, featureLabels, startIndices, allRepDicts, externalRanks)
        xFinal - an NxM matrix where N is the number of variants in our test/training set and M is the number of features; contains feature values
        yFinal - an N length array where N is the number of variants in our test/training set; 1 if the variant was reported
        catMeta - the categorical metadata 
        featureLabels - an M length array containing feature labels for our output benefit
        startIndices - a (C+1) length array containing the start indices of variants from an individual case
        allRepDicts - an N length array of None or dictionaries; if a dictionary, it contains data on a returned variant
        externalRanks - a dictionary where key is a label and value is a list of lists containing exomiser ranks on a per-case basis; missing values are -1*len(ranked variants)
    '''
    pyxisRootDir = '/Users/matt/githubProjects/VarSight/pyxis_ranks_'+PYXIS_DATE
    
    #load the feature metadata
    metaFN = '/Users/matt/githubProjects/VarSight/CODI_metadata/fields_metadata.json'
    fp = open(metaFN, 'rt')
    catMeta = json.load(fp)
    fp.close()

    #we need the altIDMap - not sure we actually need it here since the HPO query was done beforehand, but I suppose it doesn't hurt
    ontFN = '/Users/matt/githubProjects/LayeredGraph/HPO_graph_data/hp.obo'
    nodes, edges, parents, altIDMap = OntologyUtil.loadGraphStructure(ontFN)
    
    #1 - load the summary database 
    caseData = SummaryDBUtil.loadSummaryDatabase(altIDMap, args.path_only)
    
    #2 - load each CODI dump that HAS a result from the summary database and reformat the data
    allValues = []
    allCatValues = []
    allCatBreak = {}
    allClassifications = []
    allRepDicts = []
    
    #this will track the breaks between cases
    startIndices = [0]

    exomiserRanks = []
    exomiserHumanRanks = []
    pvpRanks = []
    phengenRanks = []

    for sl in sorted(caseData.keys()):
        #we only need to load those which actually have a return
        if len(caseData[sl]['primaries']) == 0:
            #I don't think we care about seeing these at all
            #print(sl, 'NO_PRIMARIES', sep=',')
            continue
        
        if len(caseData[sl]['hpoTerms']) == 0:
            #we don't care to see these either, just incomplete data :(
            #print(sl, 'NO_HPOTERMS', sep=',')
            continue

        primarySet = set([p['variant'] for p in caseData[sl]['primaries']])
        primaryDict = {p['variant'] : p for p in caseData[sl]['primaries']}
        assert(len(primaryDict) == len(caseData[sl]['primaries']))
        
        #if there is no CODICEM dump, we obviously can't do anything either
        #codiDump = CodiDumpUtil.loadCodiDump(sl, fieldsOnly, geneFieldsOnly, seqFieldsOnly, transLabelsOnly)
        #if len(codiDump) == 0:
        #    print(sl, 'NO_CODICEM_DUMP', sep=',')
        #    continue
        
        #make sure the codi dump exists
        codiDumpFN = '/Users/matt/githubProjects/VarSight/CLI_primary_V6/'+sl+'_results.json'
        if not os.path.exists(codiDumpFN):
            #print('Missing file: '+fn)
            continue

        #finally, make sure we have a PyxisMap dump we can use
        pyxisFN = pyxisRootDir+'/'+sl+'_pyxis.json'
        if not os.path.exists(pyxisFN):
            print(sl, 'NO_PYXISMAP_DUMP', sep=',')
            continue
        fp = open(pyxisFN, 'r')
        pyxisJson = json.load(fp)
        fp.close()
        
        #reformat for easy lookups
        pyxisRanks = pyxisJson['ranks']
        pyxisLen = pyxisJson['rankLen']
        
        #now do hpoUtil also
        hpoFN = pyxisRootDir+'/'+sl+'_hpoutil.json'
        if not os.path.exists(hpoFN):
            print(sl, 'NO_HPOUTIL_DUMP', sep=',')
            continue
        fp = open(hpoFN, 'r')
        hpoJson = json.load(fp)
        fp.close()
        
        #print("Loading "+sl+"...")

        hpoRanks = hpoJson['ranks']
        hpoLen = hpoJson['rankLen']
        
        (labelBlocks, values, catValues, featureLabels, catLabels, catBreak, codiDump) = CodiDumpUtil.prepDataShared(catMeta, codiDumpFN, pyxisRanks, pyxisLen, hpoRanks, hpoLen)
        
        #parse primary stuff
        classifications = []
        repDicts = []
        foundPrimaries = set([])
        fpCPRA = []
        
        for variant in codiDump:
            #determine whether this variant is primary or not
            vName = variant['variant'].replace('Chr', 'chr')
            isPrimary = (1.0 if (vName in primarySet) else 0.0)
            if isPrimary:
                foundPrimaries.add(vName)
                
                #fpCPRA is the same order as the repDicts (needed for later)
                fpCPRA.append((str(variant['chromosome']), str(variant['position']), variant['ref allele'], variant['alt allele']))
                repDicts.append(primaryDict[vName])
            else:
                repDicts.append(None)
            
            classifications.append(isPrimary)

        #occasionally primaries are missing, could be filter change OR from some targeted search by an analyst
        if sum(classifications) != len(primarySet):
            print(sl, 'SOME_MISSING', primarySet, foundPrimaries, sep=',')

        #if at least one primary is found, we will add to case to our test data
        if sum(classifications) > 0:
            allValues.append(values)
            startIndices.append(startIndices[-1]+len(values))
            allCatValues.append(catValues)
            for k in catBreak:
                if k not in allCatBreak:
                    allCatBreak[k] = []
                allCatBreak[k] += catBreak[k]
            allClassifications.append(classifications)
            allRepDicts += repDicts

            #if we are using this case, load info from exomiser for comparison later
            exRanks = ExomiserUtil.getTargetRanks(sl, fpCPRA, 'hiPhive')
            exomiserRanks.append(exRanks)

            exHum = ExomiserUtil.getTargetRanks(sl, fpCPRA, 'hiPhive_human')
            exomiserHumanRanks.append(exHum)

            pvpRank = PVPUtil.getTargetRanks(sl, fpCPRA)
            pvpRanks.append(pvpRank)

            phengenRank = PhenGenUtil.getTargetRanks(sl, fpCPRA)
            phengenRanks.append(phengenRank)

        else:
            print(sl, 'ALL_MISSING', primarySet, foundPrimaries, sep=',')
        
        #print()
    
    #determine the category mode
    OHE_MODE = 0 #categories are given *-hot encodings where * is the count of how many times that label appears
    PCA_MODE = 1 #all categories are PCA-ed together and the top X PCA components are used
    PCA_BREAK_MODE = 2 #each category has its own PCA, only X PCA components are allowed per category
    currentCatMode = OHE_MODE#PCA_BREAK_MODE

    if currentCatMode == OHE_MODE:
        #basic *-hot encoding mode, this was our default originally
        vVals = np.vstack(allValues)
        vCat = np.vstack(allCatValues)
        xFinal = np.hstack((vVals, vCat))
        featureLabels += catLabels

    elif currentCatMode == PCA_MODE:
        #PCA all categorical values together
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
        #PCA each category individually and allow PCA1 and PCA2 as inputs
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
                elif fieldType == 'float' or fieldType == 'float_reduce':
                    #TODO: does something need to be moved here for consistency if we refactor code?
                    pass
                else:
                    raise Exception('Unexpected fieldType')
        xFinal = np.hstack([vVals]+pcaBlocks)
    else:
        raise Exception('Unexpected currentCatMode')

    #return the values
    yFinal = np.array(np.hstack(allClassifications), dtype='int64')

    externalRankings = {
        'Exomiser(hiPhive)' : exomiserRanks, 
        'Exomiser(hiPhive, human only)' : exomiserHumanRanks,
        'DeepPVP' : pvpRanks,
        'Phen-Gen' : phengenRanks
    }

    return (xFinal, yFinal, catMeta, featureLabels, startIndices, allRepDicts, externalRankings)

def runClassifiers(args, values, classifications, featureLabels, startIndices, allRepDicts, externalRanks):
    '''
    @param args - any arguments from the command line argparse can be accessed here
    @param values - a matrix with R rows and C columns, where there are "C" features
    @param classifications - an array of length R corresponding to the above values
    @param featureLabels - C length array containing labels (strings) for the features
    @param startIndices - the startIndices of individual cases in the training set
    @param allRepDicts - a list of dictionaries for variants that were reported; non-reported vars are None
    @param externalRanks - a dictionary where key is a label and value is a paired exomiser ranks for each case; needs to be split up with train/test; missing values are -1*len(ranked variants)
    @return tuple
        resultsDict - a dictionary containing many results we wish to include in a paper
        trainedClassifierResults - a dictionary containing the models and features
    '''
    #TODO: I should consider ways to break this function apart; for sanity if nothing else
    resultsDict = {}
    resultsDict['FEATURE_LABELS'] = featureLabels
    print('Values:', values.shape)
    print('Classifications (and bincount):', classifications.shape, np.bincount(classifications))
    print('Num cases:', len(startIndices)-1)

    #split the data
    CASE_BASED_SPLIT = True
    TEST_SIZE = 0.5

    pDict = {
        'VARIANT_OF_UNCERTAIN_SIGNIFICANCE' : 3,
        'LIKELY_PATHOGENIC' : 4,
        'PATHOGENIC' : 5
    }
    if args.path_only:
        pList = ['LIKELY_PATHOGENIC', 'PATHOGENIC']
    else:
        pList = ['VARIANT_OF_UNCERTAIN_SIGNIFICANCE', 'LIKELY_PATHOGENIC', 'PATHOGENIC']
    

    if CASE_BASED_SPLIT:
        #we decided to split the train/test by case, going to need some custom stuff here
        #we want to make sure ~25% of the true positives are in the test set and the remaining in the training set
        #false positives will likely be out of balance but that's okay because we have a metric crap ton of them
        trainXArray = []
        trainYArray = []
        train_indices = [0]
        train_dicts = []
        tpTrain = 0
        #exomiserTrain = []
        externalTrain = {k : [] for k in externalRanks.keys()}

        testXArray = []
        testYArray = []
        test_indices = [0]
        test_dicts = []
        tpTest = 0
        #exomiserTest = []
        externalTest = {k : [] for k in externalRanks.keys()}

        #ratio of TRAIN:TEST
        invRat = 1.0/TEST_SIZE - 1

        for x in range(0, len(startIndices)-1):
            st = startIndices[x]
            et = startIndices[x+1]
            tpCount = np.sum(classifications[st:et])
            if tpTrain <= invRat*tpTest:
                #we need more in our training
                trainXArray.append(values[st:et])
                trainYArray.append(classifications[st:et])
                train_indices.append(train_indices[-1]+(et-st))
                train_dicts += allRepDicts[st:et]
                tpTrain += tpCount

                #we moved the data manipulation down
                #exomiserTrain.append(exomiserRanks[x])
                for k in externalRanks.keys():
                    externalTrain[k].append(externalRanks[k][x])

            else:
                #we need more in our test
                testXArray.append(values[st:et])
                testYArray.append(classifications[st:et])
                test_indices.append(test_indices[-1]+(et-st))
                test_dicts += allRepDicts[st:et]
                tpTest += tpCount

                #we moved the data manipulation down
                #exomiserTest.append(exomiserRanks[x])
                for k in externalRanks.keys():
                    externalTest[k].append(externalRanks[k][x])
        
        #join them all together
        train_x = np.vstack(trainXArray)
        train_y = np.hstack(trainYArray)
        test_x = np.vstack(testXArray)
        test_y = np.hstack(testYArray)

        '''
        TODO: do we care about this? seems to make most results worse overall and doesn't seem very systematic; we
        should do some pruning in the long run, but it seems like each method will benefit from a different type
        (or maybe amount) of pruning.  I think leave it out for now unless this becomes a deal-breaker with some
        ornery reviewer.
        '''
        #more details here: https://scikit-learn.org/stable/modules/feature_selection.html
        #doing this seems to help with LogisticRegression only, but doesn't bring it up to our current results
        #select the best parameters based on training only
        FS_NONE = 0
        FS_SELECT_K_BEST = 1
        #FS_RFECV = 2 #this one took forever to run, and underperformed
        FS_MODEL_SELECTION = 3
        FEATURE_SELECTION_TYPE = FS_SELECT_K_BEST

        if FEATURE_SELECTION_TYPE == FS_NONE:
            pass

        elif FEATURE_SELECTION_TYPE == FS_SELECT_K_BEST:
            if False:
                for selectorType in [SelectKBest, SelectFdr, SelectFpr, SelectFwe]:
                    #skb = SelectKBest(mutual_info_classif, k=20)
                    skb = selectorType()
                    skb.fit(train_x, train_y)
                    #for i, f in enumerate(featureLabels):
                        #print('', f, -np.log10(skb.pvalues_[i]), sep='\t')
                    #    print('', f, skb.scores_[i], sep='\t')
                    selected = skb.get_support(False)
                    selected2 = skb.get_support(True)
                    #train_x = train_x[:, selected]
                    #test_x = test_x[:, selected]
                    #featureLabels = [featureLabels[x] for x in selected2]
                    #resultsDict['FEATURE_LABELS'] = featureLabels
                    print(selectorType)
                    print('Selected features:', [featureLabels[x] for x in selected2])
                exit()
            else:
                #FWE link: https://stats.stackexchange.com/questions/328358/fpr-fdr-and-fwe-for-feature-selection
                #skb = SelectFwe()
                skb = SelectKBest(k=20)
                skb.fit(train_x, train_y)
                selected = skb.get_support(False)
                selected2 = skb.get_support(True)
                train_x = train_x[:, selected]
                test_x = test_x[:, selected]
                featureLabels = [featureLabels[x] for x in selected2]
                resultsDict['FEATURE_LABELS'] = featureLabels
                print('Selected features:', featureLabels)

        elif FEATURE_SELECTION_TYPE == FS_MODEL_SELECTION:
            print('Running SelectFromModel...')
            lsvc = LinearSVC(C=0.0001, penalty="l1", dual=False, random_state=0, class_weight="balanced", max_iter=5000)
            standardScaler = StandardScaler()
            standardScaler.fit(train_x)
            scaled_x = standardScaler.transform(train_x)
            lsvc.fit(scaled_x, train_y)
            sfModel = SelectFromModel(lsvc, prefit=True)
            #print("Optimal number of features : %d" % rfecv.n_features_)
            selected = sfModel.get_support(False)
            selected2 = sfModel.get_support(True)
            train_x = train_x[:, selected]
            test_x = test_x[:, selected]
            print("Train/test size:", train_x.shape, test_x.shape)
            featureLabels = [featureLabels[x] for x in selected2]
            resultsDict['FEATURE_LABELS'] = featureLabels
        
    else:
        #train without any regards to case labels
        train_x, test_x, train_y, test_y = train_test_split(values, classifications, test_size=TEST_SIZE, stratify=classifications)

    print('split sizes', train_x.shape, test_x.shape, train_y.shape, test_y.shape)
    
    #save a bunch of data based on what was trained/tested on
    resultsDict['TRAIN_SHAPE'] = train_x.shape
    resultsDict['TEST_SHAPE'] = test_x.shape
    resultsDict['TRAIN_TP'] = np.sum(train_y)
    resultsDict['TEST_TP'] = np.sum(test_y)
    resultsDict['TRAIN_CASES'] = len(train_indices)-1
    resultsDict['TEST_CASES'] = len(test_indices)-1

    #class_weight="balanced" is when there is a major imbalance between the number of true negatives and true positives
    currentMode=args.training_mode

    #record the training method used
    resultsDict['TRAINING_MODE'] = ['EXACT_MODE', 'GRID_MODE', 'RANDOM_MODE'][currentMode]

    classifiers = [
        ('RandomForest(sklearn)', RandomForestClassifier(random_state=0, class_weight='balanced', max_depth=3, n_estimators=100, min_samples_split=2, max_features='sqrt'), 
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
        ('LogisticRegression(sklearn)', LogisticRegression(random_state=0, class_weight='balanced', penalty='l2', C=10.0, solver='newton-cg', max_iter=200), 
        {
            'random_state' : [0],
            'class_weight' : ['balanced'],
            'penalty' : ['l2'],
            'C' : [0.01, 0.1, 1.0, 10.0, 100.0],
            'solver' : ['newton-cg', 'liblinear'],
            'max_iter' : [200]
        }),
        #('ExtraTrees(sklearn)', ExtraTreesClassifier(random_state=0, class_weight='balanced'),
        #{
        #    'random_state' : [0],
        #    'class_weight' : ['balanced'],
        #    'n_estimators' : [100, 200, 300],
        #    'max_depth' : [2, 3, 4],
        #    'min_samples_split' : [2, 3],
        #    'max_features' : ["sqrt", "log2"]
        #}),
        ('BalancedRandomForest(imblearn)', BalancedRandomForestClassifier(random_state=0, n_estimators=300, max_depth=4, min_samples_split=2, max_features='sqrt'), 
        {
            'random_state' : [0],
            'n_estimators' : [100, 200, 300],
            'max_depth' : [2, 3, 4],
            'min_samples_split' : [2, 3],
            'max_features' : ["sqrt", "log2"]
        }),
        #('imb_rus', RUSBoostClassifier(random_state=0)), #this one never seems to have a good result on the test set, I think it's overfitting due to boosting
        ('EasyEnsembleClassifier(imblearn)', EasyEnsembleClassifier(random_state=0, n_estimators=50), 
        {
            'random_state' : [0],
            'n_estimators' : [10, 20, 30, 40, 50]
        })
    ]

    #save the labels for use later
    resultsDict['CLF_LABELS'] = [l for l, r, h in classifiers]

    #things to calculate
    aucs = []
    rocs = []
    prs = []
    pr_aucs = []
    resultsDict['CLASSIFIERS'] = {}

    #uncomment to short circuit
    #classifiers = []

    trainedClassifierResults = {
        'FEATURES' : featureLabels,
        'MODELS' : {}
    }

    for clf_label, raw_clf, hyperparam in classifiers:
        #test the classifier?
        print(clf_label)
        print('\ttraining...')

        #prep this for storing any results later
        resultsDict['CLASSIFIERS'][clf_label] = {}

        scoringMode = 'f1'
        cv = StratifiedKFold(n_splits=10)
            
        #scoringMode = 'balanced_accuracy'
        if currentMode == EXACT_MODE:
            clf = raw_clf
        elif currentMode == GRID_MODE:
            #grid search over hyperparameters
            clf = GridSearchCV(raw_clf, hyperparam, cv=cv, scoring=scoringMode, n_jobs=-1)
            #clf = GridSearchCV(raw_clf, hyperparam, cv=10, scoring=make_scorer(average_precision_score), n_jobs=-1)
        elif currentMode == RANDOM_MODE:
            clf = RandomizedSearchCV(raw_clf, hyperparam, cv=cv, scoring='balanced_accuracy', n_jobs=-1)
        else:
            raise Exception("Unexpected classifier mode")
        
        #regardless of above approach, we now fit it
        clf.fit(train_x, train_y)
        if currentMode == EXACT_MODE:
            resultsDict['CLASSIFIERS'][clf_label]['TRAINED_PARAMS'] = clf.get_params()
        elif currentMode in [GRID_MODE, RANDOM_MODE]:
            print('\tBest params:', clf.best_params_)
            resultsDict['CLASSIFIERS'][clf_label]['BEST_PARAMS'] = clf.best_params_
            resultsDict['CLASSIFIERS'][clf_label]['HYPERPARAMETER_SPACE'] = hyperparam

        try:
            print('\tfeature_important:')#, clf.best_estimator_.feature_importances_)
            for j, l in enumerate(featureLabels):
                if currentMode == EXACT_MODE:
                    print('', '', l, clf.feature_importances_[j], sep='\t')
                elif currentMode in [GRID_MODE, RANDOM_MODE]:
                    print('', '', l, clf.best_estimator_.feature_importances_[j], sep='\t')
            
            if currentMode == EXACT_MODE:
                resultsDict['CLASSIFIERS'][clf_label]['FEATURE_IMPORTANCE'] = clf.feature_importances_
            elif currentMode in [GRID_MODE, RANDOM_MODE]:
                resultsDict['CLASSIFIERS'][clf_label]['FEATURE_IMPORTANCE'] = clf.best_estimator_.feature_importances_

        except Exception as e:
            #raise e
            pass
        
        #balanced accuracy - an accuracy score that is an average across the classes
        trainAcc = balanced_accuracy_score(train_y, clf.predict(train_x))
        testAcc = balanced_accuracy_score(test_y, clf.predict(test_x))
        print('\tbalanced_train_acc', trainAcc)
        print('\tbalanced_test_acc', testAcc)
        resultsDict['CLASSIFIERS'][clf_label]['TRAIN_ACCURACY'] = trainAcc
        resultsDict['CLASSIFIERS'][clf_label]['TEST_ACCURACY'] = testAcc
        
        #confusion matrix - exactly what it sounds like, 2x2 grid in our case
        trainConf = confusion_matrix(train_y, clf.predict(train_x))
        testConf = confusion_matrix(test_y, clf.predict(test_x))
        print('\tconf_matrix_train', trainConf)
        print('\tconf_matrix_test', testConf)
        resultsDict['CLASSIFIERS'][clf_label]['TRAIN_CONFUSION_MATRIX'] = trainConf
        resultsDict['CLASSIFIERS'][clf_label]['TRAIN_FPR_RATE'] = trainConf[0, 1] / np.sum(trainConf[0, :])
        resultsDict['CLASSIFIERS'][clf_label]['TRAIN_TPR_RATE'] = trainConf[1, 1] / np.sum(trainConf[1, :])
        resultsDict['CLASSIFIERS'][clf_label]['TEST_CONFUSION_MATRIX'] = testConf
        resultsDict['CLASSIFIERS'][clf_label]['TEST_FPR_RATE'] = testConf[0, 1] / np.sum(testConf[0, :])
        resultsDict['CLASSIFIERS'][clf_label]['TEST_TPR_RATE'] = testConf[1, 1] / np.sum(testConf[1, :])
        
        #roc_curve stuff - could be misleading due to imbalance
        y_pred_rf = clf.predict_proba(test_x)[:, 1]
        false_positive_rate, true_positive_rate, thresholds = roc_curve(test_y, y_pred_rf)
        roc_auc = auc(false_positive_rate, true_positive_rate)
        aucs.append(roc_auc)
        rocs.append((false_positive_rate, true_positive_rate))
        print('\troc_auc', roc_auc)
        resultsDict['CLASSIFIERS'][clf_label]['ROC_AUC'] = roc_auc
        
        #precision-recall curve stuff - should be less misleading
        precision, recall, pr_thresholds = precision_recall_curve(test_y, y_pred_rf)
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)
        prs.append((recall, precision))
        print('\tpr_auc', pr_auc)
        resultsDict['CLASSIFIERS'][clf_label]['PR_AUC'] = pr_auc
        
        if CASE_BASED_SPLIT:
            ranks = []
            rp = []
            
            rDict = {}
            
            #now do the test sets as if they were individual cases
            for x in range(0, len(test_indices)-1):
                st = test_indices[x]
                et = test_indices[x+1]

                case_x = test_x[st:et]
                case_y = test_y[st:et]
                case_d = test_dicts[st:et]

                #get the probabilities and sort them with most likely reported first
                probs = clf.predict_proba(case_x)[:, 1]
                ordered = np.argsort(probs)[::-1]

                for i, v in enumerate(ordered):
                    if case_y[v] == 1:
                        pathLevel = pDict[case_d[v]['path_level']]
                        ranks.append(i+1)
                        rp.append(pathLevel)

                        if (pathLevel not in rDict):
                            rDict[pathLevel] = []
                        rDict[pathLevel].append(i+1)
            
            print('\tNum ranked:', len(ranks))
            print('\tRanks:', ranks)
            print('\tPatho:', rp)
            print('\tmean:', np.mean(ranks))
            print('\tstdev:', np.std(ranks))
            print('\tmedian:', np.median(ranks))
            resultsDict['CLASSIFIERS'][clf_label]['TEST_RANKINGS'] = {}
            resultsDict['CLASSIFIERS'][clf_label]['TEST_RANKINGS']['OVERALL'] = ranks

            for v in pList:
                print('\t'+v)
                print('\t\tRanks:', rDict[pDict[v]])
                print('\t\tmean:', np.mean(rDict[pDict[v]]))
                print('\t\tstdev:', np.std(rDict[pDict[v]]))
                print('\t\tmedian:', np.median(rDict[pDict[v]]))
                resultsDict['CLASSIFIERS'][clf_label]['TEST_RANKINGS'][v] = rDict[pDict[v]]
            
            #this will get written multiple times but that's okay
            resultsDict['TEST_COUNTS'] = {}
            resultsDict['TEST_COUNTS']['OVERALL'] = len(ranks)
            for p in pList:
                resultsDict['TEST_COUNTS'][p] = len(rDict[pDict[p]])


        #cross validation scores
        cv = StratifiedKFold(n_splits=10)
        
        if currentMode == EXACT_MODE:
            cvClf = clf
        elif currentMode == GRID_MODE or currentMode == RANDOM_MODE:
            cvClf = clf.best_estimator_

        #save the model in this dictionary for storage
        trainedClassifierResults['MODELS'][clf_label] = cvClf

        scoringMethods = ['balanced_accuracy', 'f1']
        allScores = cross_validate(cvClf, train_x, train_y, cv=cv, scoring=scoringMethods, n_jobs=-1)
        
        bal_cvs_scores = allScores['test_balanced_accuracy']
        print('\tbal_cvs_scores', bal_cvs_scores)
        print("\tbal_Accuracy: %0.4f (+/- %0.4f)" % (bal_cvs_scores.mean(), bal_cvs_scores.std() * 2))
        resultsDict['CLASSIFIERS'][clf_label]['CV10_BALANCED_ACCURACY'] = (bal_cvs_scores.mean(), bal_cvs_scores.std() * 2)
        
        f1_cvs_scores = allScores['test_f1']
        print('\tf1_cvs_scores', f1_cvs_scores)
        print("\tf1_Accuracy: %0.4f (+/- %0.4f)" % (f1_cvs_scores.mean(), f1_cvs_scores.std() * 2))
        resultsDict['CLASSIFIERS'][clf_label]['CV10_F1'] = (f1_cvs_scores.mean(), f1_cvs_scores.std() * 2)
        
        '''
        #Might be worth doing this in the future, but not particularly useful right now due to data redundancy
        import eli5
        from eli5.sklearn.permutation_importance import PermutationImportance
        permImp = PermutationImportance(clf, scoring='balanced_accuracy', random_state=0, cv='prefit')
        permImp.fit(train_x, train_y)
        explan = eli5.explain_weights(permImp, feature_names=featureLabels)
        print(eli5.formatters.text.format_as_text(explan))
        exit()
        '''
            
        print()
    
    #before making these curves, lets add in some controls for comparison
    comparedScores = ['CADD Scaled', 'HPO-cosine']
    csRocs = []
    csAucs = []
    csPrs = []
    csPrAucs = []
    resultsDict['COMPARISON'] = {}
    resultsDict['CS_LABELS'] = comparedScores
    for csLabel in comparedScores:
        resultsDict['COMPARISON'][csLabel] = {}

        dataIndex = featureLabels.index(csLabel)
        dataColumn = test_x[:, dataIndex]
        
        #these have to be reversed because they are rank-based (i.e. smaller is better)
        if csLabel in ['HPO-cosine', 'PyxisMap']:
            dataColumn = 1-dataColumn

        false_positive_rate, true_positive_rate, thresholds = roc_curve(test_y, dataColumn)
        roc_auc = auc(false_positive_rate, true_positive_rate)
        csAucs.append(roc_auc)
        csRocs.append((false_positive_rate, true_positive_rate))
        resultsDict['COMPARISON'][csLabel]['ROC_AUC'] = roc_auc

        precision, recall, pr_thresholds = precision_recall_curve(test_y, dataColumn)
        pr_auc = auc(recall, precision)
        csPrAucs.append(pr_auc)
        csPrs.append((recall, precision))
        resultsDict ['COMPARISON'][csLabel]['PR_AUC'] = pr_auc
        
        if CASE_BASED_SPLIT:
            ranks = []
            rp = []

            rDict = {}
            
            #now do the test sets as if they were individual cases
            for x in range(0, len(test_indices)-1):
                st = test_indices[x]
                et = test_indices[x+1]

                #dataColumn is sorted such that highest is biggest rank
                case_x = dataColumn[st:et]
                case_y = test_y[st:et]
                case_d = test_dicts[st:et]

                #get the probabilities and sort them with most likely reported first
                ordered = np.argsort(case_x)[::-1]

                for i, v in enumerate(ordered):
                    if case_y[v] == 1:
                        pathLevel = pDict[case_d[v]['path_level']]
                        ranks.append(i+1)
                        rp.append(pathLevel)

                        if (pathLevel not in rDict):
                            rDict[pathLevel] = []
                        rDict[pathLevel].append(i+1)
            print(csLabel+':')
            print('\tNum ranked:', len(ranks))
            print('\tRanks:', ranks)
            print('\tPatho:', rp)
            print('\tmean:', np.mean(ranks))
            print('\tstdev:', np.std(ranks))
            print('\tmedian:', np.median(ranks))
            resultsDict['COMPARISON'][csLabel]['TEST_RANKINGS'] = {}
            resultsDict['COMPARISON'][csLabel]['TEST_RANKINGS']['OVERALL'] = ranks

            for v in pList:
                print('\t'+v)
                print('\t\tRanks:', rDict[pDict[v]])
                print('\t\tmean:', np.mean(rDict[pDict[v]]))
                print('\t\tstdev:', np.std(rDict[pDict[v]]))
                print('\t\tmedian:', np.median(rDict[pDict[v]]))
                resultsDict['COMPARISON'][csLabel]['TEST_RANKINGS'][v] = rDict[pDict[v]]
    
    #now do exomiser stuff
    if CASE_BASED_SPLIT:
        externalDatasets = [
            'Exomiser(hiPhive)', 
            'Exomiser(hiPhive, human only)',
            'Phen-Gen',
            'DeepPVP'
        ]

        resultsDict['EXTERNAL'] = {}
        resultsDict['EXT_LABELS'] = externalDatasets

        for externalLabel in externalDatasets:
            resultsDict['EXTERNAL'][externalLabel] = {}
            ranks = []
            i = 0
            rp = []

            rDict = {}
            missingDict = {}
            missingArray = []
            
            #now do the test sets as if they were individual cases
            for x in range(0, len(test_indices)-1):
                st = test_indices[x]
                et = test_indices[x+1]

                #dataColumn is sorted such that highest is biggest rank
                case_d = test_dicts[st:et]

                #ranks contain the 0-based index, so we need to add one to it for output comparison
                #ranks += [r+1 for r in exomiserArray[x]]
                #missingArray += [(r < 0) for r in exomiserTest[x]]
                missingArray += [(r < 0) for r in externalTest[externalLabel][x]]
                #TODO: remove excess code here
                if True or externalLabel == 'EXOMISER_BEST':
                    #absolute value, then 1-base it
                    #ranks += [abs(v)+1 for v in exomiserTest[x]]
                    ranks += [abs(v)+1 for v in externalTest[externalLabel][x]]
                elif externalLabel == 'EXOMISER_AVERAGE':
                    #average value if missing, then 1-base it
                    #if less than 0, the value is missing SO find the midpoint of unranked values and set the number to that
                    #midpoint calc = (# of variants - #ranked)/2.0 + (rank of highest ranked)
                    #ranks += [v+1 if v >= 0 else ((et-st+v)/2.0-v)+1 for v in exomiserTest[x]]
                    ranks += [v+1 if v >= 0 else ((et-st+v)/2.0-v)+1 for v in externalTest[externalLabel][x]]
                else:
                    raise Exception('Unimplemented exomiser output label: '+externalLabel)

                for d in case_d:
                    if d != None:
                        pathLevel = pDict[d['path_level']]
                        rp.append(pathLevel)
                        
                        if (pathLevel not in rDict):
                            rDict[pathLevel] = []
                        rDict[pathLevel].append(ranks[i])

                        if (pathLevel not in missingDict):
                            missingDict[pathLevel] = []
                        missingDict[pathLevel].append(missingArray[i])

                        i += 1

            ranks = list(filter(lambda v: v != None, ranks))
            missingCounts = sum(missingArray)

            print(externalLabel+':')
            print('\tNum ranked:', len(ranks))
            print('\tMissing:', missingCounts)
            print('\tRanks:', ranks)
            print('\tPatho:', rp)
            print('\tmean:', np.mean(ranks))
            print('\tstdev:', np.std(ranks))
            print('\tmedian:', np.median(ranks))
            resultsDict['EXTERNAL'][externalLabel]['MISSING'] = {}
            resultsDict['EXTERNAL'][externalLabel]['MISSING']['OVERALL'] = missingCounts
            resultsDict['EXTERNAL'][externalLabel]['FOUND'] = {}
            resultsDict['EXTERNAL'][externalLabel]['FOUND']['OVERALL'] = len(ranks)
            resultsDict['EXTERNAL'][externalLabel]['TEST_RANKINGS'] = {}
            resultsDict['EXTERNAL'][externalLabel]['TEST_RANKINGS']['OVERALL'] = ranks

            for v in pList:
                missingCounts = rDict[pDict[v]].count(None)
                rDict[pDict[v]] = list(filter(lambda v: v != None, rDict[pDict[v]]))
                print('\t'+v)
                print('\t\tMissing:', sum(missingDict[pDict[v]]))
                print('\t\tRanks:', rDict[pDict[v]])
                print('\t\tmean:', np.mean(rDict[pDict[v]]))
                print('\t\tstdev:', np.std(rDict[pDict[v]]))
                print('\t\tmedian:', np.median(rDict[pDict[v]]))
                resultsDict['EXTERNAL'][externalLabel]['MISSING'][v] = sum(missingDict[pDict[v]])
                resultsDict['EXTERNAL'][externalLabel]['FOUND'][v] = 'NO_IMPL'
                resultsDict['EXTERNAL'][externalLabel]['TEST_RANKINGS'][v] = rDict[pDict[v]]

    #ROC curve
    if args.path_only:
        plotFN = '/Users/matt/githubProjects/VarSight/paper/codi_rf_roc_pathOnly.png'
    else:
        plotFN = '/Users/matt/githubProjects/VarSight/paper/codi_rf_roc.png'
    plt.figure()
    plt.plot([0, 1], [0, 1], 'k--')
    
    for i, (label2, raw_clf, raw_params) in enumerate(classifiers):
        plt.plot(rocs[i][0], rocs[i][1], label=('%s (%0.4f)' % (label2, aucs[i])))
    
    #Mana pointed out these don't make much sense
    #leaving as comment if a reviewer gets angsty and we want to add back in
    #for i, label2 in enumerate(comparedScores):
    #    plt.plot(csRocs[i][0], csRocs[i][1], label=('%s (%0.4f)' % (label2, csAucs[i])))

    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.legend(loc='best')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid()
    plt.savefig(plotFN)
    plt.close()
    
    #precision recall curve
    if args.path_only:
        plotFN = '/Users/matt/githubProjects/VarSight/paper/codi_rf_pr_pathOnly.png'
    else:
        plotFN = '/Users/matt/githubProjects/VarSight/paper/codi_rf_pr.png'
    plt.figure()
    
    for i, (label2, raw_clf, raw_params) in enumerate(classifiers):
        plt.plot(prs[i][0], prs[i][1], label=('%s (%0.4f)' % (label2, pr_aucs[i])))
    
    #Mana pointed out these don't make much sense
    #leaving as comment if a reviewer gets angsty and we want to add back in
    #for i, label2 in enumerate(comparedScores):
    #    plt.plot(csPrs[i][0], csPrs[i][1], label=('%s (%0.4f)' % (label2, csPrAucs[i])))
    
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='best')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid()
    plt.savefig(plotFN)
    plt.close()

    return (resultsDict, trainedClassifierResults)

def jsonDumpFix(o):
    '''
    Patch for doing a JSON dump when numpy int64s are present
    see https://stackoverflow.com/questions/11942364/typeerror-integer-is-not-json-serializable-when-serializing-json-in-python
    '''
    if isinstance(o, np.int64):
        return int(o)
    elif isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(str(o)+' '+str(type(o)))

def generateLaTeXResult(args, d):
    '''
    This function takes a dictionary of values from our results and makes the calls to render those results through
    Jinja2 into LaTeX files for the paper.
    @param args - the command line arguments
    @param d - the data dictionary to hand to any templating
    '''
    from jinja2 import Environment, FileSystemLoader
    rootPath = '/Users/matt/githubProjects/VarSight/paper/'
    env = Environment(loader=FileSystemLoader(rootPath))
    
    if args.path_only:
        dataFN = rootPath+'/data_rendered_pathOnly.tex'
        data1FN = rootPath+'/data_rank_rendered_pathOnly.tex'
        data2FN = rootPath+'/data_top_rendered_pathOnly.tex'
        classifierFN = rootPath+'/classifier_rendered_pathOnly.tex'
        featuresFN = rootPath+'/features_rendered_pathOnly.tex'
        boxPlotFN = rootPath+'/supplements/results_boxplot_pathOnly.png'
    else:
        dataFN = rootPath+'/data_rendered.tex'
        data1FN = rootPath+'/data_rank_rendered.tex'
        data2FN = rootPath+'/data_top_rendered.tex'
        classifierFN = rootPath+'/classifier_rendered.tex'
        featuresFN = rootPath+'/features_rendered.tex'
        boxPlotFN = rootPath+'/supplements/results_boxplot.png'

    #this is independent of args
    featureDescFN = rootPath+'/supplements/feature_info_rendered.tex'
    hyperFN = rootPath+'/supplements/hyperparameter_rendered.tex'

    #this is the main data template
    template = env.get_template('data_template.tex')
    rendered = template.render(d)
    fp = open(dataFN, 'wt+')
    fp.write(rendered)
    fp.close()

    #sub-templates for smaller figures
    template = env.get_template('data_rank_template.tex')
    rendered = template.render(d)
    fp = open(data1FN, 'wt+')
    fp.write(rendered)
    fp.close()
    template = env.get_template('data_top_template.tex')
    rendered = template.render(d)
    fp = open(data2FN, 'wt+')
    fp.write(rendered)
    fp.close()

    #render the classifier only template
    template = env.get_template('classifier_template.tex')
    rendered = template.render(d)
    fp = open(classifierFN, 'wt+')
    fp.write(rendered)
    fp.close()

    #feature importance table - we need to parse our features
    orderedFeatures = []
    classWithFeat = ['RandomForest(sklearn)', 'BalancedRandomForest(imblearn)']
    shortClassNames = ['RF(sklearn)', 'BRF(imblearn)']
    CUTOFF_IMPORTANCE = 0.00
    assert(len(classWithFeat) == len(shortClassNames))
    for i, label in enumerate(d['FEATURE_LABELS']):
        orderedFeatures.append(np.mean([d['CLASSIFIERS'][c]['FEATURE_IMPORTANCE'][i] for c in classWithFeat]))
    
    argSorted = np.argsort(orderedFeatures)[::-1]
    featDict = {}
    featDict['FEATURE_CLASSIFIERS'] = shortClassNames
    featDict['TABLE'] = []
    for v in argSorted:
        if orderedFeatures[v] >= CUTOFF_IMPORTANCE:
            #print(d['FEATURE_LABELS'][v], orderedFeatures[v], [d['CLASSIFIERS'][c]['FEATURE_IMPORTANCE'][v] for c in classWithFeat], sep='\t')
            featDict['TABLE'].append([d['FEATURE_LABELS'][v]]+[d['CLASSIFIERS'][c]['FEATURE_IMPORTANCE'][v] for c in classWithFeat])
    featDict['TABLE_SUM'] = [
        np.sum([t[1] for t in featDict['TABLE']]),
        np.sum([t[2] for t in featDict['TABLE']])
    ]
    template = env.get_template('features_template.tex')
    rendered = template.render(featDict)
    fp = open(featuresFN, 'wt+')
    fp.write(rendered)
    fp.close()

    #feature description table - supplements
    metaFN = '/Users/matt/githubProjects/VarSight/CODI_metadata/fields_metadata.json'
    fp = open(metaFN, 'rt')
    catMeta = json.load(fp)
    fp.close()
    annotKeys = {
        'sequencingFields' : 'variant',
        'allelicInformation' : 'variant',
        'genes' : 'gene',
        'variantTranscripts' : 'transcript'
    }
    tempDict = {
        'DATA' : catMeta,
        'KEYS' : annotKeys
    }
    template = env.get_template('supplements/feature_info_template.tex')
    rendered = template.render(tempDict)
    fp = open(featureDescFN, 'wt+')
    fp.write(rendered)
    fp.close()

    #hyperparameter table - supplements
    d['len'] = len
    template = env.get_template('supplements/hyperparameter_template.tex')
    rendered = template.render(d)
    fp = open(hyperFN, 'wt+')
    fp.write(rendered)
    fp.close()

    #boxplots of results
    generateBoxplot(d, boxPlotFN)

    if args.path_only:
        #do nothing here
        pass
    else:
        #create a combined image
        #from: https://stackoverflow.com/questions/30227466/combine-several-images-horizontally-with-python
        from PIL import Image
        images = [Image.open(fn) for fn in ['/Users/matt/githubProjects/VarSight/paper/codi_rf_roc.png', '/Users/matt/githubProjects/VarSight/paper/codi_rf_pr.png']]
        widths, heights = zip(*(i.size for i in images))
        totalWidth = sum(widths)
        maxHeight = max(heights)
        new_im = Image.new('RGB', (totalWidth, maxHeight))
        x_offset = 0
        for im in images:
            new_im.paste(im, (x_offset, 0))
            x_offset += im.size[0]
        new_im.save('/Users/matt/githubProjects/VarSight/paper/codi_rf_combined.png')

def generateBoxplot(d, boxPlotFN):
    '''
    Generates a boxplot summary of the results that we will include in the supplements
    @param d - the dictionary of results, mostly loaded from JSON
    @param boxPlotFN - the output filename for the boxplot figure
    '''
    plt.figure()
    labels = []
    dataPoints = []

    for csLabel in d['CS_LABELS']:
        labels.append(csLabel)
        dataPoints.append(d['COMPARISON'][csLabel]['TEST_RANKINGS']['OVERALL'])
    for extLabel in d['EXT_LABELS']:
        labels.append(extLabel)
        dataPoints.append(d['EXTERNAL'][extLabel]['TEST_RANKINGS']['OVERALL'])
    for clfLabel in d['CLF_LABELS']:
        labels.append(clfLabel)
        dataPoints.append(d['CLASSIFIERS'][clfLabel]['TEST_RANKINGS']['OVERALL'])

    for i, dv in enumerate(dataPoints):
        jitter = np.random.normal(i+1, 0.02, size=len(dv))
        plt.plot(jitter, dv, 'r.', alpha=0.2)
    plt.boxplot(dataPoints)
    
    plt.xticks(range(1, len(dataPoints)+1), labels, rotation=90)
    plt.xlabel('Ranking Method')
    plt.ylabel('Rank (lower is better)')
    plt.grid()
    plt.savefig(boxPlotFN, bbox_inches="tight")

def runAnalysis(args):
    '''
    Core function for joining our pieces together
    '''
    if args.path_only:
        resultJsonFN = '/Users/matt/githubProjects/VarSight/paper/results_path_only.json'
        modelPickleFN = '/Users/matt/githubProjects/VarSight/models/VarSight_path_only.p'
    else:
        resultJsonFN = '/Users/matt/githubProjects/VarSight/paper/results.json'
        modelPickleFN = '/Users/matt/githubProjects/VarSight/models/VarSight.p'
    REGENERATE_DATA = args.regenerate

    if REGENERATE_DATA or not os.path.exists(resultJsonFN):
        (xFinal, yFinal, catMeta, featureLabels, startIndices, allRepDicts, externalRanks) = loadFormattedData(args)
        resultsDict, trainedModelDict = runClassifiers(args, xFinal, yFinal, featureLabels, startIndices, allRepDicts, externalRanks)
        trainedModelDict['CATEGORICAL_METADATA'] = catMeta

        fp = open(resultJsonFN, 'w+')
        json.dump(resultsDict, fp, default=jsonDumpFix, sort_keys=True, indent=4)
        fp.close()

        fp = open(modelPickleFN, 'wb+')
        pickle.dump(trainedModelDict, fp)
        fp.close()

    fp = open(resultJsonFN, 'rt')
    resultsDict = json.load(fp)
    fp.close()

    if args.path_only:
        print('PATH_ONLY, not integrated into main paper')
        resultsDict['np'] = np
        resultsDict['len'] = len
        #resultsDict['LEVELS'] = ['OVERALL', 'LIKELY_PATHOGENIC', 'PATHOGENIC']
        resultsDict['LEVELS'] = ['OVERALL', 'VARIANT_OF_UNCERTAIN_SIGNIFICANCE', 'LIKELY_PATHOGENIC', 'PATHOGENIC']
        generateLaTeXResult(args, resultsDict)
        
    else:
        print(json.dumps(resultsDict, indent=4, sort_keys=True))
        resultsDict['np'] = np
        resultsDict['len'] = len
        resultsDict['LEVELS'] = ['OVERALL', 'VARIANT_OF_UNCERTAIN_SIGNIFICANCE', 'LIKELY_PATHOGENIC', 'PATHOGENIC']
        generateLaTeXResult(args, resultsDict)

if __name__ == '__main__':
    p = ap.ArgumentParser(description='script for running analysis for VarSight', formatter_class=ap.RawTextHelpFormatter)

    sp = p.add_subparsers(dest='subparser')
    sp1 = sp.add_parser('gather', help='gather HPO-based results for use later')

    sp2 = sp.add_parser('analyze', help='train and test the classifiers')
    sp2.add_argument('-R', '--regenerate', dest='regenerate', action='store_true', default=False, help='regenerate all test results even if already available')
    sp2.add_argument('-P', '--path-only', dest='path_only', action='store_true', default=False, help='only uses pathogenic and likely pathogenic variants as true positives')
    ex_group = sp2.add_mutually_exclusive_group()
    ex_group.add_argument('-e', '--exact-mode', dest='training_mode', action='store_const', const=EXACT_MODE, help='use the dev-specified hyperparameters (single-execution)', default=EXACT_MODE)
    ex_group.add_argument('-g', '--grid-mode', dest='training_mode', action='store_const', const=GRID_MODE, help='perform a grid search for the best hyperparameters (long multi-execution)', default=EXACT_MODE)
    ex_group.add_argument('-r', '--random-mode', dest='training_mode', action='store_const', const=RANDOM_MODE, help='perform a random search for the best hyperparameters (short multi-execution)', default=EXACT_MODE)
    
    args = p.parse_args()
    
    if args.subparser == 'gather':
        #this will re-run the PyxisMap tests, should update dates though
        getPyxisMapResults()
    elif args.subparser == 'analyze':
        #core tests
        runAnalysis(args)
