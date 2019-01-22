'''
Primary analysis script for testing various learning methods
'''

import argparse as ap
import datetime
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import re

import CodiDumpUtil
import ExomiserUtil
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
    @return - tuple (xFinal, yFinal, featureLabels, startIndices, allRepDicts, exomiserRanks)
        xFinal - an NxM matrix where N is the number of variants in our test/training set and M is the number of features; contains feature values
        yFinal - an N length array where N is the number of variants in our test/training set; 1 if the variant was reported
        featureLabels - an M length array containing feature labels for our output benefit
        startIndices - a (C+1) length array containing the start indices of variants from an individual case
        allRepDicts - an N length array of None or dictionaries; if a dictionary, it contains data on a returned variant
        exomiserRanks - a list of lists containing exomiser ranks on a per-case basis; missing values are -1*len(ranked variants)
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
    
    seqFieldsOnly = set([l for l, d in seqValues]) | set(['chromosome', 'position', 'ref allele', 'alt allele'])
    
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
        codiDump = CodiDumpUtil.loadCodiDump(sl, fieldsOnly, geneFieldsOnly, seqFieldsOnly, transLabelsOnly)
        if len(codiDump) == 0:
            print(sl, 'NO_CODICEM_DUMP', sep=',')
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
        
        #now, we have all the pieces, put them together into tests
        values = []
        catValues = []
        catBreak = {}
        classifications = []
        repDicts = []
        
        #print('Primary set: ', primarySet)
        foundPrimaries = set([])
        fpCPRA = []
        
        #go through each variant now and format the inputs to be seen by a classifier
        for variant in codiDump:
            vData = []
            catData = []
            
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
            
            #first, add in the best PyxisMap rank
            geneList = variant['genes']
            pyxisBest = min([pyxisRanks.get(g, (pyxisLen+1, 0.0))[0] / pyxisLen for g in geneList])
            vData.append(pyxisBest)
            
            #now, add the best HPOUtil rank
            hpoBest = min([hpoRanks.get(g, (hpoLen+1, 0.0))[0] / hpoLen for g in geneList])
            vData.append(hpoBest)
            
            #now add all float values from CODICEM
            for fvl, invalidValue in floatValues:
                if variant[fvl] == 'NA':
                    vData.append(invalidValue)
                else:
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

            #add categorical information as described by "fields_metadata.json"
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
            
            #we add all values to our lists for the case
            values.append(vData)
            catValues.append(catData)
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
            exRanks = ExomiserUtil.getTargetRanks(sl, fpCPRA)
            exomiserRanks.append(exRanks)

        else:
            print(sl, 'ALL_MISSING', primarySet, foundPrimaries, sep=',')
        
        #print()
    
    #determine the category mode
    OHE_MODE = 0 #categories are given *-hot encodings where * is the count of how many times that label appears
    PCA_MODE = 1 #all categories are PCA-ed together and the top X PCA components are used
    PCA_BREAK_MODE = 2 #each category has its own PCA, only X PCA components are allowed per category
    currentCatMode = PCA_BREAK_MODE

    featureLabels = (['PyxisMap', 'HPO-cosine']+
        [l for l, d in floatValues]+
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
                else:
                    raise Exception('Unexpected fieldType')
        xFinal = np.hstack([vVals]+pcaBlocks)
    else:
        raise Exception('Unexpected currentCatMode')

    #return the values
    yFinal = np.array(np.hstack(allClassifications), dtype='int64')

    return (xFinal, yFinal, featureLabels, startIndices, allRepDicts, exomiserRanks)

def runClassifiers(args, values, classifications, featureLabels, startIndices, allRepDicts, exomiserRanks):
    '''
    @param args - any arguments from the command line argparse can be accessed here
    @param values - a matrix with R rows and C columns, where there are "C" features
    @param classifications - an array of length R corresponding to the above values
    @param featureLabels - C length array containing labels (strings) for the features
    @param startIndices - the startIndices of individual cases in the training set
    @param allRepDicts - a list of dictionaries for variants that were reported; non-reported vars are None
    @param exomiserRanks - paired exomiser ranks for each case; needs to be split up with train/test; missing values are -1*len(ranked variants)
    @return resultsDict - a dictionary containing many results we wish to include in a paper
    '''
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
        exomiserTrain = []

        testXArray = []
        testYArray = []
        test_indices = [0]
        test_dicts = []
        tpTest = 0
        exomiserTest = []

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
                exomiserTrain.append(exomiserRanks[x])

            else:
                #we need more in our test
                testXArray.append(values[st:et])
                testYArray.append(classifications[st:et])
                test_indices.append(test_indices[-1]+(et-st))
                test_dicts += allRepDicts[st:et]
                tpTest += tpCount

                #we moved the data manipulation down
                exomiserTest.append(exomiserRanks[x])
        
        #join them all together
        train_x = np.vstack(trainXArray)
        train_y = np.hstack(trainYArray)
        test_x = np.vstack(testXArray)
        test_y = np.hstack(testYArray)

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
        ('RandomForest(sklearn)', RandomForestClassifier(random_state=0, class_weight='balanced', max_depth=4, n_estimators=300, min_samples_split=2, max_features='sqrt'), 
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
        ('LogisticRegression(sklearn)', LogisticRegression(random_state=0, class_weight='balanced', penalty='l2', C=0.01, solver='liblinear', max_iter=200), 
        {
            'random_state' : [0],
            'class_weight' : ['balanced'],
            'penalty' : ['l2'],
            'C' : [0.01, 0.1, 1.0, 10.0, 100.0],
            'solver' : ['newton-cg', 'liblinear'],
            'max_iter' : [200]
        }),
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

    for clf_label, raw_clf, hyperparam in classifiers:
        #test the classifier?
        print(clf_label)
        print('\ttraining...')

        #prep this for storing any results later
        resultsDict['CLASSIFIERS'][clf_label] = {}

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
        if currentMode == EXACT_MODE:
            resultsDict['CLASSIFIERS'][clf_label]['TRAINED_PARAMS'] = clf.get_params()
        elif currentMode in [GRID_MODE, RANDOM_MODE]:
            print('\tBest params:', clf.best_params_)
            resultsDict['CLASSIFIERS'][clf_label]['BEST_PARAMS'] = clf.best_params_

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
        resultsDict['CLASSIFIERS'][clf_label]['TRAIN_FPR_RATE'] = trainConf[0, 0] / np.sum(trainConf[0, :])
        resultsDict['CLASSIFIERS'][clf_label]['TRAIN_TPR_RATE'] = trainConf[1, 1] / np.sum(trainConf[1, :])
        resultsDict['CLASSIFIERS'][clf_label]['TEST_CONFUSION_MATRIX'] = testConf
        resultsDict['CLASSIFIERS'][clf_label]['TEST_FPR_RATE'] = testConf[0, 0] / np.sum(testConf[0, :])
        resultsDict['CLASSIFIERS'][clf_label]['TEST_TPR_RATE'] = testConf[1, 1] / np.sum(testConf[1, :])
        
        #roc_curve stuff - could be misleading due to imbalance
        y_pred_rf = clf.predict_proba(test_x)[:, 1]
        false_positive_rate, true_positive_rate, thresholds = roc_curve(test_y, y_pred_rf)
        roc_auc = auc(false_positive_rate, true_positive_rate)
        aucs.append(roc_auc)
        rocs.append((false_positive_rate, true_positive_rate))
        print('\troc_auc', roc_auc)
        #resultsDict['CLASSIFIERS'][clf_label]['FALSE_POSITIVE_RATE_CURVES'] = false_positive_rate
        #resultsDict['CLASSIFIERS'][clf_label]['TRUE_POSITIVE_RATE_CURVES'] = true_positive_rate
        resultsDict['CLASSIFIERS'][clf_label]['ROC_AUC'] = roc_auc
        
        #precision-recall curve stuff - should be less misleading
        precision, recall, pr_thresholds = precision_recall_curve(test_y, y_pred_rf)
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)
        prs.append((recall, precision))
        print('\tpr_auc', pr_auc)
        #resultsDict['CLASSIFIERS'][clf_label]['PRECISION_CURVES'] = precision
        #resultsDict['CLASSIFIERS'][clf_label]['RECALL_CURVES'] = recall
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

        if currentMode == EXACT_MODE:
            #cross validation scores
            cv = StratifiedKFold(n_splits=10)
            
            #balanced accuracy - "The balanced accuracy in binary and multiclass classification problems to deal with imbalanced datasets. It is defined as the average of recall obtained on each class."
            bal_cvs_scores = cross_val_score(clf, train_x, train_y, cv=cv, scoring='balanced_accuracy')
            print('\tbal_cvs_scores', bal_cvs_scores)
            print("\tbal_Accuracy: %0.4f (+/- %0.4f)" % (bal_cvs_scores.mean(), bal_cvs_scores.std() * 2))
            resultsDict['CLASSIFIERS'][clf_label]['CV10_BALANCED_ACCURACY'] = (bal_cvs_scores.mean(), bal_cvs_scores.std() * 2)
            
            #TODO: I think we should be using F1 here, only matters if we ever go back to EXACT_MODE for reporting results
            #f1_weighted - see https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html#sklearn.metrics.f1_score
            f1w_cvs_scores = cross_val_score(clf, train_x, train_y, cv=cv, scoring='f1_weighted')
            print('\tf1w_cvs_scores', f1w_cvs_scores)
            print("\tf1w_Accuracy: %0.4f (+/- %0.4f)" % (f1w_cvs_scores.mean(), f1w_cvs_scores.std() * 2))
            resultsDict['CLASSIFIERS'][clf_label]['CV10_BALANCED_F1_WEIGHTED'] = (f1w_cvs_scores.mean(), f1w_cvs_scores.std() * 2)
            
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

        elif currentMode in [GRID_MODE, RANDOM_MODE]:
            #CV is already calculated here so we can just dump it out
            print('\tscoring systems:')
            means = clf.cv_results_['mean_test_score']
            stds = clf.cv_results_['std_test_score']
            for mean, std, params in zip(means, stds, clf.cv_results_['params']):
                print("\t\t%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))
                if params == clf.best_params_:
                    resultsDict['CLASSIFIERS'][clf_label]['CV10_BALANCED_ACCURACY'] = (mean, std*2)
            
            #TODO: do we want anything from these results or do we care?
            
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
        exomiserDatasets = [
            'EXOMISER_BEST',
            'EXOMISER_AVERAGE'
        ]

        for exomiserLabel in exomiserDatasets:
            resultsDict[exomiserLabel] = {}
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
                missingArray += [(r < 0) for r in exomiserTest[x]]
                if exomiserLabel == 'EXOMISER_BEST':
                    #absolute value, then 1-base it
                    ranks += [abs(v)+1 for v in exomiserTest[x]]
                elif exomiserLabel == 'EXOMISER_AVERAGE':
                    #average value if missing, then 1-base it
                    #if less than 0, the value is missing SO find the midpoint of unranked values and set the number to that
                    #midpoint calc = (# of variants - #ranked)/2.0 + (rank of highest ranked)
                    ranks += [v+1 if v >= 0 else ((et-st+v)/2.0-v)+1 for v in exomiserTest[x]]
                else:
                    raise Exception('Unimplemented exomiser output label: '+exomiserLabel)

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

            print(exomiserLabel+':')
            print('\tNum ranked:', len(ranks))
            print('\tMissing:', missingCounts)
            print('\tRanks:', ranks)
            print('\tPatho:', rp)
            print('\tmean:', np.mean(ranks))
            print('\tstdev:', np.std(ranks))
            print('\tmedian:', np.median(ranks))
            resultsDict[exomiserLabel]['MISSING'] = {}
            resultsDict[exomiserLabel]['MISSING']['OVERALL'] = missingCounts
            resultsDict[exomiserLabel]['FOUND'] = {}
            resultsDict[exomiserLabel]['FOUND']['OVERALL'] = len(ranks)
            resultsDict[exomiserLabel]['TEST_RANKINGS'] = {}
            resultsDict[exomiserLabel]['TEST_RANKINGS']['OVERALL'] = ranks

            for v in pList:
                missingCounts = rDict[pDict[v]].count(None)
                rDict[pDict[v]] = list(filter(lambda v: v != None, rDict[pDict[v]]))
                print('\t'+v)
                print('\t\tMissing:', sum(missingDict[pDict[v]]))
                print('\t\tRanks:', rDict[pDict[v]])
                print('\t\tmean:', np.mean(rDict[pDict[v]]))
                print('\t\tstdev:', np.std(rDict[pDict[v]]))
                print('\t\tmedian:', np.median(rDict[pDict[v]]))
                resultsDict[exomiserLabel]['MISSING'][v] = sum(missingDict[pDict[v]])
                resultsDict[exomiserLabel]['FOUND'][v] = 'NO_IMPL'
                resultsDict[exomiserLabel]['TEST_RANKINGS'][v] = rDict[pDict[v]]

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

    return resultsDict

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
        classifierFN = rootPath+'/classifier_rendered_pathOnly.tex'
        featuresFN = rootPath+'/features_rendered_pathOnly.tex'
    else:
        dataFN = rootPath+'/data_rendered.tex'
        classifierFN = rootPath+'/classifier_rendered.tex'
        featuresFN = rootPath+'/features_rendered.tex'

    #this is the main data template
    template = env.get_template('data_template.tex')
    rendered = template.render(d)
    fp = open(dataFN, 'wt+')
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
    CUTOFF_IMPORTANCE = 0.02
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

def generateViolinPlots(d):
    '''
    This function takes the results and makes some ranked violion plots for us
    TODO: doesn't seem particularly helpful right now, the plot doesn't visually stratify enough; anything we can do about that?
    '''
    #create a violin plot stacking all the types
    plotFN = '/Users/matt/githubProjects/VarSight/paper/violin_ranks.png'
    plt.figure()
    
    overall = []
    vus = []
    lp = []
    p = []
    for clf_label in d['CLF_LABELS']:
        overall.append(d['CLASSIFIERS'][clf_label]['TEST_RANKINGS']['OVERALL'])
        vus.append(d['CLASSIFIERS'][clf_label]['TEST_RANKINGS']['VARIANT_OF_UNCERTAIN_SIGNIFICANCE'])
        lp.append(d['CLASSIFIERS'][clf_label]['TEST_RANKINGS']['LIKELY_PATHOGENIC'])
        p.append(d['CLASSIFIERS'][clf_label]['TEST_RANKINGS']['PATHOGENIC'])
    
    for cs_label in d['CS_LABELS']:
        overall.append(d['COMPARISON'][cs_label]['TEST_RANKINGS']['OVERALL'])
        vus.append(d['COMPARISON'][cs_label]['TEST_RANKINGS']['VARIANT_OF_UNCERTAIN_SIGNIFICANCE'])
        lp.append(d['COMPARISON'][cs_label]['TEST_RANKINGS']['LIKELY_PATHOGENIC'])
        p.append(d['COMPARISON'][cs_label]['TEST_RANKINGS']['PATHOGENIC'])

    #for i, (label2, raw_clf, raw_params) in enumerate(classifiers):
    #    plt.plot(prs[i][0], prs[i][1], label=('%s (%0.4f)' % (label2, pr_aucs[i])))
    
    #for i, label2 in enumerate(comparedScores):
    #    plt.plot(csPrs[i][0], csPrs[i][1], label=('%s (%0.4f)' % (label2, csPrAucs[i])))
    
    #plt.xlabel('Recall')
    #plt.ylabel('Precision')
    #plt.title('PR curve - all data')
    plt.violinplot(overall)
    #plt.violinplot(vus)
    #plt.violinplot(lp)
    #plt.violinplot(p)
    #plt.legend(loc='best')
    #plt.yscale('log')
    plt.grid()
    plt.savefig(plotFN)
    plt.close()

def runAnalysis(args):
    '''
    Core function for joining our pieces together
    '''
    if args.path_only:
        resultJsonFN = '/Users/matt/githubProjects/VarSight/paper/results_path_only.json'
    else:
        resultJsonFN = '/Users/matt/githubProjects/VarSight/paper/results.json'
    REGENERATE_DATA = args.regenerate

    if REGENERATE_DATA or not os.path.exists(resultJsonFN):
        (xFinal, yFinal, featureLabels, startIndices, allRepDicts, exomiserRanks) = loadFormattedData(args)
        resultsDict = runClassifiers(args, xFinal, yFinal, featureLabels, startIndices, allRepDicts, exomiserRanks)
        
        fp = open(resultJsonFN, 'wt+')
        json.dump(resultsDict, fp, default=jsonDumpFix)
        fp.close()

    fp = open(resultJsonFN, 'rt')
    resultsDict = json.load(fp)
    fp.close()

    if args.path_only:
        print('PATH_ONLY, not integrated into main paper')
        resultsDict['np'] = np
        resultsDict['len'] = len
        resultsDict['LEVELS'] = ['OVERALL', 'LIKELY_PATHOGENIC', 'PATHOGENIC']
        generateLaTeXResult(args, resultsDict)
        
    else:
        print(json.dumps(resultsDict, indent=4, sort_keys=True))
        resultsDict['np'] = np
        resultsDict['len'] = len
        resultsDict['LEVELS'] = ['OVERALL', 'VARIANT_OF_UNCERTAIN_SIGNIFICANCE', 'LIKELY_PATHOGENIC', 'PATHOGENIC']
        generateLaTeXResult(args, resultsDict)
        generateViolinPlots(resultsDict)

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
