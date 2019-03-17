
import argparse as ap
import numpy as np
import os
import pickle

import CodiDumpUtil

def loadClassifier(args):
    '''
    @param args - the argparse values
    @return - tuple
        classifier - the trained model from the pickle file
        catMeta - the categorical metadata for all possible features
        features - the features as ordered in the model
    '''
    #pull out any args we care about
    modelsFN = args.modelsFN
    modelName = "LogisticRegression(sklearn)" #TODO: make a parameter
    
    #load the pickle file
    fp = open(modelsFN, 'rb')
    modelDict = pickle.load(fp)
    fp.close()

    allowedModels = modelDict['MODELS'].keys()
    if modelName not in allowedModels:
        raise Exception('No model named "%s" in "%s"' % (modelName, modelsFN))
    classifier = modelDict['MODELS'][modelName]
    features = modelDict['FEATURES']
    catMeta = modelDict['CATEGORICAL_METADATA']

    return (classifier, catMeta, features)

def loadHpoTerms(fn):
    '''
    @param fn - a text file with one HPO term per line in HP:####### format
    @return - a set containing all HPO terms in the file
    '''
    hpoTerms = set([])
    fp = open(fn, 'r')
    for l in fp:
        t = l.rstrip()
        hpoTerms.add(t)
    fp.close()
    return hpoTerms

def extractRelevantFeatures(xValuesFull, featureLabelsFull, features):
    '''
    @param xValuesFull - the full list of features
    @param featureLabelsFull - the labels for each column
    @param features - the features we want to extract
    @return - an array containing only the features we care about (aka, fewer columns)
    '''
    columnsKeep = []
    for f in features:
        columnsKeep.append(featureLabelsFull.index(f))
    
    colArr = np.array(columnsKeep)
    return xValuesFull[:, colArr]

def runPredictions(classifier, xFeatureVals):
    '''
    @param classifier - the classifier to run
    @param xFeatureVals - the feature values to test on
    '''
    probs = classifier.predict_proba(xFeatureVals)[:, 1]
    ordered = np.argsort(probs)[::-1]
    return (probs, ordered)

if __name__ == "__main__":
    #first set up the arg parser
    DESC = "VarSight - prioritizing Codicem filtered variants using trained classifiers"
    p = ap.ArgumentParser(description=DESC, formatter_class=ap.RawTextHelpFormatter)
    
    #optional arguments with default
    #p.add_argument('-d', metavar='depth', dest='depth', type=int, default=8, help='minimum read depth to consider a variant (default: 8)')
    #p.add_argument('-q', metavar='quality', dest='quality', type=int, default=0, help='minimum quality to consider a variant (default: 0)')
    #TODO: URL to PyxisMap as an option (with a default to our public instance)
    #TODO: model name should have a default and the option to change
    #p.add_argument('modelName', type=str, help='the model name to use in the pickle file')

    #required main arguments
    p.add_argument('codicemJson', type=str, help='a Codicem JSON containing pre-filtered variants')
    p.add_argument('hpoFN', type=str, help='a file containing one HPO term per line')
    p.add_argument('modelsFN', type=str, help='the pickle file containing trained models')
    
    #parse the arguments
    args = p.parse_args()

    #load the model and parameters
    (classifier, catMeta, features) = loadClassifier(args)
    
    #load the hpo terms
    hpoTerms = loadHpoTerms(args.hpoFN)

    #load the codi dump since it's local
    variantInfo, xValuesFull, featureLabelsFull = CodiDumpUtil.prepDataFromHpo(catMeta, args.codicemJson, hpoTerms)

    #reduce the data down to the features we care about
    xFeatureVals = extractRelevantFeatures(xValuesFull, featureLabelsFull, features)

    #run the classifier
    probs, ordered = runPredictions(classifier, xFeatureVals)

    top20 = []
    for i, o in enumerate(ordered):
        print(i, probs[o], variantInfo[o])

        if i < 20:
            top20 += variantInfo[o][1]
    print()
    
    print('List of genes for top 20 variants: '+';'.join(top20))