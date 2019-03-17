'''
This file should contain helper functions to parse and load CODICEM filter files
'''

import json
import numpy as np
import os
import re

import HPOUtil
import PyxisMapUtil

def loadCodiDump(codiDumpFN, EXTRA_FIELDS, GENE_FIELDS, SEQ_FIELDS, TRANSCRIPT_FIELDS):
    '''
    This file should open the latest CODI dump file from a pre-determined directory and load particular fields to be returned
    @return - a list of variants:
    [
        {
            "variant" : chrX:123184969A>T,
            "genes" : ["BMPR2"],
            "extraField" : "value",
            ...
        }
    ]
    '''
    #make sure the file exists
    #fn = '/Users/matt/githubProjects/VarSight/CLI_primary_V6/'+sl+'_results.json'
    if not os.path.exists(codiDumpFN):
        #print('Missing file: '+fn)
        return []
    
    #load the json file
    fp = open(codiDumpFN, 'rt')
    j = json.load(fp)
    fp.close()
    
    #parse it and store in the appropriate format
    ret = []
    for variant in j:
        '''
        #this dumps key-values if you need it
        for k in sorted(variant.keys()):
            print(k, variant[k])
        exit()
        '''
        #fill this is in with appropriate values
        vDict = {}
        
        #the higher this number is, the bigger it is in the circle pack
        #variantOrderedRank = variant['variantOrderedRank']
        
        #change this to the format we recognize from the table
        vd = variant['variantDisplayName']
        pieces = vd.split(' ')
        vdMod = pieces[0].lower()+pieces[1]+':'+pieces[2]
        vDict['variant'] = vdMod
        
        #add in all the gene symbols to the set
        #also go through fields under 'genes' during this
        genes = set([])
        for g in variant['genes']:
            #get gene centered fields
            for e in g['info']:
                if e['label'] == 'gene symbol':
                    genes.add(e['data'])
                elif (e['label'] in GENE_FIELDS):
                    if (e['label'] in vDict):
                        vDict[e['label']].append(e['data'])
                    else:
                        vDict[e['label']] = [e['data']]
            
            #go through each transcript for the gene while we're here
            for tDict in g['transcripts']:
                for vtDict in tDict['variantTranscripts']:
                    for e in vtDict['info']:
                        if (e['label'] in TRANSCRIPT_FIELDS):
                            if (e['label'] in vDict):
                                vDict[e['label']].append(e['data'])
                            else:
                                vDict[e['label']] = [e['data']]
        
        vDict['genes'] = list(genes)
        
        #go through all fields under 'allelicInformation'
        for ai in variant['allelicInformation']:
            if ai['label'] in EXTRA_FIELDS:
                assert(ai['label'] not in vDict)
                vDict[ai['label']] = ai['data']
            
            #if ai['label'] == 'segDup fracMatchIndel':
            #    print('segDup fracMatchIndel', ai['data'])
            
        #go through all fields under 'seq_info'
        for s in variant['seq_info']:
            if s['label'] in SEQ_FIELDS:
                assert(s['label'] not in vDict)
                vDict[s['label']] = s['data']

        ret.append(vDict)
        
    return ret

def prepDataFromHpo(catMeta, codiDumpFN, hpoTerms):
    '''
    This method is intended to be used by VarSight.py. As such, this is basically a wrapper around prepDataShared(...)
    such that code isn't duplicated all over the place.
    @param catMeta - the category metadata to be passed on
    @param codiDumpFN - the Codicem JSON dump
    @param hpoTerms - the set of hpo terms used to query
    '''
    #query PyxisMap
    pyxisRanks, pyxisLen, response = PyxisMapUtil.getPyxisMapResults(hpoTerms, 'http://0.0.0.0:5000')

    #query HPO
    ontFN = '/Users/matt/githubProjects/LayeredGraph/HPO_graph_data/hp.obo'
    annotFN = '/Users/matt/githubProjects/LayeredGraph/HPO_graph_data/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt'
    
    #we need the altIDMap
    hpoOnt = HPOUtil.HPOUtil(ontFN, annotFN)
    hpoRanks, hpoLen = hpoOnt.rank(hpoTerms)

    #call the shared command
    (labelBlocks, values, catValues, featureLabels, catLabels, catBreak, codiDump) = prepDataShared(catMeta, codiDumpFN, pyxisRanks, pyxisLen, hpoRanks, hpoLen)

    variantInfo = []
    for v in codiDump:
        variantInfo.append((v['variant'], v['genes']))

    #create the full data
    vVals = np.vstack(values)
    vCat = np.vstack(catValues)
    xFinal = np.hstack((vVals, vCat))
    featureLabels += catLabels

    #now return the full dataset along with featureLabels
    return (variantInfo, xFinal, featureLabels)

def prepDataShared(catMeta, codiDumpFN, pyxisRanks, pyxisLen, hpoRanks, hpoLen):
    '''
    this function should do all the work
    '''
    #build some lookup dicts
    indexDict = {}
    indexLen = {}
    ignoreDict = {}
    catLabels = []

    #treat all categorical things the same other than where it is sourced from
    labelBlocks = ['allelicInformation', 'genes', 'sequencingFields', 'variantTranscripts']
    floatDict = {k : [] for k in labelBlocks}
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
            elif d['interpret'] == 'float':
                floatDict[lg].append((d['key'], d['defaultValue']))
            elif d['interpret'] == 'float_reduce':
                if d['reduceFunction'] == 'max':
                    reduceFunc = max
                elif d['reduceFunction'] == 'min':
                    reduceFunc = min
                else:
                    raise Exception("Unknown reduce function: %s" % (d['reduceFunction'], ))
                floatDict[lg].append((d['key'], reduceFunc, d['defaultValue']))
            else:
                raise Exception("Unknown interpretation: %s" % (d['interpret'], ))
    
    #These are the fields used from CODICEM dumps
    floatValues = floatDict['allelicInformation']
    fieldsOnly = set([l for l, d in floatValues]) | set([d['key'] for d in catMeta['allelicInformation']])
    
    #(label, reduce function, default)
    geneValues = floatDict['genes']
    geneFieldsOnly = set([l for l, r, d in geneValues]) | set([d['key'] for d in catMeta['genes']])
    
    #(label, default)
    seqValues = floatDict['sequencingFields']
    seqFieldsOnly = set([l for l, d in seqValues]) | set(['chromosome', 'position', 'ref allele', 'alt allele'])
    
    #numerical variant-transcript information
    transValues = floatDict['variantTranscripts']
    transLabelsOnly = set([l for l, r, d in transValues]) | set([d['key'] for d in catMeta['variantTranscripts']])

    #load the raw codi dump in json format
    codiDump = loadCodiDump(codiDumpFN, fieldsOnly, geneFieldsOnly, seqFieldsOnly, transLabelsOnly)
    if len(codiDump) == 0:
        raise Exception('%s is empty' % (codiDumpFN, ))

    #now, we have all the pieces, put them together into tests
    values = []
    catValues = []
    catBreak = {}
    #classifications = []
    #repDicts = []
    
    #go through each variant now and format the inputs to be seen by a classifier
    for variant in codiDump:
        vData = []
        catData = []
        
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
                elif fieldType == 'float' or fieldType == 'float_reduce':
                    #these are handled elsewhere, although eventually we should maybe move them for code consistency
                    #TODO: above comment?
                    pass
                else:
                    raise Exception("Unexpected interpretation type: "+fieldType)

        #uncomment to see each variant's info
        #print(vName, isPrimary, vData, sep='\t')
        
        #we add all values to our lists for the case
        values.append(vData)
        catValues.append(catData)
        
    featureLabels = (['PyxisMap', 'HPO-cosine']+
        [l for l, d in floatValues]+
        [l for l, r, d in geneValues]+
        [l for l, d in seqValues]+
        [l for l, r, d in transValues]
    )

    #values is the variant data that is for single values
    #labelBlocks is the labels in the metadata
    #values is the numerical data
    #feature labels contains the labels for numerical data
    #catValues is the categorical data
    #catLabels contains the categorical labels
    #catBreak contains categorical data broken up by category (used for PCA version)
    #codiDump contains the raw codi information in case its needed
    return (labelBlocks, values, catValues, featureLabels, catLabels, catBreak, codiDump)