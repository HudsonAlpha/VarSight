'''
This file should contain helper functions to parse and load Exomiser ranking files
'''

import json
import os

def loadAllExomiserRanks(sl):
    '''
    Loads a .json file from Exomiser
    @param sl - the SL number to load
    '''
    fn = '/Users/matt/githubProjects/VarSight/exomiser_results/'+sl+'-hiphive-genome-PASS_ONLY.json'
    fp = open(fn, 'rt')
    j = json.load(fp)
    fp.close()

    #the file is a little weird, ranks are done by gene first, and the variants are a sub-group inside the gene
    orderedRet = []
    for gDict in j:
        #inside the gene dictionary, we want variant evaluations
        #print(gDict['geneSymbol'])
        for vEval in gDict['variantEvaluations']:
            chrom = vEval['chromosome']
            pos = vEval['position']
            ref = vEval['ref']
            alt = vEval['alt']
            #these are not the normal gDots we see typically; changing to chrom pos ref alt and figuring it out
            #gDot = vEval['transcriptAnnotations'][0]['hgvsGenomic']
            #print('\t', chrom, pos, ref, alt)
            
            #make sure there are no doubles for overlapping genes
            if (chrom, pos, ref, alt) not in orderedRet:
                orderedRet.append((str(chrom), str(pos), ref, alt))
    
    return orderedRet

def getTargetRanks(sl, foundPrimaries):
    '''
    Simple utility to load the SL file and find the rank of each primary that Exomiser actually ranked; None for any it didn't
    @param sl - the identifier for the patient
    @param foundPrimaries - list of primaries in tuple form (chrom, position, ref allele, alt allele)
    @return - a list of their 0-based rank from the Exomiser JSON output
    '''
    fullRanks = loadAllExomiserRanks(sl)
    ret = []
    for fp in foundPrimaries:
        try:
            ind = fullRanks.index(fp)
        except:
            ind = len(fullRanks)
        ret.append(ind)
    return ret

if __name__=='__main__':
    #for testing purposes
    pass