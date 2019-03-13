'''
This file should contain helper functions to parse and load Exomiser ranking files
'''

import json
import os

def loadAllExomiserRanks(sl, subtype):
    '''
    Loads a .json file from Exomiser
    @param sl - the SL number to load
    @param subtype - the particular run of exomiser
    '''
    #fn = '/Users/matt/githubProjects/VarSight/exomiser_results/'+sl+'-hiphive-genome-PASS_ONLY.json'
    fn = '/Users/matt/githubProjects/VarSight/exomiser_results_v2/%s-%s.json' % (sl, subtype)
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

def getTargetRanks(sl, foundPrimaries, subtype):
    '''
    Simple utility to load the SL file and find the rank of each primary that Exomiser actually ranked; None for any it didn't
    @param sl - the identifier for the patient
    @param foundPrimaries - list of primaries in tuple form (chrom, position, ref allele, alt allele)
    @return - a list of their 0-based rank from the Exomiser JSON output
    '''
    fullRanks = loadAllExomiserRanks(sl, subtype)
    ret = []
    for fpRaw in foundPrimaries:
        #X is changed to 23 for some reason
        if fpRaw[0] == 'X':
            fp = ('23', )+fpRaw[1:]
        else:
            fp = fpRaw
        
        try:
            ind = fullRanks.index(fp)
        except:
            assert(len(fullRanks) != 0)
            ind = -1*len(fullRanks)
            
            for fr in fullRanks:
                print('fr', fr)
            print('ERROR: missing reported variant in Exomiser output')
            print('Missing:', fp)
            raise Exception('See above')
            
        ret.append(ind)
    return ret

if __name__=='__main__':
    #for testing purposes
    pass