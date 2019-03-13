'''
This file should contain helper functions to parse and load Exomiser ranking files
'''

import os

def loadAllPVPRanks(sl):
    '''
    Loads a .res file from PVP
    @param sl - the SL number to load
    '''
    #fn = '/Users/matt/githubProjects/VarSight/exomiser_results/'+sl+'-hiphive-genome-PASS_ONLY.json'
    fn = '/Users/matt/githubProjects/VarSight/pvp_results/%s_results.txt.res' % (sl, )
    fp = open(fn, 'rt')

    #first line is the header
    header = fp.readline().rstrip()
    expectedHeader = ['Chr', 'Start', 'Ref', 'Alt', 'GT', 'Gene', 'CADD', 'GWAVA', 'DANN', 'Sim_Score', 'Prediction_Score']
    #print(header.split('\t'))
    #print(expectedHeader)
    assert(header.split('\t') == expectedHeader)

    #the file is a little weird, ranks are done by gene first, and the variants are a sub-group inside the gene
    orderedRet = []
    for l in fp:
        pieces = l.split('\t')
        chrom = pieces[0]
        pos = pieces[1]
        ref = pieces[2]
        alt = pieces[3]
        orderedRet.append((chrom, pos, ref, alt))
    
    fp.close()
    
    return orderedRet

def getTargetRanks(sl, foundPrimaries):
    '''
    Simple utility to load the SL file and find the rank of each primary that Exomiser actually ranked; None for any it didn't
    @param sl - the identifier for the patient
    @param foundPrimaries - list of primaries in tuple form (chrom, position, ref allele, alt allele)
    @return - a list of their 0-based rank from the Exomiser JSON output
    '''
    fullRanks = loadAllPVPRanks(sl)
    ret = []
    for fp in foundPrimaries:
        try:
            ind = fullRanks.index(fp)
        except:
            assert(len(fullRanks) != 0)
            ind = -1*len(fullRanks)
            
            for fr in fullRanks:
                print('fr', fr)
            print('ERROR: missing reported variant in PVP output')
            print('Missing:', fp)
            raise Exception('See above')
            
        ret.append(ind)
    return ret

if __name__=='__main__':
    #for testing purposes
    pass