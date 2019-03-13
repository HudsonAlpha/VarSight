'''
This file should contain helper functions to parse and load Exomiser ranking files
'''

import os
import vcf

def loadAllPhenGenRanks(sl):
    '''
    Loads the output from PhenGen
    @param sl - the SL number to load
    '''
    geneFN = '/Users/matt/githubProjects/VarSight/phengen_results/%s_hpo_%s_results_genescores.txt' % (sl, sl)
    vcfFN = '/Users/matt/githubProjects/VarSight/phengen_results/%s_hpo_%s_results_variants_for_topgenes.vcf' % (sl, sl)

    fp = open(geneFN, 'rt')

    #first line is the header
    header = fp.readline().rstrip()
    expectedHeader = ['#GENE_ID', 'PROBABILITY_DAMAGING']
    assert(header.split('\t') == expectedHeader)
    
    gDict = {}
    for l in fp:
        pieces = l.rstrip().split('\t')
        geneName = pieces[0]
        score = float(pieces[1])
        gDict[geneName] = score
    fp.close()

    #now go through the VCF
    vcfReader = vcf.Reader(filename=vcfFN, compressed=False)
    unsorted = []
    for variant in vcfReader:
        chrom = str(variant.CHROM)
        pos = str(variant.POS)
        ref = variant.REF
        alt = str(variant.ALT[0])
        #print(chrom, pos, ref, alt, variant.INFO)
        geneList = variant.INFO['GNID']
        for i, g in enumerate(geneList):
            if g.endswith('_neighboring'):
                geneList[i] = g[:-len('_neighboring')]
        dcod = float(variant.INFO.get('DCOD', 0.0)) #can be absent if non-coding
        dreg = float(variant.INFO['DREG'])
        
        maxGeneVal = max([gDict[g] for g in geneList])
        #print(maxGeneVal, dcod, dreg, chrom, pos, ref, alt)
        unsorted.append((maxGeneVal, dcod, dreg, chrom, pos, ref, alt))
    
    #sorting by gene, then DCOD, then DREG (bigger is more deleterious for all of them)
    #for v in unsorted:
    #    print(v)
    unsorted.sort(reverse=True)

    #now pull out the part we care about
    orderedRet = []
    for e in unsorted:
        orderedRet.append(e[3:])
    
    return orderedRet

def getTargetRanks(sl, foundPrimaries):
    '''
    Simple utility to load the SL file and find the rank of each primary that Exomiser actually ranked; None for any it didn't
    @param sl - the identifier for the patient
    @param foundPrimaries - list of primaries in tuple form (chrom, position, ref allele, alt allele)
    @return - a list of their 0-based rank from the Exomiser JSON output
    '''
    fullRanks = loadAllPhenGenRanks(sl)
    ret = []
    for fp in foundPrimaries:
        try:
            ind = fullRanks.index(fp)
        except:
            assert(len(fullRanks) != 0)
            ind = -1*len(fullRanks)
            
            #for fr in fullRanks:
            #    print('fr', fr)
            #print('WARNING: missing reported variant in PhenGen output for '+sl)
            #print('Missing:', fp)
            #raise Exception('See above')
            
        ret.append(ind)
    return ret

if __name__=='__main__':
    #for testing purposes
    pass