'''
This file should contain helper functions to load the UDN summary DB information
'''

import csv
import openpyxl

def loadSummaryDatabase(altIDMap, path_only):
    '''
    @param altIDMap - an alternate ID map for the HPO terms; dict of HPO term and the term to replace it with
    @param path_only - if True, only load pathogenic and likely pathogenic variants as true positive
    @return - a JSON dictionary with the following structure:
    {
        "SL##" : {
            "hpoTerms" : set(["HP:0000001", ...]),
            "primaries" : [
                {
                    "genes" : ["BMPR2"],
                    "variant" : "chr7:146818173G>A",
                    "path_level" : "PATHOGENIC"
                },
                ...
            ]
        },
        ...
    }
    '''
    #need to get UDN to SL##s
    TRACKING_FILENAME = '/Users/matt/githubProjects/VarSight/data/UDN Research Tracking_v1.xlsx'
    wbt = openpyxl.load_workbook(TRACKING_FILENAME)
    mainSheet = wbt['Original case order']
    mainDict = {}
    udnToSL = {}
    ignoreSet = set([None, 'SANGER ONLY', 'N/A', 'SANGER_ONLY', 'NO DNA LEFT'])
    for i, row in enumerate(mainSheet):
        if i == 0:
            continue
        elif i == 1:
            for j, c in enumerate(row):
                if c.value != None:
                    mainDict[c.value] = j
        else:
            sl_raw = row[mainDict['HA_SL#']].value
            udnID = row[mainDict['UDN#']].value

            if (sl_raw not in ignoreSet) and udnID != None:
                for slid in sl_raw.split('/'):
                    sl = slid.rstrip('R')
                    if len(sl) != 8:
                        raise Exception('Unexpected slid:'+sl+' from '+sl_raw)
                    
                    if udnID in udnToSL:
                        udnToSL[udnID].append(sl)
                    else:
                        udnToSL[udnID] = [sl]
    
    #core info to read and fill in
    SUMMARY_DB_FILENAME = '/Users/matt/githubProjects/VarSight/data/UDN Summary Database_v4.xlsx'
    wb = openpyxl.load_workbook(SUMMARY_DB_FILENAME)
    ret = {}

    #first load HPO terms for each patient
    hpoSheet = wb['Phenotips_HPO_04Jan19']
    hpoDict = {}
    for i, row in enumerate(hpoSheet):
        if i == 0:
            for j, c in enumerate(row):
                if c.value != None:
                    hpoDict[c.value] = j
        else:
            udnID = row[hpoDict['UDN#']].value
            hpo = row[hpoDict['Phenotips HPO']].value

            slids = udnToSL.get(udnID, [])
            for slid in slids:
                if hpo != None:
                    #create this entry if it doesn't exist
                    if (slid not in ret):
                        ret[slid] = {
                            "hpoTerms" : set([])
                        }
                    
                    #add each HPO term
                    for h in hpo.split(','):
                        if h[0:3] != 'HP:' or len(h) != 10:
                            print('Ignoring term:', h)
                        else:
                            #get the alternate ID if it exists, otherwise use the normal
                            ret[slid]['hpoTerms'].add(altIDMap.get(h, h))

    #second, load the primaries
    fpcsv = open('/Users/matt/githubProjects/VarSight/CODI_metadata/3Jan2019_UDN_Variants_all.csv', 'rt')
    primarySheet = csv.DictReader(fpcsv)
    revOrderPath = ['PATHOGENIC', 'LIKELY_PATHOGENIC', 'VARIANT_OF_UNCERTAIN_SIGNIFICANCE', 'BENIGN']
    if path_only:
        ALLOWED_PATHOGENICITY = set(revOrderPath[0:2])
    else:
        ALLOWED_PATHOGENICITY = set(revOrderPath[0:3])
        
    for row in primarySheet:
        slid = row['Sample']
        genes = row['Gene']
        coords = row['Genomic Location']
        pathLevel = row['Classification']

        weakestPL = max([revOrderPath.index(pl) for pl in pathLevel.split(',')])
        assert(weakestPL >= 0)

        if (slid != None and genes != None and pathLevel != None and coords != None and
            pathLevel in ALLOWED_PATHOGENICITY):
            #strip out all g.
            coordsMod = coords.replace('g.', '').replace('Chr', 'chr')

            #there are some instances with multiple reporting categories, default to the lower category
            if (slid not in ret):
                ret[slid] = {}
            if ("primaries" not in ret[slid]):
                ret[slid]['primaries'] = []

            for pDict in ret[slid]['primaries']:
                if coordsMod == pDict['variant']:
                    #there are instances where the old version was missing a gene the new had
                    pDict['genes'] = list(set(pDict['genes']) | set(genes.split(',')))
                    weakestPL = max(weakestPL, revOrderPath.index(pDict['path_level']))
                    pDict['path_level'] = revOrderPath[weakestPL]
                    break
            else:
                #print(slid, genes, pathLevel, coordsMod, sep='\t')
                ret[slid]['primaries'].append({
                    "genes" : list(set(genes.split(','))),
                    "variant" : coordsMod,
                    "path_level" : revOrderPath[weakestPL]
                })

    fpcsv.close()

    #add empty things for consistency to all keys that are missing them
    for sl in ret:
        if ('hpoTerms' not in ret[sl]):
            ret[sl]['hpoTerms'] = set([])
        if ('primaries' not in ret[sl]):
            ret[sl]['primaries'] = []
    
    return ret
    