'''
This file should contain helper functions to load the UDN summary DB information
'''

import openpyxl

def loadSummaryDatabase(altIDMap):
    '''
    @param altIDMap - an alternate ID map for the HPO terms; dict of HPO term and the term to replace it with
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
    #core info to read and fill in
    SUMMARY_DB_FILENAME = '/Users/matt/githubProjects/VarSight/data/UDN Summary Database_v3.xlsx'
    wb = openpyxl.load_workbook(SUMMARY_DB_FILENAME)
    ret = {}
    
    #first load HPO terms for each patient
    hpoSheet = wb['Patient HPO terms']
    hpoDict = {}
    for i, row in enumerate(hpoSheet):
        if i == 0:
            for j, c in enumerate(row):
                if c.value != None:
                    hpoDict[c.value] = j
        else:
            slids = row[hpoDict['SL#']].value
            hpo = row[hpoDict['Phenotips HPO']].value
            
            if slids != None:
                for slid in slids.split(','):
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
    
    #second load the primaries
    primarySheet = wb['Primary variants']
    ALLOWED_PATHOGENICITY = set(['PATHOGENIC', 'LIKELY_PATHOGENIC', 'VARIANT_OF_UNCERTAIN_SIGNIFICANCE'])
    primaryDict = {}
    for i, row in enumerate(primarySheet):
        if i == 0:
            for j, c in enumerate(row):
                if c.value != None:
                    primaryDict[c.value] = j
        else:
            slid = row[primaryDict['SL#']].value
            genes = row[primaryDict['GENE']].value
            coords = row[primaryDict['GENOMIC COORDINATES']].value
            pathLevel = row[primaryDict['REPORTED PATHOGENICITY']].value
            
            if (slid != None and genes != None and pathLevel != None and coords != None and
                pathLevel in ALLOWED_PATHOGENICITY):
                #strip out all g.
                coordsMod = coords.replace('g.', '')
                
                if (slid not in ret):
                    ret[slid] = {}
                if ("primaries" not in ret[slid]):
                    ret[slid]['primaries'] = []
                    
                print(slid, genes, pathLevel, coordsMod, sep='\t')
                ret[slid]['primaries'].append({
                    "genes" : list(set(genes.split(','))),
                    "variant" : coordsMod,
                    "path_level" : pathLevel
                })
    
    #add empty things for consistency to all keys that are missing them
    for sl in ret:
        if ('hpoTerms' not in ret[sl]):
            ret[sl]['hpoTerms'] = set([])
        if ('primaries' not in ret[sl]):
            ret[sl]['primaries'] = []
    
    return ret
    