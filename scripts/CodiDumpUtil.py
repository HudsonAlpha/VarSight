'''
This file should contain helper functions to parse and load CODICEM filter files
'''

import json
import os

def loadCodiDump(sl, EXTRA_FIELDS, GENE_FIELDS, SEQ_FIELDS, TRANSCRIPT_FIELDS):
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
    fn = '/Users/matt/githubProjects/VarSight/CLI_primary_V6/'+sl+'_results.json'
    if not os.path.exists(fn):
        #print('Missing file: '+fn)
        return []
    
    #load the json file
    fp = open(fn, 'rt')
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
        variantOrderedRank = variant['variantOrderedRank']
        
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
