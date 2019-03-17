'''
Set of functions that pull the PyxisMap data for us
'''

import json
import requests

import pronto

def getPyxisMapResults(hpoTerms, pyxisRoot='http://0.0.0.0:4999'):
    '''
    This function hits the PyxisMap /rank endpoint to get a single ranking of the genes
    @param hpoTerms - a set of core HPO terms of interest
    @return - tuple
        a dictionary where key is a gene name and value is a tuple of (0-based ranking, weight)
        number of genes ranked
    '''
    #TODO: make this URL configurable?
    pyxisURL = pyxisRoot+'/rank'
    headers = {'Content-type': 'application/json'}
    data_json = json.dumps(sorted(hpoTerms))
    
    response = requests.post(pyxisURL, data=data_json, headers=headers)
    if response.status_code == 200:
        respJson = response.json()
        ret = {d['label']: (i, d['weight']) for i, d in enumerate(respJson['rankings'])}
        return ret, len(respJson['rankings']), response
    else:
        return {}, 0, response

def getPyxisPPIResults(geneTerms):
    '''
    This function hits the PyxisMap /protrank endpoint to get a single ranking of the genes by PPI
    @param geneTerms - the set of core genes of interest
    @return - tuple
        a dictionary where key is a gene name and value is a tuple of (0-based ranking, weight)
        number of genes ranked
    '''
    #TODO: make this URL configurable?
    pyxisURL = 'http://0.0.0.0:4999/protrank'
    headers = {'Content-type': 'application/json'}
    data_json = json.dumps(sorted(geneTerms))
    
    response = requests.post(pyxisURL, data=data_json, headers=headers)
    if response.status_code == 200:
        respJson = response.json()
        ret = {d['label']: (i, d['weight']) for i, d in enumerate(respJson['rankings'])}
        return ret, len(respJson['rankings']), response
    else:
        return {}, 0, response
    