#!/usr/bin/env python3

import numpy as np

import pronto

class MicaDS:
    def __init__(self, ontologyFN):
        '''
        @param ontologyFN - the filename to load
        '''
        self.nodes, self.edges, self.parents, self.altIDMap = loadGraphStructure(ontologyFN)
        self.ancestors = deriveAncestors(self.nodes, self.parents)
        self.synonyms = getTermsAndSyns(ontologyFN)
        
        #assign a specific number to each node label for use in other calculations
        self.nToK = {n : i for i, n in enumerate(sorted(self.nodes))}

    def cosineSets(self, s1, s2):
        '''
        Calculates the full and normalized similarity between two phenotype sets
        @param s1 - set 1
        @param s2 - set 2
        @return - tuple (full, normalized) similarity between sets
        '''
        numTerms = len(self.nToK.keys())
        pv1 = np.zeros(dtype='float', shape=(numTerms, ))
        pv2 = np.zeros(dtype='float', shape=(numTerms, ))
        
        for p in s1:
            for a in self.ancestors[p]:
                pv1[self.nToK[a]] = 1.0

        for p in s2:
            for a in self.ancestors[p]:
                pv2[self.nToK[a]] = 1.0

        dp = np.dot(pv1, pv2)
        l1 = np.sqrt(np.dot(pv1, pv1))
        l2 = np.sqrt(np.dot(pv2, pv2))

        return (dp, dp/(l1*l2))

    def projectSets(self, s1, s2):
        '''
        Calculates the projected similarity of s1 onto s2
        @param s1 - set 1, any terms present here that are not in s2 are removed
        @param s2 - set 2
        @return - tuple (full, normalized) projection of s1 onto s2
        '''
        numTerms = len(self.nToK.keys())
        pv1 = np.zeros(dtype='float', shape=(numTerms, ))
        pv2 = np.zeros(dtype='float', shape=(numTerms, ))

        #get the shared list of HPO terms
        p1 = s1 & s2

        for p in p1:
            for a in self.ancestors[p]:
                pv1[self.nToK[a]] = 1.0

        for p in s2:
            for a in self.ancestors[p]:
                pv2[self.nToK[a]] = 1.0

        dp = np.dot(pv1, pv2)
        l1 = np.sqrt(np.dot(pv1, pv1))
        l2 = np.sqrt(np.dot(pv2, pv2))

        if (l1 == 0.0 or l2 == 0.0):
            return (0.0, 0.0)
        else:
            return (dp, dp/(l1*l2))

    def jaccardSets(self, s1, s2):
        '''
        Calculates the Jaccard index for two phenotype sets
        @param s1 - set 1
        @param s2 - set 2
        @return - jaccard index (float between 0.0 and 1.0)
        '''
        numTerms = len(self.nToK.keys())
        pv1 = np.zeros(dtype='float', shape=(numTerms, ))
        pv2 = np.zeros(dtype='float', shape=(numTerms, ))

        for p in s1:
            for a in self.ancestors[p]:
                pv1[self.nToK[a]] = 1.0

        for p in s2:
            for a in self.ancestors[p]:
                pv2[self.nToK[a]] = 1.0

        minVals = np.minimum(pv1, pv2)
        maxVals = np.maximum(pv1, pv2)

        return np.sum(minVals) / np.sum(maxVals)

def getTermsAndSyns(ontologyFN):
    '''
    @param ontologyFN - the ontology to load
    @return - a list of tuples (id, string) where id is the term ID and string is a display string combined the term name and synonyms
    '''
    ont = pronto.Ontology(ontologyFN)
    ret = []
    for t in ont:
        fullList = [t.name]
        for s in t.synonyms:
            fullList.append(s.desc)
        ret.append((t.id, '; '.join(fullList)))
    ret.sort()
    return ret

def deriveAncestors(nodes, parents):
    '''
    This function will calculate all parent terms for each HPO node
    @param nodes - the set of nodes in the HPO graph
    @param parents - the set of edges where key is the child node in the graph and the value is a set of parent terms
    @return - a dictionary where key is a term and value is a set of nodes that are ancestors of that term (including that term)
    '''
    ret = {}
    
    for n in nodes:
        ret[n] = set([n])
        
        processed = set([])
        todoList = [n]
        while len(todoList) > 0:
            currNode = todoList.pop(0)
            if (currNode not in processed):
                processed.add(currNode)
                ret[n] |= parents[currNode]
                todoList += list(parents[currNode])

    return ret

def loadGraphStructure(fn):
    '''
    This function parses the HPO graph structure and returns information to build the corresponding graph
    @param fn - the .obo file containing the structure
    @return - tuple (nodes, edges, altIdToMain)
        nodes - the set of HPO terms that are nodes
        edges - a dictionary where key is a parent term and value is a set of child terms
        parents - a dictionary where key is a child term and value is a set of parent terms
        altIdToMain - a dictionary where key is an old term and value is the current term (for cases where terms were consolidated)
    '''
    ont = pronto.Ontology(fn)

    nodes = set([])
    edges = {}
    parents = {}
    altIdToMain = {}

    for term in ont:
        if ('true' not in term.other.get('is_obsolete', [])):
            #term isn't obsolete, so add it
            nodes.add(term.id)
            parents[term.id] = set([p.id for p in term.parents])

            #add parent edges
            for p in term.parents:
                if p.id in edges:
                    edges[p.id].add(term.id)
                else:
                    edges[p.id] = set([term.id])

            for a in term.other.get('alt_id', []):
                altIdToMain[a] = term.id
    
    #there is an issue where some nodes have parents that aren't actually in the ontology for one reason or another, this strips those out
    for k in parents:
        missing = set([])
        for p in parents[k]:
            if p not in nodes:
                assert(p not in altIdToMain)
                missing.add(p)
        parents[k] -= missing
        
    return nodes, edges, parents, altIdToMain
