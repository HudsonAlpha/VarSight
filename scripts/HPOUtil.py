
import numpy as np

import OntologyUtil

class HPOUtil:
    '''
    This class handles all the complications with loading HPO and returned ranks to users
    '''
    def __init__(self, ontologyFN, annotFN):
        '''
        Constructor of the class, loads the ontology and preps for queries
        @param ontologyFN - the ontology filepath, expects "hp.obo"
        @param annotFN - the annotation file path, expects "ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt"
        '''
        #load the MONDO ontology
        self.ont = OntologyUtil.MicaDS(ontologyFN)
        
        #load the MONDO annotations
        self.g2p, self.p2g = loadHPOAnnot(annotFN, self.ont)
        
        self.numGenes = len(self.g2p.keys())
        self.IC = {}
        for h in self.p2g.keys():
            if len(self.p2g[h]) > 0:
                #if we have a p2g set, then this is the total number of genes divided by the number of genes in that set
                self.IC[h] = np.log2(self.numGenes / len(self.p2g[h]))
            else:
                #if we have no genes, then we can't divide by zero, so set this to a semi-low value of 1.0
                #TODO: should it actually be 0 since it doesn't have genes? not sure
                self.IC[h] = 1.0
        
    def rank(self, hpoSet):
        '''
        This function ranks genes based on their similarity to HPO terms.
        @param hpoSet - the set of HPO terms to compare against
        @return - tuple
            a dictionary where key is a gene name and value is a tuple of (0-based ranking, weight)
            number of genes ranked
        '''
        #go through each gene and score according to the similarity
        scores = []
        for g in self.g2p.keys():
            scores.append((
                self.ont.cosineSets(hpoSet, self.g2p.get(g, set([]))),
                g
            ))
        
        #sort and build the return value
        scores.sort(reverse=True, key=lambda x: (x[0][1], x[0][0]))
        ret = {}
        for i, ((rawScore, cosScore), g) in enumerate(scores):
            ret[g] = (i, cosScore)
        
        return ret, len(self.g2p.keys())

def loadHPOAnnot(fn, hpoOnt):
    '''
    This will load the HPO annotations file
    @param fn - the file to load, expected is ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
    @param hpoOnt - the MicaDS hpoOnt used for filtering out things that are absent
    @return - dictionary where key is a gene name and value is a set of associated MONDO terms
    '''
    g2p = {}
    p2g = {}
    fp = open(fn, 'rt')
    fp.readline()
    for l in fp:
        pieces = l.rstrip().split('\t')
        
        #get the gene and mondo pair
        hpoID = pieces[0]
        gpSymbol = pieces[3]
        
        assert(hpoID not in hpoOnt.altIDMap)
        
        #add this hpo ID and all of its ancestors
        if gpSymbol in g2p:
            g2p[gpSymbol].add(hpoID)
        else:
            g2p[gpSymbol] = set([hpoID])
        g2p[gpSymbol] |= hpoOnt.ancestors[hpoID]
        
        if hpoID in p2g:
            p2g[hpoID].add(gpSymbol)
        else:
            p2g[hpoID] = set([gpSymbol])
        
    fp.close()
    
    #before returning, push p2g up the list
    p2g = pushP2gUp(p2g, hpoOnt.nodes, hpoOnt.edges)
    
    return g2p, p2g

def pushP2gUp(p2g, nodes, edges):
    '''
    This function push the HPO phenotype-to-gene information completely up the "tree"-like graph structure
    @param p2g - the dictionary of HPO terms to genes
    @param nodes - the set of nodes in the HPO graph
    @param edges - the set of edges where key is the parent node in the HPO graph and the value is a set of children
    @return - a modified p2g that pushes all phenotype info full up the tree
    '''
    print('Running pushP2gUp(...)')
    analyzed = set([])
    rootNode = 'HP:0000001'
    stack = [rootNode]
    
    while len(stack) > 0:
        currNode = stack[-1]
        for end in edges.get(currNode, set([])):
            if end not in analyzed:
                stack.append(end)
        
        if currNode == stack[-1]:
            #nothing was added
            x = len(p2g.get(currNode, set([])))
            for end in edges.get(currNode, set([])):
                p2g[currNode] = p2g.get(currNode, set([])) | p2g.get(end, set([]))
            stack.pop()
            analyzed.add(currNode)
    
    return p2g