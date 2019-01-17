
import json
import os

def createFilterImages(filterFN, prefix):
    '''
    This script will generate DOT graph files and then images for a given filter.  It first draws the overall filter
    layout and then an image for each sub-filter (which may be composed of ANDs or ORs or some combo)
    '''
    fp = open(filterFN, 'rt')
    j = json.load(fp)
    fp.close()

    overallLayers = ['All Variants']
    for i, filt in enumerate(j['queries']):
        layerPrefix = prefix+'layer'+str(i)
        createLayerImage(layerPrefix, filt)
        overallLayers.append(filt['name'])
    overallLayers.append('Passed Variants')

    #draw the overall
    overallDotFN = prefix+'overall.dot'
    overallPngFN = prefix+'overall.png'

    fp = open(overallDotFN, 'wt+')
    fp.write('digraph {\n')
    prev = None
    for i, n in enumerate(overallLayers):
        nodeLabel = 'node'+str(i)
        if n in ['All Variants', 'Passed Variants']:
            shape = 'Mdiamond'
        else:
            shape='rect'
        fp.write('\t'+nodeLabel+' [label="'+n+'", shape='+shape+'];\n')

        if prev != None:
            fp.write('\t'+prev+' -> '+nodeLabel+';\n')
        prev = nodeLabel
    fp.write('}\n')
    fp.close()

    #cmd: dot -Tpng -o <filename.png> <input.dot>
    os.system('dot -Tpng -o '+overallPngFN+' '+overallDotFN)

def createLayerImage(layerPrefix, layerDict):
    '''
    Creates a DOT and image file for a given layer
    @param layerPrefix - the prefix of the filename, creates a <layerPrefix>.dot and <layerPrefix>.png file
    @param layerDict - the dictionary corresponding to the layer
    '''
    dotFN = layerPrefix+'.dot'
    pngFN = layerPrefix+'.png'
    fp = open(dotFN, 'wt+')
    fp.write('digraph {\n')
    
    print(layerDict['name'])
    nodes = [('source', 'All Variants'), ('dest', 'Pass Variants')]
    for n, l in nodes:
        fp.write('\t'+n+' [label="'+l+'" shape=Mdiamond];\n')

    #create a subgraph for labeling purposes
    fp.write('\tsubgraph cluster {\n')
    fp.write('\t\tlabel="'+layerDict['name']+'";\n')
    fNode, lNode, edgeList = recurseCriteria(fp, layerDict, [0])
    fp.write('\t}\n')
    fp.write('\t%s -> %s;\n' % (nodes[0][0], fNode))
    fp.write('\t%s -> %s;\n' % (lNode, nodes[1][0]))
    for s, d in edgeList:
        fp.write('\t'+s+' -> '+d+';\n')
    fp.write('}\n')
    fp.close()

    #cmd: dot -Tpng -o <filename.png> <input.dot>
    os.system('dot -Tpng -o '+pngFN+' '+dotFN)

REP_DICT = {
    "GREATER_THAN_OR_EQUAL" : "&ge;",
    "GREATER_THAN" : "&gt;",
    "EQUALS" : "==",
    "DOES_NOT_EQUAL" : "&ne;",
    "LESS_THAN_OR_EQUAL" : "&le;",
    "CONTAINS" : "&ni;",
    "LESS_THAN" : "&lt;"
}
def recurseCriteria(fp, currDict, nodeLabelList):
    '''
    This will recursively generate the appropriate things for parsing the filters and save them to the DOT file
    @param fp - the file pointer we are writing to 
    @param currDict - current dictionary, this is recursively handed down
    @param nodeLabelList - list of numerical identifiers for parsing things
    @return - tuple (firstNode, lastNode, edges)
        firstNode - the ID of the source node in this recursion
        lastnode - the ID of the last node in this recursion (can be the same as firstNode)
        edges - a list of (source, dest) pairs of edges
    '''
    ret = []
    for k in currDict:
        print('\t', k, currDict[k])
    if 'criteria' in currDict:
        #this is a simple filter, we won't need to recurse
        nodeLabel = 'node'+'v'.join([str(x) for x in nodeLabelList])
        assert(len(currDict['criteria']['op']) == 1)
        
        #figure out the relationship value
        relat = REP_DICT.get(currDict['criteria']['op'][0], currDict['criteria']['op'][0])
        if 'join' in currDict:
            relatStr = relat+'\\n('+currDict['join']+')'
        else:
            relatStr = relat

        critStr = '{'+'|'.join([currDict['field'], relatStr, '\\n'.join(currDict['criteria']['val'])])+'}'
        nodeText = '{'+'|'.join([currDict['logic'], critStr])+'}'
        fp.write('\t\t'+nodeLabel+' [shape=record label="'+nodeText+'"];\n')
        
        #save the results
        fRet = lRet = nodeLabel

    elif (currDict.get('join', None) in ['OR', 'AND']):
        allSingles = True
        for subFilter in currDict['filters']:
            if len(subFilter['filters']) > 0:
                print('sub', subFilter)
                allSingles = False
        
        if allSingles:
            #this is all single values comparisons, make one box
            recordLabel = '<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">'
            recordLabel += '<TR><TD COLSPAN="4">'+currDict['logic']+'</TD></TR>'
            recordLabel += '<TR><TD COLSPAN="4">BEGIN: JOIN('+currDict['join']+')</TD></TR>'
            for subFilter in currDict['filters']:
                if subFilter['field'] == 'Effects':
                    print('SUB:')
                    print(subFilter)
                    print('join' in subFilter)
                    print()
                recordLabel += '<TR>'
                assert(len(subFilter['criteria']['op']) == 1)
                
                #figure out relationship string
                relat = REP_DICT.get(subFilter['criteria']['op'][0], subFilter['criteria']['op'][0])
                if 'join' in subFilter:
                    relatStr = relat+'<br/>('+subFilter['join']+')'
                else:
                    relatStr = relat
                
                rowEntry = [subFilter['logic'], subFilter['field'], relatStr, '<br/>'.join(subFilter['criteria']['val'])]
                for rEntry in rowEntry:
                    recordLabel += '<TD>'+rEntry+'</TD>'
                recordLabel += '</TR>'
            recordLabel += '<TR><TD COLSPAN="4">END: JOIN('+currDict['join']+')</TD></TR>'
            recordLabel += '</TABLE>'
            
            nodeJoin = 'join'+'v'.join([str(x) for x in nodeLabelList])
            fp.write('\t\t'+nodeJoin+' [shape=none margin=0 label=<'+recordLabel+'>];\n')
            fRet = lRet = nodeJoin
        else:
            #this has as least one non-simple value
            nodeSplit = 'split'+'v'.join([str(x) for x in nodeLabelList])
            nodeJoin = 'join'+'v'.join([str(x) for x in nodeLabelList])
            fp.write('\t\t'+nodeSplit+' [shape=record label="{'+currDict['logic']+'|BEGIN: JOIN('+currDict['join']+')}"];\n')
            fp.write('\t\t'+nodeJoin+' [shape=record label="END: JOIN('+currDict['join']+')"];\n')
            for i, subFilter in enumerate(currDict['filters']):
                newSubList = nodeLabelList+[i]
                fSub, lSub, subEdges = recurseCriteria(fp, subFilter, newSubList)
                ret += subEdges
                ret.append((nodeSplit, fSub))
                ret.append((lSub, nodeJoin))
            fRet = nodeSplit
            lRet = nodeJoin

    else:
        fRet = lRet = 'tempBadNode'
        print('UNHANDLED')

    return (fRet, lRet, ret)

if __name__ == '__main__':
    filterFN = '../CODI_metadata/v6_clinical_filter_with_disease.json'
    filterImagePrefix = '../paper/supplements/filter_supplement_'
    createFilterImages(filterFN, filterImagePrefix)