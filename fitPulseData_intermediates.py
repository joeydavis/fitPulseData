"""
Created on Mon Apr  1 23:13:53 2013

@author: jhdavis
"""

import pylab
import vizLib
import numpy
import qMS
import glob
from scipy import optimize

#############FITTING THE 70S PEAKS####################
def poolResid1000(P, y, t):
    k = qMS.growthRate(doublingTime_1000)
    return y - qMS.poolFunc(k, t, P)
def poolResid10(P, y, t):
    k = qMS.growthRate(doublingTime_10)
    return y - qMS.poolFunc(k, t, P)

def fitPoolSizes(Proteins, medDict10, medDict1000):
    p0_10 = 0.1
    p0_1000 = 0.1
    
    protPoolDict = {}
    for prot in Proteins:
        ts10 = numpy.array([e[0] for e in medDict10[prot]])
        meas10 = numpy.array([e[1] for e in medDict10[prot]])

        ts1000 = numpy.array([i[0] for i in medDict1000[prot]])
        meas1000 = numpy.array([i[1] for i in medDict1000[prot]])

        p_10 = optimize.leastsq(poolResid10, [p0_10], args=(meas10,ts10))
        p_1000 = optimize.leastsq(poolResid1000, [p0_1000], args=(meas1000,ts1000))
        protPoolDict[prot] = {'10':p_10[0][0], '1000':p_1000[0][0]}
    return protPoolDict


#############FITTING THE INTERMEDIATES####################
def poolIntermediatesResid1000(P, y, t):
    k = qMS.growthRate(doublingTime_1000)
    return y - qMS.poolInterFunc(k, t, P)

def poolIntermediatesResid10(P, y, t):
    k = qMS.growthRate(doublingTime_10)
    return y - qMS.poolInterFunc(k, t, P)

def fitPoolIntermediateSizes(Proteins, medDict10, medDict1000):
    p0_10 = 0.01
    p0_1000 = 0.01
    
    protPoolDict = {}
    for prot in Proteins:
        ts10 = numpy.array([e[0] for e in medDict10[prot]])
        meas10 = numpy.array([e[1] for e in medDict10[prot]])

        ts1000 = numpy.array([i[0] for i in medDict1000[prot]])
        meas1000 = numpy.array([i[1] for i in medDict1000[prot]])

        p_10 = optimize.leastsq(poolIntermediatesResid10, [p0_10], args=(meas10,ts10))
        p_1000 = optimize.leastsq(poolIntermediatesResid1000, [p0_1000], args=(meas1000,ts1000))
        protPoolDict[prot] = {'10':p_10[0][0], '1000':p_1000[0][0]}
    return protPoolDict

def fracInter(flObsInt, PObsTerm, time, k):
    fl70S = qMS.poolFunc(k, time, PObsTerm)
    flInter = qMS.poolInterFunc(k, time, PObsTerm)
    X = (flObsInt - fl70S)/(flInter - fl70S)
    return X

def calcFractionIntermediate(time, proteins, datInter, datTerminal, k):
    fitInter = {}
    P = {}
    X = {}
    for p in proteins:
        X[p] = {}
        fitInter[p] = qMS.poolInterFunc(k, time, datInter[p])
        P[p] = datTerminal[p]
        X[p]['fracX'] = fracInter(fitInter[p], P[p], time, k)
    return X

def poolIntermedFracXResid(P, y, t):
    k = qMS.growthRate(doublingTime_10)
    fracX = 0.65
    return y - qMS.poolInterFracXFunc(k, t, P, fracX)

def fitPoolFixFracX(Proteins, medDict10):
    p0_10 = 0.01
    
    protPoolDict = {}
    for prot in Proteins:
        ts10 = numpy.array([e[0] for e in medDict10[prot]])
        meas10 = numpy.array([e[1] for e in medDict10[prot]])
        
        p_10 = optimize.leastsq(poolIntermedFracXResid, [p0_10], args=(meas10,ts10))
        protPoolDict[prot] = {'10':p_10[0][0], '1000':'blah'}
    return protPoolDict

#############PLOTTING THE INDIVIDUAL FITS####################
def plotPoolPage(Proteins, protPoolDict, medDict10, medDict1000, saveFile=None, 
                 figSize=(7,10), yMax=1.0, funcToFit=qMS.maxLabFunc, title=None, double=True, showPools=True):
    time = numpy.linspace(0, doublingTime_10*2, num=200)
    k10 = qMS.growthRate(doublingTime_10)
    k1000 = qMS.growthRate(doublingTime_1000)
    
    fig = pylab.figure(figsize=figSize)
    Large1Ax = []
    a = 0
    for prot in Proteins:
        p_10 = protPoolDict[prot]['10'] 
        p_10String = qMS.calcPercent(p_10, sigfig=2).split('.')[0]+'%'
        if double:
            p_1000 = protPoolDict[prot]['1000']
            p_1000String = qMS.calcPercent(p_1000, sigfig=2).split('.')[0]+'%'

        ax = fig.add_subplot(5,2,a+1)
        if showPools:
            textstr2 = '\n\n'+'10 $\mu$M : P=' + p_10String
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.025, .95, textstr2, transform=ax.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', bbox=props, color='r')
            if double:
                textstr1 = '1 mM : P=' + p_1000String
                ax.text(0.025, .925, textstr1, transform=ax.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='b')
        
        ax.scatter([i[0] for i in medDict10[prot]], [i[1] for i in medDict10[prot]], c = 'r')
        ax.plot(time, qMS.maxLabFunc(k10, time), c='r')
        if double:
            ax.scatter([i[0] for i in medDict1000[prot]], [i[1] for i in medDict1000[prot]], c = 'b')
            ax.plot(time, qMS.maxLabFunc(k1000, time), c='b')

        ax.plot(time, funcToFit(k10, time, p_10), c='r', ls='--')
        if double:
            ax.plot(time, funcToFit(k1000, time, p_1000), c='b', ls='--')
            
        ax.set_xlim([0,doublingTime_10*2])
        ax.set_ylim([0,yMax])
        ax.set_yticks([0,yMax/4, yMax/2, 3*yMax/4, yMax])
        ax.set_xticks([0,doublingTime_10*2/5, doublingTime_10*2*2/5, doublingTime_10*2*3/5, doublingTime_10*2*4/5, doublingTime_10*2])
        ax.yaxis.tick_left()
        ax.xaxis.tick_bottom()
        
        ax.text(0.5, 0.95, prot[4:], transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='center')
        Large1Ax.append(ax)
        a = a+1
    
    fig.text(0.005, 0.1165, "fraction labeled", fontsize=12, rotation='vertical', verticalalignment='center')
    fig.text(0.005, 0.317, "fraction labeled", fontsize=12, rotation='vertical', verticalalignment='center')
    fig.text(0.005, 0.505, "fraction labeled", fontsize=12, rotation='vertical', verticalalignment='center')
    fig.text(0.005, 0.7, "fraction labeled", fontsize=12, rotation='vertical', verticalalignment='center')
    fig.text(0.005, 0.9, "fraction labeled", fontsize=12, rotation='vertical', verticalalignment='center')
    fig.text(0.25, 0.0025, "time (mins)", fontsize=12, horizontalalignment='center')
    fig.text(0.75, 0.0025, "time (mins)", fontsize=12, horizontalalignment='center')
    if not (title is None):
        fig.text(0.5, 1, title, fontsize=16, horizontalalignment='center', verticalalignment='top')
    fig.tight_layout()
    if not (saveFile is None):
        pylab.savefig(saveFile)
    return fig
    
def plotPoolBar(data, keys, colors, title, figname=None, barNames=None, sortingKey=0, figSize=(10,4)):
    """plotPoolBar plots a poolSizDict datastructure as bar graphs
        A pointer to the figure is returned

    :param poolSizeDict: a dictionary of dictionaries of the values to be plotted. Example generated by fitPoolSizes
    :type poolSizeDict: a dict of dicts: Should be of the form {entry name : {valueName1 : entryValue1, valueName2 : entryValue2}
    :param keys: a list of keys to be plotted
    :type keys: list of strings
    :param colors: a list of colors to be plotted
    :type colors: list of strings
    :param title: the title of the plot
    :type title: string
    :param offset: a float offset factor if you want to alter the values uniformly
    :type offset: float
    :param figname: the full path/name to save the file
    :type figname: string
    :param barNames: a list of the keys from the dictionary to be plotted (if None, all will be plotted)
    :type barNames: a list of strings
    :param sortingKey: the index of the key to sort on (will be sorted highest to lowest in value)
    :type sortingKey: int
    :returns:  the figure generated. If figname!=None, a figure is also saved

    """
    bar = pylab.figure(figsize=figSize)
    barax = bar.add_subplot(111)
    bars = len(barNames)
    ind = numpy.arange(bars)
    width = 0.9/len(keys)
    sortedTuples = []
    if barNames is None:
        barNames = data.keys()
    barNames = qMS.sort_nicely(barNames)
    for name in barNames:
        myTuple = [name]
        for k in keys:
            myTuple.append(data[name][k])
        sortedTuples.append(myTuple)
    sortedTuples = sorted(sortedTuples, key=lambda pool:pool[sortingKey], reverse=True)
    names = [i[0] for i in sortedTuples]
    rects = []
    for b in range(0,len(keys)):
        vs = [i[b+1] for i in sortedTuples]
        xs = ind+width*b
        rects.append(barax.bar(xs, vs, width, color=colors[b]))
    
    barax.set_xticks(ind+0.5)
    barax.set_xticklabels([i[4:] for i in names], rotation='vertical', horizontalalignment='center', verticalalignment='top')    
    barax.set_xlim([0,len(ind)])
    barax.set_ylabel('fit pool size (P)')
    barax.tick_params(axis='y', direction='out')
    barax.tick_params(axis='x', direction='out')
    barax.xaxis.tick_bottom()
    barax.yaxis.tick_left()
    #barax.set_ylim([0,6.5])
    barax.set_title(title)
    bar.tight_layout()
    if not (figname is None):
        pylab.savefig(figname)
    return bar

def plotDataSets(files, names, num, den, subunits, title, yLabel, colors, saveName=None, 
                 yMax=1.25, figSize=(22,7), median=False, legendLoc='upper left',
                legendCols=3, normProtein=None):
                    
    ax = vizLib.makePlotWithFileList(files, num, den, AllProteins=subunits, yMax=yMax, 
                                     names=names, colors=colors, figSize=figSize, 
                                     median=median, normProtein=normProtein)
                                     
    pylab.legend(loc=legendLoc, ncol=legendCols)
    pylab.xticks(numpy.arange(1,len(AllSubunits)+1,1), [item[4:] for item in AllSubunits], rotation=45)
    ax.set_title(title, multialignment='center')
    ax.set_ylabel(yLabel)
    pylab.tight_layout()
    if not (saveName is None):
        pylab.savefig(saveName)
    return ax

def plotCompData(xdat, ydat, Proteins, title=None, xlabel='dat1', ylabel='dat2', xMax=1.5, yMax=1.5, figSize=(10,10), saveFile=None):
    x = [numpy.median(xdat[i]) for i in Proteins]
    y = [numpy.median(ydat[i]) for i in Proteins]    
    scat = pylab.figure(figsize=figSize)
    scatAx = scat.add_subplot(111)    
    scatAx.scatter(x,y, c='b', s=150)
    scatAx.set_title(title)
    scatAx.set_xlabel(xlabel)
    scatAx.set_ylabel(ylabel)
    scatAx.set_xlim([-0.1,xMax])
    scatAx.set_ylim([-0.1,yMax])
    scatAx.set_xticks([0,xMax/5,xMax/5*2,xMax/5*3,xMax/5*4,xMax])
    scatAx.set_yticks([0,yMax/5,yMax/5*2,yMax/5*3,yMax/5*4,yMax])
    scatAx.yaxis.tick_left()
    scatAx.xaxis.tick_bottom()
    for prot, xl, yl in zip(Proteins, x, y):
        scatAx.annotate(str(prot[4:]), xy = (float(xl), float(yl)), xytext = (15,15), textcoords = 'offset points', arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    pylab.tight_layout()
    scatAx.plot(numpy.linspace(0, 10), numpy.linspace(0,10))
    return scatAx
    
##################Initialize Varibles###########################    
if __name__ == "__main__":
    vizLib.setRcs(scale=12)
    
    ##########GLOBAL VARIABLES##########
    SmallSubunit = ['BSubS02', 'BSubS03', 'BSubS04', 'BSubS05', 'BSubS06', 'BSubS07', 'BSubS08', 'BSubS09', 'BSubS10',
                    'BSubS11', 'BSubS12', 'BSubS13', 'BSubS14', 'BSubS15', 'BSubS16', 'BSubS17', 'BSubS18', 'BSubS20']
    LargeSubunit = ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL27', 'BSubL28', 'BSubL29', 'BSubL30', 'BSubL31a', 'BSubL32',
                    'BSubL33a', 'BSubL35', 'BSubL36']    
    Inter45S =      ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL29', 'BSubL30']
    AllSubunits = LargeSubunit + SmallSubunit
    
    doublingTime_10 = 92
    doublingTime_1000 = 47
    
    path = '/home/jhdavis/McMasterMSUDataSets/PulseLabeling/2013_05_28-MSUPulse/filtered/'
    
    figWidth = 11
    
    global currentPool
    currentPool = 0.0
   
   ##########Open Files##########
    files = qMS.sort_nicely([i.split('/')[-1] for i in glob.glob(path+'*.csv')])

    reds = ['#fee5d9', '#fcbba1', '#fc9272', '#fb6a4a', '#de2d26', '#a50f15']
    blues = ['#eff3ff', '#c6dbef', '#93cae1', '#6baed6', '#3182bd', '#08519c']
    
    names70S_10 = [path+n for n in files[0:6]]
    names70S_1000 = [path+n for n in files[6:12]]

    namesInter_10 = [path+n for n in files[12:17]]
    namesInter_1000 = [path+n for n in files[18:]]
    
    times10 = [22, 49, 66, 88, 110, 132]
    times1000 = [11, 22, 33, 44, 55, 66]
    
    labels10 = ['10 $\mu$M, 22 min', '10 $\mu$M, 49 min', '10 $\mu$M, 66 min', '10 $\mu$M, 88 min', '10 $\mu$M, 110 min', '10 $\mu$M, 132 min']
    labels1000 = ['1 mM, 11 min', '1 mM, 22 min', '1 mM, 33 min', '1 mM, 44 min', '1 mM, 55 min', '1 mM, 66 min']
    
    labels1000 = labels1000[0:3]
    times1000 = times1000[0:3]
    names70S_1000 = names70S_1000[0:3]
    namesInter_1000 = namesInter_1000[0:3]

##################Extract the data###########################
    num = ['AMP_L']
    den = ['AMP_U', 'AMP_L']
    dataByProtein70S = qMS.multiStatsDict(names70S_10+names70S_1000, num, den)
    dataByProteinInter = qMS.multiStatsDict(namesInter_10+namesInter_1000, num, den)

##################Find the medians###########################
    medDict70S_10 = {}
    medDict70S_1000 = {}
    
    medDictInter_10 = {}
    medDictInter_1000 = {}
    
    ntDict70S_10 = {'names': names70S_10, 'times': times10}
    ntDict70S_1000 = {'names': names70S_1000, 'times': times1000}
    
    ntDictInter_10 = {'names': namesInter_10, 'times': times10}
    ntDictInter_1000 = {'names': namesInter_1000, 'times': times1000}

    for prot in AllSubunits:
        medDict70S_10[prot] = [[ntDict70S_10['times'][i],numpy.median(dataByProtein70S[ntDict70S_10['names'][int(i)]][prot])] for i in [0,1,2,3,4,5] if len(dataByProtein70S[ntDict70S_10['names'][int(i)]][prot]) > 0]
        medDict70S_1000[prot] = [[ntDict70S_1000['times'][i],numpy.median(dataByProtein70S[ntDict70S_1000['names'][int(i)]][prot])] for i in [0,1,2] if len(dataByProtein70S[ntDict70S_1000['names'][int(i)]][prot]) > 0]
        
        medDictInter_10[prot] = [[ntDictInter_10['times'][i],numpy.median(dataByProteinInter[ntDictInter_10['names'][int(i)]][prot])] for i in [0,1,2,3,4] if len(dataByProteinInter[ntDictInter_10['names'][int(i)]][prot]) > 0]
        medDictInter_1000[prot] = [[ntDictInter_1000['times'][i],numpy.median(dataByProteinInter[ntDictInter_1000['names'][int(i)]][prot])] for i in [0,1,2] if len(dataByProteinInter[ntDictInter_1000['names'][int(i)]][prot]) > 0]

##################Plot the datasets###########################
    pylab.close('all')    
    '''
    yMax=1.5
    figSize=(figWidth,figWidth/3.0)
    median=False    
    
    num = ['AMP_L']
    den = ['AMP_U', 'AMP_L']
    filtPlots_70S_10 = plotDataSets(names70S_10, labels10, num, den, AllSubunits, 'non-permissive conditions : 70S\nfraction labeled', 'Labeled/[Labeled+Unlabeled]', reds,
                 saveName=None, yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=3)

    filtPlots_70S_1000 = plotDataSets(names70S_1000, labels1000, num, den, AllSubunits, 'permissive conditions : 70S\nfraction labeled', 'Labeled/[Labeled+Unlabeled]', blues,
                 yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=3)

    filtPlots_Inter_10 = plotDataSets(namesInter_10, labels10, num, den, AllSubunits, 'non-permissive conditions : 45S\nfraction labeled', 'Labeled/[Labeled+Unlabeled]', reds,
                 saveName=None, yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=3)

    filtPlots_Inter_1000 = plotDataSets(namesInter_1000, labels1000, num, den, AllSubunits, 'permissive conditions : 50S\nfraction labeled', 'Labeled/[Labeled+Unlabeled]', blues,
                 yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=3)
    
    num = ['AMP_L', 'AMP_U']
    den = ['AMP_U', 'AMP_L', 'AMP_S']
    filtPlots_70S_10 = plotDataSets(names70S_10, labels10, num, den, AllSubunits, 'non-permissive conditions : 70S\nprotein levels', '[Labeled+Unlabeled/[All]', reds,
                 saveName=None, yMax=yMax, figSize=figSize, median=median, legendLoc='lower left', legendCols=3, normProtein='BSubL24')

    filtPlots_70S_1000 = plotDataSets(names70S_1000, labels1000, num, den, AllSubunits, 'permissive conditions : 70S\nprotein levels', '[Labeled+Unlabeled/[All]', blues,
                 yMax=yMax, figSize=figSize, median=median, legendLoc='lower left', legendCols=3, normProtein='BSubL24')

    filtPlots_Inter_10 = plotDataSets(namesInter_10, labels10, num, den, AllSubunits, 'non-permissive conditions : 45S\nprotein levels', '[Labeled+Unlabeled/[All]', reds,
                 saveName=None, yMax=yMax, figSize=figSize, median=median, legendLoc='lower left', legendCols=3, normProtein='BSubL24')

    filtPlots_Inter_1000 = plotDataSets(namesInter_1000, labels1000, num, den, AllSubunits, 'permissive conditions : 50S\nprotein levels', '[Labeled+Unlabeled/[All]', blues,
                 yMax=yMax, figSize=figSize, median=median, legendLoc='lower left', legendCols=3, normProtein='BSubL24')
    
    '''
##################Fit data using pool size equation###########################
    poolSize70SDict = fitPoolSizes(AllSubunits, medDict70S_10, medDict70S_1000)
    poolSizeInterDict = fitPoolIntermediateSizes(Inter45S, medDictInter_10, medDictInter_1000)
    poolSizeInterFixFracXDict = fitPoolFixFracX(Inter45S, medDictInter_10)
    
    poolSize70S_10 = {key : poolSize70SDict[key]['10'] for key in poolSize70SDict.keys()}
    poolSizeInter_10 = {key : poolSizeInterDict[key]['10'] for key in poolSizeInterDict.keys()}    
    poolSizeInterFixFracX_10 = {key : poolSizeInterFixFracXDict[key]['10'] for key in poolSizeInterFixFracXDict.keys()}    
    
    fracX = calcFractionIntermediate(i, poolSizeInter_10.keys(), poolSizeInter_10, poolSize70S_10, qMS.growthRate(doublingTime_10))
    plotPoolBar(fracX, ['fracX'], ['r'], 'calculated fraction intermediate', figname=None, barNames=poolSizeInter_10.keys(), sortingKey=1)

    #print poolSizeInterFixFracXDict
##################Plot the fits###########################
    figSize=(7,10)
    #saveFile='/home/jhdavis/Dropbox/page.png'
    saveFile=None
    '''
    for i in [1,2,3]:
        plotPoolPage(LargeSubunit[(i-1)*10:i*10], poolSize70SDict, medDict70S_10, medDict70S_1000, saveFile=saveFile, figSize=figSize, yMax=1.0, funcToFit=qMS.poolInterFracXFunc, title='45S labeling fitFracX', double=False, showPools=False)
    
    for i in [1,2]:
        plotPoolPage(SmallSubunit[(i-1)*10:i*10], poolSize70SDict, medDict70S_10, medDict70S_1000, saveFile=saveFile, figSize=figSize, yMax=1.0, title='70S labeling')
    
    for i in [1,2,3]:
        plotPoolPage(Inter45S[(i-1)*10:i*10], poolSize70SDict, medDictInter_10, medDictInter_1000, saveFile=saveFile, figSize=figSize, yMax=1.0, funcToFit=qMS.poolInterFracXFunc, title='Intermediate labeling fit FracX', double=False, showPools=False)
    '''
    for i in [1,2,3]:
        plotPoolPage(Inter45S[(i-1)*10:i*10], poolSizeInterFixFracXDict, medDictInter_10, medDictInter_1000, saveFile=saveFile, figSize=figSize, yMax=1.0, funcToFit=qMS.poolInterFracXFunc, title='Intermediate labeling fixed FracX', double=False, showPools=True)
    
##################Make bar graphs###########################
    '''
    plotPoolBar(poolSize70SDict, ['10'], ['r'], '10 $\mu$M', figname=None, barNames=AllSubunits, sortingKey=1)
    #plotPoolBar(poolSize70SDict, ['1000'], ['b'], '1 mM', figname=None, barNames=AllSubunits, sortingKey=1)
    plotPoolBar(poolSize70SDict, ['10', '1000'], ['r', 'b'], '1 mM', figname=None, barNames=LargeSubunit, sortingKey=1)
    plotPoolBar(poolSizeInterDict, ['10', '1000'], ['r', 'b'], '1 mM', figname=None, barNames=Inter45S, sortingKey=1)
    '''
##################Read protein inventory data###########################
    '''
    path = '/home/jhdavis/data/originalMSUProteinLevels/McMasterMSUCsvs/'
    McMaster45S = [path+i for i in ["McMaster45S_esi-run1_filt.csv", "McMaster45S_esi-run2.1_filt.csv", "McMaster45S_esi-run2_filt.csv", "McMaster45S_qtof_filt_filtppm.csv"]]
    McMaster50S = [path+i for i in ["McMaster50S_esi-run1_filt.csv", "McMaster50S_esi-run2.1_filt.csv", "McMaster50S_esi-run2_filt.csv", "McMaster50S_qtof_filt_filtppm.csv"]]
    
    num = ['AMP_U']
    den = ['AMP_U', 'AMP_S']
    normProtein = 'BSubL24'

    fileLists = [McMaster45S, McMaster50S]
    merged45 = qMS.mergeFiles(McMaster45S, num, den, normProtein)
    merged50 = qMS.mergeFiles(McMaster50S, num, den, normProtein)
    verifiedZero = ['BSubL16', 'BSubL28', 'BSubL36', 'BSubL31a']    
    for i in verifiedZero:
        merged45[i] = numpy.array([0.0])
    '''
    
##################Plot protein inventory data###########################
    '''
    myPlot = vizLib.plotStatsDict(merged45, name='45SMerged', proteins=LargeSubunit, offset=0.4, markerSize=12, color='#e31a1c', yMax = 1.5, median=False)
    myPlot = vizLib.addStatsDictToPlot(merged50, myPlot, name='50SMerged', offset=0.6, markerSize=12, color='#377db8', median=False)
    myPlot.set_ylabel('protein occupancy\nnormalized to L24', multialignment='center')
    myPlot.set_title('protein occupancy 45S vs. 50S')
    pylab.legend(loc='lower left', prop={'size':12})
    pylab.tight_layout()
    '''
##################Plot pool data vs. protein inventory data###########################
    '''
    poolSize70S_10 = {key : poolSize70SDict[key]['10'] for key in poolSize70SDict.keys()}
    scatAx = plotCompData(poolSize70S_10, merged45, LargeSubunit, title='test', xlabel='precursor pool size (P) from 70S measurement',
                 ylabel='relative abundance in 45S particle', xMax=1.25, yMax=1.25, saveFile=None)
    scatAx.plot([-2, 2], [0.5, 0.5], color='g', linestyle='--')
    scatAx.plot([0.18, 0.18], [-2, 2], color='g', linestyle='--')
    scatAx.text(0.14, 0.45,'I', size=20, horizontalalignment='right', verticalalignment='top', weight='bold')
    scatAx.text(0.14, 0.55,'II', size=20, horizontalalignment='right', verticalalignment='bottom', weight='bold')
    scatAx.text(0.22, 0.55,'III', horizontalalignment='left', size=20, verticalalignment='bottom', weight='bold')
    scatAx.text(0.22, 0.45,'IV', horizontalalignment='left', verticalalignment='top', size=20, weight='bold')
    
    textstr1 =  'I : absent from 45S, small pools\n' + \
                'II: present in 45S, small pools\n' + \
                'III: present in 45S, large pools\n' + \
                'IV: absent from 45S, large pools'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    scatAx.text(0.5, .4, textstr1, fontsize=15, verticalalignment='top', horizontalalignment='left', color='black', bbox=props)
    '''
##################Plot pool data vs. overLabelingPoolSizeData###########################
    
    poolSizeInter_10 = {key : poolSizeInterDict[key]['10'] for key in poolSizeInterDict.keys()}
    scatAx = plotCompData(poolSize70S_10, poolSizeInter_10, Inter45S, title='45S Pool measurements', xlabel='precursor pool size (P) from 70S measurement',
                 ylabel='precursor pool size (P) from 45S measurement', xMax=1.25, yMax=6.0, saveFile=None)
    scatAx = plotCompData(poolSize70S_10, poolSizeInterFixFracX_10, Inter45S, title='45S Pool measurements', xlabel='precursor pool size (P) from 70S measurement',
                 ylabel='precursor pool size (P) from 45S measurement_corrected', xMax=1.25, yMax=1.25, saveFile=None)
    
##################Calculate fraction intermediate###########################

    pylab.show('all')