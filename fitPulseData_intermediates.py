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
    #k = qMS.growthRate(47)
    k = qMS.growthRate(doublingTime_1000)
    return y - qMS.poolFunc(k, t, P)
def poolResid10(P, y, t):
    #k = qMS.growthRate(92)
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
    #k = qMS.growthRate(47)
    k = qMS.growthRate(doublingTime_1000)
    return y - qMS.poolOverFunc(k, t, P)

def poolIntermediatesResid10(P, y, t):
    #k = qMS.growthRate(92)
    k = qMS.growthRate(doublingTime_10)
    return y - qMS.poolOverFunc(k, t, P)

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

def plotLabelKinsPage(Proteins, protPoolDict, medDict10, medDict1000, fig, name):
    time = numpy.linspace(0, 140, 200)
    k10 = qMS.growthRate(92)
    k1000 = qMS.growthRate(47)
    
    Large1Ax = []
    a = 0
    for prot in Proteins:

        p_10 = protPoolDict[prot]['10']
        p_10String = qMS.calcPercent(p_10, sigfig=2).split('.')[0]+'%'
        #p_1000 = protPoolDict[prot]['1000']
        #p_1000String = qMS.calcPercent(p_1000, sigfig=2).split('.')[0]+'%'

        ax = fig.add_subplot(5,2,a+1)
        
        #textstr1 = '1 mM : P=' + p_1000String
        textstr2 = '\n\n'+'10 $\mu$M : P=' + p_10String
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.025, .95, textstr2, transform=ax.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', bbox=props, color='r')
        #ax.text(0.025, .925, textstr1, transform=ax.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='b')
        
        ax.scatter([i[0] for i in medDict10[prot]], [i[1] for i in medDict10[prot]], c = 'r')
        ax.scatter([i[0] for i in medDict1000[prot]], [i[1] for i in medDict1000[prot]], c = 'b')
        ax.plot(time, qMS.maxLabFunc(k10, time), c='r')
        ax.plot(time, qMS.maxLabFunc(k1000, time), c='b')
        
        ax.plot(time, qMS.poolOverFunc(k10, time, p_10), c='r', ls='--')
        #ax.plot(time, qMS.poolOverFunc(k1000, time, p_1000), c='b', ls='--')
        
        ax.set_xlim([0,140])
        ax.set_ylim([0,1.00])
        ax.set_yticks([0,.25,.50,.75,1.00])
        #if a < 7:
        #    plt.setp(ax.get_xticklabels(), visible=False)
        #if a%2==1:
        #    plt.setp(ax.get_yticklabels(), visible=False)
        ax.set_xticks([0,20,40,60,80,100,120,140])
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
    fig.tight_layout()
    #pylab.savefig(name)
    return Large1Ax

        
def plotPoolPage(Proteins, protPoolDict, medDict10, medDict1000, fig, name):
    time = numpy.linspace(0, 140, 200)
    k10 = qMS.growthRate(92)
    k1000 = qMS.growthRate(47)
    
    Large1Ax = []
    a = 0
    for prot in Proteins:
        p_10 = protPoolDict[prot]['10']
        p_10String = qMS.calcPercent(p_10, sigfig=2).split('.')[0]+'%'
        p_1000 = protPoolDict[prot]['1000']
        p_1000String = qMS.calcPercent(p_1000, sigfig=2).split('.')[0]+'%'

        ax = fig.add_subplot(5,2,a+1)
        textstr1 = '1 mM : P=' + p_1000String
        textstr2 = '\n\n'+'10 $\mu$M : P=' + p_10String
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.025, .95, textstr2, transform=ax.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', bbox=props, color='r')
        ax.text(0.025, .925, textstr1, transform=ax.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='b')
        
        ax.scatter([i[0] for i in medDict10[prot]], [i[1] for i in medDict10[prot]], c = 'r')
        ax.scatter([i[0] for i in medDict1000[prot]], [i[1] for i in medDict1000[prot]], c = 'b')
        ax.plot(time, qMS.maxLabFunc(k10, time), c='r')
        ax.plot(time, qMS.maxLabFunc(k1000, time), c='b')
        
        ax.plot(time, qMS.poolFunc(k10, time, p_10), c='r', ls='--')
        ax.plot(time, qMS.poolFunc(k1000, time, p_1000), c='b', ls='--')
        ax.set_xlim([0,140])
        ax.set_ylim([0,1.00])
        ax.set_yticks([0,.25,.50,.75,1.00])
        #if a < 7:
        #    plt.setp(ax.get_xticklabels(), visible=False)
        #if a%2==1:
        #    plt.setp(ax.get_yticklabels(), visible=False)
        ax.set_xticks([0,20,40,60,80,100,120,140])
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
    fig.tight_layout()
    #pylab.savefig(name)
    return Large1Ax
    
def plotPoolBar(data, keys, colors, title, figname=None, barNames=None, sortingKey=1):
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
    bar = pylab.figure(figsize=(20,4))
    barax = bar.add_subplot(111)
    bars = len(poolSizeDict.keys())
    ind = numpy.arange(bars)
    width = 0.9/len(keys)
    
    sortedTuples = []
    if barNames is None:
        barNames = data.keys()
    barNames.sort()
    for name in barNames:
        myTuple = [name]
        valueIds = data[name].keys()
        valueIds.sort()
        for k in valueIds:
            myTuple.append(data[name][k])
        sortedTuples.append(myTuple)
    sortedTuples = sorted(sortedTuples, key=lambda pool: pool[sortingKey], reverse=True)
    #sortedTuples = sorted([[name, poolSizeDict[name][key]] for name in poolSizeDict.keys()], key=lambda pool: pool[1], reverse=True)
    #sortedTuplesk1 = sorted([[name, poolSizeDict[name][key1], poolSizeDict[name][key2]] for name in Proteins], key=lambda pool: pool[1], reverse=True)
    names = [i[0] for i in sortedTuples]
    rects = []
    for b in len(keys):
        vs = [i[b+1] for i in sortedTuples]
        rects.append(barax.bar(ind+0.05+width, vs, width, color=colors[b]))
    
    barax.set_xticks(ind+0.5)
    barax.set_xticklabels([i[4:] for i in names], rotation='vertical', horizontalalignment='center', verticalalignment='top')    
    barax.set_ylabel('fit pool size (P)')
    barax.tick_params(axis='y', direction='out')
    barax.tick_params(axis='x', direction='out')
    barax.xaxis.tick_bottom()
    barax.yaxis.tick_left()
    barax.set_ylim([0,6.5])
    barax.set_title(title)
    barax.legend(rects, barNames)
    bar.tight_layout()
    if not (figname is None):
        pylab.savefig(figname)
    return bar

'''
def plotTwoPoolBars(poolSizeDict, Proteins, key1, color1, title1, key2, color2, title2, figname):
    bar = pylab.figure(figsize=(20,4))
    barax = bar.add_subplot(111)
    bars = len(Proteins)
    ind = numpy.arange(bars)
    width = 0.45
    
    sortedTuplesk1 = sorted([[name, poolSizeDict[name][key1], poolSizeDict[name][key2]] for name in Proteins], key=lambda pool: pool[1], reverse=True)
    names = [i[0] for i in sortedTuplesk1]
    vsk1 = [i[1] for i in sortedTuplesk1]
    vsk2 = [i[2] for i in sortedTuplesk1]
    rects1 = barax.bar(ind+0.05, vsk1, width, color=color1)
    rects2 = barax.bar(ind++0.05+width, vsk2, width, color=color2)

    barax.set_xticks(ind+0.5)
    barax.set_xticklabels([i[4:] for i in names], rotation='vertical', horizontalalignment='center', verticalalignment='top')

    barax.set_ylabel('fit pool size (P)')
    barax.tick_params(axis='y', direction='out')
    barax.tick_params(axis='x', direction='out')
    barax.xaxis.tick_bottom()
    barax.yaxis.tick_left()
    barax.legend((rects1[0], rects2[0]), (title1, title2))
    barax.set_title(figname.split('.')[0])
    bar.tight_layout()
    pylab.savefig(figname)
    return bar
'''

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
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL29', 'BSubL30']
    AllSubunits = LargeSubunit + SmallSubunit
    
    doublingTime_10 = 92
    doublingTime_1000 = 47
    
    path = '/home/jhdavis/data/2013_05_28-MSUPulse/filtered/'
   
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
    
    #labels1000 = labels1000[0:3]
    #times1000 = times1000[0:3]
    #names1000 = names1000[0:3]

##################Extract the data###########################
    num = ['AMP_L']
    den = ['AMP_U', 'AMP_L']
    dataByProtein70S = qMS.multiStatsDict(names70S_10+names70S_1000, num, den)
    dataByProteinInter = qMS.multiStatsDict(namesInter_10+namesInter_1000, num, den)
##################Plot the datasets###########################
    pylab.close('all')    
    '''
    yMax=1.5
    figSize=(22,7)
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
                 saveName=None, yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=3, normProtein='BSubL24')

    filtPlots_70S_1000 = plotDataSets(names70S_1000, labels1000, num, den, AllSubunits, 'permissive conditions : 70S\nprotein levels', '[Labeled+Unlabeled/[All]', blues,
                 yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=3, normProtein='BSubL24')

    filtPlots_Inter_10 = plotDataSets(namesInter_10, labels10, num, den, AllSubunits, 'non-permissive conditions : 45S\nprotein levels', '[Labeled+Unlabeled/[All]', reds,
                 saveName=None, yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=3, normProtein='BSubL24')

    filtPlots_Inter_1000 = plotDataSets(namesInter_1000, labels1000, num, den, AllSubunits, 'permissive conditions : 50S\nprotein levels', '[Labeled+Unlabeled/[All]', blues,
                 yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=3, normProtein='BSubL24')
    '''
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
    
##################Fit data using pool size equation###########################
    '''
    poolSizeDict = fitPoolOverSizes(Inter45S, medDict10, medDict1000)
    
    f1 = 'muspulse1p2_iso_res_filt.csv'
    f2 = 'muspulse2_iso_res_filt.csv'
    f3 = 'muspulse3_iso_res_filt.csv'
    f4 = 'muspulse4_iso_res_filt.csv'
    f5 = 'muspulse5_iso_res_filt.csv'
    f6 = 'muspulse6_iso_res_filt.csv'
    
    f7 = 'muspulse7_iso_res_filt.csv'
    f8 = 'muspulse8_iso_res_filt.csv'
    f9 = 'muspulse9_iso_res_filt.csv'
    f10 = 'muspulse10_iso_res_filt.csv'
    f11 = 'muspulse11_iso_res_filt.csv'
    f12 = 'muspulse12_iso_res_filt.csv'
    
    names10 = [path+n for n in [f1, f2, f3, f4, f5, f6]]
    names1000 = [path+n for n in [f7, f8, f9, f10, f11, f12]]
    
    csvs = {}
    dataByProtein = {}
    allProteins = set()
    dataByProtein = qMS.multiStatsDict(names10+names1000, num, den)
    
    times1000 = [11, 22, 33, 44, 55, 66]
    times10 = [22, 49, 66, 88, 110, 132]
    
    ntDict10 = {'names': names10, 'times': times10}
    ntDict1000 = {'names': names1000, 'times': times1000}

    labels10 = ['10 $\mu$M, 22 min', '10 $\mu$M, 49 min', '10 $\mu$M, 66 min', '10 $\mu$M, 88 min', '10 $\mu$M, 110 min', '10 $\mu$M, 132 min']
    labels1000 = ['1 mM, 11 min', '1 mM, 22 min', '1 mM, 33 min', '1 mM, 44 min', '1 mM, 55 min', '1 mM, 66 min']
    
    medDict10 = {}
    medDict1000 = {}
    
    
    labels1000 = labels1000[0:3]
    times1000 = times1000[0:3]
    names1000 = names1000[0:3]
    
 
    num = ['AMP_L']
    den = ['AMP_U', 'AMP_L']
    
    for prot in AllSubunits:
        medDict10[prot] = [[ntDict10['times'][i],numpy.median(dataByProtein[ntDict10['names'][int(i)]][prot])] for i in [0,1,2,3,4] if len(dataByProtein[ntDict10['names'][int(i)]][prot]) > 0]
        #medDict10[prot] = [[ntDict10['times'][i],numpy.median(dataByProtein[ntDict10['names'][int(i)]][prot])] for i in range(6,12) if len(dataByProtein[ntDict10['names'][int(i)]][prot]) > 0]
        #medDict1000[prot] = [[ntDict1000['times'][i],numpy.median(dataByProtein[ntDict1000['names'][int(i)]][prot])] for i in range(6,12) if len(dataByProtein[ntDict1000['names'][int(i)]][prot]) > 0]
        #medDict1000[prot] = [[ntDict1000['times'][i],numpy.median(dataByProtein[ntDict1000['names'][int(i)]][prot])] for i in [0,1,2,3,4,5] if len(dataByProtein[ntDict1000['names'][int(i)]][prot]) > 0]
        medDict1000[prot] = [[ntDict1000['times'][i],numpy.median(dataByProtein[ntDict1000['names'][int(i)]][prot])] for i in [0,1,2] if len(dataByProtein[ntDict1000['names'][int(i)]][prot]) > 0]
        #medDict1000[prot] = [[ntDict1000['times'][i],numpy.median(dataByProtein[ntDict1000['names'][int(i)]][prot])] for i in [0,1,2] if len(dataByProtein[ntDict1000['names'][int(i)]][prot]) > 0]
    
    
    
    poolSize70SDict = fitPoolSizes(Inter45S, medDict10, medDict1000)
    print poolSize70SDict
    '''
##################Plot the fits###########################
    '''
    ext = '.png'
    fLarge1 = pylab.figure(figsize=(7,10))
    f1axs = plotPoolPage(LargeSubunit[0:10], poolSizeDict, medDict10, medDict1000, fLarge1, '50S-1'+ext)
    
    fLarge2 = pylab.figure(figsize=(7,10))
    f1axs = plotPoolPage(LargeSubunit[10:20], poolSizeDict, medDict10, medDict1000, fLarge2, '50S-2'+ext)    
    
    fLarge3 = pylab.figure(figsize=(7,10))
    f3axs = plotPoolPage(LargeSubunit[20:30], poolSizeDict, medDict10, medDict1000, fLarge3, '50S-3'+ext)
    
    fSmall1 = pylab.figure(figsize=(7,10))
    f4axs = plotPoolPage(SmallSubunit[0:10], poolSizeDict, medDict10, medDict1000, fSmall1, '30S-1'+ext)

    fSmall2 = pylab.figure(figsize=(7,10))
    f5axs =plotPoolPage(SmallSubunit[10:], poolSizeDict, medDict10, medDict1000, fSmall2, '30S-2'+ext)
    
    ext = '.pdf'
    fLarge1 = pylab.figure(figsize=(7,10))
    f1axs = plotLabelKinsPage(Inter45S[0:10], poolSizeDict, medDict10, medDict1000, fLarge1, '50S-1'+ext)
    
    fLarge2 = pylab.figure(figsize=(7,10))
    f1axs = plotLabelKinsPage(Inter45S[10:20], poolSizeDict, medDict10, medDict1000, fLarge2, '50S-2'+ext)    
    
    fLarge3 = pylab.figure(figsize=(7,10))
    f3axs = plotLabelKinsPage(Inter45S[20:30], poolSizeDict, medDict10, medDict1000, fLarge3, '50S-3'+ext)
    
    fSmall1 = pylab.figure(figsize=(7,10))
    f4axs = plotLabelKinsPage(SmallSubunit[0:10], medDict10, medDict1000, fSmall1, '30S-1'+ext)

    fSmall2 = pylab.figure(figsize=(7,10))
    f5axs =plotLabelKinsPage(SmallSubunit[10:], medDict10, medDict1000, fSmall2, '30S-2'+ext)
    '''
    
    
    
##################Make bar graphs###########################
    
    #plotPoolBar(poolSizeDict, '10', 'r', '10 $\mu$M', '10uM.png')
    #plotPoolBar(poolSizeDict, '1000', 'b', '1 mM', '1000uM.pdf')
    #plotTwoPoolBars(poolSizeDict, LargeSubunit, '10', 'r', '10 $\mu$M', '1000', 'b', '1 mM', 'LargeSubunit.pdf')
    #plotTwoPoolBars(poolSizeDict, SmallSubunit, '10', 'r', '10 $\mu$M', '1000', 'b', '1 mM', 'SmallSubunit.pdf')
    
##################Read protein inventory data###########################
    '''
    proteinToNormalizeTo = "BSubL24"

    LargeSubunit = ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL27', 'BSubL28', 'BSubL29', 'BSubL30', 'BSubL32',
                    'BSubL33a', 'BSubL35', 'BSubL36']
    path = '/home/jhdavis/data/originalMSUProteinLevels/McMasterMSUCsvs/'
    McMaster45S = [path+i for i in ["McMaster45S_esi-run1_filt.csv", "McMaster45S_esi-run2.1_filt.csv", "McMaster45S_esi-run2_filt.csv", "McMaster45S_qtof_filt_filtppm.csv"]]
    McMaster50S = [path+i for i in ["McMaster50S_esi-run1_filt.csv", "McMaster50S_esi-run2.1_filt.csv", "McMaster50S_esi-run2_filt.csv", "McMaster50S_qtof_filt_filtppm.csv"]]
    
    num = ['AMP_U']
    den = ['AMP_U', 'AMP_S']
    normProtein = 'BSubL24'

    fileLists = [McMaster45S, McMaster50S]
    merged45 = qMS.mergeFiles(McMaster45S, num, den, normProtein)
    merged50 = qMS.mergeFiles(McMaster50S, num, den, normProtein)
    
    
    myPlot = vizLib.plotStatsDict(merged45, name='45SMerged', proteins=LargeSubunit, offset=0.4, markerSize=12, color='#e31a1c', yMax = 1.5, median=False)
    myPlot = vizLib.addStatsDictToPlot(merged50, myPlot, name='50SMerged', offset=0.6, markerSize=12, color='#377db8', median=False)
    myPlot.set_ylabel('protein occupancy\nnormalized to L24', multialignment='center')
    myPlot.set_title('protein occupancy 45S vs. 50S')
    #pylab.xticks(numpy.arange(1.5,len(LargeSubunit)+1,1), [item[4:] for item in LargeSubunit], rotation=45)
    pylab.legend(loc='lower left', prop={'size':12})
    pylab.tight_layout()
    '''
    '''
    myPlot = qMS.makePlotWithDataSets(merged, LargeSubunit, ["McMaster45S_merged", 'McMaster50S_merged'])
    for i in LargeSubunit:
        pVal = stats.ttest_ind(merged[0][i]['vals'], merged[1][i]['vals'], equal_var=False)
    print merged[0]['BSubL02']['vals']    
    print merged[1]['BSubL20']['vals']
    L30ttest = stats.ttest_ind(merged[0]['BSubL02']['vals'], merged[1]['BSubL02']['vals'])
    L12ttest = stats.ttest_ind(merged[0]['BSubL34']['vals'], merged[1]['BSubL34']['vals'])
    print L30ttest
    print L12ttest
    '''
    
##################Plot pool data vs. protein inventory data###########################
    '''
    verifiedZero = ['BSubL16', 'BSubL28', 'BSubL36']    
    for z in verifiedZero:    
        merged45[z] = numpy.array([0.0])
        
    poolSizeVsPIDict = {}
    for prot in LargeSubunit:
        poolSizeVsPIDict[prot] = {'PI':numpy.median(merged45[prot]), 'PS':poolSizeDict[prot]['10']}
    names = poolSizeVsPIDict.keys()
    names.sort()
    #x = [numpy.log2(poolSizeVsPIDict[i]['PS']+1) for i in names]
    x = [poolSizeVsPIDict[i]['PS'] for i in names]
    y = [poolSizeVsPIDict[i]['PI'] for i in names]
    scat = pylab.figure(figsize=(10,10))
    scatAx = scat.add_subplot(111)    
    scatAx.scatter(x,y, c='b', s=150)
    scatAx.set_title('protein abundance vs. pool size under RbgA-limiting conditions')
    scatAx.set_xlabel('precursor pool size (P)')
    scatAx.set_ylabel('relative abundance in 45S particle')
    scatAx.set_xlim([-0.1,1.25])
    scatAx.set_ylim([-0.1,1.25])
    scatAx.set_yticks([0,.25,.50,.75, 1.0, 1.25])
    scatAx.set_xticks([0,0.25,0.5, 0.75, 1, 1.25])
    scatAx.plot([-2, 2], [0.5, 0.5], color='g', linestyle='--')
    scatAx.plot([0.18, 0.18], [-2, 2], color='g', linestyle='--')
    scatAx.text(0.14, 0.45,'I', size=20, horizontalalignment='right', verticalalignment='top', weight='bold')
    scatAx.text(0.14, 0.55,'II', size=20, horizontalalignment='right', verticalalignment='bottom', weight='bold')
    scatAx.text(0.22, 0.55,'III', horizontalalignment='left', size=20, verticalalignment='bottom', weight='bold')
    scatAx.text(0.22, 0.45,'IV', horizontalalignment='left', verticalalignment='top', size=20, weight='bold')
    #scatAx.text('I', 0,5, 0.5)
    scatAx.yaxis.tick_left()
    scatAx.xaxis.tick_bottom()
    for prot, xl, yl in zip(names, x, y):
        scatAx.annotate(str(prot[4:]), xy = (float(xl), float(yl)), xytext = (15,15), textcoords = 'offset points', arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    
    textstr1 =  'I : absent from 45S, small pools\n' + \
                'II: present in 45S, small pools\n' + \
                'III: present in 45S, large pools\n' + \
                'IV: absent from 45S, large pools'
    #textstr2 = '\n\n'+'10 $\mu$M : P=' + p_10String
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    scatAx.text(0.5, .4, textstr1, fontsize=15, verticalalignment='top', horizontalalignment='left', color='black', bbox=props)
    pylab.tight_layout()
    '''
##################Plot pool data vs. overLabelingPoolSizeData###########################
    '''
    poolSizeVsPIDict = {}
    for prot in Inter45S:
        poolSizeVsPIDict[prot] = {'PI':poolSize70SDict[prot]['10'], 'PS':poolSizeDict[prot]['10']}
    names = poolSizeVsPIDict.keys()
    names.sort()
    #x = [numpy.log2(poolSizeVsPIDict[i]['PS']+1) for i in names]
    x = [poolSizeVsPIDict[i]['PS'] for i in names]
    y = [poolSizeVsPIDict[i]['PI'] for i in names]
    scat = pylab.figure(figsize=(10,10))
    scatAx = scat.add_subplot(111)    
    scatAx.scatter(x,y, c='b', s=150)
    scatAx.set_title('protein abundance vs. pool size under RbgA-limiting conditions')
    scatAx.set_xlabel('apparent pool using 45S overlabeling')
    scatAx.set_ylabel('apparent pool measuring 70S lag')
    scatAx.set_xlim([-0.1,7])
    scatAx.set_ylim([-0.1,1.5])
    #scatAx.set_yticks([0,.25,.50,.75, 1.0, 1.25])
    #scatAx.set_xticks([0,0.25,0.5, 0.75, 1, 1.25])
    #scatAx.plot([-2, 2], [0.5, 0.5], color='g', linestyle='--')
    #scatAx.plot([0.18, 0.18], [-2, 2], color='g', linestyle='--')
    #scatAx.text(0.14, 0.45,'I', size=20, horizontalalignment='right', verticalalignment='top', weight='bold')
    #scatAx.text(0.14, 0.55,'II', size=20, horizontalalignment='right', verticalalignment='bottom', weight='bold')
    #scatAx.text(0.22, 0.55,'III', horizontalalignment='left', size=20, verticalalignment='bottom', weight='bold')
    #scatAx.text(0.22, 0.45,'IV', horizontalalignment='left', verticalalignment='top', size=20, weight='bold')
    #scatAx.text('I', 0,5, 0.5)
    scatAx.yaxis.tick_left()
    scatAx.xaxis.tick_bottom()
    for prot, xl, yl in zip(names, x, y):
        scatAx.annotate(str(prot[4:]), xy = (float(xl), float(yl)), xytext = (15,15), textcoords = 'offset points', arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    
    #textstr1 =  'I : absent from 45S, small pools\n' + \
    #            'II: present in 45S, small pools\n' + \
    #            'III: present in 45S, large pools\n' + \
    #            'IV: absent from 45S, large pools'
    #textstr2 = '\n\n'+'10 $\mu$M : P=' + p_10String
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #scatAx.text(0.5, .4, textstr1, fontsize=15, verticalalignment='top', horizontalalignment='left', color='black', bbox=props)
    pylab.tight_layout()
    
    '''
    pylab.show('all')