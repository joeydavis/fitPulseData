# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 23:13:53 2013

@author: jhdavis
"""

import pylab
from scipy import stats
import numpy
import qMS
from scipy import optimize

def getAllCsvData(fileName, pulse=True):
    d = dict()
    [dd, prs, pel] = qMS.readCSV(fileName, pulse)
    d['data'] = dd
    d['proteinSet'] = prs
    d['peptideList'] = pel
    return d

def getMedian(dList):
    return numpy.median(dList)

def poolResid1000(P, y, t):
    k = qMS.growthRate(47)
    return y - qMS.poolFunc(k, t, P)

def poolResid10(P, y, t):
    k = qMS.growthRate(92)
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
        
def plotPoolPage(Proteins, protPoolDict, medDict10, medDict1000, fig, name):
    time = numpy.linspace(0, 80, 100)
    k10 = qMS.growthRate(80)
    k1000 = qMS.growthRate(42)
    
    Large1Ax = []
    a=0
    for prot in Proteins:
        ax = fig.add_subplot(5,2,a+1)
        ax.scatter([i[0] for i in medDict10[prot]], [i[1] for i in medDict10[prot]], c = 'r')
        ax.scatter([i[0] for i in medDict1000[prot]], [i[1] for i in medDict1000[prot]], c = 'b')
        ax.plot(time, qMS.maxLabFunc(k10, time), c='r')
        ax.plot(time, qMS.maxLabFunc(k1000, time), c='b')
        
        p_10 = protPoolDict[prot]['10']
        p_10String = qMS.calcPercent(p_10, sigfig=3)
        p_1000 = protPoolDict[prot]['1000']
        p_1000String = qMS.calcPercent(p_1000, sigfig=3)

        ax.plot(time, qMS.poolFunc(k10, time, p_10), c='r', ls='--')
        ax.plot(time, qMS.poolFunc(k1000, time, p_1000), c='b', ls='--')
        ax.set_xlim([0,80])
        ax.set_ylim([0,0.75])
        ax.set_yticks([0,.25,.50,.75])
        ax.set_xticks([0,20,40,60,80])
        ax.yaxis.tick_left()
        ax.xaxis.tick_bottom()

        textstr1 = '1 mM : P=' + p_1000String
        textstr2 = '\n\n'+'10 $\mu$M : P=' + p_10String
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.025, .95, textstr2, transform=ax.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', bbox=props, color='r')
        ax.text(0.025, .925, textstr1, transform=ax.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='b')
        
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
    #fig.savefig(name)
    return Large1Ax
        
def plotPoolBar(poolSizeDict, key, color, title):
    bar = pylab.figure(figsize=(20,4))
    barax = bar.add_subplot(111)
    bars = len(poolSizeDict.keys())
    ind = numpy.arange(bars)
    width = 0.9
    sortedTuples = sorted([[name, poolSizeDict[name][key]] for name in poolSizeDict.keys()], key=lambda pool: pool[1], reverse=True)
    names = [i[0] for i in sortedTuples]
    vs = [i[1] for i in sortedTuples]
    barax.bar(ind, vs, width, color=color)
    barax.set_xticks(ind+0.5)
    barax.set_xticklabels([i[4:] for i in names], rotation='vertical', horizontalalignment='center', verticalalignment='top')
    barax.set_ylabel('fit pool size (P)')
    barax.tick_params(axis='y', direction='out')
    barax.tick_params(axis='x', direction='out')
    barax.xaxis.tick_bottom()
    barax.yaxis.tick_left()
    barax.set_title(title)
    bar.tight_layout()
    return bar

def plotTwoPoolBars(poolSizeDict, Proteins, key1, color1, title1, key2, color2, title2):
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
#    barax.set_title(title)
    bar.tight_layout()
    return bar

##################Initialize Varibles###########################    
if __name__ == "__main__":
    SmallSubunit = ['BSubS02', 'BSubS03', 'BSubS04', 'BSubS05', 'BSubS06', 'BSubS07', 'BSubS08', 'BSubS09', 'BSubS10',
                    'BSubS11', 'BSubS12', 'BSubS13', 'BSubS14', 'BSubS15', 'BSubS16', 'BSubS17', 'BSubS18', 'BSubS20']

    LargeSubunit = ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL27', 'BSubL28', 'BSubL29', 'BSubL30', 'BSubL31a', 'BSubL32',
                    'BSubL33a', 'BSubL35', 'BSubL36']

    AllSubunits = LargeSubunit + SmallSubunit
    
    path = '/home/jhdavis/data/2013_05_28-MSUPulse/filtered/'

    reds = ['#ae2221', '#d72c2b', '#e78180']
    blues = ['#25557d', '#3170a4', '#5696cc']

    d80t26 = 'msuPulseA_iso_filt_filt.csv'
    d80t52 = 'msuPulseB_iso_filt_filt.csv'
    d80t78 = 'msuPulseC_iso_filt_filt.csv'

    d42t14 = 'msuPulseD_iso_filt.csv'
    d42t28 = 'msuPulseE_iso_filt.csv'
    d42t42 = 'msuPulseF_iso_filt.csv'

    names1000 = [d42t14, d42t28, d42t42]
    names10 = [d80t26, d80t52, d80t78]
    
    times1000 = [14, 28, 42]
    times10 = [26, 52, 78]
    
    ntDict10 = {'names': names10, 'times': times10}
    ntDict1000 = {'names': names1000, 'times': times1000}

    labels10 = ['10 uM, 26 min', '10 uM, 52 min', '10 uM, 78 min']
    labels1000 = ['1 mM 14 min', '1 mM 28 min', '1 mM 42 min']
    
    medDict10 = {}
    medDict1000 = {}

##################Extract the data###########################
    csvs = {}
    dataByProtein = {}
    allProteins = set()
    for i in names10+names1000:
        csvs[i] = getAllCsvData(path+i)
        [allProteins.add(p) for p in csvs[i]['proteinSet']]
        dataByProtein[i] = qMS.calculateStatsFile(csvs[i]['data'], csvs[i]['proteinSet'], ['ampl'], ['ampl', 'ampu'], normalization=1.0)

##################Plot the datasets###########################
    
    pylab.close('all')
    '''
    ax = qMS.plotStatsDataStruct(AllSubunits, dataByProtein[names1000[0]], labels1000[0], offset = 1.0/float(len(names1000)+1), markerSize = 5/float(len(names1000))+4, color=blues[0])
    j = 1
    for i in names1000[1:]:
        ax = qMS.addStatsDataStructToPlot(AllSubunits, dataByProtein[i], ax, labels1000[j], offset = (1.0/float(len(names1000)+1))*(j+1), markerSize = 5/float(len(names1000))+4, color=blues[j])
        j = j +1
    pylab.legend(loc='upper left', ncol=3)
    
    ax2 = qMS.plotStatsDataStruct(AllSubunits, dataByProtein[names10[0]], labels10[0], offset = 1.0/float(len(names10)+1), markerSize = 5/float(len(names10))+4, color=blues[0])
    j = 1
    for i in names10[1:]:
        ax = qMS.addStatsDataStructToPlot(AllSubunits, dataByProtein[i], ax, labels10[j], offset = (1.0/float(len(names10)+1))*(j+1), markerSize = 5/float(len(names10))+4, color=blues[j])
        j = j +1
    pylab.legend(loc='upper left', ncol=3)
    
    ax = qMS.plotStatsDataStruct(AllSubunits, dataByProtein[names1000[0]], labels1000[0], offset = 1.0/float(len(names1000)+1), markerSize = 5/float(len(names1000))+4, color=blues[0])
    pylab.legend(loc='upper left')
    pylab.tight_layout()
    pylab.savefig('1000uM_t1.png')
    
    ax = qMS.plotStatsDataStruct(AllSubunits, dataByProtein[names1000[1]], labels1000[1], offset = 1.0/float(len(names1000)+1), markerSize = 5/float(len(names1000))+4, color=blues[1])
    pylab.tight_layout()
    pylab.savefig('1000uM_t2.png')
    pylab.legend(loc='upper left')
    
    ax = qMS.plotStatsDataStruct(AllSubunits, dataByProtein[names1000[2]], labels1000[2], offset = 1.0/float(len(names1000)+1), markerSize = 5/float(len(names1000))+4, color=blues[2])
    pylab.tight_layout()
    pylab.savefig('1000uM_t3.png')
    pylab.legend(loc='upper left')
    
    ax = qMS.plotStatsDataStruct(AllSubunits, dataByProtein[names10[0]], labels10[0], offset = 1.0/float(len(names10)+1), markerSize = 5/float(len(names10))+4, color=reds[0])
    pylab.tight_layout()
    pylab.savefig('10uM_t1.png')
    pylab.legend(loc='upper left')
    
    ax = qMS.plotStatsDataStruct(AllSubunits, dataByProtein[names10[1]], labels10[1], offset = 1.0/float(len(names10)+1), markerSize = 5/float(len(names10))+4, color=reds[1])
    pylab.tight_layout()
    pylab.savefig('10uM_t2.png')
    pylab.legend(loc='upper left')
    
    ax = qMS.plotStatsDataStruct(AllSubunits, dataByProtein[names10[2]], labels10[2], offset = 1.0/float(len(names10)+1), markerSize = 5/float(len(names10))+4, color=reds[2])
    pylab.tight_layout()
    pylab.savefig('10uM_t3.png')
    pylab.legend(loc='upper left')
    '''
    
    
##################Find the medians###########################
    for prot in AllSubunits:
        medDict10[prot] = [[ntDict10['times'][i],getMedian(dataByProtein[ntDict10['names'][i]][prot]['vals'])] for i in [0,1,2]]
        medDict1000[prot] = [[ntDict1000['times'][i],getMedian(dataByProtein[ntDict1000['names'][i]][prot]['vals'])] for i in [0,1,2]]

##################Fit data using pool size equation###########################
    poolSizeDict = fitPoolSizes(AllSubunits, medDict10, medDict1000)
    
##################Plot the fits###########################
    '''
    fLarge1 = pylab.figure(figsize=(7,10))
    f1axs = plotPoolPage(LargeSubunit[0:10], poolSizeDict, medDict10, medDict1000, fLarge1, '50S-1.pdf')    
    
    fLarge2 = pylab.figure(figsize=(7,10))
    f1axs = plotPoolPage(LargeSubunit[10:20], poolSizeDict, medDict10, medDict1000, fLarge2, '50S-1.pdf')    
    
    fLarge3 = pylab.figure(figsize=(7,10))
    f3axs = plotPoolPage(LargeSubunit[20:30], poolSizeDict, medDict10, medDict1000, fLarge3, '50S-3.pdf')
    
    fSmall1 = pylab.figure(figsize=(7,10))
    f4axs = plotPoolPage(SmallSubunit[0:10], poolSizeDict, medDict10, medDict1000, fSmall1, '30S-1.pdf')

    fSmall2 = pylab.figure(figsize=(7,10))
    f5axs =plotPoolPage(SmallSubunit[10:], poolSizeDict, medDict10, medDict1000, fSmall2, '30S-2.pdf')
    '''
##################Make bar graphs###########################
    '''
    plotPoolBar(poolSizeDict, '10', 'r', '10 $\mu$M')
    plotPoolBar(poolSizeDict, '1000', 'b', '1 mM')
    plotTwoPoolBars(poolSizeDict, LargeSubunit, '10', 'r', '10 $\mu$M', '1000', 'b', '1 mM')
    plotTwoPoolBars(poolSizeDict, SmallSubunit, '10', 'r', '10 $\mu$M', '1000', 'b', '1 mM')
    '''
##################Read protein inventory data###########################
    pulse = False
    numerator = ["ampu"]
    denominator = ["ampu", "ampl"]
    path = '/home/jhdavis/scripts/python/figures/fitPuseData/msuPulseData/'
    proteinToNormalizeTo = "BSubL24"

    LargeSubunit = ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL27', 'BSubL28', 'BSubL29', 'BSubL30', 'BSubL32',
                    'BSubL33a', 'BSubL35', 'BSubL36']
    


    McMaster45S = [path+i for i in ["McMaster45S_esi-run1_filt.csv", "McMaster45S_esi-run2.1_filt.csv", "McMaster45S_esi-run2_filt.csv", "McMaster45S_qtof_filt_filtppm.csv"]]
    McMaster50S = [path+i for i in ["McMaster50S_esi-run1_filt.csv", "McMaster50S_esi-run2.1_filt.csv", "McMaster50S_esi-run2_filt.csv", "McMaster50S_qtof_filt_filtppm.csv"]]

    fileLists = [McMaster45S, McMaster50S]
    merged = []
    for listOfFiles in fileLists:
        merged.append(qMS.mergeFiles(listOfFiles, pulse, numerator, denominator, proteinToNormalizeTo, LargeSubunit))

    myPlot = qMS.makePlotWithDataSets(merged, LargeSubunit, ["McMaster45S_merged", 'McMaster50S_merged'])
    for i in LargeSubunit:
        pVal = stats.ttest_ind(merged[0][i]['vals'], merged[1][i]['vals'], equal_var=False)
    print merged[0]['BSubL02']['vals']    
    print merged[1]['BSubL20']['vals']
    L30ttest = stats.ttest_ind(merged[0]['BSubL02']['vals'], merged[1]['BSubL02']['vals'])
    L12ttest = stats.ttest_ind(merged[0]['BSubL34']['vals'], merged[1]['BSubL34']['vals'])
    print L30ttest
    print L12ttest

##################Plot pool data vs. protein inventory data###########################
    '''
    verifiedZero = ['BSubL16', 'BSubL28', 'BSubL36']    
    for z in verifiedZero:    
        merged[0][z] = {'vals':[0.05], 'nvals':1, 'loc':0, 'flab':0.05}
        
    poolSizeVsPIDict = {}
    for prot in LargeSubunit:
        poolSizeVsPIDict[prot] = {'PI':numpy.median(numpy.array(merged[0][prot]['vals'])), 'PS':poolSizeDict[prot]['10']}
    names = poolSizeVsPIDict.keys()
    names.sort()
    #x = [numpy.log2(poolSizeVsPIDict[i]['PS']+1) for i in names]
    x = [poolSizeVsPIDict[i]['PS'] for i in names]
    y = [poolSizeVsPIDict[i]['PI'] for i in names]
    scat = pylab.figure(figsize=(10,7))
    scatAx = scat.add_subplot(111)    
    scatAx.scatter(x,y, c='g', s=100)
    scatAx.set_title('protein abundance vs. pool size under RbgA-limiting conditions')
    scatAx.set_xlabel('precursor pool size (P)')
    scatAx.set_ylabel('relative abundance in 45S particle')
    scatAx.set_xlim([-0.1,2.5])
    scatAx.set_xlim([-0.1,2])
    scatAx.set_ylim([0,1.2])
    scatAx.set_yticks([0,.25,.50,.75, 1.0])
    scatAx.set_xticks([0,0.5,1, 1.5, 2, 2.5])
    #scatAx.set_xticks([0,0.5,1,1.5])
    scatAx.yaxis.tick_left()
    scatAx.xaxis.tick_bottom()
    for prot, xl, yl in zip(names, x, y):
        scatAx.annotate(str(prot[4:]), xy = (float(xl), float(yl)), xytext = (15,15), textcoords = 'offset points', arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    '''
    pylab.show('all')