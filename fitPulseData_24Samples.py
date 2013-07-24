# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 23:13:53 2013

@author: jhdavis
"""

import pylab
import matplotlib.pyplot as plt
import vizLib
import numpy
import qMS
from scipy import optimize

def poolResid1000(P, y, t):
    k = qMS.growthRate(47)
    return y - qMS.poolFunc(k, t, P)

def poolResid10(P, y, t):
    k = qMS.growthRate(92)
    return y - qMS.poolFunc(k, t, P)

def turnoverResid1000(d, y, t):
    k = qMS.growthRate(47)
    #return y - qMS.overLabelingFunc(k, t, d)
    return qMS.overLabelingFunc(k, t, d) - y

def turnoverResid10(d, y, t):
    k = qMS.growthRate(92)
    #return y - qMS.overLabelingFunc(k, t, d)
    return qMS.overLabelingFunc(k, t, d) - y

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

def fitTurnover(Proteins, medDict10, medDict1000):
    d0_10 = 0
    d0_1000 = 0
    
    protTODict = {}
    for prot in Proteins:
        print prot
        ts10 = numpy.array([e[0] for e in medDict10[prot]])
        meas10 = numpy.array([e[1] for e in medDict10[prot]])

        ts1000 = numpy.array([i[0] for i in medDict1000[prot]])
        meas1000 = numpy.array([i[1] for i in medDict1000[prot]])

        d_10 = optimize.leastsq(turnoverResid10, [d0_10], args=(meas10,ts10))
        d_1000 = optimize.leastsq(turnoverResid1000, [d0_1000], args=(meas1000,ts1000))
        protTODict[prot] = {'10':d_10[0][0], '1000':d_1000[0][0]}
        print protTODict[prot]
    return protTODict

def plotLabelKinsPage(Proteins, protTODict, medDict10, medDict1000, fig, name):
    time = numpy.linspace(0, 140, 200)
    k10 = qMS.growthRate(92)
    k1000 = qMS.growthRate(47)
    
    Large1Ax = []
    a = 0
    for prot in Proteins:
        
        d_10 = protTODict[prot]['10']
        #p_10String = qMS.calcPercent(p_10, sigfig=2).split('.')[0]+'%'
        d_1000 = protTODict[prot]['1000']
        #p_1000String = qMS.calcPercent(p_1000, sigfig=2).split('.')[0]+'%
        
        ax = fig.add_subplot(5,2,a+1)
        
        ax.scatter([i[0] for i in medDict10[prot]], [i[1] for i in medDict10[prot]], c = 'r')
        ax.scatter([i[0] for i in medDict1000[prot]], [i[1] for i in medDict1000[prot]], c = 'b')
        ax.plot(time, qMS.maxLabFunc(k10, time), c='r')
        ax.plot(time, qMS.maxLabFunc(k1000, time), c='b')
        
        ax.plot(time, qMS.overLabelingFunc(k10, time, d_10), c='r', ls='--')
        ax.plot(time, qMS.overLabelingFunc(k1000, time, d_1000), c='b', ls='--')
        
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
        
def plotPoolBar(poolSizeDict, key, color, title, figname):
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
    pylab.savefig(figname)
    return bar

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

##################Initialize Varibles###########################    
if __name__ == "__main__":
    vizLib.setRcs(scale=12)
    SmallSubunit = ['BSubS02', 'BSubS03', 'BSubS04', 'BSubS05', 'BSubS06', 'BSubS07', 'BSubS08', 'BSubS09', 'BSubS10',
                    'BSubS11', 'BSubS12', 'BSubS13', 'BSubS14', 'BSubS15', 'BSubS16', 'BSubS17', 'BSubS18', 'BSubS20']

    LargeSubunit = ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL27', 'BSubL28', 'BSubL29', 'BSubL30', 'BSubL31a', 'BSubL32',
                    'BSubL33a', 'BSubL35', 'BSubL36']

    AllSubunits = LargeSubunit + SmallSubunit
    
    path = '/home/jhdavis/data/2013_05_28-MSUPulse/filtered/'

    reds = ['#fee5d9', '#fcbba1', '#fc9272', '#fb6a4a', '#de2d26', '#a50f15']
    blues = ['#eff3ff', '#c6dbef', '#93cae1', '#6baed6', '#3182bd', '#08519c']
    '''
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
    '''
    f1 = 'muspulse13_iso_res_filt.csv'
    f2 = 'muspulse14_iso_res_filt.csv'
    f3 = 'muspulse15_iso_res_filt.csv'
    f4 = 'muspulse16_iso_res_filt.csv'
    f5 = 'muspulse17_iso_res_filt.csv'
    f6 = 'muspulse18_iso_res_filt.csv'
    
    f7 = 'muspulse19_iso_res_filt.csv'
    f8 = 'muspulse20_iso_res_filt.csv'
    f9 = 'muspulse21_iso_res_filt.csv'
    f10 = 'muspulse22_iso_res_filt.csv'
    f11 = 'muspulse23_iso_res_filt.csv'
    f12 = 'muspulse24_iso_res_filt.csv'    
    
    
    names10 = [path+n for n in [f1, f2, f3, f4, f5, f6]]
    names1000 = [path+n for n in [f7, f8, f9, f10, f11, f12]]
    
    
    
    times1000 = [11, 22, 33, 44, 55, 66]
    times10 = [22, 49, 66, 88, 110, 132]
    
    ntDict10 = {'names': names10, 'times': times10}
    ntDict1000 = {'names': names1000, 'times': times1000}

    labels10 = ['10 $\mu$M, 22 min', '10 $\mu$M, 49 min', '10 $\mu$M, 66 min', '10 $\mu$M, 88 min', '10 $\mu$M, 110 min', '10 $\mu$M, 132 min']
    labels1000 = ['1 mM, 11 min', '1 mM, 22 min', '1 mM, 33 min', '1 mM, 44 min', '1 mM, 55 min', '1 mM, 66 min']
    
    medDict10 = {}
    medDict1000 = {}
    
    '''
    labels1000 = labels1000[0:3]
    times1000 = times1000[0:3]
    names1000 = names1000[0:3]
    '''
 
    num = ['AMP_L']
    den = ['AMP_U', 'AMP_L']
    
##################Extract the data###########################
    
    csvs = {}
    dataByProtein = {}
    allProteins = set()
    dataByProtein = qMS.multiStatsDict(names10+names1000, num, den)
    
##################Plot the datasets###########################
    '''
    yMax=1.25
    pylab.close('all')
    ax = vizLib.makePlotWithFileList(names10, num, den, AllProteins=AllSubunits, yMax=yMax, names=labels10, colors=reds, figSize=(22,7),median=False)
    pylab.legend(loc='upper left', ncol=3)
    pylab.xticks(numpy.arange(1,len(AllSubunits)+1,1), [item[4:] for item in AllSubunits], rotation=45)
    ax.set_ylim([0,1])
    ax.set_title('non-permissive conditions : 45S\nfraction labeled', multialignment='center')
    ax.set_ylabel('Labeled/[Labeled+Unlabeled]')
    pylab.tight_layout()
    #pylab.savefig('10uM_fracLab_med.pdf')
    #pylab.savefig('10uM_fracLab_med.png')
    
    ax2 = vizLib.makePlotWithFileList(names1000, num, den, AllProteins=AllSubunits, yMax=yMax, names=labels1000, colors=blues, figSize=(22,7), median=False)
    pylab.xticks(numpy.arange(1,len(AllSubunits)+1,1), [item[4:] for item in AllSubunits], rotation=45)
    ax2.set_title('permissive conditions : 50S\nfraction labeled', multialignment='center')
    ax2.set_ylabel('Labeled/[Labeled+Unlabeled]')
    pylab.legend(loc='upper left', ncol=3)
    pylab.tight_layout()
    #pylab.savefig('1mM_fracLab_med.pdf')
    #pylab.savefig('1mM_fracLab_med.png')
    
    
    num = ['AMP_U', 'AMP_L']
    den = ['AMP_U', 'AMP_L', 'AMP_S']
    
    ax = vizLib.makePlotWithFileList(names10[:-1], num, den, AllProteins=AllSubunits, normProtein='BSubL24', yMax=yMax, names=labels10, colors=reds, figSize=(22,7), median=False)
    pylab.legend(loc='lower left', ncol=3, numpoints=1)
    pylab.xticks(numpy.arange(1,len(AllSubunits)+1,1), [item[4:] for item in AllSubunits], rotation=45)
    ax.set_title('non-permissive conditions : 45S\ntotal protein', multialignment='center')
    ax.set_ylabel('[Labeled+Unlabeled]/[Total]')
    pylab.tight_layout()
    #pylab.savefig('10uM_proteinLevels_med.pdf')
    #pylab.savefig('10uM_proteinLevels_med.png')
    
    ax2 = vizLib.makePlotWithFileList(names1000, num, den, AllProteins=AllSubunits, normProtein='BSubL24', yMax=yMax, names=labels1000, colors=blues, figSize=(22,7), median=False)
    ax2.set_title('permissive conditions : 50S\ntotal protein', multialignment='center')
    ax2.set_ylabel('[Labeled+Unlabeled]/[Total]')
    pylab.xticks(numpy.arange(1,len(AllSubunits)+1,1), [item[4:] for item in AllSubunits], rotation=45)
    pylab.legend(loc='lower left', ncol=3, numpoints=1)
    pylab.tight_layout()
    #pylab.savefig('1mM_proteinLevels_med.pdf')
    #pylab.savefig('1mM_proteinLevels_med.png')
    '''
    
##################Find the medians###########################
    
    for prot in AllSubunits:
        medDict10[prot] = [[ntDict10['times'][i],numpy.median(dataByProtein[ntDict10['names'][int(i)]][prot])] for i in [0,1,2,3,4] if len(dataByProtein[ntDict10['names'][int(i)]][prot]) > 0]
        #medDict10[prot] = [[ntDict10['times'][i],numpy.median(dataByProtein[ntDict10['names'][int(i)]][prot])] for i in range(6,12) if len(dataByProtein[ntDict10['names'][int(i)]][prot]) > 0]
        #medDict1000[prot] = [[ntDict1000['times'][i],numpy.median(dataByProtein[ntDict1000['names'][int(i)]][prot])] for i in range(6,12) if len(dataByProtein[ntDict1000['names'][int(i)]][prot]) > 0]
        medDict1000[prot] = [[ntDict1000['times'][i],numpy.median(dataByProtein[ntDict1000['names'][int(i)]][prot])] for i in [0,1,2,3,4,5] if len(dataByProtein[ntDict1000['names'][int(i)]][prot]) > 0]
        #medDict1000[prot] = [[ntDict1000['times'][i],numpy.median(dataByProtein[ntDict1000['names'][int(i)]][prot])] for i in [0,1,2] if len(dataByProtein[ntDict1000['names'][int(i)]][prot]) > 0]
    
##################Fit data using pool size equation###########################
    
    #poolSizeDict = fitPoolSizes(AllSubunits, medDict10, medDict1000)
    turnoverSizeDict = fitTurnover(SmallSubunit, medDict10, medDict1000)
    
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
    '''
    ext = '.pdf'
    '''
    
    fLarge1 = pylab.figure(figsize=(7,10))
    f1axs = plotLabelKinsPage(LargeSubunit[0:10], turnoverSizeDict, medDict10, medDict1000, fLarge1, '50S-1'+ext)
    
    fLarge2 = pylab.figure(figsize=(7,10))
    f1axs = plotLabelKinsPage(LargeSubunit[10:20], turnoverSizeDict, medDict10, medDict1000, fLarge2, '50S-2'+ext)    
    
    fLarge3 = pylab.figure(figsize=(7,10))
    f3axs = plotLabelKinsPage(LargeSubunit[20:30], turnoverSizeDict, medDict10, medDict1000, fLarge3, '50S-3'+ext)
    '''
    fSmall1 = pylab.figure(figsize=(7,10))
    f4axs = plotLabelKinsPage(SmallSubunit[0:10], turnoverSizeDict, medDict10, medDict1000, fSmall1, '30S-1'+ext)

    fSmall2 = pylab.figure(figsize=(7,10))
    f5axs =plotLabelKinsPage(SmallSubunit[10:], turnoverSizeDict, medDict10, medDict1000, fSmall2, '30S-2'+ext)
    
    
    
    
##################Make bar graphs###########################
    '''
    plotPoolBar(poolSizeDict, '10', 'r', '10 $\mu$M', '10uM.pdf')
    plotPoolBar(poolSizeDict, '1000', 'b', '1 mM', '1000uM.pdf')
    plotTwoPoolBars(poolSizeDict, LargeSubunit, '10', 'r', '10 $\mu$M', '1000', 'b', '1 mM', 'LargeSubunit.pdf')
    plotTwoPoolBars(poolSizeDict, SmallSubunit, '10', 'r', '10 $\mu$M', '1000', 'b', '1 mM', 'SmallSubunit.pdf')
    '''
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
    
    pylab.show('all')