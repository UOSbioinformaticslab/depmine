#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 13:56:02 2022

@author: lp217
"""

# =============================================================================
# 28/11/2022 changed Cell definition to capture disruptive 
# and non-disruptive mutations so that non-disruptives are excluded from the 
# not-mutated set for significance calculation
# 
# 31/1/2023 working on MUTS option to scatter plot dependency 
# against mutation position along gene length - just needs ploting tidied up
# 
# 1/2/2023 Major fork to include CNV data as a measure of gene disruption.
# initialy do this by adding additional to end of positionally disruptive mutations
# 
# 3/2/2023 rejig of data structures and compilation of CNV data completed - now to 
# modify significance calcualtions
# 
# 4/2/2023 - CNV fully incorported - still need to prot chnges to batch version 
# for all by all calculation
# 
# 5/2/2023 - Introduce options in ALL, to allow analysis of GeneAs showing 
# significant switching behaviour, for a given mutated GeneB 
# 
# 9/2/2023 - Identify MSI and TP53 status of cell lines using the database from 
# Chan et al Nature volume 568, pages 551–556 (2019). May add our own MSI caller later

# 10/2/2023 - Option in TYPE to plot all dependency for a given GeneA regardless of GeneB status
#
# 18/4/2023 = Added chromosome' option to ALL to analyse chromosomal location of GeneBs that give
# significant switching of a GeneA
#
# 5/7/2023 - added LINE to select cell lines on the basis of genetic profile
# 10/7/2023 - added PROF to manage genetic profiles
# NB - these eliminate the need for the hardwired MSI and TP53 definitions,
# so can move to a new profile version
#
# 13/7/2023 =  adding new function to allow specific pairwise analysis based on defined non-synonomous
#             change in a geneB - e.g. KRAS G12C 
#
# 31/7/2023 -  basics of set-based Genetic Profile's now in place. Need to create
#              functions to enable use of GPs in switch analysis
# 
# 8/8/2023 -   full Genetic Profile system now installed, including handling tools (LINE,PROF), 
#              switch plot for specified target wrt a GP (PAIR) ? may want to change that name,
#              and all gene search for targets switched by a profile. Also defined mutations
#              can have wild-cards eg V600? - ALL DONE
# 30/8/2023 -  adapt code for DepMap 2023Q2 data release - DONE
# 
# 1/9/2023  -  secondary effect analysis added to PAIR for disruptive mutations and GOFs - CNV amp/delstill
#              to add
#
# 26/9/2023 -  coalescing repetitious dependency plots routines into functions - ONGOING
# 
# 6/10/2023 -  expand switch routine to allow comparison of target dependency
#              against two defined profiles eg WRN KMT2B KMT2D
#              DONE
# 10/10/2023   cell selection by cytogenetic band and subband 
#              DONE for deletion, and amplification. NullHype also
#              modified
# 13/10/2023 - GOF added to secondary analysis in PAIR
#              TO DO - add CNV to secondary selection
#              TO DO - allow user to select profile file so different profile
#                      datasets can be used for different projects
# 20/10/2023 - add LOCUs function to analyse genes in a cytogenetic band - PENDING

# 23/10/2023 - Major fork - adding in expression data
#              stage 1 - adding set selection based on significant under- or over-expression
# 2/11/2023  - stage 2 - added over/underexpression to database build
#
# 5/11/2023  - new EXPR function to visualise expression patterns and 
#              correlations between expression and CNV 

# 14/11/2023 - need nameable profile project files that can be selectively saved and loaded - DONE
#
# 12/12/2023 - adding new HIGH function to perform various analyses on high dependency genes in given cell subset
#
# 29/01/2024 - changed expression selection so that <gene or >gene without a number 
#              uses the precalcualted 'over' and 'under' thresholds to define the sets,
#              rather than a SD-based calculation. IF PROBLEM - REVERT TO version 15

# 29/01/2024 - All gene names in muts,cnvs and exprnow made to conform to HUGO reference, so larger numbers available
#              for profile definition.

# 5/2/2024   - Added mining by tissue type, CNV amplification or deletion to PAIR function
#            - Rewrote core of ProfPairPValue to use pandas access routines giving big speed increase
#              now called AltProfPairPValue

# 22/2/2024  - added mining by under/over-expression to PAIR

# 1/3/2024   - Good to parallelise main loops in mutation, CNV and especially expression mining

# 14/3/2024  - Modified newlist to allow tissue restriction when processing external pair lists
#              Required addition of optional 'mask' parameter to PairPValue 
#
# 19/3/2024  - Need to add age and sex in secondary mining
# 22/3/2024  - Add selection based on presence of gene fusions - substantial additions required :
#              1) Load and normalise fusion data file - DONE
#              2) Build FUS primitive selection routine - DONE
#              3) incorporate into parsing = DONE
#
# 27/3/2024  - Mining and secondary factor optimisation system now in place and working.
#            - Still need to add mining for gene fusions

# =============================================================================


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
import matplotlib.patches as mpatch
from collections import OrderedDict,Counter
from scipy import stats
from scipy.signal import argrelextrema, find_peaks
from timeit import default_timer as timer
import multiprocessing as mp
# from igraph import Graph,plot
import re
import requests, sys, json
from math import log2,log10,trunc
import seaborn as sns
import random
from glob import glob


def mut_burden():
    cm = []
    for c in info['ModelID']:
        cm.append(len(muts[muts['ModelID'] == c]))
    burden = pd.DataFrame(info['ModelID'])
    burden['N'] = cm
    return burden[burden['N'] > 0]
    
def dict_merge(a,b):
    c = b.copy()
    for i in a:
        if i in c:
            c[i] = a[i]+c[i]
        else:
            c[i] = a[i]
    return dict(sorted(c.items(), key = lambda x: x[1], reverse = True))

def dict_set_redux(d,s):
    c = {}
    for i in s:
        if i in d:
            c[i] = d[i]
    return dict(sorted(c.items(), key = lambda x: x[1], reverse = True))
    

def match_keyword(u,l):
    """ matches ambiguous keywords """
    result = ''
    if u in l:# check for exact match first
        return u    
    for q in l:
        if u == q[:len(u)]:
            if result != '':
                return '*'
            else :
                result = q
    return result

def most_frequent(List):
    counter = 0
    num = List[0]

    for i in List:
        curr_frequency = List.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i

    return num

def Flistor(r):
    
    """ flattens list and replaces final comma in list with or """
    s = ', '.join(r)
    return "".join(reversed("".join(reversed(s)).replace(',','ro ',1)))

def Flistand(r):
    
    """ flattens list and replaces final comma in list with and """
    s = ', '.join(r)
    return "".join(reversed("".join(reversed(s)).replace(',','dna ',1)))


def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

global debug, muts_loaded, info_loaded, cnvs_loaded, seaborn_plot, phys_chro_loaded, expr_loaded
global oncogenes_loaded, exval_loaded, fuse_loaded, gofs_loaded, lofs_loaded, sigs_loaded, \
        paths_loaded,chromlen_loaded, ensem_loaded 
global thresh,frac_limit,p_limit,null_trim, cohen_min, sig_thresh
thresh = 0.65
frac_limit = 0.025
p_limit = 0.05
cohen_min = 0.8
sig_thresh = 0.05

deps = []
cnvs = []
muts = []
phys_chro = []
expr = []
exval = []
oncogenes = []
fuse = []
gofs = []
lofs = []
sigs = []
paths = []
ensem = []

# =============================================================================
#  ██████╗██╗      █████╗ ███████╗███████╗███████╗███████╗
# ██╔════╝██║     ██╔══██╗██╔════╝██╔════╝██╔════╝██╔════╝
# ██║     ██║     ███████║███████╗███████╗█████╗  ███████╗
# ██║     ██║     ██╔══██║╚════██║╚════██║██╔══╝  ╚════██║
# ╚██████╗███████╗██║  ██║███████║███████║███████╗███████║
#  ╚═════╝╚══════╝╚═╝  ╚═╝╚══════╝╚══════╝╚══════╝╚══════╝
#                                                         
# =============================================================================

class Gene:
    def __init__(self, gene):
        self.name = gene
        self.cells = []
        self.cnvdel = []
        self.cnvamp = []
        self.over = []
        self.under = []
        self.frac = 0.0
        self.p = 1.0
        
class Cell:
    def __init__(self, cell):
        self.name = cell
        self.tissue = 'undefined'
        self.cancer = 'undefined'
        self.dgenes = [] # genes with 'damaging' mutations
        self.mgenes = [] # genes with 'non-silent' and non-damaging mutations
        self.cnvdel = [] # genes with CNV indicative of deletion
        self.cnvamp = [] # genes with CNV indicative of amplification (and gain ?)
        self.over = [] # genes with significant overexpression relative to mean
        self.under = [] # genes with significant underexpression relative to mean
               
class Hit:
    def __init__(self,pv,fv,m,n):
        self.p = pv
        self.frac = fv
        self.mut_deps = m
        self.notmut_deps = n
           
class Pair:
    def __init__(self, gA,gB,p):
        self.geneA = gA
        self.geneB = gB
        self.p = 1.0

class Mut:
    def __init__(self,g,c,t,m,ty,p,d):
        self.gene = g # gene HUGO name
        self.cell = c # cell line
        self.trans = t # associated annotation transcript
        self.mut = m # protein change mutation string
        self.type = ty # mutational type
        self.pos = p # transcript flag position for mutation
        self.dep = d # dependency value for this gene in this cell line

# =============================================================================
# ███████╗██╗███╗   ██╗ ██████╗ ███████╗██████╗ ██████╗ ██████╗ ██╗███╗   ██╗████████╗███████╗
# ██╔════╝██║████╗  ██║██╔════╝ ██╔════╝██╔══██╗██╔══██╗██╔══██╗██║████╗  ██║╚══██╔══╝██╔════╝
# █████╗  ██║██╔██╗ ██║██║  ███╗█████╗  ██████╔╝██████╔╝██████╔╝██║██╔██╗ ██║   ██║   ███████╗
# ██╔══╝  ██║██║╚██╗██║██║   ██║██╔══╝  ██╔══██╗██╔═══╝ ██╔══██╗██║██║╚██╗██║   ██║   ╚════██║
# ██║     ██║██║ ╚████║╚██████╔╝███████╗██║  ██║██║     ██║  ██║██║██║ ╚████║   ██║   ███████║
# ╚═╝     ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚═╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝
#                                                                                             
# =============================================================================

def MakeGoFFing(d,info,muts):
    p = set(info[info['OncotreePrimaryDisease'] == d]['ModelID'])
    q = set(info['ModelID']) ^ p

# get all GoF mutations in p

    pGoFs = muts[muts['ModelID'].isin(p) & muts['LikelyGoF'] == True] \
                [['HugoSymbol','ProteinChange']]
    lg = list(pGoFs['HugoSymbol'])
    lc = list(pGoFs['ProteinChange'])
# combine into mutation
    pmts = [lg[i]+lc[i] for i in range(len(lg))]

    qGoFs = muts[muts['ModelID'].isin(q) & muts['LikelyGoF'] == True] \
                [['HugoSymbol','ProteinChange']]
    lg = list(qGoFs['HugoSymbol'])
    lc = list(qGoFs['ProteinChange'])
# combine into mutation
    qmts = [lg[i]+lc[i] for i in range(len(lg))]
    
    all_muts = sorted(set(pmts) | set(qmts)) 
    
    cont = [[0,0],[0,0]] # initialise contingency table
    
    lp = len(pmts)
    lq = len(qmts)
    r = {}
    
    for m in all_muts:
        Np = pmts.count(m)
        cont[0][0] = Np
        cont[1][0] = lp-Np
        Nq = qmts.count(m)
        cont[0][1] = Nq
        cont[1][1] = lq-Nq
        v = stats.boschloo_exact(cont,alternative="greater").pvalue
        if np.isnan(v):
            r[m] = 0.0
        else:
            r[m] = abs(-log10(v))
    return (d,dict(sorted(r.items(),key = lambda x: x[1], reverse = True)))    
    return (d,r)    

# =============================================================================
# ███████╗██╗ ██████╗ ███╗   ██╗ █████╗ ████████╗██╗   ██╗██████╗ ███████╗███████╗
# ██╔════╝██║██╔════╝ ████╗  ██║██╔══██╗╚══██╔══╝██║   ██║██╔══██╗██╔════╝██╔════╝
# ███████╗██║██║  ███╗██╔██╗ ██║███████║   ██║   ██║   ██║██████╔╝█████╗  ███████╗
# ╚════██║██║██║   ██║██║╚██╗██║██╔══██║   ██║   ██║   ██║██╔══██╗██╔══╝  ╚════██║
# ███████║██║╚██████╔╝██║ ╚████║██║  ██║   ██║   ╚██████╔╝██║  ██║███████╗███████║
# ╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚══════╝╚══════╝
# =============================================================================
"""
Signature aetiologies assigned from https://cancer.sanger.ac.uk/signatures/sbs/ 
"""

               
SigGroups = {
   '5MC':['SBS1'],
   'UNASSIGNED':['SBS5','SBS8', 'SBS28','SBS12','SBS16','SBS17A','SBS17B','SBS19','SBS23','SBS28','SBS33','SBS34','SBS37',
                 'SBS39','SBS40A','SBS40B','SBS40C','SBS41','SBS89','SBS91','SBS93','SBS94','SBS96','SBS97','SBS98'],
   'HALO':['SBS42'],
   'ROS':['SBS18'],
   'POLH_HM':['SBS9'],
   'POLED_MMR': ['SBS14','SBS20'],
   'MMR':['SBS6', 'SBS15', 'SBS21', 'SBS26', 'SBS44'],
   'POLED':['SBS10A', 'SBS10B', 'SBS10C', 'SBS10D'],
   'TMZ':['SBS11'],
   'AFLA':['SBS24'],
   'HR' : ['SBS3'],
   'BER' : ['SBS30','SBS36'],
   'IMMUNOSUPPRESSANTS' : ['SBS32'] ,
   'TREATMENT' : ['SBS25', 'SBS31', 'SBS32', 'SBS35', 'SBS86', 'SBS87', 'SBS90','SBS99'],
   'APOBEC' : ['SBS2','SBS13'],
   'TOBACCO' : ['SBS4', 'SBS29', 'SBS92'],
   'UV' : ['SBS7A', 'SBS7B', 'SBS7C', 'SBS7D', 'SBS38'],
   'A_A' : ['SBS22A','SBS22B'],
   'COLIBACTIN' : ['SBS88'],
   'LYMPHOID' : ['SBS9','SBS84','SBS85']}                                   
    

def WhichSigType(s):
    r = 'UNASSIGNED'
    for i in SigGroups:
        if s in SigGroups[i]:
            r = i
    return r

def SigType(s,p): # plots distribution of signature occurances for a set of cells
    # stypes = dict(zip(sigs.columns[1:-1],np.zeros(len(sigs.columns[1:-1]))))
# get counts for each signature in all cells in set

    sig_count =   dict(sigs[sigs['SAMPLES'].isin(s)].sum().loc[sigs.columns[1:-1]])
    norm = sum(sig_count.values())
    sname = [i for i in sig_count if sig_count[i]/norm > 0.001]
    n = [sig_count[i]/norm for i in sig_count if sig_count[i]/norm > 0.001]
    
    cmp = plt.get_cmap('tab20')
    
    sgs = {}
    for i,nm in enumerate(SigGroups):
        sgs[nm] = i

    fig, ax = plt.subplots()
    
    # ax.set_ylim(0.0,0.5)
    ax.set_ylabel('mutation fraction',rotation='vertical',fontsize=16)
    bars = ax.bar(sname, n)
    ax.set_title(f'Signatures for profile : {p}',fontsize=16)
    ax.set_xticks(range(len(sname)))
    ax.set_xticklabels(sname,rotation=90)
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=12)
    
# assign group to individual bars

    legs = {}

    for i,j in enumerate(sname):
        w = WhichSigType(j)
        if w == 'UNASSIGNED':
            bars[i].set_color([0.5,0.5,0.5,0.1]) # grey
            bars[i].set_edgecolor('grey')
            if w not in legs:
                legs[w] = [0.5,0.5,0.5,0.1]
        else:
            bars[i].set_color(list(cmp.colors[sgs[w]]))
            if w not in legs:
                legs[w] = list(cmp.colors[sgs[w]])
                
# add colour key
    ptch = []
    for i in legs:
        ptch.append(mpatch.Patch(color=legs[i],label = i ))
    plt.legend(handles=ptch, loc='upper right',fontsize = 8,ncol = 2)
    
    plt.show()
    

    
# =============================================================================
# ██████╗ ███████╗██████╗ ██████╗ ██╗      ██████╗ ████████╗
# ██╔══██╗██╔════╝██╔══██╗██╔══██╗██║     ██╔═══██╗╚══██╔══╝
# ██║  ██║█████╗  ██████╔╝██████╔╝██║     ██║   ██║   ██║   
# ██║  ██║██╔══╝  ██╔═══╝ ██╔═══╝ ██║     ██║   ██║   ██║   
# ██████╔╝███████╗██║     ██║     ███████╗╚██████╔╝   ██║   
# ╚═════╝ ╚══════╝╚═╝     ╚═╝     ╚══════╝ ╚═════╝    ╚═╝   
# =============================================================================
                                                          
def DepPlot2(d1,d2,t1,t2,p,D,l1='matched',l2='not matched'): # plot double violin plot
    if seaborn_plot:
        not_subset = random.choices(d2,k=len(d1)*2)
        data_both = pd.DataFrame({'profile':['match']*len(d1)+['notmatch']*len(not_subset),'dependency':d1+not_subset})
        sns.violinplot(data=data_both, x='profile', y='dependency',inner = None,palette=['cyan', 'orange'],scale='area')
        sns.swarmplot(data=data_both, x='profile',y="dependency", size=3,palette=['black', 'black'])
        plt.ylim(0.0,1.0)
        # plt.legend_ = None
        plt.title(f'{t1} : {t2} p = {p:.4} d = {D:.4}',fontsize = '20',pad=20)
        plt.tick_params(axis='x', labelsize=20)
        plt.ylabel('dependency',fontsize = '20')
        plt.show()
    else:
        
        plt.violinplot((d1,d2),showmeans=True, showextrema=False,widths = [0.5,0.5])
        plt.ylim(0.0,1.0)
        plt.xticks([1,2],labels = [l1,l2],fontsize = '20')
        plt.title(f'{t1} : {t2} p = {p:.3} d = {D:.3}',fontsize = '20',pad=20)
        plt.tick_params(axis='both', labelsize=12)
        plt.ylabel('dependency',fontsize = '20')
        plt.show()
                    
# =============================================================================
# ██████╗ ██╗      ██████╗ ████████╗ ██████╗███████╗██╗     ██╗  ████████╗██╗   ██╗██████╗ ███████╗
# ██╔══██╗██║     ██╔═══██╗╚══██╔══╝██╔════╝██╔════╝██║     ██║  ╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝
# ██████╔╝██║     ██║   ██║   ██║   ██║     █████╗  ██║     ██║     ██║    ╚████╔╝ ██████╔╝█████╗  
# ██╔═══╝ ██║     ██║   ██║   ██║   ██║     ██╔══╝  ██║     ██║     ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  
# ██║     ███████╗╚██████╔╝   ██║   ╚██████╗███████╗███████╗███████╗██║      ██║   ██║     ███████╗
# ╚═╝     ╚══════╝ ╚═════╝    ╚═╝    ╚═════╝╚══════╝╚══════╝╚══════╝╚═╝      ╚═╝   ╚═╝     ╚══════╝
#                                                                                                          
# =============================================================================
         
def PlotCellType(p):
    cancer = {}

    nincl = 0
    for i in p:
        c = info.loc[info['ModelID'] == i]
        if len(c) > 0:
            cname = c['OncotreeLineage'].iloc[0]
            nincl += 1
        else:
            continue
        if cname in cancer:
            cancer[cname] += 1
        else:
            cancer[cname] = 1
    print(f'{nincl} cell lines included')

    cancer = dict(sorted(cancer.items(), key = lambda x: x[1], reverse = True))

# plot bar chart
                    

    fig, ax = plt.subplots()
    c_name = [i.split()[0].lower() for i in cancer]
    c_count = [cancer[i] for i in cancer]
    ax.bar(c_name, c_count,alpha=0.3)

    ax.set_ylabel('number of cell lines  ',rotation='vertical',fontsize=16)
    ax.set_title(f'Cancer tissue distribution for profile : {m}',fontsize=16)
    ax.set_xticks(range(len(c_name)))
    ax.set_xticklabels(c_name,rotation=90)
    ax.tick_params(axis='both', labelsize=10)
    plt.show()

def GetOncoTreeSubtypes():
    return dict(sorted(Counter(list(info['OncotreeSubtype'])).items(),
                       key = lambda x: x[1], reverse = True))
def GetOncotreePrimaryDisease():
    return dict(sorted(Counter(list(info['OncotreePrimaryDisease'])).items(),
                       key = lambda x: x[1], reverse = True))    


# =============================================================================
#  ██████╗ ███████╗████████╗██████╗  █████╗ ██████╗  █████╗ ███╗   ███╗███████╗
# ██╔════╝ ██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██╔══██╗██╔══██╗████╗ ████║██╔════╝
# ██║  ███╗█████╗     ██║   ██████╔╝███████║██████╔╝███████║██╔████╔██║███████╗
# ██║   ██║██╔══╝     ██║   ██╔═══╝ ██╔══██║██╔══██╗██╔══██║██║╚██╔╝██║╚════██║
# ╚██████╔╝███████╗   ██║   ██║     ██║  ██║██║  ██║██║  ██║██║ ╚═╝ ██║███████║
#  ╚═════╝ ╚══════╝   ╚═╝   ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝
#                                                                                      
# =============================================================================
        
def GetParams():
    global thresh,frac_limit,p_limit,null_trim, cohen_min
    
    t = input(f'Dependency threshhold [{thresh}] : ').strip()
    if isFloat(t) and float(t) >=0.0 and float(t) <=1.0:
        thresh = float(t)
    print(f'Using threshold of {thresh}')
        
    t = input(f'Minimun mutational frequency [{frac_limit}] : ').strip()
    if isFloat(t) and float(t) >=0.0 and float(t) <=1.0:
        frac_limit = float(t)
    print(f'Using minimun mutational frequency of {frac_limit}')
         
    t = input(f'Maximum p value for hit [{p_limit}] : ').strip()
    if isFloat(t) and float(t) >=0.0 and float(t) <=1.0:
        p_limit = float(t)
    print(f'Using maximum p value for hit of {p_limit}')
    
    t = input(f'Minimum Cohen d value for hit [{cohen_min}] : ').strip()
    if isFloat(t) and float(t) >=0.0 and float(t) <=1.0:
        p_limit = float(t)
    print(f'Using minimum Cohen d value for hit of {cohen_min}')
    

# # Get factor to trim Null sets by 
    
#     t = input(f'Null set trim factor [{null_trim}] : ').strip()
#     if isFloat(t) and float(t) >=0.0 :
#         null_trim = float(t)
#     print(f'Using null set trim factor of {null_trim}')
    
    t = input(f'Minimum mutational signature fraction [{sig_thresh}] : ').strip()
    if isFloat(t) and float(t) >=0.0 and float(t) <=1.0:
        p_limit = float(t)
    print(f'Using minimum mutational signature fraction of {sig_thresh}')
    

# =============================================================================
# Load/Save profiles
# =============================================================================
def prof_load(fname):
    try:
        f = open(fname,'r')
    except:
        print(f'Profiles file {fname} does not exist')
        return False
    while True:
        s = f.readline()
        if s == '':
            break
        t = f.readline()
        if t == '':
            print(f'Definition missing for profile : {s}')
            break
        profiles[s[:-1].strip().upper()] = t[:-1].strip()
    f.close()
    return True

def prof_save(fname):
    f = open(fname,'w')
    for n in profiles:
        f.write(n.upper() +'\n')
        f.write(prof_flatten(profiles[n])+'\n')
    f.close()

def list_profile_files():
    prof_file_list = glob('*rofiles')
    if len(prof_file_list) > 0:
        print('Available profile files :')
        for p in prof_file_list:
            print(p)
    else:
        print('No profile files found')

# =============================================================================
# ██████╗ ██╗  ██╗██╗   ██╗███████╗         ██████╗██╗  ██╗██████╗  ██████╗ 
# ██╔══██╗██║  ██║╚██╗ ██╔╝██╔════╝        ██╔════╝██║  ██║██╔══██╗██╔═══██╗
# ██████╔╝███████║ ╚████╔╝ ███████╗        ██║     ███████║██████╔╝██║   ██║
# ██╔═══╝ ██╔══██║  ╚██╔╝  ╚════██║        ██║     ██╔══██║██╔══██╗██║   ██║
# ██║     ██║  ██║   ██║   ███████║███████╗╚██████╗██║  ██║██║  ██║╚██████╔╝
# ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ 
# routines for handling physical chromosome sets                                                                   
# =============================================================================

def chro_loc_set(locus): # set of genes at locus cleaned by hgnc names
    if ('chr'+locus) in phys_chro.columns:
        return set(phys_chro['chr'+locus]['geneSymbols']) & set(hgnc.index)
    else:
        return set()

def Locify(glist): # cluster genes into loci
    locs = {}
    for g in glist:
        l = hgnc.loc[g]["location"].split('-')[0]
        if l in locs:
            locs[l].append(g)
        else:
            locs[l] = [g]
    return locs
    
def Pathify(glist): # cluster genes into pathways
    pth = {}
    for g in glist:
        ps = set([x for x in paths if g in paths[x]])
        for j in ps:
            if j in pth:
                pth[j] += 1
            else:
                pth[j] = 1
    return dict(sorted(pth.items(), key = lambda x: x[1], reverse = True))
                
                
        
    

# =============================================================================
# ██████╗ ██████╗  ██████╗ ███████╗   ███████╗██╗      █████╗ ████████╗████████╗███████╗███╗   ██╗
# ██╔══██╗██╔══██╗██╔═══██╗██╔════╝   ██╔════╝██║     ██╔══██╗╚══██╔══╝╚══██╔══╝██╔════╝████╗  ██║
# ██████╔╝██████╔╝██║   ██║█████╗     █████╗  ██║     ███████║   ██║      ██║   █████╗  ██╔██╗ ██║
# ██╔═══╝ ██╔══██╗██║   ██║██╔══╝     ██╔══╝  ██║     ██╔══██║   ██║      ██║   ██╔══╝  ██║╚██╗██║
# ██║     ██║  ██║╚██████╔╝██║███████╗██║     ███████╗██║  ██║   ██║      ██║   ███████╗██║ ╚████║
# ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═══╝
#                                                                                                 
# flattens references to predefined profiles
# =============================================================================

def prof_flatten(s): # flattens a genetic profile
# expand operators and brackets
    for o in '()&|^#':
        s = s.replace(o,' '+o+' ')
    t = s.split()
    for i,w in enumerate(t):
        if w  in '()&|^':  # an operator
            continue
        elif w[0] not in '~*;+-=!/$#@%><{': # no prefix so must be preexisting profile
            if w.upper() in profiles: # replace with parse of stored profile
                t[i] = '(' + prof_flatten(profiles[w.upper()]) + ')'
            else:
                print(f'unrecognised token : {w}')
                return ''
        else:
            t[i] = w
    
    return ' '.join(t)

def prof_nullhype(s): # creates a null hypothesis set from a flattened genetic profile
    for o in '()&|^#':
        s = s.replace(o,' '+o+' ')
    t = s.split()
    for i,w in enumerate(t):
        if w  in '()&|^':  # an operator
            continue
# deal with Universal set    
        elif w[0] == '#': # Universal set 
            t[i] = "UNI()"
        elif w[0] == '/': # specific mutation <gene>:<change>
            w = w.upper()
            if w[-1] in '-+': # annotated missense GoF / LoF
                x = w[1:-1] # extract gene name
                if x not in hgnc.index:
                    print(f'mutation specification {w} not valid')
                    return ''
                t[i] = f"MUT('{x}')"
            else:
                x = w[1:].upper().split(':')
                if len(x) != 2 or (x[0] not in deps.columns and x[0] not in cnvs.columns): # valid gene name ?
                    print(f'mutation specification {w} not valid')
                    return ''
                t[i] = f"MUT('{x[0]}')"
        elif w[0] == '@': # load set from file
            t[i] = f"FIL('{w[1:].upper()}')"    
             
        elif w[0] == '{': # mutational signature
            s = w[1:].upper()
            if s in sigs.columns[1:-1]:
                t[i] = f"SIG('{s}')"
            else:
                m = match_keyword(s,SigGroups)
                if m == '':
                    print(f'mutational signature specification {w} not valid')
                    return
                elif m == '*':
                    print(f'mutational signature specification {w} ambiguous')
                    return
                t[i] = f"SIG('{m}')"

        elif w[0] == '%': # cancer type selection            
            c = match_keyword(w[1:].upper(),set(info['OncotreeLineage']))
            if c == '':
                print(f'cancer type specification {w} not valid')
                return
            elif c == '*':
                print(f'cancer type specification {w} ambiguous')
                return
            t[i] = f"TYP('{c}')"
        elif w[0] == '*': # sex definition
            if w[1].upper() not in 'MFU':
                print(f'cancer sex specification {w} not valid')   
                return
            if w[1].upper() == 'U':
                s = 'Unknown'
            elif w[1].upper() == 'F':
                s = 'Female'
            else:
                s = 'Male'
            t[i] = f"SEX('{s}')"
        elif w[0] == ';': # age definition
            if w[1].upper() not in 'APU':
                print(f'cancer age specification {w} not valid')   
                return
            if w[1].upper() == 'U':
                s = 'Unknown'
            elif w[1].upper() == 'P':
                s = 'Pediatric'
            else:
                s = 'Adult'
            t[i] = f"AGE('{s}')"
        else: # deal with prefix operators
            if w[0] in '<>' and w[1] in '0123':
                g = w[2:].upper() # extract gene name
            else:
                g = w[1:].upper() # extract gene name
                gg = ('chr'+g).replace('P','p').replace('Q','q')
            if g in deps.columns or g in cnvs.columns:

                if w[0] in '+-': # amplification or deletion
                    t[i] = f"(AMP('{g}') | DEL('{g}'))"
                elif w[0] == '=': # diploid eg wild type copy number
                    t[i] = f"DIP('{g}')"
                elif w[0] == '!': # gene disruption
                    t[i] = f"MUT('{g}')"
                elif w[0] == '$': # mutationally wild-type
                    t[i] = f"(WT('{g}') & DIP('{g}'))"
                elif w[0] == '~': # fusion
                    t[i] = f"FUS('{g}')"
                    
                elif w[0] == '>': # over-expression
                    if w[1] in '0123':
                        s = '0123'.index(w[1])
                        t[i] = f"OX('{g}',{s})"
                    else:
                        t[i] = f"OX('{g}',-1)" # -ve signifies use pre-calc thresholds
                elif w[0] == '<': # under-expression
                    if w[1] in '0123':
                        s = '0123'.index(w[1])
                        t[i] = f"UX('{g}',{s})"
                    else:
                        t[i] = f"UX('{g}',-1)" # -ve signifies use pre-calc thresholds             
                                       
            elif gg.split('.')[0] in phys_chro.columns:

                 if w[0] == '+': # amplification of locus
                     t[i] = f"LAMP('{gg}')"
                 elif w[0] == '-': # gene deletion
                     t[i] = f"LDEL('{gg}')"
            else:
                 print(f'function prefix {f} not valid')
                 return ''
    a = ' '.join(t)
    return f"UNI()-({a})"
            

def WithDeps(g):
    return list(set(g) & set(deps.columns))
    

# =============================================================================
# ██████╗ ██████╗  ██████╗ ███████╗   ██████╗  █████╗ ██████╗ ███████╗███████╗
# ██╔══██╗██╔══██╗██╔═══██╗██╔════╝   ██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔════╝
# ██████╔╝██████╔╝██║   ██║█████╗     ██████╔╝███████║██████╔╝███████╗█████╗  
# ██╔═══╝ ██╔══██╗██║   ██║██╔══╝     ██╔═══╝ ██╔══██║██╔══██╗╚════██║██╔══╝  
# ██║     ██║  ██║╚██████╔╝██║███████╗██║     ██║  ██║██║  ██║███████║███████╗
# ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝╚══════╝╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚══════╝
# convert a genetic profile definition into an evaluatable Python statement                                                                            
# =============================================================================
def prof_parse(s):
# expand operators and brackets
    for o in '()&|^#':
        s = s.replace(o,' '+o+' ')
    t = s.split()
    for i,w in enumerate(t):
        if w  in '()&|^':  # an operator
            continue
        elif w[0] not in '~*;+-<>=!/$#@%{': # no prefix so must be preexisting profile
            if w.upper() in profiles: # replace with parse of stored profile
                t[i] = '(' + prof_parse(profiles[w.upper()]) + ')'
            else:
                print(f'unrecognised token : {w}')
                return ''
# deal with Universal set    
        elif w[0] == '#': # Universal set
            t[i] = "UNI()"
        elif w[0] == '/': # specific mutation <gene>:<change>
            w = w.upper()
            if w[-1] == '-': # annotated missense LoF
                x = w[1:-1] # extract gene name
                if x not in hgnc.index:
                    print(f'mutation specification {w} not valid')
                    return ''
                agg = GeneAggLOF(x)
                if len(agg) > 0:
                    t[i] = f"DEM('{x}','{agg}')"
                else:
                    # print(f'No annotated LoF mutations found for {x}')
                    return ''
            elif  w[-1] == '+': # annotated missense GoF
                x = w[1:-1] # extract gene name
                if x not in hgnc.index:
                    print(f'mutation specification {w} not valid')
                    return ''
                agg = GeneAggGOF(x)
                if len(agg) > 0:
                    t[i] = f"DEM('{x}','{agg}')"
                else:
                    # print(f'No annotated GoF mutations found for {x}')
                    return ''
            else:
                x = w[1:].upper().split(':')
                if len(x) != 2 or (x[0] not in hgnc.index): # valid gene name ?
                    print(f'mutation specification {w} not valid')
                    return ''
                if x[1][-1] == 'X':
                    x[1][-1] = '?'
                if x[1][-1] not in 'ACDEFGHIKLMNPQRSTVWY?' : # invalid mutation specifier
                    print(f'mutation specification {w} not valid')
                    return ''
                t[i] = f"DEM('{x[0]}','{x[1]}')"
        elif w[0] == '@': # load set from file
            t[i] = f"FIL('{w[1:].upper()}')"    
            
        elif w[0] == '{': # mutational signature
            s = w[1:].upper()
            if s in sigs.columns[1:-1]:
                t[i] = f"SIG('{s}')"
            else:
                m = match_keyword(s,SigGroups)
                if m == '':
                    print(f'mutational signature specification {s} not valid')
                    return
                elif m == '*':
                    print(f'mutational signature specification {s} ambiguous')
                    return
                t[i] = f"SIG('{m}')"
                            
        elif w[0] == '%': # cancer type selection           
            c = match_keyword(w[1:].upper(),set(info['OncotreeLineage']))

            if c == '':
                print(f'cancer type specification {w} not valid')
                return
            elif c == '*':
                print(f'cancer type specification {w} ambiguous')
                return
            t[i] = f"TYP('{c}')"
        elif w[0] == '*': # sex definition
            if w[1].upper() not in 'MFU':
                print(f'cancer sex specification {w} not valid')   
                return
            if w[1].upper() == 'U':
                s = 'Unknown'
            elif w[1].upper() == 'F':
                s = 'Female'
            else:
                s = 'Male'
            t[i] = f"SEX('{s}')"
        elif w[0] == ';': # age definition
            if w[1].upper() not in 'APU':
                print(f'cancer age specification {w} not valid')   
                return
            if w[1].upper() == 'U':
                s = 'Unknown'
            elif w[1].upper() == 'P':
                s = 'Pediatric'
            else:
                s = 'Adult'
            t[i] = f"AGE('{s}')"
        else: # deal with prefix operators
            if w[0] in '<>' and w[1] in '0123':
                g = w[2:].upper() # extract gene name
                gg = ''
            else:
                g = w[1:].upper() # extract gene name
                gg = ('chr'+g).replace('P','p').replace('Q','q')
# 
            if g in deps.columns or g in cnvs.columns:

                if w[0] == '+': # amplification
                    t[i] = f"AMP('{g}')"
                elif w[0] == '-': # gene deletion
                    t[i] = f"DEL('{g}')"
                elif w[0] == '>': # over-expression
                    if w[1] in '0123':
                        s = '0123'.index(w[1])
                        t[i] = f"OX('{g}',{s})"
                    else:
                        t[i] = f"OX('{g}',-1)" # -ve signifies use pre-calc thresholds
                elif w[0] == '<': # under-expression
                    if w[1] in '0123':
                        s = '0123'.index(w[1])
                        t[i] = f"UX('{g}',{s})"
                    else:
                        t[i] = f"UX('{g}',-1)" # -ve signifies use pre-calc thresholds             
                elif w[0] == '=': # diploid eg wild type copy number
                    t[i] = f"DIP('{g}')"
                elif w[0] == '!': # gene disruption
                    t[i] = f"DIS('{g}')"
                elif w[0] == '$': # mutationally wild-type
                    t[i] = f"WT('{g}')"
                elif w[0] == '~': # gene fusion
                    t[i] = f"FUS('{g}')"
                    
            elif gg.split('.')[0] in phys_chro.columns:

                if w[0] == '+': # amplification of locus
                    t[i] = f"LAMP('{gg}')"
                elif w[0] == '-': # gene deletion
                    t[i] = f"LDEL('{gg}')"
# gene exists but has no mutation or CNV data recorded
            elif g in list(hgnc.index) and g not in list(muts['HugoSymbol']): 
                # print(f'gene {g} is valid but has no mutation or CNV data recorded')
                return ''
            else:
                print(f'function prefix {w[0]} and/or gene name {g} not valid')
                return ''
                
    return ' '.join(t)

# =============================================================================
#  ██████╗███████╗██╗     ██╗         ███╗   ███╗██╗   ██╗████████╗    ███████╗██╗   ██╗███╗   ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗███████╗
# ██╔════╝██╔════╝██║     ██║         ████╗ ████║██║   ██║╚══██╔══╝    ██╔════╝██║   ██║████╗  ██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║██╔════╝
# ██║     █████╗  ██║     ██║         ██╔████╔██║██║   ██║   ██║       █████╗  ██║   ██║██╔██╗ ██║██║        ██║   ██║██║   ██║██╔██╗ ██║███████╗
# ██║     ██╔══╝  ██║     ██║         ██║╚██╔╝██║██║   ██║   ██║       ██╔══╝  ██║   ██║██║╚██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║╚════██║
# ╚██████╗███████╗███████╗███████╗    ██║ ╚═╝ ██║╚██████╔╝   ██║       ██║     ╚██████╔╝██║ ╚████║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║███████║
#  ╚═════╝╚══════╝╚══════╝╚══════╝    ╚═╝     ╚═╝ ╚═════╝    ╚═╝       ╚═╝      ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝
#                                                                                                                                                
# =============================================================================

def CellMuts(c): # returns list of gene profile symbols for mutations in a specified cell line
    m = []
# get mutation records for this cell line
    mgenes = muts.loc[(muts['ModelID'] == c) & (muts['VariantInfo'] != 'SILENT')].reset_index()
    
# loop through building mutation formulae
    for i in mgenes.index:
        g = mgenes.iloc[i]['HugoSymbol']
        t = mgenes.iloc[i]['CCLEDeleterious']
        p = mgenes.iloc[i]['ProteinChange']
        # print(g,t,p)
        if t :
            m.append('!'+g)
        elif type(p) == str:
            m.append('/'+g+':'+p[2:])
    return m
# 
def CellCNVS(c):# returns list of gene profile symbols for CNV del or amp in a specified cell line
# get CNVS values for this cell line
    mc = cnvs.loc[cnvs['ModelID'] == c].drop('ModelID', axis=1)
# find genes i.e. columns where the CNV value indicates deletion
    gdel = ['-'+g for g in list(mc.columns[mc.iloc[0] < log2(0.87/2+1)])]
# find genes i.e. columns where the CNV value indicates amplification
    gamp = ['+'+g for g in list(mc.columns[mc.iloc[0] > log2(3.36/2+1)])]
    return gdel+gamp

# def LociCNVS(g) # aggregates genes in a list into cytological bands

# =============================================================================
#  ██████╗ ███████╗███╗   ██╗███████╗    ███████╗███████╗████████╗    ███████╗██╗   ██╗███╗   ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗███████╗
# ██╔════╝ ██╔════╝████╗  ██║██╔════╝    ██╔════╝██╔════╝╚══██╔══╝    ██╔════╝██║   ██║████╗  ██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║██╔════╝
# ██║  ███╗█████╗  ██╔██╗ ██║█████╗      ███████╗█████╗     ██║       █████╗  ██║   ██║██╔██╗ ██║██║        ██║   ██║██║   ██║██╔██╗ ██║███████╗
# ██║   ██║██╔══╝  ██║╚██╗██║██╔══╝      ╚════██║██╔══╝     ██║       ██╔══╝  ██║   ██║██║╚██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║╚════██║
# ╚██████╔╝███████╗██║ ╚████║███████╗    ███████║███████╗   ██║       ██║     ╚██████╔╝██║ ╚████║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║███████║
#  ╚═════╝ ╚══════╝╚═╝  ╚═══╝╚══════╝    ╚══════╝╚══════╝   ╚═╝       ╚═╝      ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝
# =============================================================================

def GeneAggGOF(g): # returns minimised profile string for missense GOFs for gene
    a = sorted(GetGOF(g))
    for i in range(len(a)):
        for j in range(i+1,len(a)):
            if a[j][:-1] == a[i][:-1]:
                a[i] = a[i][:-1]+'?'
                a[j] = a[j][:-1]+'?'
    a = list(set(a))
    if len(a) == 0:
        return ''
    # p = '/'+g+':'
    p = ''
    for i in a:
        p += i[2:]+','
    return p[:-1]

def GeneAggLOF(g): # returns minimised profile string for missense LOFs for gene
    a = sorted(GetLOF(g))
    for i in range(len(a)):
        for j in range(i+1,len(a)):
            if a[j][:-1] == a[i][:-1]:
                a[i] = a[i][:-1]+'?'
                a[j] = a[j][:-1]+'?'
    a = list(set(a))
    if len(a) == 0:
        return ''
    p = ''
    for i in a:
        p += i[2:]+','
    return p[:-1]

def GetGOF(g=''):
    if g == '':
        a = muts.loc[(muts['LikelyGoF']) &
                     (muts['VariantInfo'] == 'MISSENSE')][['HugoSymbol','ProteinChange']]
        return list(set(zip(a['HugoSymbol'],a['ProteinChange'])))
    else:
        return(list(set(muts.loc[(muts['LikelyGoF'] ) &
                                 (muts['HugoSymbol'] == g) &
                                 (muts['VariantInfo'] == 'MISSENSE'),'ProteinChange'])))
        
def GetLOF(g=''):
    if g == '':
        a = muts.loc[(muts['LikelyLoF']) &
                     (muts['VariantInfo'] == 'MISSENSE')][['HugoSymbol','ProteinChange']]
        return list(set(zip(a['HugoSymbol'],a['ProteinChange'])))
    else:
        return(list(set(muts.loc[(muts['LikelyLoF'] ) &
                                 (muts['HugoSymbol'] == g) &
                                 (muts['VariantInfo'] == 'MISSENSE'),'ProteinChange'])))
                                                                                                                                            
def UNI():
    # return set(set(info['ModelID']) | set(muts['ModelID']) | set(cnvs['ModelID']) | set(expr['ModelID'])) 
    return set(info['ModelID'])
def AMP(g):
    if g in cnvs.columns:
        return set(cnvs.loc[(cnvs[g] > log2(3.36/2+1)),'ModelID'])
    else :
        return set()    
def DIP(g):
    if g in cnvs.columns:
        return set(cnvs.loc[(cnvs[g] < log2(2.64/2+1)) &
                           (cnvs[g] > log2(1.32/2+1)),'ModelID'])
    else :
        return set()   
def DEL(g):
    if g in cnvs.columns:
        return set(cnvs.loc[(cnvs[g] < log2(0.87/2+1)),'ModelID'])
    else :
        return set()
def SIL(g):
    return set(muts.loc[(muts['HugoSymbol'] == g) & (muts['VariantInfo'] == 'SILENT'),'ModelID'])    
def DIS(g):
    return set(muts.loc[(muts['HugoSymbol'] == g) & (muts['CCLEDeleterious'] | muts['LikelyLoF']),'ModelID'])

def SIG(s): # return set of cells with non-zero occurence of specified signature
    # print(f's = {s}')
    if s in sigs.columns[1:-1]:
        return set(sigs[(sigs[s].div(sigs['Burden']) > sig_thresh)]['SAMPLES'])
    else:
        rr = set()
        # print(f's is still {s}, SigGroups[s] = {SigGroups[s]}')
              
        for i in SigGroups[s]:
            # print(f' i = {i}')
            rr |= set(sigs[(sigs[i].div(sigs['Burden']) > sig_thresh)]['SAMPLES'])
        return rr
            
def WT(g): # anything without a mutation or deletion recorded for gene g or a silent mutation
    return UNI() ^ (MUT(g) | DEL(g) | AMP(g))
# def DEML(g,c): # DEfined Mutation with wild card substitions allowed
#     s2 = []
#     m = re.compile(r'p#'+c.replace('?','[ACDEFGHIKLMNPQRSTVWY]{1}'))
#     # print(m)
#     s1 = muts.loc[(muts['HugoSymbol'] == g) & (muts['VariantInfo'] != 'silent')].reset_index()
#     for i in s1.index:
#         a = str(s1.iloc[i]['ProteinChange']).replace('.','#')
#         b = m.search(a)
#         if b != None:              
#             # print(f"matching {b.group()[2:]} in {s1.iloc[i]['ModelID']}")
#             s2.append(s1.iloc[i]['ModelID'])
#     return(set(s2))
def DEM(g,c): # DEfined Mutation list with wild card substitions allowed
# substitute, deblank and split mutation list
    # print(c)
    l = c.replace('?','[ACDEFGHIKLMNPQRSTVWY]{1}').replace(' ','').split(',')
    # print(l)
    s1 = muts.loc[(muts['HugoSymbol'] == g) & (muts['VariantInfo'] != 'SILENT')].reset_index()
    s2 = []
    for j in l: # loop over list
        # print(j)
        if j !='' :
            m = re.compile(r'p#'+j)
            # print(m)
            for i in s1.index:
                a = str(s1.iloc[i]['ProteinChange']).replace('.','#')
                b = m.search(a)
                if b != None:              
                    # print(f"matching {b.group()[2:]} in {s1.iloc[i]['ModelID']}")
                    s2.append(s1.iloc[i]['ModelID'])
    return(set(s2))

    # return set(muts.loc[(muts['HugoSymbol'] == g) & (muts['ProteinChange'] == 'p.'+c),'ModelID'])
def MUT(g): # any non-silent mutation
    return set(muts.loc[(muts['HugoSymbol'] == g) & (muts['VariantInfo'] != 'SILENT'),'ModelID'])    

def FIL(g): # make set from file
    ff = open(g,'r')
    s = ff.read()[:-1].split(',')
    ff.close()
    return set(s)

def TYP(c): # make set from cancer type
    return set(info.loc[info['OncotreeLineage'] == c]['ModelID'])

def SEX(s):
    return set(info.loc[info['Sex'] == s]['ModelID'])

def AGE(s):
    return set(info.loc[info['AgeCategory'] == s]['ModelID'])

def FUS(g): # set of cells in which gene g is involved in a fusion
    r = set(fuse[fuse['RightGene'] == g]['ModelID'])  
    l = set(fuse[fuse['LeftGene'] == g]['ModelID'])
    return r | l

def LDEL(gg): # make set from match of cnvdels to locus
    s = []
    w = gg.split('.')
    # print(f'w {w}')
    band = str(w[0])
    if len(w) > 1:
        subband = w[1]
        # print(f'BAND : {band} - SUBBAND : {subband}')
    else:
        subband = ''
    for c in Cdict:
        d = set(Cdict[c].cnvdel)
        l = set(phys_chro[band]['geneSymbols'])
        f = d & l # main band common set
        # print(f'c {c} f {f}')

        if len(f) == 0:
            continue # no band matches, ignore
# now test if subband filter is needed
        if subband == '':
            s.append(c)
        else: # look for band.subband match
            # print(f'gg {gg} subband {subband}')
            for g in f:
                # print(f'Looking for subband for gene {g}')
                if g in hgnc.index:
                    gl = hgnc.loc[g]['location']
                    # print(f'location {gl}')
                    if gl == gg[3:]: # strip 'chr'
                        s.append(c)
                        break
                    
    return set(s)

def LAMP(gg): # make set from match of cnvdels to locus
    s = []
    w = gg.split('.')
    band = w[0]
    if len(w) > 1:
        subband = w[1]
        # print(f'SUBBAND : {subband}')
    else:
        subband = ''
    for c in Cdict:
        d = set(Cdict[c].cnvamp)
        l = set(phys_chro[band]['geneSymbols'])
        f = d & l # main band common set
        # print(f'c {c} f {f}')

        if len(f) == 0:
            continue # no band matches, ignore
# now test if subband filter is needed
        if subband == '':
            s.append(c)
        else: # look for band.subband match
            # print(f'gg {gg} subband {subband}')
            for g in f:
                # print(f'Looking for subband for gene {g}')
                if g in hgnc.index:
                    gl = hgnc.loc[g]['location']
                    # print(f'location {gl}')
                    if gl == gg[3:]: # strip 'chr'
                        s.append(c)
                        break
                    
    return set(s)

def OX(g,s):
    if s < 0: # using pre-calculated high/low threshold assigments
        return(set(list(expr[expr[g] >= exval.loc[g]['High']]['ModelID'])))
    if g in expr.columns:
        h = expr[g].mean()+s*expr[g].std()
        return set(expr[expr[g] >= h]['ModelID'])
    else :
        return set()    

def UX(g,s):
    if s < 0: # using pre-calculated high/low threshold assigments
        return(set(expr[expr[g] <= exval.loc[g]['Low']]['ModelID']))
    if g in expr.columns:
        h = expr[g].mean()-s*expr[g].std()
        return set(expr[expr[g] <= h]['ModelID'])
    else :
        return set()    
   

def prof_tree(p,t):  # recursive tree walk 
    for i in profiles:
        if p in profiles[i].upper().split():
            t.append(i)
            prof_tree(i,t)
# =============================================================================
# 
# calculate Cohen D and matching sample size            
#             
# =============================================================================

def DandCount(p,q,targ,deps):
    # print(f'Profile : {prof}')
    # p = list(eval(prof_parse(prof)))
    # q = list(eval(prof_nullhype(prof_flatten(prof))))
    hype = []
    null = []
    D = 0
    
    for c in p: # loop through set of cells matching profile
        # print(f'cell line : {c}')
        dep = deps.loc[deps.index == c][targ]
        if len(dep) > 0:
            hype.append(float(dep.iloc[0]))
            
    for c in q: # loop through set of cells forming null hypothesis
        dep = deps.loc[deps.index == c][targ]
        if len(dep) > 0:
            null.append(float(dep.iloc[0]))
# catch empty slices
    if len(hype) > 0 and len(null) > 0:
        D = (np.array(hype).mean()-np.array(null).mean())/np.array(hype+null).std()
    return (D,len(hype))                        
            
            
# =============================================================================
# calculate P-value for a target gene against a genetic profile
# 
# =============================================================================

def AltProfPairPvalue(targ,prof,p,q,deps,th,null_t): # p is the set of cells with the profile, q is the null set
    cont = [[0,0],[0,0]] # initialise contingency table

# get count of p matches above thresh

    cont[0][0] = sum(deps[deps.index.isin(p)][targ]>thresh)
    cont[1][0] = len(p) - cont[0][0]
    
# get count of q matches above thresh

    cont[0][1] = sum(deps[deps.index.isin(q)][targ]>thresh)
    cont[1][1] = len(q) - cont[0][1]
    
    pv = stats.boschloo_exact(cont,alternative='greater').pvalue

    return [targ,pv]

# =============================================================================
# calculate p-value and D-valuefor a target gene against a genetic profile
# 
# =============================================================================

def AltProfPairPD(targ,p,q,deps,th): # p is the set of cells with the profile, q is the null set
    cont = [[0,0],[0,0]] # initialise contingency table

# get count of p matches above thresh

    a = deps[deps.index.isin(p)][targ]
    b = deps[deps.index.isin(q)][targ]

    cont[0][0] = sum(a > th)
    cont[1][0] = len(p) - cont[0][0]
    
# get count of q matches above thresh

    cont[0][1] = sum(b > th)
    cont[1][1] = len(q) - cont[0][1]
    
    pv = stats.boschloo_exact(cont,alternative='greater').pvalue
    with np.errstate(divide='ignore', invalid='ignore'):
        D = (a.mean()-b.mean())/(np.concatenate((a,b)).std())

    return [targ,pv,D]




def ProfPairPvalue(targ,prof,p,q,deps,th,null_t): # p is the set of cells with the profile, q is the null set
    cont = [[0,0],[0,0]] # initialise contingency table
    
    for c in p: # loop through set of cells matching profile
        dep = deps.loc[deps.index == c][targ]
        if len(dep) > 0:
            if float(dep.iloc[0]) >=th :
                cont[0][0] +=1
            else :
                cont[1][0] +=1
                
# apply trim factor to null set
    if null_t > 0:
        print(f'q before trim {len(q)}')
        nt = min(len(q),len(p)*null_t)
        q = random.choices(q,k=int(nt))
        print(f'q after trim {len(q)}')

    for c in q: # loop through set of cells forming null hypothesis
        dep = deps.loc[deps.index == c][targ]
        if len(dep) > 0:
            if float(dep.iloc[0]) >=th :
                cont[0][1] +=1
            else :
                cont[1][1] +=1
    # if cont[0][0] > 5 and cont[1][0] > 5: # minimum counts for valid ChiSq
    pv = stats.boschloo_exact(cont,alternative='greater').pvalue
# =============================================================================
#     else:
#         p = 1.0
# =============================================================================
        # print(f'Too few observations for {targ} with {prof}')
    # print(f'{targ} : {prof} p = {p:.5}')
    return [targ,pv]
    
    
# =============================================================================
# calculate P-value for a defined GeneA/GeneB pair
# 
# =============================================================================
def PairPvalue(GeneA,GeneB,CellDict,deps,th,usemuts,usecnv,useexp,mask = None):
    mutated = []
    notmutated = []
    cont = [[0,0],[0,0]]
    if mask == None:
        mask = UNI()
    
    if usemuts:
        for c in CellDict:
            if c not in mask:
                continue
            dep = deps.loc[c][GeneA]
            if dep < 0.0001: # ignore NaNs
                continue
    #        print(f'GeneA {GeneA} Cell {c} dep {dep}')
    
            if GeneB in CellDict[c].dgenes: # damaging mutation ?
                mutated.append(dep)
                if dep >=th :
                    cont[0][0] +=1
                else :
                    cont[1][0] +=1
            elif GeneB not in CellDict[c].mgenes: # not a missense
                notmutated.append(dep)
                if dep >=th :
                    cont[0][1] +=1
                else :
                    cont[1][1] +=1
                    
    if usecnv:
        for c in CellDict:
            if c not in mask:
                continue

            dep = deps.loc[c][GeneA]
            if dep < 0.0001: # ignore NaNs
                continue

    #        print(f'GeneA {GeneA} Cell {c} dep {dep}')
    
            if GeneB in CellDict[c].cnvdel: # # deep deletion ?
                mutated.append(dep)
                if dep >=th :
                    cont[0][0] +=1
                else :
                    cont[1][0] +=1
            elif GeneB not in CellDict[c].cnvamp: # not an amplification
                notmutated.append(dep)
                if dep >=th :
                    cont[0][1] +=1
                else :
                    cont[1][1] +=1
    
    if useexp:
        for c in CellDict:
            if c not in mask:
                continue

            dep = deps.loc[c][GeneA]
            if dep < 0.0001: # ignore NaNs
                continue

    #        print(f'GeneA {GeneA} Cell {c} dep {dep}')
    
            if GeneB in CellDict[c].under: # significant underexpression ?
                mutated.append(dep)
                if dep >=th :
                    cont[0][0] +=1
                else :
                    cont[1][0] +=1
            elif GeneB not in CellDict[c].over: # not an amplification
                notmutated.append(dep)
                if dep >=th :
                    cont[0][1] +=1
                else :
                    cont[1][1] +=1
    
    p = stats.boschloo_exact(cont,alternative='greater').pvalue
    # print(f'{GeneA} : {GeneB} p = {p:.5}')
    return p

# =============================================================================
# generate mutated/not-mutated dependency distributions
# for a defined target - Genetic Profile pair
# 
# =============================================================================
def ProfPairDist(targ,p,q,deps,null_t): # p is the set of cells with the profile, q is the null set:
    hype = []
    null = []
    
    for c in p: # loop through set of cells matching profile
        # print(f'cell line : {c}')
        dep = deps.loc[deps.index == c][targ]
        if len(dep) > 0:
            # print(f'dep : {dep}')
            hype.append(float(dep.iloc[0]))
            
# apply trim factor to null set
    if null_t > 0:
        print(f'q before trim {len(q)}')
        nt = min(len(q),len(p)*null_t)
        q = random.choices(q,k=int(nt))
        print(f'q after trim {len(q)}')

    for c in q: # loop through set of cells forming null hypothesis
        dep = deps.loc[deps.index == c][targ]
        if len(dep) > 0:
            null.append(float(dep.iloc[0]))
# purge NaNs

    return [[i for i in hype if i > 0.0000],[i for i in null if i > 0.0000]]

# =============================================================================
# 
#    parallelisable D values
#   
#   
# =============================================================================
def ProfPairDval(prof,deps):
    return

# =============================================================================
# generate mutated/not-mutated dependency distributions
# for a defined GeneA/GeneB pair
# 
# =============================================================================

def NewPairDist(geneB, geneA, deps):
    
    c_mut = DIS(geneB) | DEL(geneB)
    if c_mut == set():
        return []  # no mutations or deletions of geneB
    c_wt = WT(geneB)
           
    d1 = deps[geneA].reset_index()
    
    mut = d1[d1['ModelID'].isin(c_mut)]
    wt = d1[d1['ModelID'].isin(c_wt)]
    return [mut,wt]


def PairDist(GeneA,GeneB,CellDict,deps,usemuts,usecnv,useexp):
    mutated = []
    notmutated = []
    if usemuts:
        for c in CellDict:
            dep = deps.loc[c][GeneA]
            if GeneB in CellDict[c].dgenes: # damaging mutation ?
                mutated.append(dep)
            elif GeneB not in CellDict[c].mgenes: # not a missense
                notmutated.append(dep)
    if usecnv:
        for c in CellDict:
            dep = deps.loc[c][GeneA]
            if GeneB in CellDict[c].cnvdel: # deep deletion ?
                mutated.append(dep)
            elif GeneB not in CellDict[c].cnvamp: # not an amplification
                notmutated.append(dep)
                
    if useexp:
        for c in CellDict:
            dep = deps.loc[c][GeneA]
    #        print(f'GeneA {GeneA} Cell {c} dep {dep}')
    
            if GeneB in CellDict[c].under: # significant underexpression ?
                mutated.append(dep)
            elif GeneB not in CellDict[c].over: # not an amplification
                notmutated.append(dep)
    
    return [[i for i in mutated if i > 0.0000],[i for i in notmutated if i > 0.0000]]

    # return [mutated,notmutated]

# =============================================================================
# calculate P-value routine for parallel search of all mutated genes
#  ██████╗ ███████╗███╗   ██╗███████╗██╗  ██╗██╗████████╗███████╗
# ██╔════╝ ██╔════╝████╗  ██║██╔════╝██║  ██║██║╚══██╔══╝██╔════╝
# ██║  ███╗█████╗  ██╔██╗ ██║█████╗  ███████║██║   ██║   ███████╗
# ██║   ██║██╔══╝  ██║╚██╗██║██╔══╝  ██╔══██║██║   ██║   ╚════██║
# ╚██████╔╝███████╗██║ ╚████║███████╗██║  ██║██║   ██║   ███████║
#  ╚═════╝ ╚══════╝╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝╚═╝   ╚═╝   ╚══════╝
# =============================================================================

def AltGeneHits(g,targ,CellDict,deps,th,fl,pl,usemuts,usecnv,useexp):
    mutated = []
    notmutated = []
    cont = [[0,0],[0,0]]
    prob = 0.0

    if usemuts:
# get set of cells in which gene g is mutated
        p = DIS(g)
        
# get set of cells in which gene g has no non-synonymous mutations

        q = UNI() ^ MUT(g)
        
        mutated = list(deps[deps.index.isin(p)][targ])
        notmutated = list(deps[deps.index.isin(q)][targ])
        
        nmut = len(mutated)
        nnotmut = len(notmutated)
        frac = nmut/(nmut+nnotmut)
        
        cont[0][0] = len([x for x in mutated if x > thresh])
        cont[1][0] = len(mutated) - cont[0][0]
        
        cont[0][0] = len([x for x in notmutated if x > thresh])
        cont[1][0] = len(notmutated) - cont[0][0]

    # only calculate p for genes with more than the specified frequency of mutation    

        if frac >= fl:
            prob = stats.boschloo_exact(cont,alternative='greater').pvalue
            
# is it a hit ?

        if prob <= pl:
#            HitList[k] = Hit(Gdict[k].p,Gdict[k].frac,mutated,notmutated)
            return [g,prob,frac,mutated,notmutated]
        else:
            return False
    return False
            

def GeneHits(k,targ,CellDict,deps,th,fl,pl,usemuts,usecnv,useexp):
    
    
    mutated = []
    notmutated = []
    cont = [[0,0],[0,0]]
    
# for each cell line
    
    if usemuts:
        for c in CellDict:
            maxdep = 0.0
            for t in targ:
                dep = deps.loc[c][t] # change to take list of targets
                if dep > maxdep:
                    maxdep = dep
        #            if debug : print(f'Gene {g} Cell {c} dep {dep}')
            if k in CellDict[c].dgenes: # damaging mutation or ?
                mutated.append(maxdep)
                if maxdep >=th :
                    cont[0][0] +=1
                else :
                    cont[1][0] +=1
            elif k not in CellDict[c].mgenes: # not a missense
                notmutated.append(maxdep)
                if maxdep >=th :
                    cont[0][1] +=1
                else :
                    cont[1][1] +=1

    if usecnv:
        for c in CellDict:
            maxdep = 0.0
            for t in targ:
                dep = deps.loc[c][t] # change to take list of targets
                if dep > maxdep:
                    maxdep = dep
        #            if debug : print(f'Gene {g} Cell {c} dep {dep}')
            if k in CellDict[c].cnvdel: # deep deletion ?
                mutated.append(maxdep)
                if maxdep >=th :
                    cont[0][0] +=1
                else :
                    cont[1][0] +=1
            elif k not in CellDict[c].cnvamp: # not an amplification
                notmutated.append(maxdep)
                if maxdep >=th :
                    cont[0][1] +=1
                else :
                    cont[1][1] +=1

    if useexp:
        for c in CellDict:
            maxdep = 0.0
            for t in targ:
                dep = deps.loc[c][t] # change to take list of targets
                if dep > maxdep:
                    maxdep = dep
        #            if debug : print(f'Gene {g} Cell {c} dep {dep}')
            if k in CellDict[c].under: # significant underexpression ?
                mutated.append(maxdep)
                if maxdep >=th :
                    cont[0][0] +=1
                else :
                    cont[1][0] +=1
            elif k not in CellDict[c].over: # not an overexpression
                notmutated.append(maxdep)
                if maxdep >=th :
                    cont[0][1] +=1
                else :
                    cont[1][1] +=1
        
    # now do stats
    
    nmut = len(mutated)
    nnotmut = len(notmutated)
    frac = nmut/(nmut+nnotmut)
    
# only calculate p for genes with more than the specified frequency of mutation    
    
    if frac >= fl:
        p = stats.boschloo_exact(cont,alternative='greater').pvalue
#        print(f' Frequent Gene {g} contingency : {cont} : p = {Gdict[g].p}')

# is it a hit ?

        if p <= pl:
#            HitList[k] = Hit(Gdict[k].p,Gdict[k].frac,mutated,notmutated)
            return [k,p,frac,mutated,notmutated]
        else:
            return False
    return False
            
# routines for saving and loading database
# =============================================================================
# ██████╗ ██████╗     ██╗      ██████╗  █████╗ ██████╗    ██╗   ███████╗ █████╗ ██╗   ██╗███████╗
# ██╔══██╗██╔══██╗    ██║     ██╔═══██╗██╔══██╗██╔══██╗   ██║   ██╔════╝██╔══██╗██║   ██║██╔════╝
# ██║  ██║██████╔╝    ██║     ██║   ██║███████║██║  ██║████████╗███████╗███████║██║   ██║█████╗  
# ██║  ██║██╔══██╗    ██║     ██║   ██║██╔══██║██║  ██║██╔═██╔═╝╚════██║██╔══██║╚██╗ ██╔╝██╔══╝  
# ██████╔╝██████╔╝    ███████╗╚██████╔╝██║  ██║██████╔╝██████║  ███████║██║  ██║ ╚████╔╝ ███████╗
# ╚═════╝ ╚═════╝     ╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚═════╝ ╚═════╝  ╚══════╝╚═╝  ╚═╝  ╚═══╝  ╚══════╝
#                                                                                                
# =============================================================================
def DB_save(fgene,fcell):

    for i in Gdict:
        g = Gdict[i]
        fgene.write(f'{g.name}\n')
        fgene.write(f'{len(g.cells)}\n')
        for j in g.cells:
            fgene.write(f'{j}\n')
        fgene.write(f'{len(g.cnvdel)}\n')
        for j in g.cnvdel:
            fgene.write(f'{j}\n')
        fgene.write(f'{len(g.cnvamp)}\n')
        for j in g.cnvamp:
            fgene.write(f'{j}\n')
        fgene.write(f'{len(g.over)}\n')
        for j in g.over:
            fgene.write(f'{j}\n')
        fgene.write(f'{len(g.under)}\n')
        for j in g.under:
            fgene.write(f'{j}\n')


    fgene.close()
    
    for i in Cdict:
        c = Cdict[i]
        fcell.write(f'{c.name}\n')
        fcell.write(f'{c.tissue}\n')
        fcell.write(f'{c.cancer}\n')
        fcell.write(f'{len(c.dgenes)}\n')
        for j in c.dgenes:
            fcell.write(f'{j}\n')
        fcell.write(f'{len(c.mgenes)}\n')
        for j in c.mgenes:
            fcell.write(f'{j}\n')
            
        fcell.write(f'{len(c.cnvdel)}\n')
        for j in c.cnvdel:
            fcell.write(f'{j}\n')
        fcell.write(f'{len(c.cnvamp)}\n')
        for j in c.cnvamp:
            fcell.write(f'{j}\n')

        fcell.write(f'{len(c.over)}\n')
        for j in c.over:
            fcell.write(f'{j}\n')
        fcell.write(f'{len(c.under)}\n')
        for j in c.under:
            fcell.write(f'{j}\n')
            
    fcell.close()

def DB_load(fgene,fcell):
# do genes
    while True:
        gname = fgene.readline()[:-1].strip()
        if gname == '':
            break
# should be a gene name so add to dictionary
        Gdict[gname] = Gene(gname)
# get cell line count
        ncell = int(fgene.readline()[:-1].strip())
        for i in range(ncell):
            Gdict[gname].cells.append(fgene.readline()[:-1].strip())

        ncnvdel = int(fgene.readline()[:-1].strip())
        for i in range(ncnvdel):
            Gdict[gname].cnvdel.append(fgene.readline()[:-1].strip())
        ncnvamp = int(fgene.readline()[:-1].strip())
        for i in range(ncnvamp):
            Gdict[gname].cnvamp.append(fgene.readline()[:-1].strip())
            
        nover = int(fgene.readline()[:-1].strip())
        for i in range(nover):
            Gdict[gname].over.append(fgene.readline()[:-1].strip())
        nunder = int(fgene.readline()[:-1].strip())
        for i in range(nunder):
            Gdict[gname].under.append(fgene.readline()[:-1].strip())
            
    fgene.close()
    print(f'Data for {len(Gdict)} genes loaded')
# do cells
    while True:
        cname = fcell.readline()[:-1].strip()
        if cname == '':
            break
# should be a gene name so add to dictionary
        Cdict[cname] = Cell(cname)
        Cdict[cname].tissue = fcell.readline()[:-1].strip()
        Cdict[cname].cancer = fcell.readline()[:-1].strip()

# get damaged genes count
        ngene = int(fcell.readline()[:-1].strip())
        for i in range(ngene):
            Cdict[cname].dgenes.append(fcell.readline()[:-1].strip())
# get mutated genes count
        ngene = int(fcell.readline()[:-1].strip())
        for i in range(ngene):
            Cdict[cname].mgenes.append(fcell.readline()[:-1].strip())
# get CNV deleted count
        ngene = int(fcell.readline()[:-1].strip())
        for i in range(ngene):
            Cdict[cname].cnvdel.append(fcell.readline()[:-1].strip())
# get CNV amplified count
        ngene = int(fcell.readline()[:-1].strip())
        for i in range(ngene):
            Cdict[cname].cnvamp.append(fcell.readline()[:-1].strip())
# get overexpression count
        nover = int(fcell.readline()[:-1].strip())
        for i in range(nover):
            Cdict[cname].over.append(fcell.readline()[:-1].strip())
# get underexpression count
        nunder = int(fcell.readline()[:-1].strip())
        for i in range(nunder):
            Cdict[cname].under.append(fcell.readline()[:-1].strip())
        
            
    fcell.close()
    print(f'Data for {len(Cdict)} cells loaded')
    
# =============================================================================
# 
# # ██████╗██╗  ██╗██████╗  ██████╗ ██╗      ██████╗  ██████╗
# ██╔════╝██║  ██║██╔══██╗██╔═══██╗██║     ██╔═══██╗██╔════╝
# ██║     ███████║██████╔╝██║   ██║██║     ██║   ██║██║     
# ██║     ██╔══██║██╔══██╗██║   ██║██║     ██║   ██║██║     
# ╚██████╗██║  ██║██║  ██║╚██████╔╝███████╗╚██████╔╝╚██████╗
#  ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝ ╚═════╝  ╚═════╝
#                                                           
# =============================================================================

def ChroLoc(hgid,hgnc,ensem):


    if hgid in hgnc.index:
        egid = hgnc.loc[hgid]['ensembl_gene_id']
    elif hgid in list(hgnc['prev_symbol']):
        egid = hgnc.iloc[hgnc['prev_symbol'].tolist().index(hgid)]['ensembl_gene_id']
        hgid = hgnc[hgnc['prev_symbol'] == hgid].index[0]
    else:
        print(f'GeneB : {hgid} not found in HUGO symbol lists - ignored')
        return None
# print(f'HUGO : {hgid} = ENSG : {egid}')
    if egid in ensem.index:
        etid = ensem.loc[egid]['transcript_stable_id']
    else:
        print(f'GeneA : {egid} for HUGO symbol {hgid} not found in ENSEMBLE transcripts - ignored')
        return None
    # print(f'ENSG: {egid} = ENST : {etid}')

# use Ensemble API to get data                                

    server = "https://rest.ensembl.org"
    ext = f"/map/cdna/{etid}/1..1?"
    print(f'Fetching {server+ext}')
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
     
    if not r.ok:
      r.raise_for_status()
      sys.exit()
     
    d = r.json()
    chro = d['mappings'][0]['seq_region_name']
    nucl = int(d['mappings'][0]['start'])
    cyto = hgnc.loc[hgid]['location_sortable']
    return (chro,nucl,cyto)

# =============================================================================
# ███╗   ███╗     █████╗     ██╗    ███╗   ██╗
# ████╗ ████║    ██╔══██╗    ██║    ████╗  ██║
# ██╔████╔██║    ███████║    ██║    ██╔██╗ ██║
# ██║╚██╔╝██║    ██╔══██║    ██║    ██║╚██╗██║
# ██║ ╚═╝ ██║    ██║  ██║    ██║    ██║ ╚████║
# ╚═╝     ╚═╝    ╚═╝  ╚═╝    ╚═╝    ╚═╝  ╚═══╝
#                                             
# =============================================================================
if __name__ == '__main__':
    print("Number of processors: ", mp.cpu_count())

    debug = False
    
    # define data structures
    
    Gdict = {}
    Cdict = {}
    profiles = {}

    GeneDB_exists = True
    CellDB_exists = True
    
    muts_loaded = False
    cnvs_loaded = False
    info_loaded = False
    deps_loaded = False
    hgnc_loaded = False
    phys_chro_loaded = False
    expr_loaded = False
    exval_loaded = False
    prof_loaded = False
    oncogenes_loaded = False
    chromlen_loaded = False
    ensem_loaded = False
    null_trim = 0.0

    
    seaborn_plot = False
    
    spec = False
    a = input('Specify files and parameters ? ([N]/Y] : ').upper()
    if len(a) > 0 and a[0] == 'Y':
        spec = True
        
    # load all source data
    
    f = default_dep_dat_file = 'CRISPRGeneDependency.csv'
    if spec:
        f = input(f'Gene dependency file [{default_dep_dat_file}] : ')
        if f.strip() == '' :
            f = default_dep_dat_file
    print(f'Loading {f}')
    
    deps = pd.read_csv(f).fillna(0.0)
    
    default_dep_dat_file = f
    
    deps_loaded = True
      
    # make cell line index to dependency data
    
    deps.set_index('ModelID', inplace = True)
    
    # edit header names to conform to remove extraneous data
    
    gg = dict(zip(list(deps.columns),[i.split()[0] for i in list(deps.columns)]))
    deps.rename(columns = gg, inplace = True)
    # print(deps)

    f = default_cell_info_file = 'Model.csv'
    if spec:
        f = input(f'Cell sample information file [{default_cell_info_file}] : ')
        if f.strip() == '' :
            f = default_cell_info_file
    print(f'Loading {f}')

    info = pd.read_csv(f, low_memory=False).fillna('')
    info_loaded = True
    info['OncotreeLineage'] = [x.upper() for x in info['OncotreeLineage']] # make cancer types upper case
    cancer_types = set(info['OncotreeLineage'])
    # print(info)

    f = default_mut_dat_file = 'OmicsSomaticMutations_GoFLoFUpdated.csv'
    if spec:
        f = input(f'Cell line mutation file [{default_mut_dat_file}] : ')
        if f.strip() == '' :
            f = default_mut_dat_file
    print(f'Loading {f}')

    muts = pd.read_csv(f, low_memory=False).fillna('')   
    muts_loaded = True

    default_mut_dat_file = f
    # print(muts)

    f = default_phys_chro_file = 'c1.all.v2023.1.Hs.json.txt' # from GSEA
    if spec:
        f = input(f'Physical chromosome set file [{default_phys_chro_file}] : ')
        if f.strip() == '' :
            f = default_phys_chro_file
    print(f'Loading {f}')

    phys_chro = pd.read_json(f).fillna('')   
    phys_chro_loaded = True
    # print(phys_chro)

# load HUGO chromosomal location data
    f = default_HUGO_file = 'hgnc_complete_set.tsv'
    if spec:
        f = input(f'HUGO Gene file [{default_HUGO_file}] : ')
        if f.strip() == '' :
            f = default_HUGO_file
    print(f'Loading {f}')
    hgnc = pd.read_table(f, low_memory=False).fillna('')
    hgnc = hgnc[['symbol','ensembl_gene_id','prev_symbol','location','location_sortable']]
    hgnc.set_index('symbol',inplace = True)
    hgnc_loaded = True
    # print(hgnc)
# correct legacy deps gene symbols to HUGO

    bad_names = set(deps.columns) & (set(hgnc.index) ^ set(deps.columns)) 
    print(f'Updating {len(bad_names)} archaic gene names in dependency data')
    for g in bad_names:
        g2 = hgnc[hgnc['prev_symbol'].str.contains(g)].reset_index()['symbol']
        if len(g2) == 0 or (g2[0] not in hgnc.index):
            # print(f'DepMap Gene name {g} not found in HUGO - ignoring it')
            continue
        else:
            # print(f'DepMap old gene name {g} replaced by new name {g2[0]}')
            deps.rename(columns={g:g2[0]}, inplace=True)
            
    f = default_cnv_dat_file = 'OmicsCNGene.csv'
    if spec:
        f = input(f'Cell line CNV file [{default_cnv_dat_file}] : ')
        if f.strip() == '' :
            f = default_cnv_dat_file
    print(f'Loading {f}')

    cnvs = pd.read_csv(f, low_memory=False).fillna(0.0)
    cnvs_loaded = True

    default_cnv_dat_file = f

    # make cell line index to CNV data

    cnvs.rename( columns={'Unnamed: 0':'ModelID'}, inplace=True )
    cnvs.set_index('ModelID', inplace = True)

    # edit header names to conform to remove extraneous data

    gg = dict(zip(list(cnvs.columns),[i.split()[0] for i in list(cnvs.columns)]))
    cnvs.rename(columns = gg, inplace = True)
    # print(cnvs)

    f = default_fus_dat_file = 'OmicsFusionFiltered.csv'
    if spec:
        f = input(f'Cell line gene fusion file [{default_fus_dat_file}] : ')
        if f.strip() == '' :
            f = default_fus_dat_file
    print(f'Loading {f}')
    
    fuse = pd.read_csv(f, low_memory=False).fillna(0.0)
    fuse_loaded = True
    fuse['LeftGene'] = [x.split()[0] for x in list(fuse['LeftGene'])]
    fuse['RightGene'] = [x.split()[0] for x in list(fuse['RightGene'])]

    f = default_expr_dat_file = 'OmicsExpressionProteinCodingGenesTPMLogp1.csv'
    if spec:
        f = input(f'Cell line expression file [{default_expr_dat_file}] : ')
        if f.strip() == '' :
            f = default_expr_dat_file
    print(f'Loading {f}')

    expr = pd.read_csv(f, low_memory=False).fillna(0.0)
    expr_loaded = True

    default_expr_dat_file = f
    
    # make cell line index to expression data

    expr.rename( columns={'Unnamed: 0':'ModelID'}, inplace=True )
    expr.set_index('ModelID', inplace = True)


    # edit header names to conform to remove extraneous data

    gg = dict(zip(list(expr.columns),[i.split()[0] for i in list(expr.columns)]))
    expr.rename(columns = gg, inplace = True)

    f = default_file = 'ExpressionThresholds.csv'
    if spec:
        f = input(f'Expresssion Threshold File [{default_file}] : ')
        if f.strip() == '' :
            f = default_file
    print(f'Loading {f}')
    
    exval = pd.read_csv(f)
    exval.set_index('HugoSymbol', inplace = True)
    exval_loaded = True

    f = default_file = 'Final_GOF_all_table.csv'
    if spec:
        
        f = input(f'GoF File [{default_file}] : ')
        if f.strip() == '' :
            f = default_file
    print(f'Loading {f}')
    gofs = pd.read_csv(f)
    gofs_loaded = True

    if spec:
        ans = input('Update GoFs with external data ? (Y,[N]) : ').upper()
        if len(ans) > 0 and ans[0] != 'N':
            print('Annotating additional GoF mutations ...')
            gof1 = sum(muts['LikelyGoF'])
            for m in gofs.index:
                gn = gofs.iloc[m]['Gene_name']
                pc = gofs.iloc[m]['Substitution']
                if len(muts.loc[(muts['HugoSymbol'] == gn) & (muts['ProteinChange'] == 'p.'+pc),'LikelyGoF']) > 0:
                    # print(f'Assigning {gn+"p."+pc} as GoF')
                    muts.loc[(muts['HugoSymbol'] == gn) & (muts['ProteinChange'] == 'p.'+pc),'LikelyGoF'] = True
        
            gof2 = sum(muts['LikelyGoF'])
            print(f'Expanded GoF count from {gof1} to {gof2}')
    # muts.loc[(muts['HugoSymbol'] == 'BRAF') & (muts['ProteinChange'] == 'p.V600E'),'VariantInfo'] = 'fish'

    f = default_signatures_file = 'Assignment_Solution_Activities.tsv'
    if spec:
        f = input(f'Mutation signatures file [{default_signatures_file}] : ')
        if f.strip() == '' :
            f = default_signatures_file
    print(f'Loading {f}')
    sigs = pd.read_table(f, low_memory=False).fillna('')
    gg = dict(zip(list(sigs.columns),[i.upper() for i in list(sigs.columns)]))
    sigs.rename(columns = gg, inplace = True)

    # # make mutational burden column
    sigs['Burden'] = sigs.sum(axis=1, numeric_only=True)
    # # normalise signature counts
    # sigs[sigs.columns[1:-1]] = sigs[sigs.columns[1:-1]].div(sigs['Burden'],axis = 0)


    f = default_pathways_file = 'c2.cp.reactome.v2023.2.Hs.symbols.gmt'
    if spec:
        f = input(f'Pathways file [{default_pathways_file}] : ')
        if f.strip() == '' :
            f = default_pathways_file
    print(f'Loading {f}')
    paths = {}
    fgmt = open(f,'r')
    while True:
        line = fgmt.readline()
        if line == '':
            break
        w = line.split()
        paths[w[0]] = set(w[2:])


    try:
        f_gene = open('AB_GeneDatabaseExpr2023Q2','r')
    except FileNotFoundError:
        GeneDB_exists = False
    try:
        f_cell = open('AB_CellDatabaseExpr2023Q2','r')
    except FileNotFoundError:
        CellDB_exists = False

# if both databases exist, just read them in

    if GeneDB_exists and CellDB_exists:
        print('Loading pre-existing Gene and Cell databases')
        DB_load(f_gene,f_cell)
        # muts = ''
        
    else:
        print('Creating new Gene and Cell databases from DepMap source files')

# create DB from DepMap source files
# =============================================================================
#  ██████╗██████╗ ███████╗ █████╗ ████████╗███████╗    ██████╗ ██████╗ 
# ██╔════╝██╔══██╗██╔════╝██╔══██╗╚══██╔══╝██╔════╝    ██╔══██╗██╔══██╗
# ██║     ██████╔╝█████╗  ███████║   ██║   █████╗      ██║  ██║██████╔╝
# ██║     ██╔══██╗██╔══╝  ██╔══██║   ██║   ██╔══╝      ██║  ██║██╔══██╗
# ╚██████╗██║  ██║███████╗██║  ██║   ██║   ███████╗    ██████╔╝██████╔╝
#  ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝   ╚═╝   ╚══════╝    ╚═════╝ ╚═════╝ 
#                                                                                                         |         |                        
# =============================================================================
        
        #load cell mutation data
              
  
#create list of cell lines and build into cell dictionary
        
        dep_cells = deps.index.tolist()
        for c in dep_cells:
            Cdict[c] = Cell(c)
        
        print('\nDependency dataset : ' ,default_dep_dat_file, 'contains ',len(Cdict),'cell lines')
               

# loop through sample info assigning tissue and cancer info where cell name matches

        for i in range(len(info)):
            cname = info.iloc[i]['ModelID']
            if cname in Cdict:
                Cdict[cname].tissue = info.iloc[i]['SampleCollectionSite']
                Cdict[cname].cancer = info.iloc[i]['OncotreeLineage']

        # loop through mutational data acquiring list of genes that are mutated in
        # any valid cell line
        
# initialise Gdict with HUGO gene list

        for g in hgnc.index:
            Gdict[g] = Gene(g)    
        
        
        t_base = timer()
        ng = 1
        
        print('Building mutation dictionary ...')
        for i in muts.index:
            cell_line = muts.iloc[i]['ModelID']
            if cell_line in Cdict : # mutation in valid cell line 
                if muts.iloc[i]['CCLEDeleterious'] : #
                    mut_gen = muts.iloc[i]['HugoSymbol']
        # add this mutated gene to the list for this cell line if not already there
                    if not (mut_gen in Cdict[cell_line].dgenes):
                        Cdict[cell_line].dgenes.append(mut_gen)
        #            print('Damaging mutation recorded in gene :',mut_gen)
                    if mut_gen not in Gdict: # not found - try prev symbol
                        m2 = hgnc[hgnc['prev_symbol'].str.contains(mut_gen)].reset_index()['symbol']
                        if len(m2) == 0 or (m2[0] not in Gdict):
                            print(f'Gene name {mut_gen} not found in HUGO - skipping')
                            continue
                        else:
                            print(f'Old gene name {mut_gen} replaced by new name {m2[0]}')
                            mut_gen = m2[0]
                            
                    Gdict[mut_gen].cells.append(cell_line)
                    
                elif muts.iloc[i]['VariantInfo'] != 'SILENT':
                    mut_gen = muts.iloc[i]['HugoSymbol']
        # add this mutated gene to the list for this cell line if not already there
                    if not (mut_gen in Cdict[cell_line].mgenes):
                        Cdict[cell_line].mgenes.append(mut_gen)
 
                if ng %100000 == 0 :
                    print(f'{ng} muts time : {timer()-t_base:.4}s')
                ng += 1
        
        print(f'Dictionary build time : {timer()-t_base:.4}s')
                   
        print(f'{len(Gdict)} genes with damaging mutations found in dependency dataset')
                   
        # sort the gene list
        
        Gdict = OrderedDict(sorted(Gdict.items()))
        
        
        
# update CNV and EXPR data to new gene names

        for g in cnvs.columns:
            if g not in hgnc.index: # look for previous symbol
                g2 = hgnc[hgnc['prev_symbol'].str.contains(g)].reset_index()['symbol']
                if len(g2) == 0 or (g2[0] not in hgnc.index):
                    print(f'CNV Gene name {g} not found in HUGO - removing it')
                    cnvs.rename(columns={g:'IGNORE'}, inplace=True)
                    continue
                else:
                    print(f'CNV old gene name {g} replaced by new name {g2[0]}')
                    cnvs.rename(columns={g:g2[0]}, inplace=True)
                    
        if 'IGNORE' in cnvs.columns:
            cnvs.drop(columns = ['IGNORE'], inplace = True)
            
# remove duplicates

        cnvs.drop(columns = [cnvs.columns[i] for i,t in enumerate(cnvs.columns.duplicated()) if t],inplace = True)

        for g in expr.columns:
            if g not in hgnc.index: # look for previous symbol
                g2 = hgnc[hgnc['prev_symbol'].str.contains(g)].reset_index()['symbol']
                if len(g2) == 0 or (g2[0] not in hgnc.index):
                    print(f'EXPR Gene name {g} not found in HUGO - removing it')
                    expr.rename(columns={g:'IGNORE'}, inplace=True)
                    continue
                else:
                    print(f'EXPR old gene name {g} replaced by new name {g2[0]}')
                    expr.rename(columns={g:g2[0]}, inplace=True)

        if 'IGNORE' in expr.columns:
            expr.drop(columns = ['IGNORE'], inplace = True)
            
        expr.drop(columns = [expr.columns[i] for i,t in enumerate(expr.columns.duplicated()) if t],inplace = True)


# Purge CNV data of cell lines without dependency data
        l1 = len(cnvs.index)
        c1 = set(cnvs.index)
        c2 = set(deps.index)
        c3 = c1 & (c1 ^ c2)
        cnvs.drop(list(c3),inplace = True)
        l2 = len(cnvs.index)
        print(f'{l1-l2} cell lines without dependency data removed from CNV data')
        
# Purge EXPR data of cell lines without dependency data
        l1 = len(expr.index)
        c1 = set(expr.index)
        c2 = set(deps.index)
        c3 = c1 & (c1 ^ c2)
        expr.drop(list(c3),inplace = True)
        l2 = len(expr.index)
        print(f'{l1-l2} cell lines without dependency data removed from EXPR data')

# =============================================================================
# now add in CNV data 
# =============================================================================
        
        print('Adding Copy Number Variation data ...')
        
# =============================================================================
# set acceptance parameters for CNV : based on 
# deep del < log2(0.87/ploidy+1) < het loss < log2(1.32/ploidy+1) < diploid 
# < log2(2.64/ploidy+1) < gain < log2(3.36/ploidy+1) < amp
#                                                                                           
# =============================================================================
        ploidy = 2
        CNVval = np.zeros(4)
        CNVval[0] = log2(0.87/ploidy+1)
        CNVval[1] = log2(1.32/ploidy+1)
        CNVval[2] = log2(2.64/ploidy+1)
        CNVval[3] = log2(3.36/ploidy+1)
        
        CNVtype = ('Del','Het','Dip','Gain','Amp')
                
        CNVhist = np.zeros(5)


# add cell lines where genes are deleted and amplified
        
        ng = 1
        t_base = timer()

        for g in cnvs.columns: #loop over genes in CNV data
            Gdict[g].cnvamp = list(cnvs[cnvs[g] > CNVval[3]].index)
            Gdict[g].cnvdel = list(cnvs[cnvs[g] < CNVval[0]].index)
            if ng %1000 == 0 :
                print(f'{ng} genes time : {timer()-t_base:.4}s')
            ng += 1
            
        cnvtemp = cnvs.transpose()

        ng = 1
        t_base = timer()

# add deleted/amplified genes to cell lines

        for c in cnvtemp.columns[1:]:
            Cdict[c].cnvamp = list(cnvtemp[cnvtemp[c] > CNVval[3]].index)
            Cdict[c].cnvdel = list(cnvtemp[cnvtemp[c] < CNVval[0]].index)
            if ng %1000 == 0 :
                print(f'{ng} cells time : {timer()-t_base:.4}s')
            ng += 1
            
        del cnvtemp

    #     for g in cnvs.columns: #loop over genes in CNV data
    #         for c in cnvs.index:
    #             if c in Cdict : # CNV in valid cell line
    #                 cnv = cnvs.loc[c][g]
                    
    #                 if isFloat(cnv):  # ignore garbage
                    
    # # accumulate in cnv array 
                    
    #                     if cnv < CNVval[0]:
    #                         # print(f'{g} deleted in {c} - CNV : {cnv}')
    #                         Cdict[c].cnvdel.append(g)
    #                         Gdict[g].cnvdel.append(c)
    #                         CNVhist[0] +=1
    #                     elif cnv > CNVval[3]:
    #                         Cdict[c].cnvamp.append(g)
    #                         Gdict[g].cnvamp.append(c)
    #                         # print(f'{g} amplified in {c} - CNV : {cnv}')
    #                         CNVhist[4] +=1
    #                     elif cnv > CNVval[0] and cnv < CNVval[1]:
    #                         CNVhist[1] += 1
    #                     elif cnv > CNVval[1] and cnv < CNVval[2]:
    #                         CNVhist[2] += 1
    #                     elif cnv > CNVval[2] and cnv < CNVval[3]:
    #                         CNVhist[3] += 1
                
    #             if ng %1000 == 0 :
    #                 print(f'{ng} CNVs time : {timer()-t_base:.4}s')
    #             ng += 1
                    


# # optional histogram

#         p = input('Plot a histogram of Copy Number variations ? : [Y]')
#         if p == '' or p.lower()[0] == 'y':
            
            
            
#             plt.stairs(CNVhist)
#             plt.xticks(list(range(5)),labels = [' '*20+i for i in CNVtype ])
#             plt.title('Distribution of CNV categories')
#             plt.ylabel('counts')
#             plt.show()
            
#         print(f'{int(CNVhist[0])} CNVs consistent with loss of function were found')
                                                                                                  
# =============================================================================
# now add in CNV data 
# =============================================================================
        
        print('Adding Expression data ...')
        default_file = 'ExpressionThresholds.csv'
        
        f = input(f'Expression Threshold File [{default_file}] : ')
        if f.strip() == '' :
            f = default_file
        print(f'Loading {f}')
        
        exval = pd.read_csv(f)
        exval.set_index('HugoSymbol', inplace = True)
        
        print(exval)

        ng = 1
        t_base = timer()
        elist = list(exval.index)
        
        for g in expr.columns: #loop over genes in EXPR data
            if g in elist: # only genes with expression threshold data
                low_val = exval.loc[g]['Low']
                high_val  = exval.loc[g]['High']
                Gdict[g].over = list(expr[expr[g] > high_val].index)
                Gdict[g].under = list(expr[expr[g] < low_val].index)
                for c in Gdict[g].over:
                    Cdict[c].over.append(g)
                for c in Gdict[g].under:
                    Cdict[c].under.append(g)
            if ng %1000 == 0 :
                print(f'{ng} genes time : {timer()-t_base:.4}s')
            ng += 1
       


        # for g in expr.columns: #loop over genes in expr data
        #     if g in elist: # only genes with dependency and expression threshold data
        #         low_val = exval.loc[g]['Low']
        #         high_val  = exval.loc[g]['High']
        #         # gm = expr[g].mean()
        #         # gs = 2*expr[g].std()
        #         for c in expr.index:
        #             if c in Cdict : # expression in valid cell line
        #                 exp = expr.loc[c][g]
        #                 if isFloat(exp):  # ignore garbage
        #                     if high_val != None and exp > high_val : # overexpressed
        #                         Cdict[c].over.append(g)
        #                         Gdict[g].over.append(c)
        #                     elif low_val != None and exp < low_val : # underexpressed
        #                         Cdict[c].under.append(g)
        #                         Gdict[g].under.append(c)
    
                
        #         if ng %1000 == 0 :
        #             print(f'{ng} expressions time : {timer()-t_base:.4}s')
        #         ng += 1


# now save the compiled database

        f_gene = open('AB_GeneDatabaseExpr2023Q2','w')
        f_cell = open('AB_CellDatabaseExpr2023Q2','w')
        DB_save(f_gene,f_cell)
        print('AB_GeneDatabaseExpr2023Q2 written')
        print('AB_CellDatabaseExpr2023Q2 written')

# reset index struture of cnvs
    cnvs.reset_index(inplace = True)
# reset index struture of expr
    expr.reset_index(inplace = True)
    
    if spec:
        GetParams()    
 
# select task - i.e. QUIT,Find synthetic realtionships, pair analysis by mutation/cancer type

    while True:
        print('This program can :')
        print('find Genes B whose mutional status affects target Gene A dependency (FIND)')
        print('analyse Target gene dependency by cancer tissue type (TYPE)')
        print('Analyse genes in a defined locus (LOCU)')
        print('analyse pair list from an AllByAll run (ALL)')
        print('calculate p values for a list of gene pairs (LIST)')
        print('analyse mutations in GeneB in cells with GeneA dependency on GeneB (MUTS)')
        # print('analyse highly represented GenesA and GenesB (HUBS)')
        print('plot and analyse dependency of a Target on a defined Genetic Profile (MINE)')
        print('analyse chromosomal locations of GeneBs that switch GeneA dependency (CHRO)')
        print('define genetic profiles for cell line selection (DEFI)')
        print('manage cell line genetic profiles (PROF)')
        print('chromosomally map mutations/deletions in a cell line (CYTO)')
        print('find genes whose dependency is significantly switched by a defined Genetic Profile (SWITCH)')
        print('find genes whose dependency is significantly switched by know Gain of Function mutations (GOF)')
        print('change paramaters (PARA)')
        print('analyse gene expression data (EXPR)')
        print('analyse mutational signatures (SIG)')
        print('output cell lines mtching profile with specified dependent gene')
        print('(QUIT)\n')
        
        tasks = ['ALL',
         'CHRO',
         'CYTO',
         'DEFI',
         'DEPS',
         'EXPR',
         'FIND',
         'GOF',
         'LIST',
         'LOCU',
         'MINE',
         'MUTS',
         'PARA',
         # 'PLOT',
         'PROF',
         'QUIT',
         'SIG',
         'SWITCH',
         'TISS',
         'TYPE']
        
        task = input(f"select {'/'.join(tasks)} : ")
        if task.lower() == 'quit' :
            break
        if task.lower() == 'para' :
            GetParams()
            

# =============================================================================
# ██████╗ ██████╗ ██╗███╗   ███╗███████╗
# ██╔══██╗██╔══██╗██║████╗ ████║██╔════╝
# ██████╔╝██████╔╝██║██╔████╔██║█████╗  
# ██╔═══╝ ██╔══██╗██║██║╚██╔╝██║██╔══╝  
# ██║     ██║  ██║██║██║ ╚═╝ ██║███████╗
# ╚═╝     ╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚══════╝
# =============================================================================
                                      

# primary cancer type factors - may become free standing later

        if task.lower() == 'prime':
            dis = GetOncotreePrimaryDisease()
            
# loop through diseases            

            pool = mp.Pool(mp.cpu_count())
            results = pool.starmap_async(MakeGoFFing,
                        [(d,info,muts) for d in dis]).get()

                # results = pool.starmap_async(ProfPairPvalue,
                #         [(g,prof,p,q,deps,thresh) for g in deps.columns]).get()
            pool.close()

            print('GoF Fingerprints created')
            fing = {i[0].replace(',','').replace(' ','_'):i[1] for i in results}
            
            ans = input('Save GoF fingerprints to file ? ([Y],N) : ').upper()
            if len(ans) > 0 and ans[0] == 'N':
                continue
            fdat = pd.DataFrame(columns = ['GoF']+list(fing))
            fdat['GoF'] = list(results[0][1])
            for d in fdat.columns[1:]:
                fdat[d] = [fing[d][i] for i in fing[d]]
            print(fdat)   
            fdat.to_csv('Fingerprints_GoF.csv')
            print('Fingerprints_GoF.csv saved')
            # for d in dis:
                
# p is set of cells with this disease type
#                 p = set(info[info['OncotreePrimaryDisease'] == d]['ModelID'])
#                 q = UNI() ^ p
            
# # get all GoF mutations in p

#                 pGoFs = muts[muts['ModelID'].isin(p) & muts['LikelyGoF'] == True] \
#                             [['HugoSymbol','ProteinChange']]
#                 lg = list(pGoFs['HugoSymbol'])
#                 lc = list(pGoFs['ProteinChange'])
# # combine into mutation
#                 pmts = [lg[i]+lc[i] for i in range(len(lg))]

#                 qGoFs = muts[muts['ModelID'].isin(q) & muts['LikelyGoF'] == True] \
#                             [['HugoSymbol','ProteinChange']]
#                 lg = list(qGoFs['HugoSymbol'])
#                 lc = list(qGoFs['ProteinChange'])
# # combine into mutation
#                 qmts = [lg[i]+lc[i] for i in range(len(lg))]
                
#                 all_muts = sorted(set(pmts) | set(qmts)) 
                
#                 cont = [[0,0],[0,0]] # initialise contingency table
                
#                 lp = len(pmts)
#                 lq = len(qmts)
#                 r = {}
                
#                 for m in all_muts:
#                     Np = pmts.count(m)
#                     cont[0][0] = Np
#                     cont[1][0] = lp-Np
#                     Nq = qmts.count(m)
#                     cont[0][1] = Nq
#                     cont[1][0] = lq-Nq
#                     v = stats.boschloo_exact(cont,alternative="less").pvalue
#                     if np.isnan(v):
#                         r[m] = 0.0
#                     else:
#                         r[m] = -log10(v)
#                 r = dict(sorted(r.items(),key = lambda x: x[1], reverse = True))
#                 dpattern[d] = r
                
                
# MSI one-off function

        if task.lower() == 'msi' : 
# get temp slice of hgnc
# load HUGO chromosomal location data
            f = default_HUGO_file = 'hgnc_complete_set.tsv'
            if spec:
                f = input(f'HUGO Gene file [{default_HUGO_file}] : ')
                if f.strip() == '' :
                    f = default_HUGO_file
            print(f'Loading temp slice of {f}')
            u2h = pd.read_table(f, usecols=['symbol','ucsc_id'],low_memory=False).fillna('')
            print('Loading human microsatelite table')
            misat = pd.read_table('MSDBG000001_perf.tsv',low_memory=False).fillna('')
            misat = list(misat[misat['Location'] == 'Exon']['GeneName'])
# load WRN GeneB hit list
            wrns = pd.read_csv('WRN-DepHits.csv')

            msig = []
            for i in wrns.index:
                gb = wrns.iloc[i]['GeneB']
                l = u2h[u2h['symbol'] == gb]
                if len(l) > 0:
                    ugb = list(l)[0]
                else:
                    print(f'No uscs name for {gb} - skipping')
                    continue
                if ugb in misat:
                    msig.append(ugb)

        elif task.lower() == 'tiss' :
            if not info_loaded:
                default_cell_info_file = 'Model.csv'
                f = input(f'Cell sample information file [{default_cell_info_file}] : ')
                if f.strip() == '' :
                    f = default_cell_info_file
                print(f'Loading {f}')
            
                info = pd.read_csv(f)
                info_loaded = True
                info['OncotreeLineage'] = [x.upper() for x in info['OncotreeLineage']] # make cancer types upper case
                cancer_types = set(info['OncotreeLineage'])
            print('Known tissues : \n',sorted(list(cancer_types)))
            print('\n')
        elif task.lower() == 'plot' :
            seaborn_plot = not seaborn_plot
            
            
            
# =============================================================================
# ███████╗██╗ ██████╗ 
# ██╔════╝██║██╔════╝ 
# ███████╗██║██║  ███╗
# ╚════██║██║██║   ██║
# ███████║██║╚██████╔╝
# ╚══════╝╚═╝ ╚═════╝ 
# =============================================================================
        
        elif task.lower() == 'sig' : # analyse signature usage in a cell selection
             
            first_call = True
            while not prof_loaded:
                profiles = {}
                list_profile_files()
                default_prof_file = 'AB_Dep_Profiles'
                f = input(f'Default Profile file [{default_prof_file}] : ')
                if f.strip() == '' :
                    f = default_prof_file
                print(f'Loading {f}')
                if prof_load(f):
                    print(f'{f} loaded')
                    prof_loaded = True
                    break
                else:
                    continue
            if len(profiles) < 1:
                print('\nNo genetic profiles have been defined - do this in DEFI\n')
                continue

# get genetic profile to use

            while True:
                if first_call and len(profiles) > 1:
                    print('\nAvailable predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    first_call = False

                gp = input('\nDefined genetic profile to use (QUIT): ')
                if len(gp) == 0:
                    continue
                if gp.lower() == 'quit' :
                    break
                if gp[0] == '?':
                    print('\nAvailable predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    
# check that profile exists
                m = match_keyword(gp.upper(),profiles)
                if  m == '*':
                    print(f'{gp.upper()} is ambiguous')
                    continue
                if m not in profiles:
                    print(f'Genetic profile : {gp} is not defined\nDo this in LINE\n')
                    continue

                prof = profiles[m]
                
                p = list(eval(prof_parse(prof)))
                
                if len(p) > 0:
                    SigType(p,m)
                else:
                    print(f'profile {m} gives emopty set')

# =============================================================================
# ██╗  ██╗██╗ ██████╗ ██╗  ██╗
# ██║  ██║██║██╔════╝ ██║  ██║
# ███████║██║██║  ███╗███████║
# ██╔══██║██║██║   ██║██╔══██║
# ██║  ██║██║╚██████╔╝██║  ██║
# ╚═╝  ╚═╝╚═╝ ╚═════╝ ╚═╝  ╚═╝
# =============================================================================



        elif task.lower() == 'high' : # find genes that have high dependency in a defined cell selection
            
            first_call = True
            while not prof_loaded:
                profiles = {}
                list_profile_files()
                default_prof_file = 'AB_Dep_Profiles'
                f = input(f'Default Profile file [{default_prof_file}] : ')
                if f.strip() == '' :
                    f = default_prof_file
                print(f'Loading {f}')
                if prof_load(f):
                    print(f'{f} loaded')
                    prof_loaded = True
                    break
                else:
                    continue
            if len(profiles) < 1:
                print('\nNo genetic profiles have been defined - do this in DEFI\n')
                continue

# get genetic profile to use

            while True:
                if first_call and len(profiles) > 1:
                    print('\nAvailable predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    first_call = False

                gp = input('\nDefined genetic profile to use (QUIT): ')
                if len(gp) == 0:
                    continue
                if gp.lower() == 'quit' :
                    break
                if gp[0] == '?':
                    print('\nAvailable predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    
# check that profile exists
                m = match_keyword(gp.upper(),profiles)
                if  m == '*':
                    print(f'{gp.upper()} is ambiguous')
                    continue
                if m not in profiles:
                    print(f'Genetic profile : {gp} is not defined\nDo this in LINE\n')
                    continue

                prof = profiles[m]

# generate hypothesis set
                p = list(eval(prof_parse(prof)))
                
# generate null set
                q = list(eval(prof_nullhype(prof_flatten(prof))))
                
                print(f'{len(p)} cell lines in matching profile set')
                print(f'{len(q)} cell lines in null set')
                
                ans = input('Plot cell type distribution for this profile ? ([Y],N) : ')
                if len(ans) == 0 or ans[0].upper() == 'Y':
                    PlotCellType(p)

                print(f'Calculating dependency for all targets in {len(p)} cell lines with genetic profile {m}')


# create deps subset from p

                depsub = deps[deps.index.isin(list(p))]
                dmeans = {}
                
                for g in depsub.columns:
                    d = depsub[g].mean()
                    if d >thresh:
                        dmeans[g] = d
                
                               
                high_deps = dict(sorted(dmeans.items(),key = lambda x: x[1], reverse = True))

                print(f'Cells matching genetic profile {m} have {len(high_deps)} genes with |dep| > {thresh} ')
                
# make dataframe from these

                high_genes = pd.DataFrame(high_deps.items())
                high_genes.rename(columns = {high_genes.columns[0]:'HugoSymbol',high_genes.columns[1]:'Mean Dependency'},inplace = True)
                
# look for known oncogenes in this set

                if not oncogenes_loaded:
                    default_oncogene_file = 'oncogene_only.csv'
                    f = input(f'Oncogene file [{default_oncogene_file}] : ')
                    if f.strip() == '' :
                        f = default_oncogene_file
                    print(f'Loading {f}')
    
                    oncogenes = pd.read_csv(f)
                    oncogenes_loaded = True

# look for known GOFs in filtered list
                
                high_genes['Oncogene'] = high_genes['HugoSymbol'].isin(oncogenes['HugoSymbol'])  
                

                pgof = set([x[0] for x in sorted(GetGOF())])

                high_genes['GoF'] = high_genes['HugoSymbol'].isin(pgof)
            
                high_oncs = high_genes[high_genes['Oncogene']]
                high_gofs = high_genes[high_genes['GoF']]
                print(f'Of which {len(high_oncs)} are known oncogenes')
                print(f'Of which {len(high_gofs)} have known gain-of-function mutations')
                
# look for oncogenes with high expression in matching set

 
                
                print(high_oncs)
                high_oncs.to_csv(f'{m}-High.csv',index = False)
            
# =============================================================================
#  ██████╗  ██████╗ ███████╗
# ██╔════╝ ██╔═══██╗██╔════╝
# ██║  ███╗██║   ██║█████╗  
# ██║   ██║██║   ██║██╔══╝  
# ╚██████╔╝╚██████╔╝██║     
#  ╚═════╝  ╚═════╝ ╚═╝     
#  find genes switched by Gain of Function mutations                                     
# =============================================================================
        elif task.lower() == 'gof':
            
            a = input('Use built in GoF list or read external file ([B],E) ? :').upper()
            if a != '' and a[0] == 'E':
                gofs = pd.read_csv('Final_GOF_all_table.csv')
                gof = [(gofs.iloc[i]['Gene_name'],'p.'+gofs.iloc[i]['Substitution']) for i in gofs.index]
            
# get GOF list            
            else:
                gof = sorted(GetGOF())



            lgof = []
            for i in gof:
                if not pd.isna(i[1]):
                    tt = [i[0],i[1]]
                    lgof.append(tt)
            
            print(f'GoF list has {len(lgof)} individual mutations')
            
            agg_res = agg_gene = False
            a = input('Aggregate by genes, residues, or not aggregate (G,R,[N]) ? : ')
            if len(a) > 0 and a[0].upper() =='R':
                agg_res = True
            elif len(a) > 0 and a[0].upper() =='G':
                agg_gene = agg_res = True # perfor residue aggregation before gene aggregation
            
            if agg_gene:
                fhitname = 'GOF_AggGene.tsv'
            elif agg_res:
                fhitname = 'GOF_AggResi.tsv'
            else:
                fhitname = 'GOF_AggMuts.tsv'


            big_t_base = timer()    
            
            if agg_res:
                tgof = []
                for i in range(len(lgof)):
                    mul = False
                    p1 = lgof[i][0]
                    m1 = lgof[i][1][:-1]
                    for j in range(i+1,len(lgof)):
                        if lgof[j][0] == p1 and lgof[j][1][:-1] == m1:
                            lgof[i][1] = m1+'?'
                            lgof[j][1] = m1+'?'
                                           
            for i in range(len(lgof)):
                lgof[i] = (lgof[i][0],lgof[i][1].replace('*','?'))
            gof = sorted(list(set(lgof)))
            print(f'Residue aggregated GoF list has {len(gof)} residues')

# open output file

            fhit = open(fhitname,'w')
            fhit.write('GOF\tGeneA\tp\tD\n')
                
# loop through GOFs - for individual or residue aggregated

            if not agg_gene:
                for gx in gof:
    
    # make  profile from individual gof
    
                    print(f'generating profile for gene : {gx[0]} change : {gx[1]}')
                    prof = f'/{gx[0]}:{gx[1][2:]}'
                    print(prof)
    # generate hypothesis set
                    p = list(eval(prof_parse(prof)))
    # generate null set
                    q = list(eval(prof_nullhype(prof_flatten(prof))))
                    print(f'{len(p)} cell lines in matching profile set')
                    print(f'{len(q)} cell lines in null set')
    # check for null sets
    
                    if len(p) < 5: # should be minimum of 5 in each observation in contingency table
                        print(f'Genetic profile : {prof} returns < 5 observations')
                        continue
                    if len(q) == 0:
                        print(f'Null set for Genetic profile : {prof} is empty')
                        continue
                    
    # now parallelise over all genes in deps
    
                    results = []
                    t_base = timer()
    
                    pool = mp.Pool(mp.cpu_count())
    #                glist = list(set(deps.columns) | set(cnvs.columns[1:]))
                    glist = list(set(deps.columns))
    
                    qs = random.choices(q,k=int(min(len(q)/2,10*len(p)))) # set null set sample size
                    print(f'null set using sample size = {len(qs)}')
                    results = pool.starmap_async(AltProfPairPvalue,
                            [(g,prof,p,qs,deps,thresh,null_trim) for g in glist]).get()
    
                    # results = pool.starmap_async(ProfPairPvalue,
                    #         [(g,prof,p,q,deps,thresh) for g in deps.columns]).get()
                    pool.close()
    
                    print(f'Gene hits total time : {timer()-t_base:.5}')
                    
                    Plist = {} #initialise hitlist
                    HitList = {}
            
                    for r in results:
                        if r :
                            Plist[r[0]] = r[1]    
    # sort on p values
                    
                    Plist = dict(sorted(Plist.items(), key = lambda x: x[1], reverse = False))
                    
                    for i in Plist:
                        if Plist[i] < p_limit: # hypethesis is True
                            HitList[i] = Plist[i]
                    
                    if len(HitList) > 0 :
                        print(f'{len(HitList)} dependency switches found for genetic profile {prof} with p<={p_limit}\n')                     
                            
                for h in HitList:
                    p = eval(prof_parse(prof))
                    q = eval(prof_nullhype(prof_flatten(prof)))
                    match,notmatch = ProfPairDist(h,p,q,deps,null_trim)                   
    # calculate Cohen's D
                    D = (np.array(match).mean()-np.array(notmatch).mean())/np.array(match+notmatch).std()
                    if abs(D)> cohen_min:
                        fhit.write(f'{prof}\t{h}\t{HitList[h]:.5}\t{D}\n')
                                       
                fhit.close()
                print(f'Total GOF time : {timer()-big_t_base:.5}')

                print(f'{fhitname} written')
                continue    
                
# alternative for gene aggregated GOFs        
                
            pgof = sorted(list(set([x[0] for x in gof])))
            print(f'Gene aggregated GoF list has {len(pgof)} proteins')
            
            gof_prof = {}

# build /GENE:mut1,mut2,...mutn string

            for gg in pgof:
                prof = f'/{gg}:'
                
                for gx in gof:
                    if gx[0] == gg: # matches unique gene name
                        prof += f'{gx[1][2:]},'
                prof = prof[:-1]
                print(f'Gene aggregated profile {prof}')
                gof_prof[gg+'_GOF'] = prof
                
                p = eval(prof_parse(prof))
                print(f'{len(p)} cell lines in matching profile set')
                q = eval(prof_nullhype(prof_flatten(prof)))
                print(f'{len(q)} cell lines in null set')
                if len(p) < 5: # should be minimum of 5 in each observation in contingency table
                    print(f'Aggregated matching set for gene : {gg} returns < 5 observations - skipping\n')
                    continue
                if len(q) == 0:
                    print(f'Aggregated null set for gene : {gg} is empty')
                    continue
                
            # for gg in pgof:
            #     pp = set()
            #     qq = set()
                
            #     for gx in gof:
            #         if gx[0] == gg: # matches unique gene name
            #             print(f'generating profile for gene : {gx[0]} change : {gx[1]}')
            #             prof = f'/{gx[0]}:{gx[1][2:]}'
            #             print(prof)
        # generate hypothesis set
                        # p = eval(prof_parse(prof))
                        # print(f'{len(p)} cell lines in matching profile set')
                        # pp = pp | p
                        # print(f'{len(pp)} cell lines in aggregating matching profile set')
                        # q = eval(prof_nullhype(prof_flatten(prof)))
                        # print(f'{len(q)} cell lines in null set')
                        # qq = qq | q
                        # print(f'{len(qq)} cell lines in aggregating null set')
#                 if len(pp) < 5: # should be minimum of 5 in each observation in contingency table
#                     print(f'Aggregated matching set for gene : {gg} returns < 5 observations')
# # =============================================================================
# #                 if len(pp) == 0: # 
# #                     print(f'Aggregated matching set for gene : {gg} is empty')
# # =============================================================================

# #                     continue
#                 if len(qq) == 0:
#                     print(f'Aggregated null set for gene : {gg} is empty')
#                     continue

# now parallelise over all genes in deps

                results = []
                t_base = timer()

                pool = mp.Pool(mp.cpu_count())
#                glist = list(set(deps.columns) | set(cnvs.columns[1:]))
                glist = list(set(deps.columns))

                qs = random.choices(list(q),k = int(min(len(q)/2,10*len(p)))) # set null set sample size
                print(f'null set using sample size = {len(qs)}')
                results = pool.starmap_async(AltProfPairPD,
                        [(g,list(p),qs,deps,thresh) for g in glist]).get()
                
# =============================================================================
#                 results = pool.starmap_async(ProfPairPvalue,
#                         [(g,prof,p,q,deps,thresh) for g in glist]).get()
# 
# =============================================================================
                # results = pool.starmap_async(ProfPairPvalue,
                #         [(g,prof,p,q,deps,thresh) for g in deps.columns]).get()
                pool.close()

                print(f'Gene hits total time : {timer()-t_base:.5}')
                
                Plist = {} #initialise hitlist
                Dlist = {}
                HitList = {}
        
                for r in results:
                    if r :
                        Plist[r[0]] = r[1] # returned p value
                        Dlist[r[0]] = r[2] # returned D value
# sort on p values
                
                Plist = dict(sorted(Plist.items(), key = lambda x: x[1], reverse = False))
                
                for i in Plist:
                    if Plist[i] < p_limit and Dlist[i] > cohen_min: # hypothesis is True
                        HitList[i] = Plist[i]

                if len(HitList) > 0 :
                    print(f'{len(HitList)} dependency switches found for aggregated GOF for gene {gg} with p<={p_limit:.2} and D>{cohen_min:.2} \n')                     
                    for h in HitList:
                            fhit.write(f'{prof}\t{h}\t{HitList[h]:.5}\t{Dlist[h]}\n')
                
                                   
            fhit.close()
            print(f'Total GOF time : {timer()-big_t_base:.5}')

            print(f'{fhitname} written')
            ans = input('Save gene aggregated GoF profiles ? ([Y],N) : ').upper()
            if len(ans) > 0 and ans[0] != 'N':
                fgof = open('AB_GeneGOF.profiles','w')
                for p in gof_prof:
                    fgof.write(p.upper()+'\n')
                    fgof.write(gof_prof[p]+'\n')
                fgof.close()
                       
        elif task.lower()[:4] == 'pgof': # plot GOF switch results
        
        
            default_gof_file = 'GOF_AggGene.tsv'
            f = input(f'Default Profile file [{default_gof_file}] : ')
            if f.strip() == '' :
                f = default_gof_file
            default_gof_file = f
                
            gofs = pd.read_table(f)
            gofs.sort_values(by = 'p', inplace = True)
            gofs.reset_index(inplace = True)
            del gofs['index']
            
            print(gofs)
# loop through gofs calculating D
            for i in gofs.index:
                gB = gofs.iloc[i]['GOF']
                gA = gofs.iloc[i]['GeneA']
                pval = gofs.iloc[i]['p']
                D = gofs.iloc[i]['D']
                if abs(D) > cohen_min:
                
                    print(f'{gA:<8} : {gB:<8} p = {pval:.5f} : d = {D:.5f}')
                    t = input('Plot dependency distributions ([Y],N,QUIT) : ')
                    if t != '' and t[0].upper() == 'Q':
                        break
                    if t != '' and t[0].upper() == 'N':
                        continue
# trim gB name if gene aggregated
                    if 'ene' in default_gof_file:
                        t = gB[:gB.index(':')]+' gene aggregated'
                    else:
                        t = gB
                        
                    p = eval(prof_parse(gB))
                    q = eval(prof_nullhype(prof_flatten(gB)))
                    match,notmatch = ProfPairDist(gA,p,q,deps,null_trim)                   

                    plt.violinplot((match,notmatch),showmeans=True, showextrema=False)
                    plt.xticks([1,2],labels = ['matched','not matched'],fontsize = '20')
                    plt.title(f'{gA} : {t} p = {pval:.4g} d = {D:.4g}',fontsize = '20',pad=20)
                    plt.tick_params(axis='both', labelsize=12)
                    plt.ylabel('dependency',fontsize = '20')
                    plt.show()
                else:
                    print(f'ignoring {gA} : {gB} d = {D}')

                
                            
        elif task.lower()[:4] == 'swit':
        
# =============================================================================
# ███████╗██╗    ██╗██╗████████╗ ██████╗██╗  ██╗
# ██╔════╝██║    ██║██║╚══██╔══╝██╔════╝██║  ██║
# ███████╗██║ █╗ ██║██║   ██║   ██║     ███████║
# ╚════██║██║███╗██║██║   ██║   ██║     ██╔══██║
# ███████║╚███╔███╔╝██║   ██║   ╚██████╗██║  ██║
# ╚══════╝ ╚══╝╚══╝ ╚═╝   ╚═╝    ╚═════╝╚═╝  ╚═╝
# find target genes whose dependency is switched by a genetic profile                                                  
# =============================================================================
            first_call = True
            while not prof_loaded:
                profiles = {}
                list_profile_files()
                default_prof_file = 'AB_Dep_Profiles'
                f = input(f'Default Profile file [{default_prof_file}] : ')
                if f.strip() == '' :
                    f = default_prof_file
                print(f'Loading {f}')
                if prof_load(f):
                    print(f'{f} loaded')
                    prof_loaded = True
                    break
                else:
                    continue
            if len(profiles) < 1:
                print('\nNo genetic profiles have been defined - do this in DEFI\n')
                continue

# get genetic profile to use

            while True:
                if first_call and len(profiles) > 1:
                    print('\nAvailable predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    first_call = False

                gp = input('\nDefined genetic profile to use (QUIT): ')
                if len(gp) == 0:
                    continue
                if gp.lower() == 'quit' :
                    break
                if gp[0] == '?':
                    print('\nAvailable predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    
# check that profile exists
                m = match_keyword(gp.upper(),profiles)
                if  m == '*':
                    print(f'{gp.upper()} is ambiguous')
                    continue
                if m not in profiles:
                    print(f'Genetic profile : {gp} is not defined\nDo this in LINE\n')
                    continue

                prof = profiles[m]
                print(f'Calculating dependency for all targets in cells with genetic profile {m}')

# generate hypothesis set - cleaned up by having deps data
                p = list(eval(prof_parse(prof)) & set(deps.index))
# generate null set - cleaned up by having deps data
                q = list(eval(prof_nullhype(prof_flatten(prof))) & set(deps.index))
                print(f'{len(p)} cell lines in matching profile set')
                print(f'{len(q)} cell lines in null set')
# check for null sets

                if len(p) == 0:
                    print(f'Genetic profile : {gp} returns empty set')
                    continue
                if len(q) == 0:
                    print(f'Null set for Genetic profile : {gp} is empty')
                    continue

# trim p and q to 


# now parallelise over all genes in deps

                results = []
                t_base = timer()

                pool = mp.Pool(mp.cpu_count())
#                glist = list(set(deps.columns) | set(cnvs.columns[1:]))
                glist = list(set(deps.columns))

                
                results = pool.starmap_async(AltProfPairPvalue,
                        [(g,prof,p,q,deps,thresh,null_trim) for g in glist]).get()

                # results = pool.starmap_async(ProfPairPvalue,
                #         [(g,prof,p,q,deps,thresh) for g in deps.columns]).get()
                pool.close()

                print(f'Gene hits total time : {timer()-t_base:.5}')
                
                Plist = {} #initialise hitlist
                HitList = {}
        
                for r in results:
                    if r :
                        Plist[r[0]] = r[1]    
# sort on p values
                
                Plist = dict(sorted(Plist.items(), key = lambda x: x[1], reverse = False))
                
# =============================================================================
# # correct for multiple testing
# 
#                 pvals = list(Plist.values())
#                 pnames = list(Plist)
#                 bh = multi.multipletests(pvals,alpha = p_limit,method='fdr_tsbh',returnsorted=False)   
#                 for i in range(len(Plist)):
#                     if bh[0][i]: # hypethesis is True
#                         HitList[pnames[i]] = bh[1][i]
#                 
# =============================================================================
                for i in Plist:
                    if Plist[i] < p_limit: # hypethesis is True
                        HitList[i] = Plist[i]
# make boxplots for HitList targets        
                
                if len(HitList) > 0 :
                    print(f'{len(HitList)} dependency switches found for genetic profile {prof} with p<={p_limit}\n')
                
                
                
                fig, ax = plt.subplots()
                ax.tick_params(axis='both', labelsize=10)
                
                y = [-log10(HitList[i]) for i in HitList]
                maxy = max(y)*1.1
                
                plt.bar([i for i,x in enumerate(HitList)],y,alpha=0.3)
                plt.tick_params(bottom = False)
                plt.ylim(0.0,maxy)
                # plt.xticks(list(range(len([i for i in HitList]))),labels = [i for i in HitList], rotation = 'vertical',
                #            horizontalalignment = 'center',fontsize = 6)
                plt.ylabel('-log(p)',fontsize = '20')
                plt.title(f'Switching targets for {m} with p < {p_limit}',fontsize = '20',pad=20)

                plt.show()

                
                cohen = {}    
                for h in HitList:
                    match,notmatch = ProfPairDist(h,p,q,deps,null_trim)                   
    # calculate Cohen's D
                    cohen[h] = (np.array(match).mean()-np.array(notmatch).mean())/np.array(match+notmatch).std()

                for h in HitList:
                    match,notmatch = ProfPairDist(h,p,q,deps,null_trim)                   
                     
                    print(f'{h:<8} : p = {HitList[h]:.5f} : d = {cohen[h]:.5f}')
                    t = input('Plot dependency distributions ([Y],N,QUIT) : ')
                    if t != '' and t[0].upper() == 'Q':
                        break
                    if t != '' and t[0].upper() == 'N':
                        continue

                    plt.violinplot((match,notmatch),showmeans=True, showextrema=False)
                    plt.xticks([1,2],labels = ['matched','not matched'],fontsize = '20')
                    plt.title(f'{h} : {m} p = {HitList[h]:.4g} d = {cohen[h]:.4g}',fontsize = '20',pad=20)
                    plt.tick_params(axis='both', labelsize=12)
                    plt.ylabel('dependency',fontsize = '20')
                    plt.show()
    
                ans = input('Do you want to save the Hit List : ? ')
                if len(ans) > 0 and ans[0].lower() == 'y':
                    fhit = open(f'{m}-DepProfHits.csv','w')
                    fhit.write('TARGET,p,D\n')
                    for h in HitList:
                        fhit.write(f'{h},{HitList[h]:.5g},{cohen[h]:.5g}\n')
                    fhit.close()
                    
                    print(f'{m}-DepProfHits.csv written')

# ProfPairPvalue(target[0],prof,p,q,deps,thresh)


        elif task.lower() == 'prof':
# =============================================================================
#  ██████╗ ██████╗  ██████╗ ███████╗
#  ██╔══██╗██╔══██╗██╔═══██╗██╔════╝
#  ██████╔╝██████╔╝██║   ██║█████╗  
#  ██╔═══╝ ██╔══██╗██║   ██║██╔══╝  
#  ██║     ██║  ██║╚██████╔╝██║     
#  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     
# manage cell line profiles
# =============================================================================
            first_call = True
            while not prof_loaded:
                profiles = {}
                list_profile_files()
                default_prof_file = 'AB_Dep_Profiles'
                f = input(f'Select profile file [{default_prof_file}] : ')
                if f.strip() == '' :
                    f = default_prof_file
                print(f'Loading {f}')
                if prof_load(f):
                    print(f'{f} loaded')
                    prof_loaded = True
                    break
                else:
                    continue
            if len(profiles) < 1:
                print('\nNo genetic profiles have been defined - do this in DEFI\n')
                continue
            while True:
                if first_call and len(profiles) > 1:
                    print('\nAvailable predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    first_call = False
                    
                a = input('Show, Delete, Name, Flatten, Read, Write, Quit, ? (S,D,N,F,R,W,[Q]) : ').upper()
                if len(a) == 0 or a[0] == 'Q':
                    break
                
                elif a[0] == '?': # list profiles
                    print('\nAvailable predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
    
                elif a[0] == 'R': # read profile file
                    add = input('Add new profiles, Replace loaded profiles, Quit ([A],R,Q) : ').upper()
                    if len(add) > 0 and add[0] == 'Q':
                        continue
                    if len(add) > 0 and add[0] == 'R':
                        profiles = {}
                        prof_loaded = False
                    list_profile_files()
                    default_prof_file = 'AB_Dep_Profiles'
                    f = input(f'Select profile file [{default_prof_file}] : ')
                    if f.strip() == '' :
                        f = default_prof_file
                    print(f'Loading {f}')
                    if prof_load(f):
                        print(f'{len(profiles)} loaded from {f}')
                        prof_loaded = True
                    else:
                        continue
        
                
                elif a[0] == 'W': # write profile file
                    def_prof_file = 'AB_Dep_Profiles'
                    add = input('Create new profile file, Overwrite default, Quit ([C],O,Q) : ').upper()
                    if len(add) > 0 and add[0] == 'Q':
                        continue
                    if len(add) > 0 and add[0] == 'O':
                        prof_save(def_prof_file)
                        print(f'{len(profiles)} profiles written to default file {def_prof_file}')
                        continue
                    while True:
                        b = input('New profile file name ? (xxxx.profiles) : ' ).strip()
                        if '.profiles' in b:
                            b =  b.split('.')[0]   
                        if len(b) > 0:
                            pname = b+'.profiles'
                            if pname in glob('*rofiles'): # danger of overwriting
                                b = input(f'{pname} exists. Overwrite or Quit ? (O,[Q]) : ').upper()
                                if len(b) > 0 and b[0] != 'O':
                                    break
                            prof_save(pname)    
                            print(f'{len(profiles)} profiles written to {pname}')
                            break

                elif a[0] == 'F': # FLATTEN or Flatten <name>
                    w = a.split()
                    if len(w) == 1:
                        b = input('Flatten which profile ? : ')
                        if len(b) <1 :
                            continue
                        else:
                            w.append(b.upper())
                    
                    m = match_keyword(w[1],profiles)
                    if  m == '*':
                        print(f'{w[1]} is ambiguous')
                    elif m in profiles:
                        print(f'\n{m} : {prof_flatten(profiles[m])}')
                    else:
                        print(f'{w[1]} is not a known profile')
    
                elif a[0] == 'S': # SHOW or SHOW <name>
                    w = a.split()
                    if len(w) == 1:
                        b = input('Show which profile ? : ')
                        if len(b) <1 :
                            continue
                        else:
                            w.append(b.upper())
                    m = match_keyword(w[1],profiles)
                    if  m == '*':
                        print(f'{w[1]} is ambiguous')
                    elif m in profiles:
                        print(f'\n{m} : {profiles[m]}')
                    else:
                        print(f'{w[1]} is not a known profile')
                        
                elif a[0] == 'D': # DELETE or DELETE <name>
                    w = a.split()
                    if len(w) == 1:
                        b = input('Delete which profile ? : ')
                        if len(b) <1 :
                            continue
                        else:
                            w.append(b.upper())
                    m = match_keyword(w[1],profiles)
                    if m.upper() == 'ALL':
                        print('ALL is a protected profile and cannot be deleted')
                        continue
                    if  m == '*':
                        print(f'{w[1]} is ambiguous')
                    elif m in profiles:
# build list of profiles that depend on specified profile
                        p=[]                    
                        prof_tree(m,p)
                        if len(p) > 0:
                            print(f'\nDeleting {m} will also remove {Flistand(p)}')
                        b = input(f'Are you sure you want to delete {m} ? (Y/[N] : ')
                        if len(b) > 0 and b[0].lower() == 'y':
                            del profiles[m]
                            for i in p:
                                del profiles[i]
                                
                            add = input('Create new profile file, Overwrite default, Quit ([C],O,Q) : ').upper()
                            if len(add) > 0 and add[0] == 'Q':
                                continue
                            if len(add) > 0 and add[0] == 'O':
                                prof_save('AB_Dep_Profiles')
                                print(f'{len(profiles)} profiles written to default file')
                                continue
                            while True:
                                b = input('New profile file name ? (xxxx.profiles) : ' ).strip()
                                if '.profiles' in b:
                                    b =  b.split('.')[0]   
                                if len(b) > 0:
                                    pname = b+'.profiles'
                                    if pname in glob('*rofiles'): # danger of overwriting
                                        b = input(f'{pname} exists. Overwrite or Quit ? (O,[Q]) : ').upper()
                                        if len(b) > 0 and b[0] != 'O':
                                            break
                                    prof_save(pname)    
                                    print(f'{len(profiles)} profiles written to {pname}')
                                    break
                                                                
                    else:
                        print(f'{w[1]} is not a known profile')
                        
                elif a[0] == 'N': # NAME or NAME <src> or NAME <src> <dest>
                    w = a.split()
                    if len(w) == 1:
                        b = input('Name which profile ? : ')
                        if len(b) < 1 :
                            continue
                        elif len(b.split()) == 1 :
                            w.append(b.upper())
                            c = input(f'Name {b} as what ? : ')
                            if len(c) < 1 :
                                continue
                            else:
                                w.append(c.upper())
                    elif len(w) == 2:
                        c = input(f'Name {w[1]} as what ? : ')
                        if len(c) < 1 :
                            continue
                        else:
                            w.append(c.upper())
                    m = match_keyword(w[1],profiles)
                    if  m == '*':
                        print(f'{w[1]} is ambiguous')
                        continue
                    if m not in profiles:
                        print(f'{w[1]} is not a known profile')
                        continue
                    if w[2] in profiles:
                        print(f'{w[2]} is already a profile name')
                        continue
                    for i in profiles:
                        if m in profiles[i].upper():
                            print(f'Updating reference to {m} in {i} to {w[2]}')
                            profiles[i] = profiles[i].replace(m.lower(),w[2].lower())
                    profiles[w[2]] = profiles[m]
                    del profiles[m]
        
        elif task.lower() == 'defi':
# =============================================================================
# ██████╗ ███████╗███████╗██╗
# ██╔══██╗██╔════╝██╔════╝██║
# ██║  ██║█████╗  █████╗  ██║
# ██║  ██║██╔══╝  ██╔══╝  ██║
# ██████╔╝███████╗██║     ██║
# ╚═════╝ ╚══════╝╚═╝     ╚═╝
# find cell lines with specific mutation/deletion profiles                                         
# =============================================================================
            
            first_call = True
            while not prof_loaded:
                profiles = {}
                list_profile_files()
                default_prof_file = 'AB_Dep_Profiles'
                f = input(f'Profile file [{default_prof_file}] : ')
                if f.strip() == '' :
                    f = default_prof_file
                print(f'Loading {f}')
                if prof_load(f):
                    print(f'{f} loaded')
                    prof_loaded = True
                    break
                else:
                    continue
            if len(profiles) < 1:
                print('\nNo genetic profiles have been defined - do this in DEFI\n')
                continue

            while True:
                if first_call and len(profiles) > 0:
                    print('Available predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    first_call = False
               
                s = input('\nGenetic profile (QUIT) : ').lower()
                
                if s == 'quit':
                    break
                elif len(s) < 1:
                    continue
                elif s[0] == '?':
                    print('Available predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                    continue

# is it an 'in-line' definition i.e. MYC_UP = +MYC | >MYC
                ss = s.split()
                if len(ss) > 2 and ss[1] == '=':
                    u_name = s.split()[0].upper()
                    s = ' '.join(s.split()[2:])
                else:
                    u_name = ''

# parse SET based genetic profile definition


                u = prof_parse(s)
                if u == '':
                    # print(f'profile definition {s} parses to null')
                    continue
# loop through cell lines acquiring those that match condition string
        
                try:
                    r = list(eval(u))
                except Exception as err:
                    print(f'{err} in profile definition : \n{u}')
                    continue
                    
                hits = []    
                if len(r) > 0:
                    hits = r        
                    print(f'Found {len(hits)} cell lines matching : {s}')
                else:
                    print('profile definition returns empty set')
                    continue
                if len(hits) > 0:
                    a = input('Save list ? (Y/[N]) : ')
                    if len(a) > 0 and a.lower()[0] != 'n':
                        default_cell_save_file = 'AB_Dep_CellLine_Subset'
                        f = input(f'List file name [{default_cell_save_file}] : ')
                        if f.strip() == '' :
                            f = default_cell_save_file
                        fl = open(f,'w')
                        fl.write(','.join(hits)+'\n')
                        fl.close()
                        for h in hits:
                            f = info.loc[info['ModelID'] == h]
                            if len(f) > 0:
                                cname = f['CCLEName'].iloc[0]    
                                print(f'{h}  : {cname}')
                            else:
                                print(f'{h}  : no CCLE name found')
                    a = input(f'\nSave profile {u_name} (Y/[N]) : ')
                    if len(a) > 0 and a.lower()[0] != 'n':
                        while True:
                            if u_name == '':
                                u_name = input('profile name : (QUIT) ').upper()
                            if len(u_name) < 1:
                                continue
                            if u_name.upper() == 'QUIT':
                                break
# check name contains no reserved characters

                            if len(set(u_name) & set('~*;+-=!/$#@%><{&^|')) > 0:
                                print(f'profile name {u_name} contains reserved character(s) - choose a new name for this profile')                                                       
                                break
# check if name already exists in profiles
                            if u_name in profiles:
                                a = input(f'{u_name} is already defined, overwrite with new definition ? (Y/[N]) : ')
                                if len(a) == 0 or a[0].lower() != 'y':
                                    break
# check for self-reference
                                if u_name.lower() in s.lower():
                                    print(f'Redefining {u_name} creates a self-reference - choose a new name for this profile')
                                    continue
                                                                                            
# add to profiles
                            profiles[u_name] = s
                            prof_save('TEMP_Profiles')
                            print('working profiles saved in TEMP_Profiles')
                            break                    
        
        elif task.lower() == 'cyto':
                
# =============================================================================
#  ██████╗██╗   ██╗████████╗ ██████╗ 
# ██╔════╝╚██╗ ██╔╝╚══██╔══╝██╔═══██╗
# ██║      ╚████╔╝    ██║   ██║   ██║
# ██║       ╚██╔╝     ██║   ██║   ██║
# ╚██████╗   ██║      ██║   ╚██████╔╝
#  ╚═════╝   ╚═╝      ╚═╝    ╚═════╝ 
# show chromosomal locations of amplifications and deletions in a specified cell line                        
# =============================================================================
            
# load chromosome length data
                
            if not chromlen_loaded:    
                default_chromolen_file = 'ChromosomeLengths.csv'
                f = input(f'Chromosome length file [{default_chromolen_file}] : ')
                if f.strip() == '' :
                    f = default_chromolen_file
                print(f'Loading {f}')
                chrom = pd.read_csv(f)
                chromlen_loaded = True

# load Ensemble refseq gene-to-transcript mapping file
    
            if not ensem_loaded:
                default_ENSEMBLE_file = 'Homo_sapiens.GRCh38.110.refseq.tsv'
                f = input(f'ENSEMBLE reference sequence file [{default_ENSEMBLE_file}] : ')
                if f.strip() == '' :
                    f = default_ENSEMBLE_file
                print(f'Loading {f}')
                ensem = pd.read_table(f)
                ensem.drop_duplicates(subset = ['gene_stable_id'],inplace = True)
                ensem = ensem[['gene_stable_id','transcript_stable_id']]
                ensem.set_index('gene_stable_id',inplace = True)
                ensem_loaded = True
        
            while True:
                c = input('Cell line name ? (QUIT) : ').upper()
                if len(c) == 0:
                    continue
                if len(c) > 0 and c[0].lower() == 'q':
                    break
# try and match as ModelID
                if c not in list(info['ModelID']): # not model ID so try CellLineName
                    cn = info[info['StrippedCellLineName'].str.contains(c)]['ModelID']        
                    if len(cn) < 1:
                        print(f'Cell line name {c} not recognised')
                        continue
                    if len(cn) > 1:
                        print(f'Cell line name {c} is ambiguous')
                        print('Choose from :')
                        for i in cn:
                            clong = list(info[info['ModelID'] == i]['StrippedCellLineName'])[0]
                            print(f'{i} - {clong}')
                        continue
                    c = cn.iloc[0]
                    clong = list(info[info['ModelID'] == c]['StrippedCellLineName'])[0]
# draw chromosomes

# define Matplotlib figure and axis
                    fig, ax = plt.subplots(figsize = (25,12))
                        # plt.figure(figsize = (25,12))
                    chromwidth = 0.01
                    cyscale = 1/ (chrom.loc[0]['length']*1.1)
        #                print(cyscale)
                    cwidth = 0.005
                    centsize = 0.5
                    xstart = 0.04
                    ystart = 0.07
                    # plt.rcParams.update({'font.size': 30})
                    plt.tick_params(left = False, right = False , labelleft = False ,
                                        labelbottom = False, bottom = False)

# draw chromosomes
                
                    for i in range(len(chrom)):
                        xpos = xstart+i*8*cwidth
                        pstart = 0
                        plen = chrom.iloc[i]['centromere']
                        qstart = chrom.iloc[i]['centromere']
                        qlen = chrom.iloc[i]['length']-chrom.iloc[i]['centromere']
                        
                        ax.plot(xpos, 1-(ystart+pstart*cyscale),xpos,1-(ystart+pstart*cyscale)-plen*cyscale,lw=50,solid_capstyle='round', color='k')
                        ax.add_patch(Rectangle((xpos, 1-(ystart+pstart*cyscale)), 0.02, -plen*cyscale,fc = '0.8', ec='w'))
                        ax.add_patch(Rectangle((xpos, 1-(ystart+qstart*cyscale)), 0.02, -qlen*cyscale,fc = '0.7', ec='w'))
                    
                        ax.add_patch(Ellipse((xpos+(centsize/50)+0.001,1-(ystart+qstart*cyscale)), width=centsize/20, height=centsize/25,fc = 'w', ec='k'))
                        ax.text(xpos+0.01,0.94,chrom.iloc[i]['chromosome'],horizontalalignment='center',fontsize='30')
                
                    plt.rcParams.update({'font.size': 40})
                    ax.text(xstart-0.025,0.75,'p',horizontalalignment='center',style='italic',fontsize='30')
                    ax.text(xstart-0.025,0.25,'q',horizontalalignment='center',style='italic',fontsize='30')
                
# loop through amplified and deleted genes

                    for t in Cdict[c].cnvamp:
                        d = ChroLoc(t,hgnc,ensem)
                        if d == None:
                            print(f'Failed to get location data for gene : {t}')
                            continue
                        
        #                    print(f'{t} - {d[0]},{d[1]},{d[2]}')
# add lines for genes                                
                        icro = chrom['chromosome'].tolist().index(d[0])
                        xpos = xstart+icro*8*cwidth
                        ypos = 1-(ystart+d[1]*cyscale)
                        ax.plot([xpos,xpos+0.02],[ypos,ypos],color='blue')
                        
                    for t in Cdict[c].cnvdel:
                        d = ChroLoc(t,hgnc,ensem)
                        if d == None:
                            print(f'Failed to get location data for gene : {t}')
                            continue
                        
        #                    print(f'{t} - {d[0]},{d[1]},{d[2]}')
# add lines for genes                                
                        icro = chrom['chromosome'].tolist().index(d[0])
                        xpos = xstart+icro*8*cwidth
                        ypos = 1-(ystart+d[1]*cyscale)
                        ax.plot([xpos,xpos+0.02],[ypos,ypos],color='red')
                    
                        # ax.text(xpos-0.02,ypos,t,fontsize='10')
                        # plt.plot([xpos,xpos+0.02],[ypos,ypos],color='red')
                    # avloc = avloc/nloc
                    # ypos = 1-(ystart+avloc*cyscale)
                    # if 'p' in s:
                    #     ss = s.split('p')[1]
                    # else:
                    #     ss = s.split('q')[1]
                    # ax.text(xpos-0.02,ypos,ss,fontsize='12')
                    # ax.text(xstart+21*8*cwidth+0.01,0.1,s,horizontalalignment='center',color='blue')
                    ax.text(xstart+21*8*cwidth+0.01,0.2,'amplification',horizontalalignment='center',color='blue')
                    ax.text(xstart+21*8*cwidth+0.01,0.1,'deletion',horizontalalignment='center',color='red')
                    ax.text(xstart+8*8*cwidth+0.01,0.1,f'{clong} - {c}',horizontalalignment='center',color='black')
                    plt.show()


        elif task.lower() == 'locu' :
            
# =============================================================================
# ██╗      ██████╗  ██████╗██╗   ██╗
# ██║     ██╔═══██╗██╔════╝██║   ██║
# ██║     ██║   ██║██║     ██║   ██║
# ██║     ██║   ██║██║     ██║   ██║
# ███████╗╚██████╔╝╚██████╗╚██████╔╝
# ╚══════╝ ╚═════╝  ╚═════╝ ╚═════╝ 
# scan GeneA v. all genes in a cytogenetic locus looking for LoF SSL
# =============================================================================
                                              
# load chromosome length data
            
            if not chromlen_loaded:    
                default_chromolen_file = 'ChromosomeLengths.csv'
                f = input(f'Chromosome length file [{default_chromolen_file}] : ')
                if f.strip() == '' :
                    f = default_chromolen_file
                print(f'Loading {f}')
                chrom = pd.read_csv(f)
                chromlen_loaded = True

# load Ensemble refseq gene-to-transcript mapping file

            if not ensem_loaded:
                default_ENSEMBLE_file = 'Homo_sapiens.GRCh38.110.refseq.tsv'
                f = input(f'ENSEMBLE reference sequence file [{default_ENSEMBLE_file}] : ')
                if f.strip() == '' :
                    f = default_ENSEMBLE_file
                print(f'Loading {f}')
                ensem = pd.read_table(f)
                ensem.drop_duplicates(subset = ['gene_stable_id'],inplace = True)
                ensem = ensem[['gene_stable_id','transcript_stable_id']]
                ensem.set_index('gene_stable_id',inplace = True)
                ensem_loaded = True
                
# get locus
            while True:
                
                s = input('GeneA and chromosomal locus (QUIT) : ').lower()
                if s == 'quit':
                    break
                elif len(s) < 1:
                    continue
                
# check number of parameters given        
        
                w = s.replace(',',' ').split()
                if len(w) != 2:
                    print(f'Wrong number of parameters - 2 needed {len(w)} given')
                    continue

#  check that first parameter is valid gene name

                if w[0].upper() in hgnc.index:
                    geneA = w[0].upper()
                else:
                    print(f'{w[0]} is not a valid gene name')
                    continue
                
                

                if '.' in w[1] and w[1] in set(hgnc['location']): # is it a sub-band
                    gloc = list(hgnc[hgnc['location'] == w[1]].index)    
                    
                else: # major band
                    l1 = pd.DataFrame(list(hgnc.index),[x.split('.')[0] for x in hgnc['location']]).reset_index()
                    l1.rename(columns = {'index':'location',l1.columns[1]:'symbol'},inplace = True)
                    l1.set_index('symbol',inplace =True)
                    
                    if w[1] in set(l1['location']):
                        gloc = list(l1[l1['location'] == w[1]].index)
                    else:
                        print(f'{w[1]} is not a recognised genetic locus')
                        continue
                    
# print gene list
        
                print(f'Locus {w[1]} contains genes : {gloc}')

# now loop through gene list calculating p and d but only for LoF mutations

                ngloc = len(gloc)
                dvals = {}
                pvals = {}
                ps = {}
                qs = {}
                for g in gloc:
                    p = set()
                    pf1 = prof_parse(f'!{g}')
                    if pf1 != '':
                        p = eval(pf1)
                    pf2 = prof_parse(f'/{g}-')
                    if pf2 != '':
                        p |= eval(pf2)
                        
                    if p == set():
                        continue

                    p = list(p)
                    
                    q = list(eval(prof_parse(f'${g}')))
                        
                    r = AltProfPairPD(geneA,p,q,deps,thresh) # calc p and d
                    pvals[g] = r[1]
                    dvals[g] = r[2]
                    ps[g] = p
                    qs[g] = q
# sort pvals
                pvals = dict(sorted(pvals.items(), key = lambda x: x[1]))
                temp = dvals.copy()
                dvals = {}
                for i in temp:
                    if not np.isnan(temp[i]):
                        dvals[i] = temp[i]
                        
                dvals = dict(sorted(dvals.items(), key = lambda x: x[1],reverse=True))
                print(f'{len(pvals)} genes in {w[1]} have LoF')
                
                
# plot waterfall plot of effect size
                
                fig, ax = plt.subplots(figsize=[len(dvals)/4,4])
                ax.tick_params(axis='both', labelsize=10)
                
                y = [dvals[i] for i in dvals]
                # maxy = max(y)*1.1
                bars = plt.bar([i for i in dvals],y)
                
                for n,i in enumerate(dvals):
                    if pvals[i] < p_limit:
                        bars[n].set_alpha(0.5)
                    else:
                        bars[n].set_alpha(0.2)

                ax.set_xticklabels([i for i in dvals],rotation=90)

                plt.tick_params(bottom = False)
                # plt.ylim(0.0,maxy)
                # plt.xticks(list(range(len([i for i in HitList]))),labels = [i for i in HitList], rotation = 'vertical',
                #            horizontalalignment = 'center',fontsize = 6)
                plt.ylabel('d',fontsize = '20')
                plt.title(f'Mutation only {w[1]} switch effects for {w[0].upper()}',fontsize = '20',pad=20)

                ax.plot([-0.5,len(dvals)],[0.8,0.8],color='black',linewidth='0.5',linestyle='dashed')

                plt.show()

                

                for g in pvals:
                    if pvals[g] < p_limit:
                        ans = input(f'{geneA}:{g} - p = {pvals[g]:.4g} - d = {dvals[g]:.4g} - Plot ([Y],N,QUIT) : ').upper()
                        if len(ans) > 0 and ans[0] == 'Q':
                            break
                        elif len(ans) > 0 and ans[0] == 'N':
                            continue
                        fig, ax = plt.subplots()
                        ax.tick_params(axis='both', labelsize=10)
                        plt.ylim(0.0,1.0)
                        match,notmatch = ProfPairDist(geneA,ps[g],qs[g],deps,null_trim)                   
                         
                        plt.violinplot((match,notmatch),showmeans=True, showextrema=False)
                        plt.xticks([1,2],labels = [f'matched ({len(match)})','not matched'],fontsize = '20')
                        plt.title(f'{geneA} : {g} p = {pvals[g]:.4g} d = {dvals[g]:.4g}',fontsize = '20',pad=20)
                        plt.tick_params(axis='both', labelsize=12)
                        plt.ylabel('dependency',fontsize = '20')
                        plt.show()

                
# =============================================================================
#                 t = input('Plot genes on chromosome ([Y]/N)' ).upper()
#                 if len(t) == 0 or t[0] != 'N':
#             
#                 
# # define Matplotlib figure and axis
#                     fig, ax = plt.subplots(figsize = (25,12))
#                         # plt.figure(figsize = (25,12))
#                     chromwidth = 0.01
#                     cyscale = 1/ (chrom.loc[0]['length']*1.1)
#         #                print(cyscale)
#                     cwidth = 0.005
#                     centsize = 0.5
#                     xstart = 0.04uscs

#                     ystart = 0.07
#                     # plt.rcParams.update({'font.size': 30})
#                     plt.tick_params(left = False, right = False , labelleft = False ,
#                                         labelbottom = False, bottom = False)
# 
# # draw chromosomes
#                 
#                     for i in range(len(chrom)):
#                         xpos = xstart+i*8*cwidth
#                         pstart = 0
#                         plen = chrom.iloc[i]['centromere']
#                         qstart = chrom.iloc[i]['centromere']
#                         qlen = chrom.iloc[i]['length']-chrom.iloc[i]['centromere']
#                         
#                         ax.plot(xpos, 1-(ystart+pstart*cyscale),xpos,1-(ystart+pstart*cyscale)-plen*cyscale,lw=50,solid_capstyle='round', color='k')
#                         ax.add_patch(Rectangle((xpos, 1-(ystart+pstart*cyscale)), 0.02, -plen*cyscale,fc = '0.8', ec='w'))
#                         ax.add_patch(Rectangle((xpos, 1-(ystart+qstart*cyscale)), 0.02, -qlen*cyscale,fc = '0.7', ec='w'))
#                     
#                         ax.add_patch(Ellipse((xpos+(centsize/50)+0.001,1-(ystart+qstart*cyscale)), width=centsize/20, height=centsize/25,fc = 'w', ec='k'))
#                         ax.text(xpos+0.01,0.94,chrom.iloc[i]['chromosome'],horizontalalignment='center',fontsize='30')
#                 
#                     plt.rcParams.update({'font.size': 40})
#                     ax.text(xstart-0.025,0.75,'p',horizontalalignment='center',style='italic',fontsize='30')
#                     ax.text(xstart-0.025,0.25,'q',horizontalalignment='center',style='italic',fontsize='30')
#                    
#                     loc_list={}
#                     avloc = 0
#                     nloc = 0
#                         
#                     for t in g:
#                         d = ChroLoc(t,hgnc,ensem)
#                         if d == None:
#                             print(f'Failed to get location data for gene : {t}')
#                             continue
#                         loc_list[t] = d[2]
#                         avloc += d[1]
#                         nloc += 1
#                         
#         #                    print(f'{t} - {d[0]},{d[1]},{d[2]}')
# # add lines for genes                                
#                         icro = chrom['chromosome'].tolist().index(d[0])
#                         xpos = xstart+icro*8*cwidth
#                         ypos = 1-(ystart+d[1]*cyscale)
#                         ax.plot([xpos,xpos+0.02],[ypos,ypos],color='red')
#                         # ax.text(xpos-0.02,ypos,t,fontsize='10')
#                         # plt.plot([xpos,xpos+0.02],[ypos,ypos],color='red')
#                     avloc = avloc/nloc
#                     ypos = 1-(ystart+avloc*cyscale)
#                     if 'p' in s:
#                         ss = s.split('p')[1]
#                     else:
#                         ss = s.split('q')[1]
#                     ax.text(xpos-0.02,ypos,ss,fontsize='12')
#                     ax.text(xstart+21*8*cwidth+0.01,0.1,s,horizontalalignment='center',color='green')
# 
#                     plt.show()
# 
#             
# =============================================================================
        elif task.lower() == 'chro' :
            
# =============================================================================
#  ██████╗██╗  ██╗██████╗  ██████╗ 
# ██╔════╝██║  ██║██╔══██╗██╔═══██╗
# ██║     ███████║██████╔╝██║   ██║
# ██║     ██╔══██║██╔══██╗██║   ██║
# ╚██████╗██║  ██║██║  ██║╚██████╔╝
#  ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ 
# report chromosomal locations of Gene(s)                                          
# =============================================================================
        

# load chromosome length data
    
            if not chromlen_loaded:    
                default_chromolen_file = 'ChromosomeLengths.csv'
                f = input(f'Chromosome length file [{default_chromolen_file}] : ')
                if f.strip() == '' :
                    f = default_chromolen_file
                print(f'Loading {f}')
                chrom = pd.read_csv(f)
                chromlen_loaded = True
    
# load HUGO chromosomal location data
            
            default_HUGO_file = 'hgnc_complete_set.tsv'
            f = input(f'HUGO Gene file [{default_HUGO_file}] : ')
            if f.strip() == '' :
                f = default_HUGO_file
            print(f'Loading {f}')
            hgnc = pd.read_table(f)
            hgnc = hgnc[['symbol','ensembl_gene_id','prev_symbol','location','location_sortable']]
            hgnc.set_index('symbol',inplace = True)
            
# load Ensemble refseq gene-to-transcript mapping file
    
            if not ensem_loaded:
                default_ENSEMBLE_file = 'Homo_sapiens.GRCh38.110.refseq.tsv'
                f = input(f'ENSEMBLE reference sequence file [{default_ENSEMBLE_file}] : ')
                if f.strip() == '' :
                    f = default_ENSEMBLE_file
                print(f'Loading {f}')
                ensem = pd.read_table(f)
                ensem.drop_duplicates(subset = ['gene_stable_id'],inplace = True)
                ensem = ensem[['gene_stable_id','transcript_stable_id']]
                ensem.set_index('gene_stable_id',inplace = True)
                ensem_loaded = True
            
            default_list_file = 'DepProfHits.csv'
            
            while True:
                
                f = input(f'Gene list file [{default_list_file}] : ')
                if f.strip() == '' :
                    f = default_list_file
                print(f'Loading {f}')
                
                try:
                    genes = pd.read_csv(f)
                    print(f'{len(genes)} genes loaded')
                    break
                except FileNotFoundError:
                    print(f'file {f} not found\n')

            target = list(genes['GeneB'])
            
# define Matplotlib figure and axis
            fig, ax = plt.subplots(figsize = (25,12))
                # plt.figure(figsize = (25,12))
            chromwidth = 0.01
            cyscale = 1/ (chrom.loc[0]['length']*1.1)
#                print(cyscale)
            cwidth = 0.005
            centsize = 0.5
            xstart = 0.04
            ystart = 0.07
            # plt.rcParams.update({'font.size': 30})
            plt.tick_params(left = False, right = False , labelleft = False ,
                                labelbottom = False, bottom = False)

# draw chromosomes
            
            for i in range(len(chrom)):
                xpos = xstart+i*8*cwidth
                pstart = 0
                plen = chrom.iloc[i]['centromere']
                qstart = chrom.iloc[i]['centromere']
                qlen = chrom.iloc[i]['length']-chrom.iloc[i]['centromere']
                
                ax.plot(xpos, 1-(ystart+pstart*cyscale),xpos,1-(ystart+pstart*cyscale)-plen*cyscale,lw=50,solid_capstyle='round', color='k')
                ax.add_patch(Rectangle((xpos, 1-(ystart+pstart*cyscale)), 0.02, -plen*cyscale,fc = '0.8', ec='w'))
                ax.add_patch(Rectangle((xpos, 1-(ystart+qstart*cyscale)), 0.02, -qlen*cyscale,fc = '0.7', ec='w'))
            
                ax.add_patch(Ellipse((xpos+(centsize/50)+0.001,1-(ystart+qstart*cyscale)), width=centsize/20, height=centsize/25,fc = 'w', ec='k'))
                ax.text(xpos+0.01,0.94,chrom.iloc[i]['chromosome'],horizontalalignment='center',fontsize='30')
        
            plt.rcParams.update({'font.size': 40})
            ax.text(xstart-0.025,0.75,'p',horizontalalignment='center',style='italic',fontsize='30')
            ax.text(xstart-0.025,0.25,'q',horizontalalignment='center',style='italic',fontsize='30')
           
            loc_list={}
                
            for t in target:
                d = ChroLoc(t,hgnc,ensem)
                if d == None:
                    continue
                loc_list[t] = d[2]
                
#                    print(f'{t} - {d[0]},{d[1]},{d[2]}')
# add lines for genes                                
                icro = chrom['chromosome'].tolist().index(d[0])
                xpos = xstart+icro*8*cwidth
                ypos = 1-(ystart+d[1]*cyscale)
                ax.plot([xpos,xpos+0.02],[ypos,ypos],color='red')
                ax.text(xpos-0.02,ypos,t,fontsize='10')
                # plt.plot([xpos,xpos+0.02],[ypos,ypos],color='red')

            plt.show()
                
# sort list on the basis of cytogenetic location
                
            loc_list = dict(sorted(loc_list.items(), key=lambda x:x[1]))
            for i in loc_list:
                print(f'Gene {i}-{loc_list[i]}')
                
            
        
        elif task.lower() == 'hubs':
# =============================================================================
# cluster and rank multiply occuring GeneAs and GeneBs
# =============================================================================
            default_dep_dat_file = 'AB_Dep_AllByAll_pairs'
            
            f = input(f'Gene AB pair file [{default_dep_dat_file}] : ')
            if f.strip() == '' :
                f = default_dep_dat_file
            print(f'Loading {f}')
            
            
            pairs = pd.read_csv(f)
            pairs.sort_values(by = ['p'], inplace = True)
            print(pairs)
            
# =============================================================================
# make unique GeneA and GeneB lists
# 
# =============================================================================
            uGa = list(set([pairs['GeneA'].iloc[i] for i in range(len(pairs))]))
            uGb = list(set([pairs['GeneB'].iloc[i] for i in range(len(pairs))])) 

            GA = {}
            GB = {}
    
# =============================================================================
# loop through pairs finding GeneAs and GeneBs and build lists of their interactions
#             
# =============================================================================
            
            for i in range(len(pairs)):
                gA = pairs['GeneA'].iloc[i]
                gB = pairs['GeneB'].iloc[i]
                
                if gA in GA: # GeneA already there so add its pair to list
                    GA[gA].append(gB)
                else:
                    GA[gA] = [gB]
                    
                if gB in GB: # GeneA already there so add its pair to list
                    GB[gB].append(gA)
                else:
                    GB[gB] = [gA]
                       
# =============================================================================
# sort GA and GB to put largest hub at top
# 
# =============================================================================
            GA = dict(sorted(GA.items(), key = lambda x: len(x[1]), reverse = True))
            GB = dict(sorted(GB.items(), key = lambda x: len(x[1]), reverse = True))

# =============================================================================
# for each hub make graph from ordinal numbers of gene in gene pairs
# 
# =============================================================================

# =============================================================================
# ADD IN GRAPH STUFF            
# =============================================================================
            

        elif task.lower() == 'all' :
# =============================================================================
# plot hits from all by all run
#  █████╗ ██╗     ██╗     
# ██╔══██╗██║     ██║     
# ███████║██║     ██║     
# ██╔══██║██║     ██║     
# ██║  ██║███████╗███████╗
# ╚═╝  ╚═╝╚══════╝╚══════╝
# =============================================================================
            while True:
                
                try:
                    default_dep_dat_file = 'AB_Dep_AllByAll_pairs'
                    
                    f = input(f'Gene AB pair file [{default_dep_dat_file}] : ')
                    if f.strip() == '' :
                        f = default_dep_dat_file
                    pairs = pd.read_csv(f)
                    print(f'Loading {f}')

                    pairs.sort_values(by = 'p', inplace = True)
        
                    print(pairs)
                    break
                except:
                    print(f'File {f} not found')
            

            
            usemuts = usecnv = useexp = False
            while True:
                s = []
                w = input('Use Mutations (M), Copy Number Deletion (D), Underexpression (U) :').lower()
                if len(w) == 0:
                    continue
                if 'm' in w:
                    usemuts = True
                    s.append('mutations')
                if 'd' in w:
                    usecnv = True
                    s.append('deletions')
                if 'u' in w:
                    useexp = True
                    s.append('under-expression')
                    
                if usemuts or usecnv or useexp:
                    print(f'Using {Flistand(s)} ')
                    break
                else:
                    print(f'No valid choices in input : {w}')
# =============================================================================
# 
# choose whether to plot GeneA dependency switch plots, chromosomal location, or
# ranked GeneA plots
# 
# =============================================================================
            while True:
                t = input('GeneA Switch plot (S),GeneB Target plot (T), GeneB Chromosomal location (C), QUIT (Q): [S]').lower()
                if t != '' and t[0] =='q' :
                    break
                elif t == '' or t == 's': #SWITCH PLOT
                    
                    i = 0
                    while True:
                        gA = pairs['GeneA'].iloc[i]
                        gB = pairs['GeneB'].iloc[i]
                        pAB = pairs['p'].iloc[i]
                        mut,notmut = PairDist(gA,gB,Cdict,deps,usemuts,usecnv,useexp)
                        if len(mut) == 0 or len(notmut) == 0:
                            print(f'Zero distribution for {gA} with {gB} : {len(mut)},{len(notmut)}')
                            i += 1
                            continue
                        r = input(f'{gA} : {gB} p = {pAB:.5} - Plot or Quit ([P]/Q)')
                        
                        if r == '' or r[0].lower() != 'q':    
                            if useexp :
                                labs = ['under','not under']
                            else:
                                labs = ['mutated','not mutated']

                            
                            
                            plt.violinplot((mut,notmut),showmeans=True, showextrema=False)
                            plt.xticks([1,2],labels = labs,fontsize = '20')
                            plt.yticks(fontsize = '20')
                            plt.tick_params(axis='both', labelsize=12)
                            plt.title(f'{gA} : {gB} p = {pAB:.5}',fontsize = '20',pad=20)
                            plt.ylabel('dependency', fontsize='20')
                            plt.show()
                            i += 1
                        
                        else:
                            break
                        
                elif t == 'c': # CHROMOSOME LOCATION
                
# load chromosome length data

                    if not chromlen_loaded:    
                        default_chromolen_file = 'ChromosomeLengths.csv'
                        f = input(f'Chromosome length file [{default_chromolen_file}] : ')
                        if f.strip() == '' :
                            f = default_chromolen_file
                        print(f'Loading {f}')
                        chrom = pd.read_csv(f)
                        chromlen_loaded = True

# load HUGO chromosomal location data
                    
                    if not hgnc_loaded:
                        default_HUGO_file = 'hgnc_complete_set.tsv'
                        f = input(f'HUGO Gene file [{default_HUGO_file}] : ')
                        if f.strip() == '' :
                            f = default_HUGO_file
                        print(f'Loading {f}')
                        hgnc = pd.read_table(f)
                        hgnc = hgnc[['symbol','ensembl_gene_id','location','location_sortable','prev_symbol']]
                        hgnc.set_index('symbol',inplace = True)
                        hgnc_loaded = True
                        
# load Ensemble refseq gene-to-transcript mapping file

                    if not ensem_loaded:
                        default_ENSEMBLE_file = 'Homo_sapiens.GRCh38.110.refseq.tsv'
                        f = input(f'ENSEMBLE reference sequence file [{default_ENSEMBLE_file}] : ')
                        if f.strip() == '' :
                            f = default_ENSEMBLE_file
                        print(f'Loading {f}')
                        ensem = pd.read_table(f)
                        ensem.drop_duplicates(subset = ['gene_stable_id'],inplace = True)
                        ensem = ensem[['gene_stable_id','transcript_stable_id']]
                        ensem.set_index('gene_stable_id',inplace = True)
                        ensem_loaded = True
                       
# make unique list of GenesA in p value order                    
                    

                    uGene = []
                    for i in pairs['GeneA']:
                        if i not in uGene:
                            uGene.append(i)
                    print(f'{len(uGene)} unique GeneAs found')

# choose whether to specify GeneAs or go through in order

                    choose = True
                    ans = input('Specify GeneAs to analyse or List in p order ? ([S],L) : ').upper()
                    if len(ans) == 0 or ans[0] == 'S':
                        cGene = []
                        while True:
                            gg = input('Gene names to analyse or QUIT ? : ').upper()
                            if gg == 'QUIT':
                                break
                            cgs = gg.replace(',',' ').split()
                            for i in cgs:
                                if i in uGene:
                                    cGene.append(i)
                                else:
                                    print(f'Gene {i} not found')
                        if len(cGene) > 0:
                            uGene = cGene
                        else:
                            print('No genes specified - reverting to ordered list')
                            
# go through unique GeneAs accumulating GeneBs that swithed them and their chromosomal locations                    

                    for gA in uGene:
                        glist = {}
                        
                        print(f'finding GenesB for {gA}')
                        for i in range(len(pairs)):
                            if pairs.iloc[i]['GeneA'] == gA :
                                hgid = pairs['GeneB'].iloc[i]

# obtain Ensemble transcript from HUGO name

                                if hgid in hgnc.index:
                                    egid = hgnc.loc[hgid]['ensembl_gene_id']
                                elif hgid in list(hgnc['prev_symbol']):
                                    egid = hgnc.iloc[hgnc['prev_symbol'].tolist().index(hgid)]['ensembl_gene_id']
                                else:
                                    print(f'GeneB : {hgid} not found in HUGO symbol lists - ignored')
                                    continue
                                # print(f'HUGO : {hgid} = ENSG : {egid}')
                                if egid in ensem.index:
                                    etid = ensem.loc[egid]['transcript_stable_id']
                                else:
                                    print(f'GeneA : {egid} for HUGO symbol {hgid} not found in ENSEMBLE transcripts - ignored')
                                    continue
                                # print(f'ENSG: {egid} = ENST : {etid}')
# use Ensemble API to get data                                

                                server = "https://rest.ensembl.org"
                                ext = f"/map/cdna/{etid}/1..1?"

                                r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                                 
                                if not r.ok:
                                  r.raise_for_status()
                                  sys.exit()
                                 
                                d = r.json()
                                chro = d['mappings'][0]['seq_region_name']
                                nucl = int(d['mappings'][0]['start'])
                                glist[hgid] = (chro,nucl)              
                                # print(f'{hgid} {egid} {etid} {chro} {nucl}')
                                print(f"{hgid} - {hgnc.loc[hgid]['location']}")

# get chromosomal info for GeneA as well

                        if gA in hgnc.index:
                            egid = hgnc.loc[gA]['ensembl_gene_id']
                        elif gA in list(hgnc['prev_symbol']):
                            egid = hgnc.iloc[hgnc['prev_symbol'].tolist().index(gA)]['ensembl_gene_id']
                        else:
                            print(f'GeneB : {hgid} not found in HUGO symbol lists - ignored')
                            continue
                        # print(f'HUGO : {hgid} = ENSG : {egid}')
                        if egid in ensem.index:
                            etid = ensem.loc[egid]['transcript_stable_id']
                        else:
                            print(f'GeneA : {egid} for HUGO symbol {gA} not found in ENSEMBLE transcripts - ignored')
                            continue
                        # print(f'ENSG: {egid} = ENST : {etid}')

# use Ensemble API to get data                                

                        server = "https://rest.ensembl.org"
                        ext = f"/map/cdna/{etid}/1..1?"

                        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                         
                        if not r.ok:
                          r.raise_for_status()
                          sys.exit()
                         
                        d = r.json()
                        chro = d['mappings'][0]['seq_region_name']
                        nucl = int(d['mappings'][0]['start'])
                        gA_info = (chro,nucl)              
                        print(f'for GeneA : {gA} {egid} {etid} {chro} {nucl}')

# print list of GeneAs and their chromosomal locations
                        r = input(f'{len(glist)} targets for {gA} - Plot or Quit ([P]/Q)')
                        if r == '' or r[0].lower() != 'q':    

#define Matplotlib figure and axis
                            fig, ax = plt.subplots(figsize = (25,12))
                            
                            chromwidth = 0.01
                            cyscale = 1/ (chrom.loc[0]['length']*1.1)
                            print(cyscale)
                            cwidth = 0.005
                            centsize = 0.5
                            xstart = 0.04
                            ystart = 0.07
                            plt.rcParams.update({'font.size': 30})
                            plt.tick_params(left = False, right = False , labelleft = False ,
                                            labelbottom = False, bottom = False)
    
    # draw chromosomes
                            
                            for i in range(len(chrom)):
                                xpos = xstart+i*8*cwidth
                                pstart = 0
                                plen = chrom.iloc[i]['centromere']
                                qstart = chrom.iloc[i]['centromere']
                                qlen = chrom.iloc[i]['length']-chrom.iloc[i]['centromere']
                                
                                # ax.plot(xpos, 1-(ystart+pstart*cyscale),xpos,1-(ystart+pstart*cyscale)-plen*cyscale,lw=50,solid_capstyle='round', color='k')
                                ax.add_patch(Rectangle((xpos, 1-(ystart+pstart*cyscale)), 0.02, -plen*cyscale,fc = '0.8', ec='w'))
                                ax.add_patch(Rectangle((xpos, 1-(ystart+qstart*cyscale)), 0.02, -qlen*cyscale,fc = '0.7', ec='w'))
                            
                                ax.add_patch(Ellipse((xpos+(centsize/50)+0.001,1-(ystart+qstart*cyscale)), width=centsize/20, height=centsize/25,fc = 'w', ec='k'))
                                ax.text(xpos+0.01,0.94,chrom.iloc[i]['chromosome'],horizontalalignment='center')
                        
                            plt.rcParams.update({'font.size': 40})
                            ax.text(xstart-0.025,0.75,'p',horizontalalignment='center',style='italic')
                            ax.text(xstart-0.025,0.25,'q',horizontalalignment='center',style='italic')
                            

                            for g in glist:
                                print(f'{g} - {glist[g][0]},{glist[g][1]}')
# add lines for GeneAs                                
                                icro = chrom['chromosome'].tolist().index(glist[g][0])
                                xpos = xstart+icro*8*cwidth
                                ypos = 1-(ystart+glist[g][1]*cyscale)
                                ax.plot([xpos,xpos+0.02],[ypos,ypos],color='red')
                                
# add in Genea location
                            icro = chrom['chromosome'].tolist().index(gA_info[0])
                            xpos = xstart+icro*8*cwidth
                            ypos = 1-(ystart+gA_info[1]*cyscale)
                            ax.plot([xpos,xpos+0.02],[ypos,ypos],color='green',linewidth='3')
                            ax.text(xstart+21*8*cwidth+0.01,0.1,gA,horizontalalignment='center',color='green')

                            
 
                            plt.show()
                        else:
                            break

                elif t == 't':

# make unique list of GenesB in p value order                    
                    
                    uGene = []
                    for i in pairs['GeneB']:
                        if i not in uGene:
                            uGene.append(i)
                    print(f'{len(uGene)} unique GeneBs found')
                    
# go through unique GeneBs accumulating GeneAs and p values                    

                    maxy = -log10(pairs['p'].iloc[0])
                    plim = 0.001
                    for gB in uGene:
                        alist = []
                        plist = []
                        for i in range(len(pairs)):
                            if pairs.iloc[i]['GeneB'] == gB and pairs['p'].iloc[i]<=plim:
                                alist.append(pairs['GeneA'].iloc[i])
                                plist.append(-log10(pairs['p'].iloc[i]))
                        r = input(f'{len(alist)} targets for {gB} - Plot or Quit ([P]/Q)')
                        if r == '' or r[0].lower() != 'q':    
                            
                            plt.bar(alist, plist)
                            plt.tick_params(bottom = False)
                            plt.ylim(0.0,maxy)
                            plt.xticks(list(range(len(alist))),labels = alist, rotation = 'vertical',
                                       horizontalalignment = 'center')
                            plt.ylabel('-log(p)')
                            plt.title(f'Switching targets for {gB} with p < {plim}')

                            plt.show()
                        else:
                            break
                        
# =============================================================================
# ██████╗ ███████╗██████╗ ███████╗
# ██╔══██╗██╔════╝██╔══██╗██╔════╝
# ██║  ██║█████╗  ██████╔╝███████╗
# ██║  ██║██╔══╝  ██╔═══╝ ╚════██║
# ██████╔╝███████╗██║     ███████║
# ╚═════╝ ╚══════╝╚═╝     ╚══════╝
#                                                         
# =============================================================================
                        
        elif task.lower() == 'deps' :
            
            first_call = True

# check that profiles are loaded
            while not prof_loaded:
                profiles = {}
                list_profile_files()
                default_prof_file = 'AB_Dep_Profiles'
                f = input(f'Profile file [{default_prof_file}] : ')
                if f.strip() == '' :
                    f = default_prof_file
                print(f'Loading {f}')
                if prof_load(f):
                    print(f'{f} loaded')
                    prof_loaded = True
                    break
                else:
                    continue
            if first_call and len(profiles) > 0:
                print('Available predefined genetic profiles :')
                for n in sorted(profiles):
                    print(n)
                    first_call = False

            while True:
                t = input('\nTarget gene and defined genetic profile(s) (QUIT): ')
                if len(t) == 0:
                    continue
                if t.lower() == 'quit' :
                    break
                if t[0] == '?':
                    print('Available predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                        first_call = False
                    
# break target names into list
                target = t.upper().replace(',',' ').split()
# check that target is valid gene name
                if target[0] not in deps.columns:
                    print(f'Target gene : {target[0]} is not found')
                    continue
                if len(target) < 2:
                    print('No genetic profile provided')
                    continue
# check that profile(s) exists
                m = match_keyword(target[1],profiles)
                if  m == '*':
                    print(f'{target[1]} is ambiguous')
                    continue
                if m not in profiles:
                    print(f'Genetic profile : {target[1]} is not defined\nDo this in LINE\n')
                    continue
                target[1] = m
                
                prof = profiles[target[1]]
# generate cell line set
                p = list(eval(prof_parse(prof)))
                if len(p) > 0:
                    f = f'{target[1]}-{target[0]}_deps.tsv'
                    print(f'Saving to {f}')
                    fl = open(f,'w')
                    fl.write(f'DepMapID\tCCLEName\t{target[0]}_dependency\n')
                
                for h in p:
                    f = info.loc[info['ModelID'] == h]
                    if len(f) > 0:
                        cname = f['CCLEName'].iloc[0]    
                        print(f'{h}  : {cname}')
                    else:
                        print(f'{h}  : no CCLE name found - skipping {h}')
                        continue
# get dependency on target[0]
                
                    if h in deps.index:
                        d = deps[target[0]].loc[h]
                        fl.write(f'{h}\t{cname}\t{d}\n')
                    else:
                        fl.write(f'{h}\t{cname}\tNaN\n')
                
                fl.close()
                    
            

        elif task.lower() == 'mine' :
# =============================================================================
# ███╗   ███╗██╗███╗   ██╗███████╗
# ████╗ ████║██║████╗  ██║██╔════╝
# ██╔████╔██║██║██╔██╗ ██║█████╗  
# ██║╚██╔╝██║██║██║╚██╗██║██╔══╝  
# ██║ ╚═╝ ██║██║██║ ╚████║███████╗
# ╚═╝     ╚═╝╚═╝╚═╝  ╚═══╝╚══════╝
#                                 #                            
# Plot and analyse target dependency on genetic profile match/nomatch or
# comparison of two genetic profiles
# =============================================================================
   
            first_call = True

# check that profiles are loaded
            while not prof_loaded:
                profiles = {}
                list_profile_files()
                default_prof_file = 'AB_Dep_Profiles'
                f = input(f'Profile file [{default_prof_file}] : ')
                if f.strip() == '' :
                    f = default_prof_file
                print(f'Loading {f}')
                if prof_load(f):
                    print(f'{f} loaded')
                    prof_loaded = True
                    break
                else:
                    continue
            if first_call and len(profiles) > 0:
                print('Available predefined genetic profiles :')
                for n in sorted(profiles):
                    print(n)
                    first_call = False

            while True:
                t = input('\nTarget gene and defined genetic profile(s) (QUIT): ')
                if len(t) == 0:
                    continue
                if t.lower() == 'quit' :
                    break
                if t[0] == '?':
                    print('Available predefined genetic profiles :')
                    for n in sorted(profiles):
                        print(n)
                        first_call = False
                    
# break target names into list
                target = t.upper().replace(',',' ').split()
# check that target is valid gene name
                if target[0] not in deps.columns:
                    print(f'Target gene : {target[0]} is not found')
                    continue
                if len(target) < 2:
                    print('No genetic profile provided')
                    continue
# check that profile(s) exists
                m = match_keyword(target[1],profiles)
                if  m == '*':
                    print(f'{target[1]} is ambiguous')
                    continue
                if m not in profiles:
                    print(f'Genetic profile : {target[1]} is not defined\nDo this in LINE\n')
                    continue
                target[1] = m
                
                if len(target) == 3:
                    m = match_keyword(target[2],profiles)
                    if  m == '*':
                        print(f'{target[2]} is ambiguous')
                        continue
                    if m not in profiles:
                        print(f'Genetic profile : {target[2]} is not defined\nDo this in LINE\n')
                        continue
                    target[2] = m
    
                print(f'Plotting dependency on target gene {target[0]}')
                print(f'in cells with genetic profile {target[1]}')
                if len(target) == 3:
                    print(f'compared to cells with genetic profile {target[2]}')

                prof = profiles[target[1]]
# generate hypothesis set
                p = list(eval(prof_parse(prof)))
                
# if only one profile given generate null set
                if len(target) == 2:
                    q = list(eval(prof_nullhype(prof_flatten(prof))))
                    pd = set(p) & set(deps.index)
                    qd = set(q) & set(deps.index)
                    print(f'{len(pd)} cell lines with dependency data in matching profile set')
                    print(f'{len(qd)} cell lines with dependency data in null set')
                elif len(target) == 3: # two profiles given
                    prof2 = profiles[target[2]]
                    q = list(eval(prof_parse(prof2)))
# check for null sets

                if len(p) == 0:
                    print('Genetic profile : {[target[1]} returns empty set')
                    continue
                if len(q) == 0:
                    if len(target) == 2:
                        print('Null set for Genetic profile : {[target[1]} is empty')
                    elif len(target) == 3:
                        print('Genetic profile : {[target[2]} returns empty set')
    
                    continue

                match,notmatch = ProfPairDist(target[0],p,q,deps,null_trim)
                r = ProfPairPvalue(target[0],prof,p,q,deps,thresh,null_trim) # returns [target name,p]
                pAB = r[1]
                
# calculate Cohen's D

                with np.errstate(divide = 'ignore'):
                    D = (np.array(match).mean()-np.array(notmatch).mean())/np.array(match+notmatch).std()
                print(f'Cohen d value - {D}')
                
                if len(target) == 2:
                    DepPlot2(match,notmatch,target[0],target[1].replace('_','-'),pAB,D)
                elif len(target) == 3:
                    DepPlot2(match,notmatch,target[0],target[1].replace('_','-'),pAB,D,l1=target[1].replace('_','-'),l2=target[2].replace('_','-'))

# initialise factor list

                factors = []
                froot = '('+prof_flatten(prof)+')'
                nfact = 0
                d1 = deps[target[0]].reset_index()

                
# option to analyse sex distribution

                t = input('Analyse sex distribution (Y,[N]) : ')
                if t != '' and t[0].upper() == 'Y':
# get sex type count for the hypothesis set
                    sex = {}
                    nincl = 0
                    for i in p:
                        c = info.loc[info['ModelID'] == i]
                        if len(c) > 0:
                            cname = c['Sex'].iloc[0]
                            nincl += 1
                        else:
                            continue
                        if cname in sex:
                            sex[cname] += 1
                        else:
                            sex[cname] = 1
                    print(f'{nincl} cell lines included')
                    sex = dict(sorted(sex.items(), key = lambda x: x[1], reverse = True))                    

                    csig = {}
                    cD = {}
                    for c in sex:
                        if c == 'Unknown':
                            continue
                        pdis = set(p) & SEX(c)
                        qdis = set(p) - pdis
                        # d1 = deps[target[0]].reset_index()

                        dmut = list(d1[d1['ModelID'].isin(pdis)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(qdis)][target[0]])
                        
                        with np.errstate(divide = 'ignore', invalid = 'ignore'):
                            Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                            if len(dmut) > 0 and len(dwt) > 0:
                                Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()
                            else:
                                Dd = 0
                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            cD[c] = Dd
                            if np.average(dmut) > np.average(dwt):
                                csig[c] = Pval  # mutation of m is overrepresented in high dependency cell lines
                            else:
                                csig[c] = -Pval  # -ve signifies m is underrepresented
                    if len(csig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        csig = dict(sorted(csig.items(), key = lambda x: x[1], reverse = False))
                        for c in csig:
                            pdis = set(p) & SEX(c)
                            qdis = set(p) - pdis
                            print(f'{len(pdis)} cells derive from {c.lower()}')
                            if csig[c] < 0:
                                tt = 'low'
                            else:
                                tt = 'high'
                            print(f'Sex origin {c} is highly significantly over-represented - p = {abs(csig[c])} d = {cD[c]}\n',
                                  f'in profile matching cells with {tt} dependency on {target[0]}')
                            
                            
                            
                            t = input('Plot affected cells ([Y],N,QUIT) : ')
                  
                            if t != '' and t[0].upper() == 'Q':
                                break
                            if t != '' and t[0].upper() == 'N':
                                continue
# encode additional factor
                            if tt == 'high': # add this as a positive factor
                                factors.append(f' *{c} ')
                            else :
                                factors.append(f' (# ^ *{c}) ')
                            nfact += 1
                            
                            gdeps,nomatch = ProfPairDist(target[0],pdis,qdis,deps,null_trim) # get deps for these cells
                            fig = plt.figure(figsize=(3,4))
                            width = 0.02
                            plt.violinplot((match),showmeans=True, showextrema=False)
    
                            plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                            data = np.array(gdeps)
                            x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                            plt.scatter(x, data, s=10,c='r')
                            
                            plt.text(1.1,0.5,f'{c.lower()}({len(gdeps)})\np = {abs(csig[c]):.4} d = {cD[c]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                            # plt.title(f'{h} : {m} p = {HitList[h]:.4g} d = {cohen[h]:.4g}',fontsize = '20',pad=20)
                            gmean = np.array(gdeps).mean()
                            plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')
                            
                            plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                            plt.ylim(0.0,1.0)
                            plt.tick_params(axis='both', labelsize=12)
                            plt.ylabel('dependency',fontsize = '20')
                            plt.text(0.1,0.1,'sex',horizontalalignment='left',style='italic',fontsize='15',color='r')

                            plt.show()
                            

# option to analyse age distribution

# option to analyse sex distribution

                t = input('Analyse age distribution (Y,[N]) : ')
                if t != '' and t[0].upper() == 'Y':
# get sex type count for the hypothesis set
                    age = {}
                    nincl = 0
                    for i in p:
                        c = info.loc[info['ModelID'] == i]
                        if len(c) > 0:
                            cname = c['AgeCategory'].iloc[0]
                            nincl += 1
                        else:
                            continue
                        if cname in age:
                            age[cname] += 1
                        else:
                            age[cname] = 1

                    print(f'{nincl} cell lines included')
                    age = dict(sorted(age.items(), key = lambda x: x[1], reverse = True))                    

                    csig = {}
                    cD = {}
                    for c in age:
                        if c == 'Unknown':
                            continue
# get subset of hypothesis set that have this cancer type
                        pdis = set(p) & AGE(c)
                        qdis = set(p) - pdis
                        # d1 = deps[target[0]].reset_index()

                        dmut = list(d1[d1['ModelID'].isin(pdis)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(qdis)][target[0]])
                        
                        with np.errstate(divide = 'ignore', invalid = 'ignore'):
                            Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                            if len(dmut) > 0 and len(dwt) > 0:
                                Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()
                            else:
                                Dd = 0
                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            cD[c] = Dd
                            if np.average(dmut) > np.average(dwt):
                                csig[c] = Pval  # mutation of m is overrepresented in high dependency cell lines
                            else:
                                csig[c] = -Pval  # -ve signifies m is underrepresented
                    if len(csig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        csig = dict(sorted(csig.items(), key = lambda x: x[1], reverse = False))
                        for c in csig:
                            pdis = set(p) & AGE(c)
                            qdis = set(p) - pdis
                            print(f'{len(pdis)} cells derive from {c.lower()}')
                            if csig[c] < 0:
                                tt = 'low'
                            else:
                                tt = 'high'
                            print(f'Age origin {c} is highly significantly over-represented - p = {abs(csig[c])} d = {cD[c]}\n',
                                  f'in profile matching cells with {tt} dependency on {target[0]}')

                            t = input('Plot affected cells ([Y],N,QUIT) : ')
                  
                            if t != '' and t[0].upper() == 'Q':
                                break
                            if t != '' and t[0].upper() == 'N':
                                continue

                            if tt == 'high': # add this as a positive factor
                                factors.append(f'  ;{c} ')
                            else :
                                factors.append(f' (# ^ ;{c}) ')
                            nfact += 1
                                                        
                            gdeps,nomatch = ProfPairDist(target[0],pdis,qdis,deps,null_trim) # get deps for these cells
                            fig = plt.figure(figsize=(3,4))
                            width = 0.02
                            plt.violinplot((match),showmeans=True, showextrema=False)
    
                            plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                            data = np.array(gdeps)
                            x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                            plt.scatter(x, data, s=10,c='r')
                            
                            plt.text(1.1,0.5,f'{c.lower()}({len(gdeps)})\np = {abs(csig[c]):.4} d = {cD[c]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                            # plt.title(f'{h} : {m} p = {HitList[h]:.4g} d = {cohen[h]:.4g}',fontsize = '20',pad=20)
                            gmean = np.array(gdeps).mean()
                            plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')
                            
                            plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                            plt.ylim(0.0,1.0)
                            plt.tick_params(axis='both', labelsize=12)
                            plt.ylabel('dependency',fontsize = '20')
                            plt.text(0.1,0.1,'age',horizontalalignment='left',style='italic',fontsize='15',color='r')

                            plt.show()
    
# now give option to analyse tissue and/or mutational patterns of cell lines above and below mean of matched distribution

                t = input('Analyse tissue distribution (Y,[N]) : ')
                if t != '' and t[0].upper() == 'Y':
# get cancer type count for the hypothesis set
                    cancer = {}
                    nincl = 0
                    for i in p:
                        c = info.loc[info['ModelID'] == i]
                        if len(c) > 0:
                            cname = c['OncotreeLineage'].iloc[0]
                            nincl += 1
                        else:
                            continue
                        if cname in cancer:
                            cancer[cname] += 1
                        else:
                            cancer[cname] = 1
                    print(f'{nincl} cell lines included')
                    cancer = dict(sorted(cancer.items(), key = lambda x: x[1], reverse = True))                    
                
                    csig = {}
                    cD = {}
                    for c in cancer:
# get subset of hypothesis set that have this cancer type
                        pdis = set(p) & TYP(c)
                        qdis = set(p) - pdis
                        # d1 = deps[target[0]].reset_index()

                        dmut = list(d1[d1['ModelID'].isin(pdis)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(qdis)][target[0]])
                        
                        with np.errstate(divide = 'ignore', invalid = 'ignore'):
                            Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                            if len(dmut) > 0 and len(dwt) > 0:
                                Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()
                            else:
                                Dd = 0
                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            cD[c] = Dd
                            if np.average(dmut) > np.average(dwt):
                                csig[c] = Pval  # mutation of m is overrepresented in high dependency cell lines
                            else:
                                csig[c] = -Pval  # -ve signifies m is underrepresented
                    if len(csig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        csig = dict(sorted(csig.items(), key = lambda x: x[1], reverse = False))
                        for c in csig:
                            pdis = set(p) & TYP(c)
                            qdis = set(p) - pdis
                            print(f'{len(pdis)} cells derive from {c.lower()}')
                            if csig[c] < 0:
                                tt = 'low'
                            else:
                                tt = 'high'
                            print(f'Tissue origin {c} is highly significantly over-represented - p = {abs(csig[c])} d = {cD[c]}\n',
                                  f'in profile matching cells with {tt} dependency on {target[0]}')
                            t = input('Plot affected cells ([Y],N,QUIT) : ')
                  
                            if t != '' and t[0].upper() == 'Q':
                                break
                            if t != '' and t[0].upper() == 'N':
                                continue
                            
                            if tt == 'high': # add this as a positive factor
                                factors.append(f' %{c.split()[0]} ')
                            else :
                                factors.append(f' (# ^  %{c.split()[0]}) ')
                            nfact += 1

                            gdeps,nomatch = ProfPairDist(target[0],pdis,qdis,deps,null_trim) # get deps for these cells
                            fig = plt.figure(figsize=(3,4))
                            width = 0.02
                            plt.violinplot((match),showmeans=True, showextrema=False)
    
                            plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                            data = np.array(gdeps)
                            x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                            plt.scatter(x, data, s=10,c='r')
                            
                            plt.text(1.1,0.5,f'{c.lower()}({len(gdeps)})\np = {abs(csig[c]):.4} d = {cD[c]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                            # plt.title(f'{h} : {m} p = {HitList[h]:.4g} d = {cohen[h]:.4g}',fontsize = '20',pad=20)
                            gmean = np.array(gdeps).mean()
                            plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')
                            
                            plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                            plt.ylim(0.0,1.0)
                            plt.tick_params(axis='both', labelsize=12)
                            plt.ylabel('dependency',fontsize = '20')
                            plt.text(0.1,0.1,'tissue',horizontalalignment='left',style='italic',fontsize='15',color='r')

                            plt.show()


                t = input('Analyse mutational signature distribution (Y,[N]) : ')
                # sig_thresh = 0.1
                if t != '' and t[0].upper() == 'Y':
                    f = input(f'Minimum fraction of cell mutational burden to consider a signature [{sig_thresh}] : ')
                    if f.strip() != '' :
                        try:
                            sig_thresh = float(f)
                        except ValueError:
                            print(f'Invalid value : using {sig_thresh}')
                    if sig_thresh < 0.0 or sig_thresh > 1.0:
                        sig_thresh = 0.05
                        print(f'Invalid value : resetting {sig_thresh}')

# get signature type count above threshold for the hypothesis set
                    signatures = {}
                    psigs = sigs.loc[sigs['SAMPLES'].isin(p)] # data for cells in p
                    pd.set_option('mode.chained_assignment', None)
                    psigs[psigs.columns[1:-1]] = psigs[psigs.columns[1:-1]].div(sigs['Burden'],axis = 0)
                    pd.set_option('mode.chained_assignment', 'warn')
                    print(f'{len(psigs)} cell lines included')
                    
                    for s in psigs.columns[1:-1]:
                        sl = list(psigs[psigs[s] > sig_thresh]['SAMPLES'])
                        # print(s,sl)
                        if len(sl) > 0:
                            signatures[s] = len(sl)
                            
                    signatures = dict(sorted(signatures.items(), key = lambda x: x[1], reverse = True))                    
                
                    csig = {}
                    cD = {}
                    for c in signatures:
# get subset of hypothesis set that have this signature type
                        pdis = set(p) & SIG(c)
                        qdis = set(p) - pdis
                        # d1 = deps[target[0]].reset_index()

                        dmut = list(d1[d1['ModelID'].isin(pdis)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(qdis)][target[0]])
                        
                        with np.errstate(divide = 'ignore', invalid = 'ignore'):
                            Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                            if len(dmut) > 0 and len(dwt) > 0:
                                Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()
                            else:
                                Dd = 0
                        # print(f'c = {c}\ndmut = {dmut}\ndwt = {dwt}\np = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            cD[c] = Dd
                            if np.average(dmut) > np.average(dwt):
                                csig[c] = Pval  # signature c is overrepresented in high dependency cell lines
                            else:
                                csig[c] = -Pval  # -ve signifies c is underrepresented
                    if len(csig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        csig = dict(sorted(csig.items(), key = lambda x: x[1], reverse = False))
                        for c in csig:
                            pdis = set(p) & SIG(c)
                            qdis = set(p) - pdis
                            print(f'{len(pdis)} cells derive from {c.lower()}')
                            if csig[c] < 0:
                                tt = 'low'
                            else:
                                tt = 'high'
                            print(f'Mutational signature {c} ({WhichSigType(c)}) is highly significantly over-represented - p = {abs(csig[c])} d = {cD[c]}\n',
                                  f'in profile matching cells with {tt} dependency on {target[0]}')
                            t = input('Plot affected cells ([Y],N,QUIT) : ')
                  
                            if t != '' and t[0].upper() == 'Q':
                                break
                            if t != '' and t[0].upper() == 'N':
                                continue
                            
                            if tt == 'high': # add this as a positive factor
                                factors.append(f' {{{c.split()[0]} ') # {{{ is needed to escape { for signature prefix
                            else :
                                factors.append(f' (# ^  {{{c.split()[0]}) ')
                            nfact += 1

                            gdeps,nomatch = ProfPairDist(target[0],pdis,qdis,deps,null_trim) # get deps for these cells
                            fig = plt.figure(figsize=(3,4))
                            width = 0.02
                            plt.violinplot((match),showmeans=True, showextrema=False)
    
                            plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                            data = np.array(gdeps)
                            x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                            plt.scatter(x, data, s=10,c='r')
                            
                            plt.text(1.1,0.5,f'{c.lower()}({len(gdeps)})\np = {abs(csig[c]):.4} d = {cD[c]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                            # plt.title(f'{h} : {m} p = {HitList[h]:.4g} d = {cohen[h]:.4g}',fontsize = '20',pad=20)
                            gmean = np.array(gdeps).mean()
                            plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')
                            
                            plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                            plt.ylim(0.0,1.0)
                            plt.tick_params(axis='both', labelsize=12)
                            plt.ylabel('dependency',fontsize = '20')
                            plt.text(0.1,0.1,'signature',horizontalalignment='left',style='italic',fontsize='15',color='r')

                            plt.show()
                        

                frac_limit = 0.1    
                t = input('Analyse disruptive / LoF missense mutational distribution (Y,[N]) : ')
                if t != '' and t[0].upper() == 'Y':
                    f = input(f'Minimum fraction of cells to consider a gene effect [{frac_limit}] : ')
                    if f.strip() != '' :
                        try:
                            frac_limit = float(f)
                        except ValueError:
                            print(f'Invalid value : using {frac_limit}')
                    if frac_limit < 0.0 or frac_limit > 1.0:
                        frac_limit = 0.05
                        print(f'Invalid value : resetting {frac_limit}')
                    
                    sp = set(p)
                    sq = set(q)
                    n_limit = len(sp)*frac_limit
                    
# gather all disruptive mutations in profile matching cells
                    mall = []                    
                    for c in p:
                        if c in Cdict:
                            mall += Cdict[c].dgenes
                    muniq = set(mall)
                    
                    gcount = {} 
                    for m in muniq:
                        nm = mall.count(m)
                        if nm > n_limit:
                            gcount[m] = nm
                    gcount =   dict(sorted(gcount.items(), key = lambda x: x[1], reverse = True))
                    print(f'{len(gcount)} genes are mutated in > {(frac_limit)*100}% of the matching cell lines')
# go through genes calcualting t-test
                    gsig = {}
                    gD = {}
                    print('Mining for significant deviations ...')
                    for g in gcount:
                        # print(f'Gene {g} is disrupted in {gcount[g]}/{len(sp)} matching cells')
                        cmut = (DIS(g) | DEM(g,GeneAggLOF(g))) & sp # cells in matched set with m gene disruption
                        if cmut == set():
                            print(f'No disruptive /LoF mutations found for gene : {g} in matched set')
                            continue
                        cwt = sp - cmut 
                        if cwt == set():
                            print(f'No undisrupted cells found for gene : {g} in matched set')
                            continue
                        
                        # d1 = deps[target[0]].reset_index()
       
                        dmut = list(d1[d1['ModelID'].isin(cmut)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(cwt)][target[0]])
                        # print(f'{len(dmut)} mutant cell lines found')
                        with np.errstate(divide = 'ignore', invalid = 'ignore'):
                            Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                            if len(dmut) > 0 and len(dwt) > 0:
                                Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()
                            else:
                                Dd = 0

                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            gD[g] = Dd
                            if np.average(dmut) > np.average(dwt):
                                gsig[g] = Pval  # mutation of m is overrepresented in high dependency cell lines
                            else:
                                gsig[g] = -Pval  # -ve signifies m is underrepresented
                                
                    if len(gsig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        gsig = dict(sorted(gsig.items(), key = lambda x: x[1], reverse = False))
                        for g in gsig:
                            if gsig[g] < 0:
                                tt = 'low'
                            else:
                                tt = 'high'
                            print(f'Mutational disruption / LoF of {g} is highly significantly over-represented - p = {abs(gsig[g]):.4} d = {gD[g]:.4}\n',
                                  f'in profile matching cells with {tt} dependency on {target[0]}')
                            
                            t = input('Plot affected cells ([Y],N,QUIT) : ')
                            if t != '' and t[0].upper() == 'Q':
                                break
                            if t != '' and t[0].upper() == 'N':
                                continue
                                                        
                            if tt == 'high': # add this as a positive factor
                                factors.append(f' !{g} ')
                            else :
                                factors.append(f' (# ^  !{g}) ')
                            nfact += 1
                            
                            pdis = (DIS(g) | DEM(g,GeneAggLOF(g))) & sp # get set of cells from match set with 2ndry disrupted
                            gdeps,nomatch = ProfPairDist(target[0],pdis,q,deps,null_trim) # get deps for these cells
                            fig = plt.figure(figsize=(3,4))
                            width = 0.02
                            plt.violinplot((match),showmeans=True, showextrema=False)

                            plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                            data = np.array(gdeps)
                            x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                            plt.scatter(x, data, s=10,c='r')

                            plt.text(1.1,0.5,f'{g}({len(gdeps)})\np = {abs(gsig[g]):.4} d = {gD[g]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                            gmean = np.array(gdeps).mean()
                            plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')
                            plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                            plt.ylim(0.0,1.0)
                            plt.tick_params(axis='both', labelsize=12)
                            plt.ylabel('dependency',fontsize = '20')
                            plt.text(0.1,0.1,'loss of function',horizontalalignment='left',style='italic',fontsize='15',color='r')

                            plt.show()
                        print('\n')
                        t = input('Save the Secondary Hit List (Y,[N]) : ')
                        if t != '' and t[0].upper() == 'Y':    
                            fhit = open(f'{target[0]}-2ndryDisruptionHits.csv','w')
                            fhit.write('TARGET, P\n')
                            for g in gsig:
                                fhit.write(f'{g},{gsig[g]:.5}\n')
                            fhit.close()
                            print(f'{target[0]}-2ndryDisruptionHits.csv written')


                t = input('Analyse Gain of Function mutations (Y,[N]) : ')
                if t != '' and t[0].upper() == 'Y':
                    f = input(f'Minimum fraction of cells to consider a gene effect [{frac_limit}] : ')
                    if f.strip() != '' :
                        try:
                            frac_limit = float(f)
                        except ValueError:
                            print(f'Invalid value : using {frac_limit}')
                    if frac_limit < 0.0 or frac_limit > 1.0:
                        frac_limit = 0.05
                        print(f'Invalid value : resetting {frac_limit}')
                    
                    sp = set(p)
                    sq = set(q)
                    n_limit = len(sp)*frac_limit
# get all GOF mutations

                    GOF = muts.loc[(muts['LikelyGoF']) &
                                 (muts['VariantInfo'] == 'MISSENSE')][['ModelID','HugoSymbol','ProteinChange']]
                    GOF['mutation'] = GOF['HugoSymbol']+GOF['ProteinChange']
                    del GOF['HugoSymbol']
                    del GOF['ProteinChange']
                    
# get those cells in profile match that harbour GOF muts

                    pGOF = GOF[GOF['ModelID'].isin(sp)]
                    
                    
                    print(f'{len(pGOF)} cell lines in matched set harbour likely GOF mutations')
                    
# make unique list of which GOFs are present

                    uGOF = set(pGOF['mutation'])

                    print(f'{len(uGOF)} different GOF mutations in matched set')
 
# count occurences of the unique 
                    gcount = {}
                    
                    for m in uGOF:
                        lm = list(pGOF['mutation'])
                        nm = lm.count(m)
                        if nm > n_limit:
                            gcount[m] = nm
                    gcount =   dict(sorted(gcount.items(), key = lambda x: x[1], reverse = True))
                    print(f'{len(gcount)} genes have GOF mutation in > {(frac_limit)*100}% of the matching cell lines')

# go through genes calcualting t-test
                    gsig = {}
                    gD = {}
                    
                    for g in gcount:
                        GOFgene = g.split('p.')[0]
                        GOFchange = g.split('p.')[1]
                        
                        # print(f'Gene {g} is disrupted in {gcount[g]}/{len(sp)} matching cells')
                        cmut = DEM(GOFgene,GeneAggGOF(GOFgene)) & sp # cells in matched set with m defined GOF mutation
                        if cmut == set():
                            print(f'No GOF mutations found for gene : {g} in matched set')
                            continue
                        cwt = sp - cmut 
                        if cwt == set():
                            print(f'No unmutated cells found for gene : {g} in matched set')
                            continue
                        
                        # d1 = deps[target[0]].reset_index()
       
                        dmut = list(d1[d1['ModelID'].isin(cmut)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(cwt)][target[0]])
                        # print(f'{len(dmut)} mutant cell lines found')
                        with np.errstate(divide = 'ignore', invalid = 'ignore'):
                            Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                            if len(dmut) > 0 and len(dwt) > 0:
                                Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()
                            else:
                                Dd = 0

                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            gD[GOFgene] = Dd
                            if np.average(dmut) > np.average(dwt):
                                gsig[GOFgene] = Pval  # mutation of m is overrepresented in high dependency cell lines
                            else:
                                gsig[GOFgene] = -Pval  # -ve signifies m is underrepresented
                                
                    if len(gsig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        gsig = dict(sorted(gsig.items(), key = lambda x: x[1], reverse = False))
                        for g in gsig:
                            if gsig[g] < 0:
                                tt = 'low'
                            else:
                                tt = 'high'
                            print(f'Gain of Function mutation of {g} is highly significantly over-represented - p = {abs(gsig[g])} d = {gD[g]}\n',
                                  f'in profile matching cells with {tt} dependency on {target[0]}')
                            
                            t = input('Plot affected cells ([Y],N,QUIT) : ')
                            if t != '' and t[0].upper() == 'Q':
                                break
                            if t != '' and t[0].upper() == 'N':
                                continue
                                                        
                            if tt == 'high': # add this as a positive factor
                                factors.append(f' /{g}:{GeneAggGOF(g)} ')
                            else :
                                factors.append(f' (# ^  /{g}:{GeneAggGOF(g)}) ')
                            nfact += 1
                            
                            pdis = sp & DEM(g,GeneAggGOF(g))# get set of cells from match set with 2ndry disrupted
                            gdeps,nomatch = ProfPairDist(target[0],pdis,q,deps,null_trim) # get deps for these cells
                            fig = plt.figure(figsize=(3,4))
                            width = 0.02
                            plt.violinplot((match),showmeans=True, showextrema=False)
                            plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                            data = np.array(gdeps)
                            x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                            plt.scatter(x, data, s=10,c='r')
                            
                            plt.text(1.1,0.5,f'{g}({len(gdeps)})\np = {abs(gsig[g]):.4} d = {gD[g]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                            gmean = np.array(gdeps).mean()
                            plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')

                            plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                            plt.ylim(0.0,1.0)
                            plt.tick_params(axis='both', labelsize=12)
                            plt.ylabel('dependency',fontsize = '20')
                            plt.text(0.1,0.1,'gain of function',horizontalalignment='left',style='italic',fontsize='15',color='r')

                            plt.show()
                        print('\n')
                        t = input('Save the Secondary Hit List (Y,[N]) : ')
                        if t != '' and t[0].upper() == 'Y':    
                            fhit = open(f'{target[0]}-2ndryGOFHits.csv','w')
                            fhit.write('TARGET, P\n')
                            for g in gsig:
                                fhit.write(f'{g},{gsig[g]:.5}\n')
                            fhit.close()
                            print(f'{target[0]}-2ndryGOFHits.csv written')



                t = input('Analyse Copy Number Variation (Y,[N]) : ')
                if t != '' and t[0].upper() == 'Y':
                    f = input(f'Minimum fraction of cells to consider a gene effect [{frac_limit}] : ')
                    if f.strip() != '' :
                        try:
                            frac_limit = float(f)
                        except ValueError:
                            print(f'Invalid value : using {frac_limit}')
                    if frac_limit < 0.0 or frac_limit > 1.0:
                        frac_limit = 0.05
                        print(f'Invalid value : resetting {frac_limit}')
                    
                    sp = set(p)
                    sq = set(q)
                    n_limit = len(sp)*frac_limit
# gather all CNV deletions in profile matching cells
                    del_all = []
                    del_loc_all = []
                    for c in p:
                        if c in Cdict:
                            del_all += Cdict[c].cnvdel
                            del_loc_all += list(set(hgnc[hgnc.index.isin(Cdict[c].cnvdel)]['location_sortable']))
                    # del_uniq = set(del_all)
                    # del_loc_uniq = set(del_loc_all)
                    
# find frequent deleted loci

                    locdel_count = {key:val for key, val in 
                                    Counter(del_loc_all).items() if val > n_limit}
                    locdel_count =   dict(sorted(locdel_count.items(), key = lambda x: x[1], reverse = True))
                    print(f'{len(locdel_count)} loci contain genes that are deleted in > {(frac_limit)*100}% of the matching cell lines')
                    
                    gdel_count = {key:val for key, val in 
                                    Counter(del_all).items() if val > n_limit}                   
                    gdel_count =   dict(sorted(gdel_count.items(), key = lambda x: x[1], reverse = True))
                    print(f'{len(gdel_count)} individual genes are deleted in > {(frac_limit)*100}% of the matching cell lines')

# gather all CNV amplifications in profile matching cells

                    amp_all = [] 
                    amp_loc_all = []
                    for c in p:
                        if c in Cdict:
                            amp_all += Cdict[c].cnvamp
                            amp_loc_all += list(set(hgnc[hgnc.index.isin(Cdict[c].cnvamp)]['location_sortable']))
                    amp_uniq = set(amp_all)
                    amp_loc_uniq = set(amp_loc_all)
                    
# find frequent amplified loci

                    locamp_count = {key:val for key, val in 
                                    Counter(amp_loc_all).items() if val > n_limit}
                    locamp_count =   dict(sorted(locamp_count.items(), key = lambda x: x[1], reverse = True))
                    print(f'{len(locamp_count)} loci contain genes that are amplified in > {(frac_limit)*100}% of the matching cell lines')
                    
                    gamp_count = {key:val for key, val in 
                                    Counter(amp_all).items() if val > n_limit}                   

                    gamp_count =   dict(sorted(gamp_count.items(), key = lambda x: x[1], reverse = True))
                    print(f'{len(gamp_count)} genes are amplified in > {(frac_limit)*100}% of the matching cell lines')

# go through deletions calculating t-test

                    gsig = {}
                    gD = {}
                    
                    for g in gdel_count:
                        # print(f'Gene {g} is disrupted in {gcount[g]}/{len(sp)} matching cells')
                        cmut = DEL(g) & sp # cells in matched set with m gene deletion
                        if cmut == set():
                            print(f'No deletions found for gene : {g} in matched set')
                            continue
                        cwt = sp - cmut 
                        if cwt == set():
                            print(f'No non-deleted cells found for gene : {g} in matched set')
                            continue
                        
                        # d1 = deps[target[0]].reset_index()
       
                        dmut = list(d1[d1['ModelID'].isin(cmut)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(cwt)][target[0]])
                        # print(f'{len(dmut)} mutant cell lines found')
                        with np.errstate(divide = 'ignore', invalid = 'ignore'):
                            Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                            if len(dmut) > 0 and len(dwt) > 0:
                                Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()
                            else:
                                Dd = 0

                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            gD[g] = Dd
                            if np.average(dmut) > np.average(dwt):
                                gsig[g] = Pval  # mutation of m is overrepresented in high dependency cell lines
                            else:
                                gsig[g] = -Pval  # -ve signifies m is underrepresented
                                
                    if len(gsig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        gsig = dict(sorted(gsig.items(), key = lambda x: x[1], reverse = False))
                        
# IF MULTIPLE SIGNIFICANT GENES SHARE LOCUS - MAKE LOCUS FACTOR NOT GENES          

                        loci = Locify(list(gsig)) # cluster genes into loci
                        
                        for l in loci:
                            if len(loci[l]) > 1: # cluster of deletions
                                print(f'Multiple significantly deleted genes identified at {l}')
                                for g in loci[l]:
                                    if gsig[g] < 0:
                                        tt = 'low'
                                    else:
                                        tt = 'high'    

                                    print(f'Deletion of {g} at {hgnc.loc[g]["location"]} is highly significantly over-represented - p = {abs(gsig[g]):.4} d = {gD[g]:.4}\n',
                                          f'in profile matching cells with {tt} dependency on {target[0]}')
 
                                t = input('Plot affected cells ([Y],N,QUIT) : ')
                                if t != '' and t[0].upper() == 'Q':
                                    break
                                if t != '' and t[0].upper() == 'N':
                                    continue
                                    
                                if tt == 'high': # add this as a positive locus factor
                                    factors.append(f' -{l} ')
                                else :
                                    factors.append(f' (# ^  -{l}) ')
                                nfact += 1
                                
                                Dmean = np.array([gD[g] for g in loci[l]]).mean()
                                pmean = np.array([gsig[g] for g in loci[l]]).mean()

                                pdis = sp & LDEL('chr'+l)# get set of cells from match set with 2ndry locus deletion
                                gdeps,nomatch = ProfPairDist(target[0],pdis,q,deps,null_trim) # get deps for these cells
                                fig = plt.figure(figsize=(3,4))
                                width = 0.02
                                plt.violinplot((match),showmeans=True, showextrema=False)

                                plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                                data = np.array(gdeps)
                                x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                                plt.scatter(x, data, s=10,c='r')

                                plt.text(1.1,0.5,f'{l}({len(gdeps)})\np = {pmean:.4} d = {Dmean:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                                gmean = np.array(gdeps).mean()
                                plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')

                                plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                                plt.ylim(0.0,1.0)
                                plt.tick_params(axis='both', labelsize=12)
                                plt.ylabel('dependency',fontsize = '20')
                                plt.text(0.1,0.1,'locus deletion',horizontalalignment='left',style='italic',fontsize='15',color='r')

                                plt.show()

                                
                            else: # singleton deletion
                                g = loci[l][0]
                                if gsig[g] < 0:
                                    tt = 'low'
                                else:
                                    tt = 'high'    

                                print(f'Singleton deletion of {g} at {l} is highly significantly over-represented - p = {abs(gsig[g]):.4} d = {gD[g]:.4}\n',
                                      f'in profile matching cells with {tt} dependency on {target[0]}')
                                t = input('Plot affected cells ([Y],N,QUIT) : ')
                                if t != '' and t[0].upper() == 'Q':
                                    break
                                if t != '' and t[0].upper() == 'N':
                                    continue      
                                if tt == 'high': # add this as a positive gene factor
                                    factors.append(f' -{g} ')
                                else :
                                    factors.append(f' (# ^  -{g}) ')
                                nfact += 1

                                pdis = sp & DEL(g)# get set of cells from match set with 2ndry disrupted
                                gdeps,nomatch = ProfPairDist(target[0],pdis,q,deps,null_trim) # get deps for these cells
                                fig = plt.figure(figsize=(3,4))
                                width = 0.02
                                plt.violinplot((match),showmeans=True, showextrema=False)
    
                                plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                                data = np.array(gdeps)
                                x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                                plt.scatter(x, data, s=10,c='r')
    
                                plt.text(1.1,0.5,f'{g}({len(gdeps)})\np = {abs(gsig[g]):.4} d = {gD[g]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                                gmean = np.array(gdeps).mean()
                                plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')
    
                                plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                                plt.ylim(0.0,1.0)
                                plt.tick_params(axis='both', labelsize=12)
                                plt.ylabel('dependency',fontsize = '20')
                                plt.text(0.1,0.1,'deletion',horizontalalignment='left',style='italic',fontsize='15',color='r')
    
                                plt.show()
                        print('\n')
                        
                        t = input('Save the Secondary Hit List (Y,[N]) : ')
                        if t != '' and t[0].upper() == 'Y':    
                            fhit = open(f'{target[0]}-2ndryDeletionHits.csv','w')
                            fhit.write('TARGET, P\n')
                            for g in gsig:
                                fhit.write(f'{g},{gsig[g]:.5}\n')
                            fhit.close()
                            print(f'{target[0]}-2ndryDeletionHits.csv written')

# go through amplifications calculating t-test

# go through deletions calculating t-test

                    gsig = {}
                    gD = {}
                    
                    for g in gamp_count:
                        # print(f'Gene {g} is disrupted in {gcount[g]}/{len(sp)} matching cells')
                        cmut = AMP(g) & sp # cells in matched set with m gene deletion
                        if cmut == set():
                            print(f'No amplifications found for gene : {g} in matched set')
                            continue
                        cwt = sp - cmut 
                        if cwt == set():
                            print(f'No non-amplified cells found for gene : {g} in matched set')
                            continue
                        
                        # d1 = deps[target[0]].reset_index()
       
                        dmut = list(d1[d1['ModelID'].isin(cmut)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(cwt)][target[0]])
                        # print(f'{len(dmut)} mutant cell lines found')
                        with np.errstate(divide = 'ignore', invalid = 'ignore'):
                            Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                            if len(dmut) > 0 and len(dwt) > 0:
                                Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()
                            else:
                                Dd = 0

                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            gD[g] = Dd
                            if np.average(dmut) > np.average(dwt):
                                gsig[g] = Pval  # mutation of m is overrepresented in high dependency cell lines
                            else:
                                gsig[g] = -Pval  # -ve signifies m is underrepresented
                                
                    if len(gsig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        gsig = dict(sorted(gsig.items(), key = lambda x: x[1], reverse = False))
                        
# IF MULTIPLE SIGNIFICANT GENES SHARE LOCUS - MAKE LOCUS FACTOR NOT GENES          

                        loci = Locify(list(gsig)) # cluster genes into loci
                        
                        for l in loci:
                            if len(loci[l]) > 1: # cluster of amplifications
                                print(f'Multiple significantly amplified genes identified at {l}')
                                for g in loci[l]:
                                    if gsig[g] < 0:
                                        tt = 'low'
                                    else:
                                        tt = 'high'    

                                    print(f'Amplification of {g} at {hgnc.loc[g]["location"]} is highly significantly over-represented - p = {abs(gsig[g]):.4} d = {gD[g]:.4}\n',
                                          f'in profile matching cells with {tt} dependency on {target[0]}')
 
                                t = input('Plot affected cells ([Y],N,QUIT) : ')
                                if t != '' and t[0].upper() == 'Q':
                                    break
                                if t != '' and t[0].upper() == 'N':
                                    continue
                                    
                                if tt == 'high': # add this as a positive locus factor
                                    factors.append(f' -{l} ')
                                else :
                                    factors.append(f' (# ^  -{l}) ')
                                nfact += 1
                                
                                Dmean = np.array([gD[g] for g in loci[l]]).mean()
                                pmean = np.array([gsig[g] for g in loci[l]]).mean()

                                pdis = sp & LAMP('chr'+l)# get set of cells from match set with 2ndry locus amplification
                                gdeps,nomatch = ProfPairDist(target[0],pdis,q,deps,null_trim) # get deps for these cells
                                fig = plt.figure(figsize=(3,4))
                                width = 0.02
                                plt.violinplot((match),showmeans=True, showextrema=False)

                                plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                                data = np.array(gdeps)
                                x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                                plt.scatter(x, data, s=10,c='r')

                                plt.text(1.1,0.5,f'{l}({len(gdeps)})\np = {pmean:.4} d = {Dmean:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                                gmean = np.array(gdeps).mean()
                                plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')

                                plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                                plt.ylim(0.0,1.0)
                                plt.tick_params(axis='both', labelsize=12)
                                plt.ylabel('dependency',fontsize = '20')
                                plt.text(0.1,0.1,'locus amplification',horizontalalignment='left',style='italic',fontsize='15',color='r')

                                plt.show()

                                
                            else: # singleton amplification
                                g = loci[l][0]
                                if gsig[g] < 0:
                                    tt = 'low'
                                else:
                                    tt = 'high'    

                                print(f'Singleton amplification of {g} at {l} is highly significantly over-represented - p = {abs(gsig[g]):.4} d = {gD[g]:.4}\n',
                                      f'in profile matching cells with {tt} dependency on {target[0]}')
                                t = input('Plot affected cells ([Y],N,QUIT) : ')
                                if t != '' and t[0].upper() == 'Q':
                                    break
                                if t != '' and t[0].upper() == 'N':
                                    continue      
                                if tt == 'high': # add this as a positive gene factor
                                    factors.append(f' -{g} ')
                                else :
                                    factors.append(f' (# ^  -{g}) ')
                                nfact += 1

                                pdis = sp & DEL(g)# get set of cells from match set with 2ndry disrupted
                                gdeps,nomatch = ProfPairDist(target[0],pdis,q,deps,null_trim) # get deps for these cells
                                fig = plt.figure(figsize=(3,4))
                                width = 0.02
                                plt.violinplot((match),showmeans=True, showextrema=False)
    
                                plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                                data = np.array(gdeps)
                                x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                                plt.scatter(x, data, s=10,c='r')
    
                                plt.text(1.1,0.5,f'{g}({len(gdeps)})\np = {abs(gsig[g]):.4} d = {gD[g]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                                gmean = np.array(gdeps).mean()
                                plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')
    
                                plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                                plt.ylim(0.0,1.0)
                                plt.tick_params(axis='both', labelsize=12)
                                plt.ylabel('dependency',fontsize = '20')
                                plt.text(0.1,0.1,'amplification',horizontalalignment='left',style='italic',fontsize='15',color='r')
    
                                plt.show()

                        print('\n')
                        t = input('Save the Secondary Hit List (Y,[N]) : ')
                        if t != '' and t[0].upper() == 'Y':    
                            fhit = open(f'{target[0]}-2ndryAmplificationHits.csv','w')
                            fhit.write('TARGET, P\n')
                            for g in gsig:
                                fhit.write(f'{g},{gsig[g]:.5}\n')
                            fhit.close()
                            print(f'{target[0]}-2ndryAmplificationHits.csv written')


                t = input('Analyse Expression Levels (slow) (Y,[N]) : ')
                if t != '' and t[0].upper() == 'Y':
                    
                    f = input(f'Minimum fraction of cells to consider a gene effect [{frac_limit}] : ')
                    if f.strip() != '' :
                        try:
                            frac_limit = float(f)
                        except ValueError:
                            print(f'Invalid value : using {frac_limit}')
                    if frac_limit < 0.0 or frac_limit > 1.0:
                        frac_limit = 0.05
                        print(f'Invalid value : resetting {frac_limit}')
                    
                    sp = set(p)
                    sq = set(q)

                    n_limit = len(sp)*frac_limit # min no of cells for acceptance

                            
# gather all OVER-EXPRESSED genes  in profile matching cells

                    t_base = timer()

                    over_all = []
                    for c in p:
                        if c in Cdict:
                            over_all += Cdict[c].over
                    # over_uniq = set(over_all)

                    # print(f'Gathering OVER-EXPRESSED genes took : {timer()-t_base:.5}')


                    t_base = timer()
                    
                    over_count = {key:val for key, val in 
                                  Counter(over_all).items() if val > n_limit}
                    over_count =   dict(sorted(over_count.items(), key = lambda x: x[1], reverse = True))
                    print(f'{len(over_count)} genes are over-expressed in > {(frac_limit)*100}% of the matching cell lines')
                    # print(f'Counting OVER-EXPRESSED genes took : {timer()-t_base:.5}')

# gather all UNDER-EXPRESSED genes  in profile matching cells
                    under_all = []
                    for c in p:
                        if c in Cdict:
                            under_all += Cdict[c].under
                    # under_uniq = set(under_all)

                    under_count = {key:val for key, val in 
                                  Counter(under_all).items() if val > n_limit}
                    under_count =   dict(sorted(under_count.items(), key = lambda x: x[1], reverse = True))
                    print(f'{len(under_count)} genes are under-expressed in > {(frac_limit)*100}% of the matching cell lines')

                    path_thresh = 0.1

# find pathways that are highly represented amongst the over/under-expressed genes
                    op = Pathify(over_count)
                    opaths = [x for x in op if op[x]/len(over_count) > path_thresh]
                    if len(opaths) > 0:
                        print(f'{len(opaths)} pathways are enriched in > {(path_thresh)*100}% of the over-expressed genes')
                    else:
                        print(f'No pathways are enriched in > {(path_thresh)*100}% of the over-expressed genes')

                    up = Pathify(under_count)
                    upaths = [x for x in op if op[x]/len(under_count) > path_thresh]
                    if len(upaths) > 0:
                        print(f'{len(upaths)} pathways are represented in > {(path_thresh)*100}% of the under-expressed genes')
                    else:
                        print(f'No pathways are enriched in > {(path_thresh)*100}% of the under-expressed genes')

# recreate gene list for over-expressed genes in enriched pathways

                    over_enrich = set()
                    for i in opaths:
                        over_enrich |= (set(paths[i]) & set(over_count))
                    
                    if len(over_enrich) > 0:
                        print(f'Testing significance for {len(over_enrich)} pathway-enriched over-expressed genes')
# go through over-expressions calculating t-test

                    gsig = {}
                    gD = {}
                     
                    for g in over_enrich:
                        cmut = OX(g,-1) & sp # cells in matched set with precalculated over-expression
                        if cmut == set():
                            print(f'No overexpressing cells found for gene : {g} in matched set')
                            continue
                        cwt = sp - cmut 
                        if cwt == set():
                            print(f'No non-overexpressing cells found for gene : {g} in matched set')
                            continue
                        
                        # d1 = deps[target[0]].reset_index()
       
                        dmut = list(d1[d1['ModelID'].isin(cmut)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(cwt)][target[0]])
                        # print(f'{len(dmut)} mutant cell lines found')
                        Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                        with np.errstate(divide = 'ignore'):
                            Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()

                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            gD[g] = Dd
                            if np.average(dmut) > np.average(dwt):
                                gsig[g] = Pval  # overexpression of g is overrepresented in high dependency cell lines
                            else:
                                gsig[g] = -Pval  # -ve signifies g is underrepresented
                                
                    if len(gsig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        gsig = dict(sorted(gsig.items(), key = lambda x: x[1], reverse = False))
                        for g in gsig:
                            if gsig[g] < 0:
                                tt = 'low'
                            else:
                                tt = 'high'
                            print(f'Overexpression of {g} is highly significantly over-represented - p = {abs(gsig[g])} d = {gD[g]}\n',
                                  f'in profile matching cells with {tt} dependency on {target[0]}')
                            t = input('Plot affected cells ([Y],N,QUIT) : ')
                            
                            if t != '' and t[0].upper() == 'Q':
                                break
                            if t != '' and t[0].upper() == 'N':
                                continue
                                                        
                            if tt == 'high': # add this as a positive factor
                                factors.append(f' >{g} ')
                            else :
                                factors.append(f' (# ^  >{g}) ')
                            nfact += 1

                            pdis = sp & OX(g,-1)# get set of cells from match set with 2ndry over-expressed
                            gdeps,nomatch = ProfPairDist(target[0],pdis,q,deps,null_trim) # get deps for these cells
                            fig = plt.figure(figsize=(3,4))
                            width = 0.02
                            plt.violinplot((match),showmeans=True, showextrema=False)

                            plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                            data = np.array(gdeps)
                            x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                            plt.scatter(x, data, s=10,c='r')

                            plt.text(1.1,0.5,f'{g}({len(gdeps)})\np = {abs(gsig[g]):.4} d = {gD[g]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                            gmean = np.array(gdeps).mean()
                            plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')

                            plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                            plt.ylim(0.0,1.0)
                            plt.tick_params(axis='both', labelsize=12)
                            plt.ylabel('dependency',fontsize = '20')
                            plt.text(0.1,0.1,'over-expression',horizontalalignment='left',style='italic',fontsize='15',color='r')

                            plt.show()
                        print('\n')
                        t = input('Save the Secondary Hit List (Y,[N]) : ')
                        if t != '' and t[0].upper() == 'Y':    
                            fhit = open(f'{target[0]}-2ndryOverexpressionHits.csv','w')
                            fhit.write('TARGET, P\n')
                            for g in gsig:
                                fhit.write(f'{g},{gsig[g]:.5}\n')
                            fhit.close()
                            print(f'{target[0]}-2ndryOverexpressionHits.csv written')

# recreate gene list for over-expressed genes in enriched pathways

                    under_enrich = set()
                    for i in opaths:
                        under_enrich |= (set(paths[i]) & set(under_count))
                    if len(under_enrich) > 0:
                        print(f'Testing significance for {len(under_enrich)} pathway-enriched under-expressed genes')
# go through under-expressions calculating t-test

                    gsig = {}
                    gD = {}
                    
                    for g in under_enrich:
                        cmut = UX(g,-1) & sp # cells in matched set with precalculated under-expression
                        if cmut == set():
                            print(f'No underexpressing cells found for gene : {g} in matched set')
                            continue
                        cwt = sp - cmut 
                        if cwt == set():
                
                            print(f'No non-underexpressing cells found for gene : {g} in matched set')
                            continue
                        
                        # d1 = deps[target[0]].reset_index()
       
                        dmut = list(d1[d1['ModelID'].isin(cmut)][target[0]])
                        dwt = list(d1[d1['ModelID'].isin(cwt)][target[0]])
                        # print(f'{len(dmut)} mutant cell lines found')
                        Tval,Pval = stats.ttest_ind(a=dmut, b=dwt, equal_var=False)
                        with np.errstate(divide = 'ignore'):
                            Dd = (np.array(dmut).mean()-np.array(dwt).mean())/np.array(dmut+dwt).std()

                        # print(f'p = {Pval}')
                        if Pval < 0.05 and abs(Dd) >= cohen_min:
                            gD[g] = Dd
                            if np.average(dmut) > np.average(dwt):
                                gsig[g] = Pval  # overexpression of g is overrepresented in high dependency cell lines
                            else:
                                gsig[g] = -Pval  # -ve signifies g is underrepresented
                                
                    if len(gsig) == 0:
                        print('No significant secondary effects were identified within the matching set')
                    else:
                        gsig = dict(sorted(gsig.items(), key = lambda x: x[1], reverse = False))
                        for g in gsig:
                            if gsig[g] < 0:
                                tt = 'low'
                            else:
                                tt = 'high'
                            print(f'Underexpression of {g} is highly significantly over-represented - p = {abs(gsig[g])} d = {gD[g]}\n',
                                  f'in profile matching cells with {tt} dependency on {target[0]}')

                            t = input('Plot affected cells ([Y],N,QUIT) : ')
                            if t != '' and t[0].upper() == 'Q':
                                break
                            if t != '' and t[0].upper() == 'N':
                                continue
                                                        
                            if tt == 'high': # add this as a positive factor
                                factors.append(f' <{g} ')
                            else :
                                factors.append(f' (# ^  <{g}) ')
                            nfact += 1

                            pdis = sp & UX(g,-1)# get set of cells from match set with 2ndry under-expression
                            gdeps,nomatch = ProfPairDist(target[0],pdis,q,deps,null_trim) # get deps for these cells

                            fig = plt.figure(figsize=(3,4))
                            width = 0.02
                            plt.violinplot((match),showmeans=True, showextrema=False)

                            plt.xticks([1],labels = [f'matched ({len(set(p))})'],fontsize = '20')
                            data = np.array(gdeps)
                            x = np.ones(data.shape[0]) + (np.random.rand(data.shape[0])*width-width/2.)
                            plt.scatter(x, data, s=10,c='r')

                            plt.text(1.1,0.5,f'{g}({len(gdeps)})\np = {abs(gsig[g]):.4} d = {gD[g]:.4g}',horizontalalignment='left',style='italic',fontsize='15',color='r')
                            gmean = np.array(gdeps).mean()
                            plt.plot([1-2*width,1+2*width],[gmean,gmean],color='r')

                            plt.title(f"{target[0]} : {target[1]} p = {pAB:.4} d = {D:.4g}",fontsize = '20',pad=20)
                            plt.ylim(0.0,1.0)
                            plt.tick_params(axis='both', labelsize=12)
                            plt.ylabel('dependency',fontsize = '20')
                            plt.text(0.1,0.1,'under-expression',horizontalalignment='left',style='italic',fontsize='15',color='r')

                            plt.show()
                            
                        print('\n')
                        t = input('Save the Secondary Hit List (Y,[N]) : ')
                        if t != '' and t[0].upper() == 'Y':    
                            fhit = open(f'{target[0]}-2ndryUnderexpressionHits.csv','w')
                            fhit.write('TARGET, P\n')
                            for g in gsig:
                                fhit.write(f'{g},{gsig[g]:.5}\n')
                            fhit.close()
                            print(f'{target[0]}-2ndryUnderexpressionHits.csv written')
                
                if nfact < 1:
                    continue
                ftemp = factors.copy()
                factors = []
                # fkeep = np.ones(len(nfact),dtype=bool)
                print(f'\n{nfact} individual secondary factors have been selected : ')
                for f in ftemp:
                    ans = input(f'{f}\tkeep ? ([Y],N) : ').upper()
                    if len(ans) == 0 or ans[0] != 'N':
                        factors.append(f)
                    else:
                        print(f'Factor {f} eliminated')
                
                        
                
                t = input('Identify Pareto Optimal combinations ? (Y,[N]) : ').upper()
                if len(t) > 0 and t[0] != 'N':
                    
# make base sets from factors 


                    pbase = []
                    # qbase = []
                    
                    for f in factors:
                        pbase.append(eval(prof_parse(f)))
                        # qbase.append(eval(prof_nullhype(prof_flatten(f))))

# generate combinations 
                    pfs = []
                    psets = []
                    # qsets = []
                    pfs.append(froot) # base profiles on original
                    psets.append(eval(prof_parse(froot)))
                    # qsets.append(eval(prof_nullhype(prof_flatten(froot))))
                    
                    for j in range(len(factors)): # add factors in turn
                        q = []
                        r = []
                        s = []
                        for i in range(len(pfs)):
                            q.append(pfs[i]+' & '+factors[j])
                            r.append(psets[i] & pbase[j])
                            # s.append(qsets[i] & qbase[j])
                        pfs = pfs + q
                        psets = psets + r
                        # qsets = qsets + s
                    print(f'{len(pfs)} profiles generated')
                    print(f'{len(psets)} p-sets generated')

                    print('Calculating scores ...')
                    us = UNI()
                    score = []
    
                    t_base = timer()

                    with mp.Pool(mp.cpu_count()) as pool:
                        score = pool.starmap_async(DandCount,
                            [(i,us-i,target[0],deps) for i in psets]).get()

                    print(f'Calculating scores took : {timer()-t_base:.5}')
                    
                    fig, ax = plt.subplots()

                    ax.set_title(f"{target[0]} : {target[1]} profile refinement",fontsize=14)
                    # ax1.title.set_text(f' r = {cc1:.4f}')
                    ax.tick_params(axis='both', labelsize=10)
                    ax.set_ylabel('cell line count',rotation='vertical',fontsize=14)
                    ax.set_xlabel('effect size',rotation='horizontal',fontsize=14)
# purge 0 scores
                    
                    xx = [i[0] for i in score if i[0] != 0 and i[1] != 0]
                    yy = [i[1] for i in score if i[0] != 0 and i[1] != 0]
                    ax.scatter(xx,yy,s=50)


# find Pareto optima

                    dom = np.ones(len(score),dtype = bool)
                    
                    for i in range(len(score)):
                        if score[i][1] < 1:
                            dom[i] = False  
                            continue
                        for j in range(len(score)):
                            if (score[j][0] <= score[i][0] and \
                                score[j][1] <= score[i][1]) and \
                                (score[j][0] < score[i][0] or \
                                 score[j][1] < score[i][1]):
                                    dom[j] = False
                                    
                    xx = [score[i][0] for i in range(len(score)) if dom[i]]
                    yy = [score[i][1] for i in range(len(score)) if dom[i]]
                    dompfs = [pfs[i] for i in range(len(score)) if dom[i]]
                    ax.scatter(xx,yy,s=15,color='red')
                    
                    xydom = sorted(list(zip(xx,yy,dompfs)), key = lambda x:x[0])
# look for duplicate scores and eliminate longer profile

                    dom = np.ones(len(xydom),dtype = bool)
                    for i in range(len(xydom)):
                        for j in range(len(xydom)):
                            if  xydom[i][0] == xydom[j][0] and \
                                xydom[i][1] == xydom[j][1] and \
                                len(xydom[i][2]) > len(xydom[j][2]):
                                    dom[j] = False
                    xydom_filt = [xydom[i] for i in range(len(xydom)) if dom[i]]
                    
                    for j,i in enumerate(xydom_filt):
                        print(f'{j} {i[0]:.4} {i[1]}\n{i[2]}')
                        ax.text(i[0],i[1],f' {j}\n',fontsize='8',color='k')
                        
                    
                    xx,yy,dompfs = zip(*xydom_filt)
                    ax.plot(xx, yy, color = 'red',linestyle='-')
                    
                    ax.text(max(xx),max(yy)*0.9,'Pareto front',horizontalalignment='right',style='italic',fontsize='15',color='r')
                    
                    plt.show()
                                        
                    ans = input('Save refined profiles to a new profile collection ? ([Y],N) : ').upper()
                    if len(ans) == 0 or ans[0] != 'N':
                        fname = target[0]+'_'+target[1]+'_refined.profiles'
                        f = open(fname,'w')
                        for i in range(len(dompfs)):
                            s = input(f"Save profile {i} ? ([Y],N) : ").upper()
                            if len(s) == 0 or s[0] != 'N':
                                f.write(f'{target[0]}_{target[1]}_ref_{i}\n')
                                f.write(prof_flatten(dompfs[i]).replace(" ",'')+'\n')
                        f.close()    

# 
#                     
# =============================================================================
# =============================================================================
# p values for a list of pairs, both ways round
# 
# ██╗     ██╗███████╗████████╗
# ██║     ██║██╔════╝╚══██╔══╝
# ██║     ██║███████╗   ██║   
# ██║     ██║╚════██║   ██║   
# ███████╗██║███████║   ██║   
# ╚══════╝╚═╝╚══════╝   ╚═╝   
# =============================================================================
            
        elif task.lower() == 'list' :
            
            default_file = 'query_list.csv'
            
            f = input(f'Gene AB pair file [{default_file}] : ')
            if f.strip() == '' :
                f = default_file
            print(f'Loading {f}')
            
            pairs = pd.read_csv(f)
            print(pairs)


            GeneA = 'GeneA'
            GeneB = 'GeneB'
            other = ''
            while True:
                ga = input(f'Column name for target gene [{GeneA}] ? : ')
                if ga  != '':
                    GeneA = ga
                gb = input(f'Column name for mutated gene [{GeneB}] ? : ')
                if gb  != '':
                    GeneB = gb
                
                if GeneA in pairs.columns and GeneB in pairs.columns:
                    break
                print(f'One or both of {GeneA} and {GeneB} are invalid column names')

            while True:
                other = input('Column name for ancillary data (<return> for none) : ')
                if other == '': # no ancillary data
                    print('No ancillary data used')
                    break
                else:
                    if other in pairs.columns:
                        break
                    print(f'{other} is not a valid column name')
                    
            while True:
                tiss = input('Tissue to restrict to (<return> for none) : ').upper()
                if tiss == '':
                    print('No tissue restriction applied')
                    tmask = UNI()
                    break
                else:
                    c = match_keyword(tiss,set(info['OncotreeLineage']))
                    if c == '':
                        print(f'cancer type specification {w} not valid')
                        continue
                    elif c == '*':
                        print(f'cancer type specification {w} ambiguous')
                        continue
                    tmask = TYP(c)
                    print(f'restricting calculations to cell lines deriving from {c}')
                    break
                        
                                           
            f_name = f.split('.')[0]+'_p_mut_cnv_exp.'+f.split('.')[1]
            
            fout = open(f_name,'w')
            if other == '':
                fout.write('Target,Affected,pmut,pcnv,pexp\n')
            else:
                fout.write(f'Target,Affected,{other},pmut,pcnv,pexp\n')
                
            for i in range(len(pairs)):
                gA = pairs[GeneA].iloc[i]
                gB = pairs[GeneB].iloc[i]
                if gA not in deps.columns:
                    print(f'No dependency data for {gA} - skipping')
                    continue
                if gB not in deps.columns:
                    print(f'No dependency data for {gB} - skipping')
                    continue

                p_mut_1 = PairPvalue(gA,gB,Cdict,deps,thresh,True,False,False,tmask)
                p_cnv_1 = PairPvalue(gA,gB,Cdict,deps,thresh,False,True,False,tmask)
                p_exp_1 = PairPvalue(gA,gB,Cdict,deps,thresh,False,False,True,tmask)
                               
                p_mut_2 = PairPvalue(gB,gA,Cdict,deps,thresh,True,False,False,tmask)
                p_cnv_2 = PairPvalue(gB,gA,Cdict,deps,thresh,False,True,False,tmask)
                p_exp_2 = PairPvalue(gB,gA,Cdict,deps,thresh,False,False,True,tmask)

                pmut = min(p_mut_1,p_mut_2)
                pcnv = min(p_cnv_1,p_cnv_2)
                pexp = min(p_exp_1,p_exp_2)

                fout.write(f'{gA},{gB},{pairs[other].iloc[i]},{p_mut_1},{p_cnv_1},{p_exp_1}\n')
                fout.write(f'{gB},{gA},{pairs[other].iloc[i]},{p_mut_2},{p_cnv_2},{p_exp_2}\n')

                fout.flush()
            
            fout.close()
            print(f'{f_name} written\n\n')            
            
        elif task.lower() == 'oldlist' :
            
            default_file = 'query_list.csv'
            
            f = input(f'Gene AB pair file [{default_file}] : ')
            if f.strip() == '' :
                f = default_file
            print(f'Loading {f}')
            
            pairs = pd.read_csv(f)
            print(pairs)


            GeneA = 'GeneA'
            GeneB = 'GeneB'
            other = ''
            while True:
                ga = input(f'Column name for target gene [{GeneA}] ? : ')
                if ga  != '':
                    GeneA = ga
                gb = input(f'Column name for mutated gene [{GeneB}] ? : ')
                if gb  != '':
                    GeneB = gb
                
                if GeneA in pairs.columns and GeneB in pairs.columns:
                    break
                print(f'One or both of {GeneA} and {GeneB} are invalid column names')
            
            usemuts = usecnv = useexp = False
            while True:
                s = []
                w = input('Use Mutations (M), Copy Number Deletion (D), Underexpression (U) :').lower()
                if len(w) == 0:
                    continue
                if 'm' in w:
                    usemuts = True
                    s.append('mutations')
                if 'd' in w:
                    usecnv = True
                    s.append('deletions')
                if 'u' in w:
                    useexp = True
                    s.append('under-expression')
                    
                if usemuts or usecnv or useexp:
                    print(f'Using {Flistand(s)} ')
                    break
                else:
                    print(f'No valid choices in input : {w}')
            
            bb = input('Calculate p both ways round ? (N,[Y]) ? :').lower()
            
            if len(bb) > 0 and bb[0] == 'n':
                both = False
            else:
                both = True
            
            
            f_name = f.split('.')[0]+'_p.'+f.split('.')[1]
            
            fout = open(f_name,'w')
            fout.write('GeneA,GeneB,p\n')
                
            for i in range(len(pairs)):
                gA = pairs[GeneA].iloc[i]
                gB = pairs[GeneB].iloc[i]
                if gA not in deps.columns:
                    print(f'No dependency data for {gA} - skipping')
                    continue
                p1 = PairPvalue(gA,gB,Cdict,deps,thresh,usemuts,usecnv,useexp)    
                fout.write(f'{gA},{gB},{p}\n')
                if both:
                    if gB not in deps.columns:
                        print(f'No dependency data for {gB} - skipping')
                        continue
                    
                    p = PairPvalue(gB,gA,Cdict,deps,thresh,usemuts,usecnv,useexp)    
                    fout.write(f'{gB},{gA},{p}\n')
            
            fout.close()
# =============================================================================
# search of GeneA for hits with given parameters            
# 
# ███████╗██╗███╗   ██╗██████╗ 
# ██╔════╝██║████╗  ██║██╔══██╗
# █████╗  ██║██╔██╗ ██║██║  ██║
# ██╔══╝  ██║██║╚██╗██║██║  ██║
# ██║     ██║██║ ╚████║██████╔╝
# ╚═╝     ╚═╝╚═╝  ╚═══╝╚═════╝ 
# =============================================================================


        elif task.lower() == 'find' :

        # get targets
            
            while True:
                gene_found = True
                
                t = input('Target gene(s)A (QUIT): ')
                
                if t.lower() == 'quit' :
                    break
        # break target names into list
        
                target = t.upper().replace(',',' ').split()
        # loop through targets and check that target name is in list
                for t in target:
                    
                    if t not in list(deps.columns) :
                        gene_found = False
                        print(f'Target gene : {t} is not found')
                        continue
                if not gene_found:
                    continue
        # now find genes
        
                usemuts = usecnv = useexp = False
                while True:
                    s = []
                    w = input('Use Mutations (M), Copy Number Deletion (D), Underexpression (U) :').lower()
                    if len(w) == 0:
                        continue
                    if 'm' in w:
                        usemuts = True
                        s.append('mutations')
                    if 'd' in w:
                        usecnv = True
                        s.append('deletions')
                    if 'u' in w:
                        useexp = True
                        s.append('under-expression')
                    if usemuts or usecnv or useexp:
                        print(f'Using {Flistand(s)} ')
                        break
                    else:
                        print(f'No valid choices in input : {w}')

                
                ng = 1
                results = []
                t_base = timer()
                if len(target) == 1:
                    print(f'Finding hits for {target[0]}')
                elif len(target) == 2:
                    print(f'Finding hits for either {Flistor(target)}')
                else:
                    print(f'Finding hits for any of {Flistor(target)}')
            # set up parallelism
            
                pool = mp.Pool(mp.cpu_count())
                glist = list(set(deps.columns))
                # glist = list(set(deps.columns) | set(cnvs.columns[1:]))

                
                results = pool.starmap_async(GeneHits,
                        [(g,target,Cdict,deps,thresh,frac_limit,p_limit,usemuts,usecnv,useexp) for g in glist]).get()
                pool.close()
            
            # =============================================================================
            #     
            #     for g in Gdict:
            #         results[g] = GeneHits(g)
            #         
            #         if ng %1000 == 0 :
            #             print(f'{ng} genes time : {time.process_time()-t_base}')
            #         ng += 1
            #     
            #     
            # =============================================================================
                print(f'Gene hits total time : {timer()-t_base:.5}')
                
                HitList = {} #initialise hitlist
        
                for r in results:
                    if r :
                        HitList[r[0]] = Hit(r[1],r[2],r[3],r[4])    
                # sort on p values
                
                HitList = dict(sorted(HitList.items(), key = lambda x: x[1].p, reverse = False))
                    
                # make boxplots for HitList targets        
                
                if len(HitList) > 0 :
                #    HitList = sorted(HitList)   
                    print(f'{len(HitList)} hits found for target {Flistor(target)} with p<={p_limit} and mutated in > {100*frac_limit}% of tested cell lines')     
                    for h in HitList:
                        print(f'{h:<8} : p = {HitList[h].p:5} : occurence = {HitList[h].frac:5}')
                #    print(str(HitList).replace("'","").replace(" ","").strip('[]'))
                    
                a = input('Do you want to plot the Hit List ? ([Y],N) : ').lower()
                if  len(a) == 0 or (len(a) > 0 and a[0] != 'n'):    
                    for i in HitList:
                        h = HitList[i]
                        D = (np.array(h.mut_deps).mean()-np.array(h.notmut_deps).mean())/np.array(h.mut_deps+h.notmut_deps).std()

                        # print(f'TEST {i} p = {h.p:5} f = {h.frac:5} Nmut = {len(h.mut_deps)},{len(h.notmut_deps)}')
                        t1 = Flistor(target)
                        t2 = f'{i} f = {h.frac:.2f}'
                        DepPlot2(h.mut_deps,h.notmut_deps,t1,t2,h.p,D)                       
                        # plt.rcParams.update({'font.size': 30})
                        # plt.violinplot((h.mut_deps,h.notmut_deps),showmeans=True, showextrema=False)
                        # plt.xticks([1,2],labels = ['mutated','not mutated'])
                        # plt.title(f'{Flistor(target)} : GeneB - {i} p = {h.p:.3} f = {h.frac:.3}')
                        # plt.ylabel('dependency')
                        # plt.show()
                        
                a = input('Do you want to save the Hit List ? ([Y],N) : ').lower()
                if  len(a) == 0 or (len(a) > 0 and a[0] != 'n'):    
                    fhit = open(f'{"_".join(target)}-DepHits.csv','w')
                    fhit.write('GeneA,GeneB,p,f\n')
                    for i in HitList:
                        h = HitList[i]
                        fhit.write(f'{target[0]},{i},{h.p:.5},{h.frac:.5}\n')
                    fhit.close()
                    print(f'{"_".join(target)}-DepHits.csv written')

# =============================================================================
# 
# MUTS  - map mutations of a GeneB onto its sequences against dependencies of a GeneA
#
# ███╗   ███╗██╗   ██╗████████╗███████╗
# ████╗ ████║██║   ██║╚══██╔══╝██╔════╝
# ██╔████╔██║██║   ██║   ██║   ███████╗
# ██║╚██╔╝██║██║   ██║   ██║   ╚════██║
# ██║ ╚═╝ ██║╚██████╔╝   ██║   ███████║
# ╚═╝     ╚═╝ ╚═════╝    ╚═╝   ╚══════╝
#                              
# =============================================================================

        elif task.lower() == 'muts' : 
# need GeneA and GeneB from user

            LOF_Muts={'NONSENSE':'r*','NONSTOP':'r!','FRAME_SHIFT_INS':'bv',
                      'START_CODON_INS':'gx','FRAME_SHIFT_DEL':'b^'}

            while True:
                gene_found = True
 
                t = input('Target GeneA and mutated GeneB (QUIT): ')
                 
                if t.lower() == 'quit' :
                    break
         # break target names into list
         
                target = t.upper().replace(',',' ').split()
                
         # loop through targets and check that target name is in list
                for t in target:
                     
                    if t not in list(deps.columns) :
                        gene_found = False
                        print(f'Gene : {t} is not found')
                        continue
                if not gene_found:
                     continue
                if len(target) != 2 :
                    print('Two valid gene names must be specified')
                    continue
        # now find genes
                        
                print(f'Analysing dependency on GeneA {target[0]}')
                print(f' in cells with mutations in GeneB {target[1]}')
                geneA = target[0]
                geneB = target[1]

# for all cell lines with mutations in GeneB, get the dependency value for GeneA

                if (geneB in Gdict) and len(Gdict[geneB].cells) > 0:
                       
                    print(f'GeneB {geneB} has disruptive mutations in {len(Gdict[geneB].cells)} cell lines')
                    
# =============================================================================
# load source mutation data if not already loaded
# =============================================================================
                    if len(muts) == 0:
                        default_mut_dat_file = 'OmicsSomaticMutations_GoFLoFUpdated.csv'
                        f = input(f'Cell line mutation file [{default_mut_dat_file}] : ')
                        if f.strip() == '' :
                            f = default_mut_dat_file
                        print(f'Loading {f}')
                        
                        muts = pd.read_csv(f)      
                        
                        default_mut_dat_file = f
                    
# =============================================================================
# get data for all  mutations of this Gene
# =============================================================================
                    mPosGeneB = []
                    for i in muts.index:
                        if muts.iloc[i]['HugoSymbol'] == geneB:
                            if muts.iloc[i]['VariantInfo'] in LOF_Muts: # a loss of function mutation
                                p = muts.iloc[i]['ProteinChange']
                                c = muts.iloc[i]['ModelID']
                                if c not in deps.index:  # cant do anything if mutation is in cell line lacking dependency data
                                    continue
                                t = muts.iloc[i]['Transcript']
                                if p != '' and p !='NA' and not pd.isna(p):  # catch poorly annotated splice disruptions
                                    pos = int(re.findall(r'\d+',p)[0])
                                    d = deps.loc[c][geneA]
                                    v = muts.iloc[i]['VariantInfo']
                                    mPosGeneB.append(Mut(geneB,c,t,p,v,pos,d))
                                    print(f'{geneB} {v} at {pos} with {geneA} dependency {d} in tumour {c}')
                                    
                    print(f'{geneB} has {len(mPosGeneB)} mutations in cells with dependency data for {geneA}')
                    
                    trans = most_frequent([j.trans for j in mPosGeneB])
                    print(f'using transcript {trans} as reference')
                    
                    server = "https://grch37.rest.ensembl.org"
                    ext = f"/sequence/id/{trans}?multiple_sequences=1;type=protein"
                     
                    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                     
                    if not r.ok:
                      r.raise_for_status()
                      sys.exit()
                     
                    d = json.loads(r.text[1:-1])
                    print(len(d['seq']))

# loop through the different types of LOF mutations and plot differently

                    fig,ax = plt.subplots()                    
                    leg = []
                    for m in LOF_Muts:
                        px = []
                        py = []
                        for i in mPosGeneB:
                            if i.type == m:
                                px.append(i.pos)
                                py.append(i.dep)
                        
                        l, = ax.plot(px,py,LOF_Muts[m])
                        leg.append(l)

                    ax.plot([0,len(d['seq'])],[-0.1,-0.1],'k-',lw=10)
                    ax.text(len(d['seq'])/2,-0.1,geneB,fontsize=10,horizontalalignment='center',
                            verticalalignment='center',color='w', weight='bold')
                    plt.ylim(ymax=1.0)
                    plt.title(f'{geneA} dependency v. position of {geneB} disrupting mutation')
                    plt.ylabel('dependency')
                    plt.xlabel(f'{geneB} codon number')
                    plt.legend(leg,list(LOF_Muts),loc='center left',bbox_to_anchor=(1,0.5))
                    plt.show()                        
                        
                            
                else:
                    print(f'GeneB {geneB} is not mutated in any cell line - analysis cannot be performed')

# =============================================================================
# show tissue type distributions of cells matching defined genetic profile
# 
# ████████╗██╗   ██╗██████╗ ███████╗
# ╚══██╔══╝╚██╗ ██╔╝██╔══██╗██╔════╝
#    ██║    ╚████╔╝ ██████╔╝█████╗  
#    ██║     ╚██╔╝  ██╔═══╝ ██╔══╝  
#    ██║      ██║   ██║     ███████╗
#    ╚═╝      ╚═╝   ╚═╝     ╚══════╝
# =============================================================================

        
        elif task.lower() == 'type' : # do the type analysis
# need genetic profile 
                

            first_call = True

# check that profiles are loaded
            while not prof_loaded:
                profiles = {}
                list_profile_files()
                default_prof_file = 'AB_Dep_Profiles'
                f = input(f'Profile file [{default_prof_file}] : ')
                if f.strip() == '' :
                    f = default_prof_file
                print(f'Loading {f}')
                if prof_load(f):
                    print(f'{f} loaded')
                    prof_loaded = True
                    break
                else:
                    continue
            if first_call and len(profiles) > 0:
                print('Available predefined genetic profiles :')
                for n in sorted(profiles):
                    print(n)
                    first_call = False

# get genetic profile to use

            while True:
                print('\nAvailable predefined genetic profiles :')
                for gp in sorted(profiles):
                    print(gp)

                gp = input('\nDefined genetic profile to use (QUIT): ')
                if len(gp) == 0:
                    continue
                if gp.lower() == 'quit' :
                    break
# check that profile exists
                m = match_keyword(gp.upper(),profiles)
                if  m == '*':
                    print(f'{gp.upper()} is ambiguous')
                    continue
                if m not in profiles:
                    print(f'Genetic profile : {gp} is not defined\nDo this in LINE\n')
                    continue

                prof = profiles[m]
# generate hypothesis set
                p = list(eval(prof_parse(prof)))
                
# generate normalisation set

                pnorm = list(eval(prof_parse('#')))
                
                cancer = {}
                nincl = 0
                for i in p:
                    c = info.loc[info['ModelID'] == i]
                    if len(c) > 0:
                        cname = c['OncotreeLineage'].iloc[0]
                        nincl += 1
                    else:
                        continue
                    if cname in cancer:
                        cancer[cname] += 1
                    else:
                        cancer[cname] = 1
                print(f'{nincl} cell lines included')

                cancer = dict(sorted(cancer.items(), key = lambda x: x[1], reverse = True))

# repeat for all cells
                cnorm = {}
                nincl = 0
                for i in pnorm:
                    c = info.loc[info['ModelID'] == i]
                    if len(c) > 0:
                        cname = c['OncotreeLineage'].iloc[0]
                        nincl += 1
                    else:
                        continue
                    if cname in cnorm:
                        cnorm[cname] += 1
                    else:
                        cnorm[cname] = 1
                print(f'{nincl} cell lines included in normalisation')

                cancer = dict(sorted(cancer.items(), key = lambda x: x[1], reverse = True))


# plot bar chart
                                

                fig, ax = plt.subplots()
                c_name = [i.split()[0].lower() for i in cancer]
                c_count = [cancer[i] for i in cancer]
                ax.bar(c_name, c_count,alpha=0.3)

                ax.set_ylabel('number of cell lines  ',rotation='vertical',fontsize=16)
                ax.set_title(f'Cancer tissue distribution for profile : {m}',fontsize=16)
                ax.set_xticks(range(len(c_name)))
                ax.set_xticklabels(c_name,rotation=90)
                ax.tick_params(axis='both', labelsize=10)
                plt.show()

                fig, ax = plt.subplots()
                c_name = [i.split()[0].lower() for i in cancer]
                c_count = [cancer[i]/cnorm[i] for i in cancer]
                ax.bar(c_name, c_count,alpha=0.3)

                ax.set_ylabel('Ratio profile / all ',rotation='vertical',fontsize=16)
                ax.set_title(f'Normalised tissue distribution for profile : {m}',fontsize=16)
                ax.set_xticks(range(len(c_name)))
                ax.set_xticklabels(c_name,rotation=90)
                ax.tick_params(axis='both', labelsize=10)
                plt.show()



# =============================================================================
# ███████╗██╗  ██╗██████╗ ██████╗ 
# ██╔════╝╚██╗██╔╝██╔══██╗██╔══██╗
# █████╗   ╚███╔╝ ██████╔╝██████╔╝
# ██╔══╝   ██╔██╗ ██╔═══╝ ██╔══██╗
# ███████╗██╔╝ ██╗██║     ██║  ██║
# ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝
# Plot expression and CNV variation                
# =============================================================================

        elif task.lower() == 'expr' : # do the expression analysis 

            while True:
                gene_found = True
                
                t = input('Gene to analyse expression of : (QUIT): ').strip().upper()
                
                if len(t) == 0:
                    continue
                if t[0] == 'Q' :
                    break
                if t[0] == '*': # do full analysis of expression dataset
                
# scatter plot of expression for all genes v. dependency
                    continue
                
                if t not in expr.columns[1:] :
                    print(f'{t} not valid gene name')
                    continue
                exp = expr[t]
                mexp = exp.mean()
                emod = exp.mode()[0]
                emed = exp.median()
                emin = exp.min()
                emax = exp.max()
                esd = exp.std()
                
# make distribution

                v = np.zeros(20)
                bn = (exp.max()-emin)/20
                for e in exp:
                    i = min(trunc((e-emin)/bn),19)
                    v[i] += 1
                    

                x = np.array(exp)
                kde = stats.gaussian_kde(x)
                ksamp = 100
                xx = np.linspace(x.min(), x.max(), ksamp)
                kdif = np.diff(xx)[0]

                fig, ax = plt.subplots()

                ax.tick_params(axis='both', labelsize=10)
                
                vals = ax.hist(x, bins=20, alpha=0.3)
                vmax = vals[0].max()
                tol = vmax * 0.1
                
                plt.ylim(0.0,vmax*1.1)

                wid = np.diff(vals[1])[0]
                ksc = wid*len(x)
                kk = kde(xx)*ksc
                
                ax.plot(xx,kk)
                
                peaks,_ = find_peaks(kk,height=tol)
                # peaks = argrelextrema(kk, np.greater,order=2)[0]
                troughs = argrelextrema(kk, np.less,order=2)[0]
                print(peaks,troughs)

                split_pop = False

                if abs(mexp-emed)>wid/2:
                    print(f'mean: {mexp:.2f}, median: {emed:.2f} - not normal')
                if len(peaks) == 0:
                    print('No peaks found')
                if len(peaks) > 1:
                    print('multiple maxima detected')
                if len(peaks) > 0 and (peaks[0] * kdif) <= bn:   # in zeroth bin
                    print('with silent component')
                if len(peaks) == 2 and len(troughs) > 0:
                    split_pop = True
                    tro = troughs[0]*kdif+emin
                    print(f'expression has two major populations - split at expression value of {tro:.4}')
                    ax.plot([tro,tro],[0,vmax],color='black',linewidth='1',linestyle='dashed')    

                for p in peaks:
                    ax.text(p*kdif+emin,kk[p],'+',fontsize=12,horizontalalignment='center')

                ax.plot([mexp,mexp],[0,vmax],color='magenta',linewidth='2')
                ax.text(mexp,vmax+vmax/50,'μ',fontsize=12,horizontalalignment='center')
                ax.plot([emed,emed],[0,vmax],color='blue',linewidth='2')
                ax.text(emed,vmax+vmax/50,'x̃',fontsize=12,horizontalalignment='center')
                # ax.plot([mmod,mmod],[0,vmax],color='orange',linewidth='2')
                ax.plot([mexp+3*esd,mexp+3*esd],[0,vmax/8],color='magenta',linewidth='3')
                ax.plot([mexp-3*esd,mexp-3*esd],[0,vmax/8],color='magenta',linewidth='3')
                ax.text(mexp+3*esd,vmax/8+vmax/50,'3',fontsize=10,horizontalalignment='center')
                ax.text(mexp-3*esd,vmax/8+vmax/50,'3',fontsize=10,horizontalalignment='center')
                ax.plot([mexp+2*esd,mexp+2*esd],[0,vmax/4],color='magenta',linewidth='2')
                ax.plot([mexp-2*esd,mexp-2*esd],[0,vmax/4],color='magenta',linewidth='2')
                ax.text(mexp+2*esd,vmax/4+vmax/50,'2',fontsize=10,horizontalalignment='center')
                ax.text(mexp-2*esd,vmax/4+vmax/50,'2',fontsize=10,horizontalalignment='center')
                ax.plot([mexp+1*esd,mexp+1*esd],[0,vmax/2],color='magenta',linewidth='2')
                ax.plot([mexp-1*esd,mexp-1*esd],[0,vmax/2],color='magenta',linewidth='2')
                ax.text(mexp+1*esd,vmax/2+vmax/50,'1',fontsize=10,horizontalalignment='center')
                ax.text(mexp-1*esd,vmax/2+vmax/50,'1',fontsize=10,horizontalalignment='center')
                ax.tick_params(axis='both', labelsize=10)
                ax.set_ylabel('cell lines',rotation='vertical',fontsize=14)
                ax.set_xlabel('expression',rotation='horizontal',fontsize=14)
                plt.title(f'{t}',fontsize=14)
                plt.show()
                
                
                s = input('Plot expression/CNV correlation ? ([Y],N) ? : ')
                if len(s) > 0 and s[0].upper() == 'N':
                    continue
                fig, ax = plt.subplots(1,2, figsize =(12,4))
                ax1 = ax[0]
                ax2 = ax[1]

                exp = expr[['ModelID',t]]
                cnv = cnvs[['ModelID',t]]
                exp_cnv = pd.merge(exp,cnv, on = 'ModelID')
                deps.reset_index(inplace = True)
                dep = deps[['ModelID',t]]
                exp_cnv_dep = pd.merge(exp_cnv,dep,on = 'ModelID')
                deps.set_index('ModelID', inplace = True)
                
                
                x = exp_cnv_dep[exp_cnv_dep.columns[1]].to_numpy()
                y = exp_cnv_dep[exp_cnv_dep.columns[2]].to_numpy()
                z = exp_cnv_dep[exp_cnv_dep.columns[3]].to_numpy()
                rain = plt.get_cmap('rainbow')

                m, b = np.polyfit(x, y, 1)
                cc1 = exp_cnv_dep[exp_cnv_dep.columns[1]].corr(exp_cnv_dep[exp_cnv_dep.columns[2]])
                cc2 = exp_cnv_dep[exp_cnv_dep.columns[1]].corr(exp_cnv_dep[exp_cnv_dep.columns[3]])
                cc3 = exp_cnv_dep[exp_cnv_dep.columns[2]].corr(exp_cnv_dep[exp_cnv_dep.columns[3]])

# plot CNV x Expr with dep colouring

                # rx = x.max()-x.min()
                # ry = y.max()-y.min()
                # xm = x.max()*(ry/rx)
                # ym = m*xm
                # theta = degrees(atan2(ym,xm))
                # print(f'plotting angle {theta}, x = {x.max()}, y = {m*x.max()}')

                ax1.plot(x, m*x + b, color = 'grey')
                
                ax1.set_title(f' r = {cc1:.4f}',fontsize=14)
                # ax1.title.set_text(f' r = {cc1:.4f}')
                ax1.tick_params(axis='both', labelsize=10)
                ax1.set_ylabel('relative copy number',rotation='vertical',fontsize=14)
                ax1.set_xlabel('expression',rotation='horizontal',fontsize=14)
                
                if split_pop:
                    ax1.plot([tro,tro],[0,cnv.max()[1]],color='black',linewidth='1',linestyle='dashed')

                
                sc = ax1.scatter(x,y,c=z,s=2,cmap = rain,vmin=0.0,vmax=1.0)
                cb = fig.colorbar(sc,ax=ax1)
                cb.set_label('dependency',fontsize=14,rotation='vertical')
                cb.set_ticks([0.0,0.2,0.4,0.6,0.8,1.0])
                cb.ax.tick_params(labelsize=10)

# plot Dep x Expr with CNV colouring

                m, b = np.polyfit(x, z, 1)
                w = np.array([p for p in x if m*p+b <=1.0 and m*p+b >= 0.0]) # clip ccnd2best fit line to Y=0.0 - 1.0
                
                
                # rx = x.max()-x.min()
                # ry = z.max()-z.min()
                # xm = x.max()*(ry/rx)
                # ym = m*xm
                # theta = degrees(atan2(ym,xm))
                # print(f'plotting angle {theta}')
                
                ax2.set_ylim(0.0,1.0)
                ax2.set_title(f' r = {cc2:.4f}',fontsize=14)

                ax2.plot(w, m*w + b, color = 'grey')
                ax2.tick_params(axis='both', labelsize=10)
                ax2.set_ylabel('dependency',rotation='vertical',fontsize=14)
                ax2.set_xlabel('expression',rotation='horizontal',fontsize=14)
                sc = ax2.scatter(x,z,c=y,s=2,cmap = rain,vmin=0.0,vmax = cnv.max().iloc[1])
                cb = fig.colorbar(sc,ax=ax2)
                cb.set_label('relative copy number',fontsize=14,rotation='vertical')
                cb.set_ticks([int(i) for i in range(int(cnv.max().iloc[1]+1.0))])
                cb.ax.tick_params(labelsize=10)

                plt.suptitle(f'{t}',fontsize=16)              
                plt.show()
                                
        else :
            print(f'{task} not understood')    
            
            
        



                