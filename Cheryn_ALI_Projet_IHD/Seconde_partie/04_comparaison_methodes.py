#!/usr/bin/env python
#Cheryn ALI
#Octobre 2019
#./04_comparaison_methodes.py -m binomial -c viridis
from scipy.stats import binom, hypergeom, chi2_contingency
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import argparse

#PARAM
parser = argparse.ArgumentParser(description='Measures comparison')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage. ')
parser.add_argument('-c', '--color', required=False, default='viridis', help='3D plots Colors viridis or rainbow')
param = parser.parse_args()

# FUNCTION
def compute_pval(g,q,t,c, measure):
    if measure=='binomial': # binom.cdf(>=success, attempts, proba)
        pval = binom.cdf( q - c, q, 1 - float(t)/g)
    elif measure=='hypergeometric': # hypergeom.sf(common-1, population, target, query) = 1-p( X <= x-1 ) = p( X >= x )
        pval = hypergeom.sf(c-1, g, t, q)
    elif measure == 'coverage':
        pval = 1-((c*c)/(q*t))
    elif measure=='chi2':
        contingency_table = np.array([[c, q-c], [t-c, g-q-t+c]])
        pval = chi2_contingency(contingency_table)[1]
    else:
        print('sorry, %s not implemented' % ( measure ))
        exit(1)
    return pval


if param.measure=='binomial':
    #BINOMIAL
    tabQ_binomial = []
    tabT_binomial = []
    tabC_binomial = []
    tabpval_binomial = []

    for q in range(2,51):
        for t in range(2,51):
            for c in range(2,min(q,t)):
                pval = compute_pval(250,q,t,c, 'binomial')
                tabQ_binomial.append(q)
                tabT_binomial.append(t)
                tabC_binomial.append(c)
                tabpval_binomial.append(pval)

    #BINOMIAL Visualisation
    Q_binomial = np.array(tabQ_binomial)
    T_binomial = np.array(tabT_binomial)
    C_binomial = np.array(tabC_binomial)
    pval_binomial = (-1)*(np.log(np.array(tabpval_binomial)))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(Q_binomial, T_binomial, C_binomial, c=pval_binomial, cmap=param.color, linewidth=0.5)
    ax.set_title('BINOMIAL')
    ax.set_xlabel('Query')
    ax.set_ylabel('Target')
    ax.set_zlabel('Commun')
    plt.show()

if param.measure=='hypergeometric':
    #HYPERGEOMETRIC
    tabQ_hypergeometric = []
    tabT_hypergeometric = []
    tabC_hypergeometric = []
    tabpval_hypergeometric = []

    for q in range(2,51):
        for t in range(2,51):
            for c in range(2,min(q,t)):
                pval = compute_pval(250,q,t,c, 'hypergeometric')
                tabQ_hypergeometric.append(q)
                tabT_hypergeometric.append(t)
                tabC_hypergeometric.append(c)
                tabpval_hypergeometric.append(pval)

    #HYPERGEOMETRIC Visualisation
    Q_hypergeometric = np.array(tabQ_hypergeometric)
    T_hypergeometric = np.array(tabT_hypergeometric)
    C_hypergeometric = np.array(tabC_hypergeometric)
    pval_hypergeometric = (-1)*np.log(np.array(tabpval_hypergeometric))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(Q_hypergeometric, T_hypergeometric, C_hypergeometric, c=pval_hypergeometric, cmap=param.color, linewidth=0.5)
    ax.set_title('HYPERGEOMETRIC')
    ax.set_xlabel('Query')
    ax.set_ylabel('Target')
    ax.set_zlabel('Commun')
    plt.show()

if param.measure=='coverage':
    #COVERAGE
    tabQ_coverage = []
    tabT_coverage = []
    tabC_coverage = []
    tabpval_coverage = []

    for q in range(2,51):
        for t in range(2,51):
            for c in range(2,min(q,t)):
                pval = compute_pval(250,q,t,c, 'coverage')
                tabQ_coverage.append(q)
                tabT_coverage.append(t)
                tabC_coverage.append(c)
                tabpval_coverage.append(pval)

    #COVERAGE Visualisation
    Q_coverage = np.array(tabQ_coverage)
    T_coverage = np.array(tabT_coverage)
    C_coverage = np.array(tabC_coverage)
    pval_coverage = (-1)*np.log(np.array(tabpval_coverage))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(Q_coverage, T_coverage, C_coverage, c=pval_coverage, cmap=param.color, linewidth=0.5)
    ax.set_title('COVERAGE')
    ax.set_xlabel('Query')
    ax.set_ylabel('Target')
    ax.set_zlabel('Commun')
    plt.show()

if param.measure=='chi2':
    #CHI2
    tabQ_chi2 = []
    tabT_chi2 = []
    tabC_chi2 = []
    tabpval_chi2 = []

    for q in range(2,51):
        for t in range(2,51):
            for c in range(2,min(q,t)):
                pval = compute_pval(250,q,t,c, 'chi2')
                tabQ_chi2.append(q)
                tabT_chi2.append(t)
                tabC_chi2.append(c)
                tabpval_chi2.append(pval)

    #CHI2 Visualisation
    Q_chi2 = np.array(tabQ_chi2)
    T_chi2 = np.array(tabT_chi2)
    C_chi2 = np.array(tabC_chi2)
    pval_chi2 = (-1)*np.log(np.array(tabpval_chi2))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(Q_chi2, T_chi2, C_chi2, c=pval_chi2, cmap=param.color, linewidth=0.5)
    ax.set_title('CHI2')
    ax.set_xlabel('Query')
    ax.set_ylabel('Target')
    ax.set_zlabel('Commun')
    plt.show()


