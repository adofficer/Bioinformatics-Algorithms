# Final Project

import profile
import numpy as np
import pandas as pd
import urllib2
import json
import sys
from scipy.stats import f_oneway
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.cluster.vq import *
from sklearn import decomposition
from random import random, randint


BASE_URL = 'http://decon.fas.harvard.edu/pyramidal/gene/' 

###### DO NOT RUN WITH DEBUG = TRUE, WILL MAKE TONS OF PCA PLOTS ########
debug = False


def get_json(gene, BASE_URL):

    """
    Grabs expression level JSON of gene from the DeCoN API
    
    Sensitive to API connection errors and JSON errors
    """
    
    try:
        api = urllib2.urlopen(BASE_URL + gene + '/expression/')
    except urllib2.URLError:
        print "Error connecting to API"
        return
    try:
        js = json.load(api)
    except:
        print "Error reading JSON file"
        return
    
    api.close()
    if debug:
        print 'For gene ' + gene + ' the JSON begins with the following:'
        print js[0]
    
    return js
    
def read_json(js):
    """
    Parses the fpkm values from the gene JSONs and merges them into the growing
    dataframe
    """

    if js:
        gene_id = js[0]['gene_id']
        fpkm = dict()
        for i in xrange(len(js)):
            if js[i]['quant_status'] == 'OK':
                fpkm[js[i]['sample_name']] = float(js[i]['fpkm'])
            else:
                fpkm[js[i]['sample_name']] = 'NaN'
        
        if debug:
            print "fpkm dictionary for gene " + gene_id + " looks like:"
            print fpkm
        
        tdf = pd.DataFrame.from_dict(fpkm.items()).T    
        tdf.index = ['sample_name', gene_id]
        tdf.columns = tdf.ix[0]
        tdf = tdf[tdf.index != 'sample_name']
        
        if debug:
            print "tdf for " + gene_id + " has the following column names:"
            print tdf.columns.values.tolist()
            print "tdf for " + gene_id + " has the following row names:"
            print tdf.index.values.tolist()
            
        return tdf
    
    if debug:
        print "Empty JSON file"
    return

# Following code snippet used for pulling the data from the DeCoN API
"""
list of MGI numbers from Loyal
f = open('decon_ref_seqs/decon_master_gene_list')
all_rna = pickle.load(f)
f.close()

# tests the whole JSON import functionality and the interfaces between it
t = pd.concat([read_json(get_json(all_rna[i], BASE_URL)) for i in xrange(5)])
assert len(t) == 3, "Final merging and function connecting going wrong, check all_rna"

# tests the get_json function
assert get_json('Tle4', BASE_URL)[0] == {u'quant_status': u'OK', u'conf_hi': 263.983, u'_state': None, u'conf_lo': 160.512, u'fpkm': 212.247, u'celltype': u'corticothal', u'gene_id': u'Tle4', u'sample_name': u'E15_corticothal', u'timepoint': u'E15'}, "Error in the get_json function"

# tests the read_json function
assert read_json(get_json('Tle4', BASE_URL)).E16_cpn.tolist() == [3.16099], "Error in the read_json function"

# for testing
t = pd.concat([read_json(get_json(all_rna[i], BASE_URL)) for i in xrange(200)])')

# writes the data to a file, only done once
# involves several tens of thousand queries to the DeCoN API
# takes about two hours, NOT ADVISED TO RUN
result = pd.concat([read_json(get_json(i, BASE_URL)) for i in all_rna])
f = open('all_decon_data.pkl', 'w')
pickle.dump(result, f)
f.close()

# for reading the data back in from the .pkl file
f = open('all_decon_data.pkl')
result = pickle.load(f)
f.close

"""

def p_filt(df, min_p):
    """
    Filters the given dataframe by p-value, all lower than min_p are retained
    
    Not optimal workflow, ideally just add a column to subset by, but it works!
    
    Uses the f-test or one-way ANOVA to determine significance between 3 cell
    types in the dataset
    """

    filt_df = df[1:2]
    sub = [x for x in list(df) if 'subcereb' in x]
    cpn = [x for x in list(df) if 'cpn' in x]
    cort = [x for x in list(df) if 'corticothal' in x]
    
    for i in xrange(2, len(df)):
        row = df[i - 1 : i]        

        if max(row.stack().value_counts()) < 11 and not ('NaN' in row.as_matrix()):
            
            if f_oneway(row[sub].as_matrix()[0], 
                        row[cpn].as_matrix()[0],
                        row[cort].as_matrix()[0])[1] < min_p:
                filt_df = pd.concat([filt_df, row])
                if debug:
                    print row.T
    
    return filt_df

# used temporarily to merge the samples by common cell type
# data munging from the format the API returns
def merge_type(df):
    # initialize growing dataframe
    merge = pd.DataFrame( columns = ['subcereb', 'cpn', 'corticothal'])

    # subset by cell type
    sub = [x for x in list(df) if 'subcereb' in x]
    cpn = [x for x in list(df) if 'cpn' in x]
    cort = [x for x in list(df) if 'corticothal' in x]
    
    sub_mat = df[sub]
    cpn_mat = df[cpn]
    cort_mat = df[cort]
    
    # takes average of all reads and merges into one row of data in the growing
    # data frame
    for i in xrange(len(df)):
        t = pd.DataFrame([sub_mat[i:i+1].as_matrix().mean(),
                          cpn_mat[i:i+1].as_matrix().mean(),
                          cort_mat[i:i+1].as_matrix().mean()],
                          columns = ['subcereb', 'cpn', 'corticothal'])
        
        merge = pd.concat([merge, t])
    return merge

"""
temp to load data in to test p_filt
f = open('all_decon_data.pkl')
results = pickle.load(f)
f.close()

this is the 0.05 filtered list on my local
f = open("p0pt5.pkl")
filt = pickle.load(f)
f.close()
t = filt[1:50]

Genetic Algorithm:

Randomly choose a starting set from the rows in df
Evaluate the objective function for all of the entries in the set
  Quantify how well the 12 samples are clustered into the 3 depths within the brain
  PCA then K-means and save which group they are in
  Score the groups, want penalties for unevenness and samples in wrong groups
Choose the top n
Genetically change the n
Repeat
"""

def genet_alg(df, set_size, targ_max, top_n, mut_rate = 0.1, max_loop = 500):
    """
    poss is a dictionary of the following form:
    str: [list of 3 str, int]
    
    str is a comma-merged list of the subset of rows in df
    list of 3 is the labeled k-means clusters
    int is the score of the list of 3
    """
    
    # initialize the variables and starting set size
    count = 0
    last = 0
    this = 0
    con = 0
    is_converged = False
    poss = dict()
    best = dict()
    
    # random initialization of the first sets
    while len(poss) < top_n * 2:
        t = [randint(0, len(df) - 1) for x in xrange(set_size)]
        if len(set(t)) == len(t):
            poss[','.join([str(x) for x in t])] = [None, None]
            
    if debug:
        print "Initial set is defined as:"
        print poss.keys()
    
    while this < targ_max and count <= max_loop and not is_converged:
    # for i in xrange(10):
        # change up the dictionary, now because of the while statement
        if len(best) != 0:
            poss = genet_switch(df, best, mut_rate, top_n)
        
        # evaluate the scoring function for all entries in the set
        poss = eval_set(df, poss)
        
        last = max([x[1] for x in poss.values()])        

        # keep only the top n that are good enough and repeat
        best = filt_poss(df, poss, top_n)
        this = max([x[1] for x in poss.values()])
        
        # tests for convergence, no sense in looping more than we need to
        if last - this == 0:
            con += 1
            if con > (100.0 / mut_rate):
                is_converged = True
        else:
            con = 0
        count += 1
    
    if is_converged:
        print "Sequence converged"
    elif count > max_loop:
        print "Maximum number of iterations performed"
    else:
        print "Hit target score!"
    
    return best


def eval_set(df, poss):
    """
    Takes the full df and the dict of possibilities and returns the full list
    of clusters with scores
    
    Input: dict: {str: None}
    Output: dict: {str: [list, int]}
    """
    
    for entry in poss.keys():
        # initialize the temp array with correct rows
        rows = [df[int(i):int(i) + 1] for i in entry.split(',')]
        t = pd.concat(rows)
    
        if debug:
            print t
            print t.T
        
        # perform the PCA and the kmeans clustering
        forpca = np.array(t.T, dtype = 'float_')
        pca = decomposition.PCA(n_components = 3)
        pca.fit(forpca)
        forpca = pca.transform(forpca)        
        res, lab = kmeans2(forpca, 3, minit = 'points')
        
        # option to visualize the clustering
        if debug:
            colors = ([([0.1,1,0.1],[1,0.1,0.1],[0.1,0.1,1])[i] for i in lab])
            
            # actually plots the figure
            fig, ax = plt.subplots()
            ax.scatter(forpca[:,0],forpca[:,1], c=colors)
             
            # mark centroids as (X)
            ax.scatter(res[:,0],res[:,1], marker='o', s = 500, linewidths=2, c='none')
            ax.scatter(res[:,0],res[:,1], marker='x', s = 500, linewidths=2)
             
            # adds labels with explained variance
            ax.set_xlabel('PC1, Explained variance: ' +
                          str(pca.explained_variance_ratio_[0] * 100)[0:4])
            ax.set_ylabel('PC2, Explained variance: ' +
                          str(pca.explained_variance_ratio_[1] * 100)[0:4])
            
            
        
        # grabs the clusters from the k-means results and scores them
        names = list(df)
        c0 = [names[i] for i in xrange(12) if lab[i] == 0]
        c1 = [names[i] for i in xrange(12) if lab[i] == 1]
        c2 = [names[i] for i in xrange(12) if lab[i] == 2]
        
        clust = [c0, c1, c2]
        
        poss[entry] = [clust, score_clust(clust)]
        
    return poss

def filt_poss(df, poss, top_n):
    """
    Takes the list of possibilities, determines which are the top N
    and then returns only those that pass the threshold
    input: dict: {str: [list, int]}
    output: same as input
    """
    
    best = dict()
    top_scores = [x[1] for x in sorted(poss.values())][0:top_n + 1]
    
    if debug:
        print "Top scores are:"
        print top_scores
    
    for s in poss.keys():
        
        # check if score of this set is in the top scores and there is room
        if poss[s][1] in top_scores and len(best) < top_n:
            best[s] = poss[s]

    if debug:
        print "Best sets are"
        print best
    
    return best

def score_clust(clust):
    score = 100
    
    s_dict = {'011': 15, '012': 10, '111': 20, '013': 10,
              '022': 15, '112': 20, '122': 20, '113': 20,
              '014': 10, '024': 25, '000': 30, '003': 5,
              '033': 30, '114': 40, '024': 25, '034': 40,
              '004': 0, '001': 25, '044': 40, '002': 10,
              '124': 50, '023': 20, '133': 55, '233': 60,
              '144': 65, '134': 55, '123': 20, '244': 70,
              '234': 65, '224': 55, '223': 25, '222': 35,
              '334': 80}
    
    for s in clust: 
        sub = len([1 for x in s if x[-9:] == '_subcereb'])
        cpn = len([1 for x in s if x[-4:] == '_cpn'])
        cort = len([1 for x in s if x[-12:] == '_corticothal'])
        summ = "".join([str(x) for x in sorted([sub, cpn, cort])])
        score -= s_dict[summ]
        
        if debug:
            print score
            print "The code is",
            print summ
            print "Just subtracted"
            print s_dict[summ]
    return score

def genet_switch(df, best, mut_rate, top_n):
    """
    performs the genetic switching, takes a dict of length top_n
    and returns a dict of length top_n * 2
    """
    
    out = dict()

    best_list = [[int(x) for x in s.split(',')] for s in best.keys()]
    bl = len(best_list)
    set_size = 2 * top_n

    if debug:
        print best_list
    
    while len(out) < set_size:
        if mut_rate > random():
            
            # choose a random item from best_list and a random index to mutate
            i = randint(0, bl - 1)
            t = list(best_list[i])
            i = randint(0, len(t) - 1)
            t[i] = randint(0, len(df) - 1)
            
            # check if all entries are unique and we aren't repeating
            if len(set(t)) == len(t) and t not in best_list:
                out[','.join([str(x) for x in t])] = [None, None]
        else:
            # crossover at a random point
            t1 = list(best_list[randint(0, bl - 1)])
            t2 = list(best_list[randint(0, bl - 1)])
            i = randint(0, len(t1) - 1)
            c1 = list(t1[:i] + t2[i:])
            c2 = list(t2[:i] + t1[i:])
            
            # check if all entries are unique and we aren't repeating
            if len(set(c1)) == len(c1) and len(set(c2)) == len(c2) and c1 not in best_list and c2 not in best_list:
                out[",".join([str(x) for x in c1])] = [None, None]
                out[",".join([str(x) for x in c2])] = [None, None]
    
    return out

# test the scoring function
assert(score_clust([[u'E18_corticothal', u'E16_corticothal', u'E15_corticothal', u'P1_corticothal'], [u'E16_subcereb', u'E15_cpn', u'E15_subcereb', u'P1_subcereb', u'E18_subcereb'], [u'E18_cpn', u'P1_cpn', u'E16_cpn']]) == 85)



####### visualization methods from here on #######

def heatmap(df,
            edgecolors='w',
            cmap=mpl.cm.RdBu,
            log=False,vmin=0,vmax=500):    
    width = len(df.columns)/4
    height = len(df.index)/4
    
    fig, ax = plt.subplots(figsize=(width,height))
      
    heatmap = ax.pcolor(df,
                        edgecolors=edgecolors,  # put white lines between squares in heatmap
                        cmap=cmap,
                        vmin=vmin, # defaults to 0
                        vmax=vmax, # defaults to 500
                        norm=mpl.colors.LogNorm() if log else None)
    
    ax.autoscale(tight=True)  # get rid of whitespace in margins of heatmap
    ax.set_aspect('equal')  # ensure heatmap cells are square
    ax.xaxis.set_ticks_position('top')  # put column labels at the top
    ax.tick_params(bottom='off', top='off', left='off', right='off')  # turn off ticks
    
    plt.yticks(np.arange(len(df.index)) + 0.5, df.index)
    plt.xticks(np.arange(len(df.columns)) + 0.5, df.columns, rotation=90)
    
    # ugliness from http://matplotlib.org/users/tight_layout_guide.html
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="1%")
    plt.colorbar(heatmap, cax=cax)
    return fig

def pca_vis_3d(df):
    """
    Takes as input an array where rows are probes and columns are samples
    Performs PCA analysis and then kmeans clusters the results down to two dims
    
    Much of this code taken from:
    http://blog.mpacula.com/2011/04/27/k-means-clustering-example-python/
    
    """
    
    # data transformation and PCA computation, transforms in-place
    forpca = np.array(df.T, dtype = 'float_')
    pca = decomposition.PCA(n_components = 3)
    pca.fit(forpca)
    forpca = pca.transform(forpca)

    # performs the kmeans clustering, 3 set as constant for now
    res, lab = kmeans2(forpca, 3, minit = 'points')
    
    colors = ([([0.4,1,0.4],[1,0.4,0.4],[0.1,0.8,1])[i] for i in lab])
    
    # actually plots the figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(forpca[:, 0], forpca[:, 1], forpca[:, 2], c = colors)

    # shows the locations of each of the centroids
    ax.scatter(res[:,0],res[:,1],res[:,2], marker='o', s = 500, linewidths=2, c='none')
    ax.scatter(res[:,0],res[:,1],res[:,2], marker='x', s = 500, linewidths=2)

    # adds labels with explained variance
    ax.set_xlabel('PC1, Explained variance: ' + str(pca.explained_variance_ratio_[0] * 100)[0:4])
    ax.set_ylabel('PC2, Explained variance: ' + str(pca.explained_variance_ratio_[1] * 100)[0:4])
    ax.set_zlabel('PC3, Explained variance: ' + str(pca.explained_variance_ratio_[2] * 100)[0:4])    
    
    plt.show()
    
    pass


def pca_vis_2d(df):
    """
    Takes as input an array where rows are probes and columns are samples
    Performs PCA analysis and then kmeans clusters the results down to two dims
    
    Much of this code taken from:
    http://blog.mpacula.com/2011/04/27/k-means-clustering-example-python/
    
    """
    
    # data transformation and PCA computation, transforms in-place
    forpca = np.array(df.T, dtype = 'float_')
    pca = decomposition.PCA(n_components = 2)
    pca.fit(forpca)
    forpca = pca.transform(forpca)

    # performs the kmeans clustering, 3 set as constant for now
    res, lab = kmeans2(forpca, 3, minit = 'points')
    
    colors = ([([0.1,1,0.1],[1,0.1,0.1],[0.1,0.1,1])[i] for i in lab])
    
    # actually plots the figure
    fig, ax = plt.subplots()
    ax.scatter(forpca[:,0],forpca[:,1], c=colors)
     
    # mark centroids as (X)
    ax.scatter(res[:,0],res[:,1], marker='o', s = 500, linewidths=2, c='none')
    ax.scatter(res[:,0],res[:,1], marker='x', s = 500, linewidths=2)
     
    # adds labels with explained variance
    ax.set_xlabel('PC1, Explained variance: ' + str(pca.explained_variance_ratio_[0] * 100)[0:4])
    ax.set_ylabel('PC2, Explained variance: ' + str(pca.explained_variance_ratio_[1] * 100)[0:4])
    
    # adds the data annotations
    for i, l in enumerate(list(df)):
        if debug:
            print (forpca[:, 0][i], forpca[:, 1][i])
        
        ax.annotate(l, (forpca[:, 0][i], forpca[:, 1][i]))
    
    pass
