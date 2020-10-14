import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os.path
from sklearn.decomposition import PCA

def load_data(samples   = ['1', '2', '3', '4','5'], 
              types     = ['3E2H','4E1H','4E2H','4H'],
              labels    = ['b', 'g', 'c'],
              mutations = ['1', '2', '3', '4'],
              dropna    = True,
              scale_id  = 0.2,
             ):
    """
    Load the output data from post-processing.
    Make sure the file path is at the same diretory of this toolbox.
    Make sure the naming of file is in the same fashion: r_[sameple]_[type]_[label]_[metation].[rgyr/prmsd/ermsd/rmsf/nc].
    Make sure the file format is in time series with the same stepsize (except rmsf).
    
    samples: an array contains all the index for samples(replica).
    types: an array contains all the names of types.
    labels: an array contains all the labels for good/bad/crystal.
    mutations: an array contains all the index for different mutaions.
    dropna: clip all timeseries to align with shorest one (except rmsf).
    scale_id: scaling factor to transform from the frame to time. default: 0.2 (5 fs -> 1 ns)
    
    output: 5 pd.DataFrame objects: rgyr, prmsd, ermsd, rmsf, nc
    """
    samples   = ['1', '2', '3', '4']
    types     = ['3E2H','4E1H','4E2H','4H']
    labels    = ['b', 'g', 'c']
    mutations = ['1', '2', '3', '4']

    df_rgyr  = []
    df_prmsd = []
    df_ermsd = []
    df_rmsf  = []
    df_nc    = []

    for s in samples:
        for t in types:
            for l in labels:
                for m in mutations:
                    protein_name = 'r_' + s + '_' + t + '_' + l + '_' + m
                    if os.path.isfile( protein_name + '.rgyr'):
                        df_rgyr.append(  pd.read_csv(protein_name + '.rgyr',  sep=" ", index_col=0, header=None, names=[protein_name]) )
                        df_prmsd.append( pd.read_csv(protein_name + '.prmsd', sep=" ", index_col=0, header=None, names=[protein_name]) )
                        df_ermsd.append( pd.read_csv(protein_name + '.ermsd', sep=" ", index_col=0, header=None, names=[protein_name]) )
                        df_rmsf.append(  pd.read_csv(protein_name + '.rmsf',  sep=" ", index_col=0, header=None, names=[protein_name]) )
                        df_nc.append(    pd.read_csv(protein_name + '.nc'  ,  sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )

    if dropna:
        rgyr  = pd.concat(df_rgyr , axis=1, sort=False).dropna()
        prmsd = pd.concat(df_prmsd, axis=1, sort=False).dropna()
        ermsd = pd.concat(df_ermsd, axis=1, sort=False).dropna()
        nc    = pd.concat(df_nc   , axis=1, sort=False).dropna()
        rmsf  = pd.concat(df_rmsf , axis=1, sort=False)

    rgyr.index  = rgyr.index  * scale_id
    prmsd.index = prmsd.index * scale_id
    ermsd.index = ermsd.index * scale_id
    nc.index    = nc.index    * scale_id
    
    return rgyr, prmsd, ermsd, rmsf, nc

def ACF(x, exclude_steps=4):
    """
    calculate the autocorrelation funtion (ACF) given a time series.
    
    x: an array or pd.DataFrame/pd.Series contains the time series of property.
    exclude_steps: the initial steps to exclude if the time series contains equilibration.
    ouput: an array of a time series of ACF
    """
    x = x[exclude_steps:]
    r = np.correlate(x-x.mean(),x-x.mean(), mode='full')[len(x):]
    
    return r/r[0]

def sample(trj, delta=500, skip_initial=4):
    """
    denoise the trj by calculating the mean for every time interval given interval size delta.
    
    trj: pd.DataFrame of trjactories with columns indicating different protein.
    delta: the size of time interval to calculate the mean. make sure it is within [1,total length of trj].
    skip_initial: the initial steps to be excluded if trj contains unwanted steps such as equilibration.
    
    output: 1 pd.DataFrame objects of trj sampled by mean.
    """
    trj_sample = pd.DataFrame({})
    for column in trj.columns:
        temp = []
        for i in range(len(trj[column])//delta+1):
            temp.append(trj[column].iloc[skip_initial+delta*i:skip_initial+delta*(i+1)].mean())
        trj_sample[column] = temp 
    trj_sample.index = np.arange(0,len(trj_sample[column]),1)*delta*0.2
    
    return trj_sample


def plot_type_trj(trj, ptype, stepsize=1, threshold=0, ls='-', marker=None, output_violate=False):
    """
    plot the (time) series of a set of trajetories given certain protein type.
    
    trj: pd.DataFrame of trjactories with columns indicating different protein.
    ptype: protein type (ex. 3E2H, 4H, 4E1H, 4E2H ... etc). Make sure the string of ptype is within the coloumn index of trj.
    stepsize: the stepsize to plot the array
    threshold: the value to plot a horizontal line.
    output_violate: output an array contains all the proteins which exceeds the threshold.
    ls: line style for matplotlib.pyplot.plot
    marker: marker style for matplotlib.pyplot.plot
    
    output: (if output_violate != 0) an array of strings contains all the proteins which exceeds the threshold.
    """
    
    violates = []
    for column in trj.columns:
        if ptype in column:
            if threshold != 0:
                if np.max(trj[column]) >= threshold:
                    violates.append(column)
                    plt.plot( trj[column].dropna().iloc[0:-1:stepsize].index, trj[column].dropna().iloc[0:-1:stepsize], label=column, ls=ls, marker=marker)
                else:
                    if 'b' in column: 
                        plt.plot( trj[column].dropna().iloc[0:-1:stepsize].index, trj[column].dropna().iloc[0:-1:stepsize], c='red', ls=ls, marker=marker)
                    elif 'g' in column: 
                        plt.plot( trj[column].dropna().iloc[0:-1:stepsize].index, trj[column].dropna().iloc[0:-1:stepsize], c='blue', ls=ls, marker=marker)
                    else: 
                        plt.plot( trj[column].dropna().iloc[0:-1:stepsize].index, trj[column].dropna().iloc[0:-1:stepsize], c='green', ls=ls, marker=marker)
            else: 
                if 'b' in column: 
                    plt.plot( trj[column].dropna().iloc[0:-1:stepsize].index, trj[column].dropna().iloc[0:-1:stepsize], c='red', ls=ls, marker=marker)
                elif 'g' in column: 
                    plt.plot( trj[column].dropna().iloc[0:-1:stepsize].index, trj[column].dropna().iloc[0:-1:stepsize], c='blue', ls=ls, marker=marker)
                else: 
                    plt.plot( trj[column].dropna().iloc[0:-1:stepsize].index, trj[column].dropna().iloc[0:-1:stepsize], c='green', ls=ls, marker=marker)              
                
    if threshold != 0:
        plt.axhline(y=threshold, ls='--', c='r')
        plt.legend(fontsize=15)
      
    plt.title(ptype, fontsize=15)
    plt.xlabel("time ($ns$)", fontsize=15)
    
    return violates if output_violate else None

def plot_hist_for_good_bad(trj, ptype, vertical=True, xlim=None):
    """
    plot the histogram of certain type of protein according to their labels (g/b).
    
    trj: pd.DataFrame of trjactories with columns indicating different protein.
    ptype: protein type (ex. 3E2H, 4H, 4E1H, 4E2H ... etc). Make sure the string of ptype is within the coloumn index of trj.
    xlim: a tuple to adjust the xlim for plotting
    """
    good, bad = [], []
    
    for column in trj.columns:
        if ptype in column and 'b' in column:
            bad = np.concatenate((bad, trj[column].values), axis=-1)
        if ptype in column and 'g' in column:
            good = np.concatenate((good, trj[column].values), axis=-1)

    hist = pd.DataFrame({ "Good":good, "Bad" :bad})

    for column in hist.columns:
        sns.distplot(hist[column], label=column, vertical=vertical)

    plt.legend(fontsize=15)
    plt.title(ptype, fontsize=15)
    plt.xlim(xlim)
    
    return None


def pca_rmsf(ptype, rmsf, clip_from_start= 5, clip_from_tail = 10, output_pca_components= False):
    """
    plot the pca analysis for certain type of protein.
    
    ptype: protein type (ex. 3E2H, 4H, 4E1H, 4E2H ... etc). Make sure the string of ptype is within the coloumn index of trj.
    rmsf: pd.DataFrame of rmsf data with columns indicating different protein.
    clip_from_start: exclude the rmsf from the starting  
    clip_from_tail: exclude the rmsf from the ending
    output_pca_components: boolean value decide to get the output or not.
    
    output: an array of 2 principle components (first and second)
    """
    epitope = {
        "3E2H":[12,17],
        "4E2H":[45,50],
        "4E1H":[11,16],
        "4H":[[11,39],[57,71]]
    }
    
    mask_ptype = rmsf.columns.str.contains(ptype) * ~rmsf.columns.str.contains("c")
    colors = np.where( [("b" in column) & ("c" not in column) for column in rmsf.columns[mask_ptype]],  'r', 'b')
    data = rmsf[rmsf.columns[mask_ptype]].dropna().iloc[clip_from_start:-clip_from_tail].transpose()

    pca = PCA(n_components=2)
    pca.fit(data)
    proj = np.dot( pca.components_, rmsf[rmsf.columns[mask_ptype]].dropna().iloc[clip_from_start:-clip_from_tail])
    cmap = np.where( ["b" in column for column in rmsf.columns[mask_ptype]],  'r', 'b')

    fig = plt.figure(figsize=(15,6))
    plt.subplots_adjust(wspace = .1)
    plt.subplots_adjust(hspace = .25)

    fig.add_subplot(1,2,1)
    plt.scatter(proj[0], proj[1], c=cmap)
    plt.xlabel("PCA-1", fontsize=15)
    plt.ylabel("PCA-2", fontsize=15)
    plt.xlim(0,25)

    fig.add_subplot(1,2,2)
    plt.plot(pca.components_[0], label="PCA-1")
    plt.plot(pca.components_[1], label="PCA-2")
    plt.xlabel("indices", fontsize=15)
    if ptype != "4H":
        plt.axvspan(epitope[ptype][0], epitope[ptype][1], alpha=0.3, color='red')
    else:
        plt.axvspan(epitope[ptype][0][0], epitope[ptype][0][1], alpha=0.3, color='red')
        plt.axvspan(epitope[ptype][1][0], epitope[ptype][1][1], alpha=0.3, color='red')
    plt.legend(fontsize=15)
    
    return pca.components_ if output_pca_components else None