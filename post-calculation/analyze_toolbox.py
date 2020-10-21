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
    
    output: 1 pd.Deries contains 9 pd.DataFrame objects: rgyr, pprmsd, eprmsd, sprmsd, eermsd, sermsd, pnc, enc, rmsf
    """
    
    df_rgyr   = []
    df_pprmsd = []
    df_eprmsd = []
    df_sprmsd = []
    df_eermsd = []
    df_sermsd = []
    df_pnc    = []
    df_enc    = []
    df_rmsf   = []

    for s in samples:
        for t in types:
            for l in labels:
                for m in mutations:
                    protein_name = 'r_' + s + '_' + t + '_' + l + '_' + m
                    if os.path.isfile( protein_name + '.rgyr'):
                        df_rgyr.append(   pd.read_csv(protein_name + '.rgyr',   sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )
                        df_pprmsd.append( pd.read_csv(protein_name + '.pprmsd', sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )
                        df_eprmsd.append( pd.read_csv(protein_name + '.eprmsd', sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )
                        df_sprmsd.append( pd.read_csv(protein_name + '.sprmsd', sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )
                        df_eermsd.append( pd.read_csv(protein_name + '.eermsd', sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )
                        df_sermsd.append( pd.read_csv(protein_name + '.sermsd', sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )
                        df_rmsf.append(   pd.read_csv(protein_name + '.rmsf',   sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )
                        df_pnc.append(    pd.read_csv(protein_name + '.enc'  ,  sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )
                        df_enc.append(    pd.read_csv(protein_name + '.pnc'  ,  sep=" ", index_col=0, header=None, names=[protein_name], skiprows=1) )

    if dropna:
        rgyr   = pd.concat(df_rgyr  , axis=1, sort=False).dropna()
        pprmsd = pd.concat(df_pprmsd, axis=1, sort=False).dropna()
        eprmsd = pd.concat(df_eprmsd, axis=1, sort=False).dropna()
        sprmsd = pd.concat(df_sprmsd, axis=1, sort=False).dropna()
        eermsd = pd.concat(df_eermsd, axis=1, sort=False).dropna()
        sermsd = pd.concat(df_sermsd, axis=1, sort=False).dropna()
        pnc    = pd.concat(df_pnc   , axis=1, sort=False).dropna()
        enc    = pd.concat(df_enc   , axis=1, sort=False).dropna()
        rmsf   = pd.concat(df_rmsf , axis=1, sort=False)

    rgyr.index   = rgyr.index   * scale_id
    pprmsd.index = pprmsd.index * scale_id
    eprmsd.index = pprmsd.index * scale_id
    sprmsd.index = pprmsd.index * scale_id
    eermsd.index = eermsd.index * scale_id
    sermsd.index = sermsd.index * scale_id
    pnc.index     = pnc.index   * scale_id
    enc.index     = enc.index   * scale_id
    
    data = pd.Series([rgyr, pprmsd, eprmsd, sprmsd, eermsd, sermsd, pnc, enc, rmsf])
    data.index = ["rgyr", "pprmsd", "eprmsd", "sprmsd", "eermsd", "sermsd", "pnc", "enc", "rmsf"]
    
    return data

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

def sample(data, parameter, delta=500, skip_initial=4, scale_id=0.2):
    """
    denoise the trj by calculating the mean for every time interval given interval size delta.
    For parameter nc, the scale_id is automatically times 10 assuming the output time frame is 10 times larger.
    
    trj: pd.DataFrame of trjactories with columns indicating different protein.
    delta: the size of time interval to calculate the mean. make sure it is within [1,total length of trj].
    skip_initial: the initial steps to be excluded if trj contains unwanted steps such as equilibration.
    scale_id: scaling factor to transform from the frame to time. default: 0.2 (5 fs -> 1 ns)
    
    output: 1 pd.DataFrame objects of trj sampled by mean.
    """
    trj_sample = pd.DataFrame({})
    trj = data[parameter]
    
    for column in trj.columns:
        temp = []
        for i in range(len(trj[column])//delta+1):
            temp.append(trj[column].iloc[skip_initial+delta*i:skip_initial+delta*(i+1)].mean())
        trj_sample[column] = temp 
    if "nc" not in parameter:
        trj_sample.index = np.arange(0,len(trj_sample[column]),1)*delta*scale_id
    else:
        trj_sample.index = np.arange(0,len(trj_sample[column]),1)*delta*scale_id*10
    
    return trj_sample

def get_sample_points(data, sample_data, standardize=True):
    
    mask_data = data.index[ data.index != 'rmsf' ][1:]
    sample_points = pd.Series({"b":pd.Series([], dtype="float64"), "g":pd.Series([], dtype="float64")})

    for column in data["rgyr"].columns:

        df = pd.concat( [sample_data[parameter][column] for parameter in mask_data], axis=1 )
        df.index = np.arange(len(df))
        df.columns = mask_data

        if "b" in column:
            sample_points["b"] = pd.concat( [sample_points["b"], df], axis=0 )
        elif "g" in column:
            sample_points["g"] = pd.concat( [sample_points["g"], df], axis=0 )

    for index in sample_points.index:
        sample_points[index] = sample_points[index].drop(0, axis=1)
        sample_points[index].index = np.arange(len(sample_points[index]))
        if standardize:
            sample_points[index] = (sample_points[index] - sample_points[index].mean()) / sample_points[index].std()
            
    return sample_points

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

def plot_all_parameter(data_trj, data_hist, ptype, outputfile=False):
    
    fig = plt.figure(figsize=(20,18))
    plt.subplots_adjust(wspace = .01)
    plt.subplots_adjust(hspace = .25)

    ptype = ptype
    mask = (data_trj.index != "rmsf") *  (data_trj.index != "rgyr")
    mask_data = data_trj.index[mask]

    assert len(mask_data) == 8

    for i in range(len(mask_data)):

        parameter = mask_data[i]
        fig.add_subplot(4,8, (1+4*i, 3+4*i))
        if "nc" in parameter:
            plot_type_trj(data_trj[parameter], ptype, stepsize=2, threshold=0)
        else:
            plot_type_trj(data_trj[parameter], ptype, stepsize=25, threshold=0)
        plt.ylabel(parameter, fontsize=15)

        fig.add_subplot(4,8, 4+4*i)
        plot_hist_for_good_bad(data_hist[parameter], ptype)
        plt.xlabel('relatvie counts', fontsize=15)
        plt.ylabel(None)

    fig.tight_layout()
    if outputfile: plt.savefig(ptype+'.jpg')

def pca_rmsf(ptype, rmsf, clip_from_start= 0, clip_from_tail = 0, output_pca_components= False):
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
    select = rmsf[rmsf.columns[mask_ptype]].dropna()
    data = select.iloc[clip_from_start:select.shape[0]-clip_from_tail].transpose()

    pca = PCA(n_components=3)
    pca.fit(data)
    proj = np.dot( pca.components_, select.iloc[clip_from_start:select.shape[0]-clip_from_tail])
    cmap = np.where( ["b" in column for column in rmsf.columns[mask_ptype]],  'r', 'b')

    fig = plt.figure(figsize=(24,6))
    plt.subplots_adjust(wspace = .1)
    plt.subplots_adjust(hspace = .25)

    fig.add_subplot(1,3,1)
    plt.scatter(proj[0], proj[1], c=cmap)
    plt.xlabel("PCA-1", fontsize=15)
    plt.ylabel("PCA-2", fontsize=15)
    plt.xlim(0,25)

    ax = fig.add_subplot(1,3,2, projection='3d')
    ax.view_init(35, 40)
    ax.scatter3D(proj[0], proj[1], proj[2], c=cmap);   
    
    fig.add_subplot(1,3,3)
    plt.plot(np.arange(len(pca.components_[0]))+clip_from_start, pca.components_[0], label="PCA-1")
    plt.plot(np.arange(len(pca.components_[0]))+clip_from_start, pca.components_[1], label="PCA-2")
    plt.xlabel("indices", fontsize=15)
    if ptype != "4H":
        plt.axvspan(epitope[ptype][0], epitope[ptype][1], alpha=0.3, color='red')
    else:
        plt.axvspan(epitope[ptype][0][0]+clip_from_start, epitope[ptype][0][1]+clip_from_start, alpha=0.3, color='red')
        plt.axvspan(epitope[ptype][1][0]+clip_from_start, epitope[ptype][1][1]+clip_from_start, alpha=0.3, color='red')
    plt.legend(fontsize=15)
    
    return pca.components_ if output_pca_components else None
