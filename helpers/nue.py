import statsmodels as sm
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
import matplotlib as mpl
import pdg
api = pdg.connect()

# dictionary mapping signal to ints. Signal == 0 is assumed to be the desired topology. 
signal_dict = {"nueCC":0,"numuCC":1,"nuNC":2,"nuOther":3,"dirt":4,"cosmic":5, "intime":6}
signal_labels = list(signal_dict.keys())
# dictionary mapping particle to pdg code

pdg_dict = {
    "e-": {"pdg":11, "mass":0.000511},
    "mu-": {"pdg":13, "mass":0.105658},
    "gamma": {"pdg":22, "mass":0},
    "p": {"pdg":2212, "mass":0.938272},
    "pi0": {"pdg":111, "mass":0.134976},
    "pi+": {"pdg":211, "mass":0.139570},
    "n": {"pdg":2112, "mass":0.939565},
    "other": {"pdg":0, "mass":0}
}

colors = ["C0","C1","C2","C3","C4","C5","C6"]

nue_scale = 1.0 # 1.0 is the default value
int_scale = 1.0 # scale for intime cosmics

def get_hist_weights(series: pd.Series):
    return np.ones_like(series)/float(len(series))

def calc_err(passed: list,total: list):
    """
    Calculates the Binomial Error for a pass/fail selection
    
    Parameters
    ----------
    passed: list of ints
        number of events in a bin that passed the selection
    total: list of ints
        original number of events in that bin
    
    Returns
    -------
    err: list of lists
        list of lower and upper errors for each bin
    """
    err = [[],[]]
    eff = passed/total
    for i in range(len(passed)):
        val_passed = passed[i]
        val_tot = total[i]
        interval = sm.stats.proportion.proportion_confint(val_passed,val_tot,method='wilson')
        err[0].append(abs(eff[i]-interval[0]))
        err[1].append(abs(eff[i]-interval[1]))
    return err

def flatten_df(df: pd.DataFrame): 
    """
    Flattens a multi-index dataframe. Changes column names from "col.sub0.sub1" to "col_sub0_sub1"
    
    Parameters
    ----------
    df: pandas dataframe
        dataframe to be flattened
    
    Returns
    -------
    flat_df: flattened pandas dataframe
    """
    flat_df = df.copy()
    flat_df.reset_index(inplace=True)
    flat_df.columns = ['_'.join(col) for col in flat_df.columns.values]
    flat_df.columns = [col.strip('_') for col in flat_df.columns]
    return flat_df

def get_slc(df: pd.DataFrame):
    """
    Returns a dataframe containing only slices (duplicate slices are not dropped, pfps are dropped).
    "Duplicate" slices references slices corresponding to the same neutrino interaction. In other words,
    allows "slice double counting." 
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.slc__index"])
    return nodup_df

def get_evt(df: pd.DataFrame):
    """
    Returns a dataframe containing only events (duplicate slices are dropped).
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.mc.nu__index"])
    return nodup_df

def get_signal_evt(df: pd.DataFrame):
    """
    Returns a dataframe containing only signal **events** (duplicate slices are dropped).
    Assumes the input dataframe has a "signal" column, where signal=0. 
    """
    # only run on flattened dataframes 
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.mc.nu__index"])
    return nodup_df[nodup_df.signal == 0]

def get_backgr_evt(df: pd.DataFrame):
    """
    Returns a dataframe containing only background events (duplicate slices are dropped).
    
    Parameters
    ----------
    df: input dataframe
    
    Returns
    -------
    backgr_df: dataframe containing only background events
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.mc.nu__index"])
    slices_df = nodup_df[nodup_df["signal"]!=0]
    return slices_df

def get_signal_slc(df: pd.DataFrame):
    """
    Returns a dataframe containing only signal slices (Duplicate slices are not dropped, pfps are dropped)
    "Duplicate" slices references slices corresponding to the same neutrino interaction. In other words,
    allows "slice double counting."  
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.slc__index"])
    return nodup_df[nodup_df.signal==0]

def get_backgr_slc(df: pd.DataFrame):
    """
    Returns a dataframe containing only background slices (Duplicate slices are not dropped, pfps are dropped)
    "Duplicate" slices references slices corresponding to the same neutrino interaction. In other words,
    allows "slice double counting." 
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.slc__index"])
    return nodup_df[nodup_df.signal!=0]

def get_slices(df: pd.DataFrame,int_type):
    """
    Returns a dataframe containing slices of a certain interaction type (Duplicate slices are not dropped.)
    "Duplicate" slices references slices corresponding to the same neutrino interaction. In other words,
    allows "slice double counting." 
    
    Parameters
    ----------
    df: input dataframe
    int_type: int
        interaction type
        
    Returns
    -------
    slices_df: dataframe containing only slices of the specified interaction type
    """
    # nodup_df = df.drop_duplicates(subset=["ntuple","nu_index","entry"])
    # slices_df = nodup_df[nodup_df["signal"]==int_type]
    slices_df = df[df["signal"]==int_type]
    return slices_df

### PLOTTING ###

def plot_var(df, var: str, bins: np.ndarray,
             label: str = "", title: str = "",
             stacked: bool = True, 
             scale: list = None,
             mult_factor: float = None,
             cut_val: list = None):
    """
    Plots a variable for each interaction type in a histogram.
    
    Parameters
    ----------
    df: input dataframe or list of dataframes
    var: str
        variable to be plotted, must be a column in the dataframe
    bins: np.ndarray
        histogram binning
    label: str
        x-axis label. If empty, uses the variable name
    title: str
        title of the plot
    stacked: bool
        if True, plots a stacked histogram
    scale: list of floats
        if specified, scales the histograms from df by the specified values 
    mult_factor: float
        if specified, multiplies the signal histogram by this factor
    cut_val: list
        if specified, plots vertical lines at the specified values
    """
    if (type(df) is not list): df = [df]
    if scale == None: scale = list(np.ones(len(df)))
    if (len(scale) != len(df)): 
        print("Error: scale must be the same length as df")
        return 
    if mult_factor==None: mult=1.0
    else: mult=mult_factor
    
    bin_steps   = bins[1:]
    bin_centers = (bins[:-1] + bins[1:])/2
    bin_spacing = bins[1] - bins[0]
    
    hists   = np.zeros((len(signal_dict),len(bin_steps)))
    steps   = np.zeros((len(signal_dict),len(bin_steps))) # for the step plot
    bottom  = np.zeros((len(signal_dict),len(bin_steps))) # for the bar plot 
    order   = np.arange(len(signal_dict),0,step=-1)

    for this_df, this_scale in zip(df,scale):
        for i, entry in enumerate(signal_dict):
            hist, _____ = np.histogram(get_slices(this_df,i)[var],bins=bins)
            hists[i] = hists[i] + this_scale*hist

    hists[0] = mult*hists[0]
    
    for i, entry in enumerate(signal_dict):
        plot_label = signal_labels[i]
        
        if ((mult_factor!=None) & (i==0)): plot_label = signal_labels[i] + f" [x{mult_factor}]"
        if stacked==False:
            plt.step(bins, np.insert(hists[i],obj=0,values=hists[i][0]), label=plot_label,zorder=order[i],color=colors[i])
        if stacked==True: 
            if i==0: steps[i] = hists[i]; 
            else:    steps[i] = hists[i] + steps[i-1]; bottom[i] = steps[i-1]
            # need to append a value at the beginning so that the step will be plotted across the first bin
            fix_step = np.insert(steps[i],obj=0,values=steps[i][0])
            plt.step(bins,fix_step, label=plot_label, zorder=order[i],color=colors[i])
    
    # fix for strange ylim behavior 
    ymin, ymax = plt.ylim()
    if stacked==True:
        for i, entry in enumerate(signal_dict):
            plt.bar( bin_centers, hists[i], bottom=bottom[i], width=bin_spacing, alpha=0.25, color=colors[i])
       
    # plot vertical lines
    if cut_val != None:
        for i in range(len(cut_val)):
            plt.vlines(cut_val[i],ymin,ymax,lw=2,color="gray",linestyles="--",zorder=6)
    plt.ylim(0,ymax)
    
    if label == "": plt.xlabel(var)
    else: plt.xlabel(label)
    if title == "": plt.title(var+" distribution (for each interaction type)")
    else: plt.title(title)
    
    plt.legend()
    plt.show()
    
def plot_var_pdg(df, var: str, bins: np.ndarray,
                 label: str = "", title: str = "",
                 scale: list = None,
                 cut_val: list = None):
    """
    Plots a variable for each pdg in a stacked histogram.
    
    Parameters
    ----------
    df: input dataframe or list of dataframes
    var: str
        variable to be plotted, must be a column in the dataframe
    bins: np.ndarray
        histogram binning
    label: str
        x-axis label
    title: str
        title of the plot
    scale: list of floats
        if specified, scales the histograms from df by the specified values 
    cut_val: list
        if specified, plots vertical lines at the specified values
    """
    if (type(df) is not list): df = [df]
    if scale == None: scale = list(np.ones(len(df)))
    if (len(scale) != len(df)): 
        print("Error: scale must be the same length as df")
        return 
    
    bin_steps   = bins[1:]
    bin_centers = (bins[:-1] + bins[1:])/2
    bin_spacing = bins[1] - bins[0]
    
    hists   = np.zeros((len(pdg_dict)-1,len(bin_steps)))
    steps   = np.zeros((len(pdg_dict)-1,len(bin_steps))) # for the step plot
    bottom  = np.zeros((len(pdg_dict)-1,len(bin_steps))) # for the bar plot 
    order   = np.arange(len(pdg_dict)-1,0,step=-1)
    
    other_df = []
    for idx in range(len(df)):
        this_df = df[idx]
        this_scale = scale[idx]
        this_other = df[idx].copy() 
        for i, key in enumerate(list(pdg_dict.keys())): 
            pdg_value = pdg_dict[key]["pdg"]
            # ignore neutrons and other (take care of "other" more explicitly at the end)
            if (key == "n" or key == "other"):
                continue
            pdg_df = this_df[abs(this_df.pfp_shw_truth_p_pdg)==pdg_value]
            hist, _____ = np.histogram(pdg_df[var],bins=bins)
            hists[i] = hists[i] + this_scale*hist
            this_other = this_other[this_other.pfp_shw_truth_p_pdg!=pdg_value]
        other_df.append(this_other)
    
    for idx in range(len(other_df)):
        this_scale = scale[idx]
        this_other = other_df[idx]
        if len(this_other)!=0: 
            hist, _____ = np.histogram(this_other[var],bins=bins)
            hists[-1] = hists[-1] + this_scale*hist        

    for i in range(len(hists)): 
        plot_label = list(pdg_dict.keys())[i]
        if i == len(hists)-1: plot_label = "other"
        if i==0: steps[i] = hists[i]; 
        else:    steps[i] = hists[i] + steps[i-1]; bottom[i] = steps[i-1]
        # need to append a value at the beginning so that the step will be plotted across the first bin
        fix_step = np.insert(steps[i],obj=0,values=steps[i][0])
        plt.step(bins,fix_step, label=plot_label, zorder=order[i])

    ymin, ymax = plt.ylim()
    for i in range(len(hists)): 
        plt.bar( bin_centers, hists[i], bottom=bottom[i], width=bin_spacing, alpha=0.25)
            
    if cut_val != None:
        for i in range(len(cut_val)):
            plt.vlines(cut_val[i],ymin,ymax,lw=2,color="gray",linestyles="--",zorder=6)
    plt.ylim(0,ymax)
    
    if label == "": plt.xlabel(var)
    else: plt.xlabel(label)
    if title == "": plt.title(var+" distribution (separated by particle type)")      
    else: plt.title(title)

    plt.legend()
    plt.show()

### FIXING ###

def shw_energy_fix(df: pd.DataFrame):
    """
    Fixes the shower energy column in the dataframe.
    The best plane is the plane with the greatest number of hits and positive reconstructed energy. 
    
    Parameters
    ----------
    df: input dataframe
    
    Returns
    -------
    shw_df: dataframe with fixed shower energy column
    """
    
    shw_df = df.copy()
    col_nhits2  = 'shw_plane_I2_nHits'
    col_nhits1  = 'shw_plane_I1_nHits'
    col_nhits0  = 'shw_plane_I0_nHits'
    col_energy2 = 'shw_plane_I2_energy'
    col_energy1 = 'shw_plane_I1_energy'
    col_energy0 = 'shw_plane_I0_energy'
    
    if "pfp_shw_plane_I2_nHits" in df.columns:
        col_nhits2 = 'pfp_' + col_nhits2
        col_nhits1 = 'pfp_' + col_nhits1
        col_nhits0 = 'pfp_' + col_nhits0
        col_energy2 = 'pfp_' + col_energy2
        col_energy1 = 'pfp_' + col_energy1
        col_energy0 = 'pfp_' + col_energy0        
    elif "shw_plane_I2_nHits" not in df.columns:
        print("Error: Dataframe does not contain shower energy columns")
        return df

    nhits2 = ((shw_df[col_nhits2] >= shw_df[col_nhits1]) & (shw_df[col_nhits2]>= shw_df[col_nhits0]))
    nhits1 = ((shw_df[col_nhits1] >= shw_df[col_nhits2]) & (shw_df[col_nhits1]>= shw_df[col_nhits0]))
    nhits0 = ((shw_df[col_nhits0] >= shw_df[col_nhits2]) & (shw_df[col_nhits0]>= shw_df[col_nhits1]))

    # if energy[plane] is positive
    energy2 = (shw_df[col_energy2] > 0 )
    energy1 = (shw_df[col_energy1] > 0 )
    energy0 = (shw_df[col_energy0] > 0 )

    conditions = [(nhits2 & energy2),
                (nhits1 & energy1),
                (nhits0 & energy0),
                (((nhits2 & energy2)== False) & (energy1) & (shw_df[col_nhits1]>= shw_df[col_nhits0])), # if 2 is invalid, and 1 is positive and 1>0, go with 1 
                (((nhits2 & energy2)== False) & (energy0) & (shw_df[col_nhits0]>= shw_df[col_nhits1])), # if 2 is invalid, and 0 is positive and 0>1, go with 0
                (((nhits1 & energy1)== False) & (energy2) & (shw_df[col_nhits2]>= shw_df[col_nhits0])), # if 1 is invalid, and 2 is positive and 2>0, go with 2 
                (((nhits1 & energy1)== False) & (energy0) & (shw_df[col_nhits0]>= shw_df[col_nhits2])), # if 1 is invalid, and 0 is positive and 0>2, go with 0
                (((nhits0 & energy0)== False) & (energy2) & (shw_df[col_nhits2]>= shw_df[col_nhits1])), # if 0 is invalid, and 2 is positive and 2>1, go with 2              
                (((nhits0 & energy0)== False) & (energy1) & (shw_df[col_nhits1]>= shw_df[col_nhits2])), # if 0 is invalid, and 1 is positive and 1>2, go with 1 
                ((shw_df[col_nhits2]==-5) & (shw_df[col_nhits1]==-5) & (shw_df[col_nhits0]==-5))]
    shw_energy_choices = [ shw_df[col_energy2],
                    shw_df[col_energy1],
                    shw_df[col_energy0],
                    shw_df[col_energy1],
                    shw_df[col_energy0],
                    shw_df[col_energy2],
                    shw_df[col_energy0],
                    shw_df[col_energy2],
                    shw_df[col_energy1],
                    -1]
    shw_plane_choices = [2,1,0,1,0,2,0,2,1,-1]

    shw_df['shw_energy'] = np.select(conditions, shw_energy_choices, default = -1)
    shw_df["shw_plane"] = np.select(conditions, shw_plane_choices, default = -1)
    shw_df.drop(columns=[col_nhits1, col_nhits2, col_nhits0, col_energy0, col_energy1, col_energy2],inplace=True)
    return shw_df

### SELECTION CUTS ### 

def cutRecoVtxFV(df: pd.DataFrame, col:str="slc_vertex_"):
    """
    Selects slices with reconstructed vertices inside the FV: x[5,180], y[-180,180], z[20,470]
    
    Parameters
    ----------
    df: input dataframe
    col: str
        column name of the reconstructed vertex position
    """
    maskRecoFV = ((abs(df[col+"x"]) > 5)& (abs(df[col+"x"]) < 180)
                    & (df[col+"y"] > -180)  & (df[col+"y"] < 180) 
                    & (df[col+"z"] > 20)    & (df[col+"z"] < 470))
    return df[maskRecoFV]

def cutNotCosmic(df: pd.DataFrame, col:str="slc_is_clear_cosmic"):
    """
    Selects slices that are not clear cosmics.
    
    Parameters
    ----------
    df: input dataframe
    col: str
        column name of the clear cosmic bool
    """
    maskNotCosmic = (df[col]==False)
    return df[maskNotCosmic]

def cutShower(df: pd.DataFrame, col:str="pfp_trackScore"):
    """
    Selects slices with at least one pfp with trackScore < 0.6. 
    This cut can be computationally intensive.
    Ideally perform this cut **after** cutRecoVtxFv and cutNotCosmic.
    
    Parameters
    ----------
    df: input dataframe
    col: str
        column name of the track score
    """
    score_df = df.query(col+" > 0") # require that we look at pfps with a sensible track score
    score_df = score_df[(score_df[col] == score_df.groupby(["ntuple","entry","rec.slc__index"])[col].transform(min))]
    score_df = score_df.query("pfp_trackScore < 0.6")[["ntuple","entry","rec.slc__index"]] # require that the minimum track score is < 0.6
    return df.merge(score_df,how="inner")

def cutShowerEnergy(df: pd.DataFrame, 
                    shw_col:str="shw_energy",
                    score_col:str="pfp_trackScore",
                    shw_cut_val:float=0.1,
                    score_cut_val:float=0.6):
    """
    Selects slices that have at least one pfp that fulfills the conditions: has a trackScore < 0.6 and reco
    shower energy greater than `cut_val` (default 0.1 Gev, or 100 MeV).
    
    Parameters
    ----------
    df: input dataframe
    shw_col: str
        column name of the shower energy
    score_col: str
        column name of the track score
    shw_cut_val: float
        cut value for the shower energy
    score_cut_val: float
        cut value for the track score
    """
    shw_df = df.query("(" + score_col+" < @score_cut_val ) & ("+score_col+" > 0)")
    max_df = shw_df[(shw_df[shw_col] == shw_df.groupby(["ntuple","entry","rec.slc__index"])[shw_col].transform(max))]
    max_df = max_df.query(shw_col+" > @shw_cut_val")
    slc_max_df = max_df[['ntuple','entry','rec.slc__index']].drop_duplicates()
    return df.merge(slc_max_df,how="inner",on=['ntuple','entry','rec.slc__index'])

def cutPreselection(df: pd.DataFrame, 
                    whereRecoVtxFV: bool=True, 
                    whereNotCosmic: bool=True, 
                    whereShower: bool=True):
    """
    Performs all presection cuts: reco vertex in FV, not clear cosmic, and at least one pfp with trackScore < 0.6.
    
    Parameters
    ----------
    df: input dataframe
    whereRecoVtxFV: bool
        if True, performs the reco vertex in FV cut
    whereNotCosmic: bool
        if True, performs the not clear cosmic cut
    whereShower: bool
        if True, performs the at least one pfp with trackScore < 0.6 cut
    
    Returns
    -------
    df: dataframe with preselection cuts applied
    """
    if whereRecoVtxFV:
        df = cutRecoVtxFV(df)
    if whereNotCosmic: 
        df = cutNotCosmic(df)
    if whereShower:
        df = cutShower(df)
    return df

def cutCRUMBS(df: pd.DataFrame, col:str="slc_crumbs_result_score", cut_val:float =-0.15):
    """
    Selects slices with CRUMBS score > -0.15 [or other cut value].
    
    Parameters
    ----------
    df: input dataframe
    col: str
        column name of the CRUMBS score
    cut_val: float
        cut value for the CRUMBS score
    """
    maskCRUMBS = (df[col] > cut_val)
    return df[maskCRUMBS]

def cutContainment(df: pd.DataFrame, whereTrkCont: bool=True, whereShwCont: bool=True, 
                   trk_col: str = "pfp_trk_end_", shw_col: str = "pfp_shw_end_"):
    """
    Selects slices with all pfps contained in the detector volume.
    To decide whether to use the pfp_trk or pfp_shw end position, pfps with trackScore >= 0.5 are considered tracks, and pfps with trackScore < 0.5 are considered showers.
    
    Parameters
    ----------
    df: input dataframe
    whereTrkCont: bool
        if True, performs the track containment cut (uses pfp_trk_end)
    whereShwCont: bool
        if True, performs the shower containment cut (uses pfp_shw_end)
    trk_col: str
        column name of the track end position
    shw_col: str
        column name of the shower end position
    """
    df["shw_exit"] = 0 
    df["trk_exit"] = 0

    df["trk_exit"] = np.where(((df.pfp_trackScore >= 0.5) & 
                                ((abs(df[trk_col+"x"]) > 195) |
                                    (df[trk_col+"y"] < -195)  | (df[trk_col+"y"] > 195) &
                                    (df[trk_col+"z"] < 5)    | (df[trk_col+"z"] > 495))),
                                    1,df["trk_exit"])
    df["shw_exit"] = np.where(((df.pfp_trackScore < 0.5) & 
                                ((abs(df[shw_col+"x"]) > 195) |
                                    (df[shw_col+"y"] < -195)  | (df[shw_col+"y"] > 195) &
                                    (df[shw_col+"z"] < 5)    | (df[shw_col+"z"] > 495))),
                                    1,df["shw_exit"])
    # sum the number of exiting trks/shws 
    pfp_exit_df = df.groupby(["ntuple","entry","rec.slc__index"]).agg(ntrk_exit = ('trk_exit','sum'),
                                                                      nshw_exit = ('shw_exit','sum')).reset_index()
    # require that there are no exiting trks/shw if specified
    if (whereTrkCont): pfp_exit_df = pfp_exit_df.query("ntrk_exit==0")
    if (whereShwCont): pfp_exit_df = pfp_exit_df.query("nshw_exit==0")

    df = df.merge(pfp_exit_df[["ntuple","entry","rec.slc__index"]],how="right")
    df.drop(columns=["trk_exit","shw_exit"],inplace=True)
    return df

### MASKS ### 

def maskTrueVtxFv(nuprim_df: pd.DataFrame, 
                  xmin: float = 5, xmax: float = 180, 
                  ymin: float = -180, ymax: float = 180, 
                  zmin: float = 20, zmax: float = 470):
    """
    Returns a mask for true vertex inside the FV: x[5,180], y[-180,180], z[20,470].
    Intended to be used on the multi-index nuprim dataframe.
    """
    whereFV = ((abs(nuprim_df.position.x) > xmin)& (abs(nuprim_df.position.x) < xmax)
                & (nuprim_df.position.y > ymin)  & (nuprim_df.position.y < ymax)
                & (nuprim_df.position.z > zmin)  & (nuprim_df.position.z < zmax))
    return whereFV

### DATAFRAMES ## 

def defineBackground(nuprim_df: pd.DataFrame):
    """
    Returns the input dataframe with a column "signal" with entries that correspond to different kinds of background.
    Intended to be used on the multi-index nuprim dataframe, thus cosmic backgrounds must be added later.

    The background is defined as [numbers can be overriden by redefining object `signal_dict`]:
    - numu cc FV: 1
    - nc FV: 2
    - other nu, FV: 3
    - dirt: 4
    """
    whereFV = maskTrueVtxFv(nuprim_df)
    if "signal" not in nuprim_df.columns:
        nuprim_df["signal"] = -1
    nuprim_df["signal"] = np.where(whereFV & (nuprim_df.iscc==1) & (abs(nuprim_df.pdg)==14), signal_dict["numuCC"], nuprim_df["signal"]) # numu cc FV 
    nuprim_df["signal"] = np.where(whereFV & (nuprim_df.iscc==0), signal_dict["nuNC"], nuprim_df["signal"]) # nc FV
    nuprim_df["signal"] = np.where(whereFV & (nuprim_df["signal"]<0), signal_dict["nuOther"], nuprim_df['signal']) # other nu, FV 
    nuprim_df["signal"] = np.where(whereFV == False, signal_dict["dirt"], nuprim_df["signal"]) # outside FV
    return nuprim_df

def getPDGCounts(nuprim_df: pd.DataFrame):
    """
    Returns a dataframe with the number of events for each pdg code, as well as the total deposited energy in the AV. 
    Intended to be used on the multi-index nuprim dataframe. Needs to be merged to the neutrino event df after.
    """
    flat_nuprim_df = flatten_df(nuprim_df)
    flat_nuprim_df["prim_abs_pdg"] = abs(flat_nuprim_df["prim_pdg"])
    flat_nuprim_df["prim_depE"] = flat_nuprim_df.prim_startE - flat_nuprim_df.prim_endE

    # if the pdg is > 1e9, it is a nucleus 
    flat_nuprim_df["prim_abs_pdg"] = np.where(flat_nuprim_df.prim_abs_pdg>1e9,
                                            1e9,
                                            flat_nuprim_df.prim_abs_pdg)

    # if the pdg is -1, it's "invisible"
    # the pi0 is always invisible! explicitly specify to keep it
    flat_nuprim_df["prim_abs_pdg"] = np.where(flat_nuprim_df.eval("(prim_depE < 0.05) & (prim_abs_pdg!=111)"),
                                            -1,
                                            flat_nuprim_df["prim_abs_pdg"])

    flat_nuprim_df["prim_exit"] =  np.where((((abs(flat_nuprim_df["prim_end_x"]) > 195) |
                                                (abs(flat_nuprim_df["prim_end_y"]) > 195) &
                                                (flat_nuprim_df["prim_end_z"] < 5)    | (flat_nuprim_df["prim_end_z"] > 495))),
                                            True,False)

    # # get total deposited energy for each neutrino event
    depE_sum =  flat_nuprim_df.groupby(["ntuple","entry","rec.mc.nu__index"])["prim_depE"].sum().reset_index().rename(columns={"prim_depE":"total_prim_depE"})
    # # get whether the primaries are contained
    prim_cont = (flat_nuprim_df.groupby(["ntuple","entry","rec.mc.nu__index"])["prim_exit"]).sum().reset_index().rename(columns={"prim_exit":"prim_exit_count"})
    prim_cont["prim_cont"] = np.where(prim_cont["prim_exit_count"] == 0,True,False)
    prim_cont.drop(columns=["prim_exit_count"],inplace=True)
    prim_cont = prim_cont.merge(depE_sum,on=["ntuple","entry","rec.mc.nu__index"])

    # get counts of each PDG type and rename the columns to the PDG names 
    pdg_counts = (flat_nuprim_df.groupby(["ntuple","entry","rec.mc.nu__index"])["prim_abs_pdg"]).value_counts().unstack(fill_value=0).reset_index()
    for col in pdg_counts.columns:
        if ((type(col) is float)):
            if ((col > 0) & (col < 1e9)):
                name = api.get_particle_by_mcid(col).name
                pdg_counts.rename(columns={col:name},inplace=True)
            if (col >=1e9):
                pdg_counts.rename(columns={col:"nucleus"},inplace=True)
            if (col ==-1):
                pdg_counts.rename(columns={col:"invisible"},inplace=True)

    pdg_counts = pdg_counts.merge(prim_cont,on=["ntuple","entry","rec.mc.nu__index"])
    return pdg_counts

def getPFP(slcshw_df: pd.DataFrame, 
           slctrk_df: pd.DataFrame,
           cheat: bool = False):
    """
    Returns a dataframe with the combined (shw+trk) pfp information for each slice.
    
    Corrects the faulty track length column with new column "fix_trk_len".
    
    Corrects the shower energy and bestplane columns with new columns "shw_energy" and "shw_plane".
    
    By default will drop the pfp with the same index as the slice (pfp_id = slc_self). When using cheat Pandora, we **do not** drop this pfp. 
    """
    
    slcshw_df = flatten_df(slcshw_df)
    slcshw_df["rec.mc.nu__index"] = slcshw_df["slc_tmatch_idx"]
    slcshw_df = shw_energy_fix(slcshw_df)
    # primary pfp corresponds to the neutrino, drop it
    if cheat==False:
        slcshw_df = slcshw_df[slcshw_df.pfp_id != slcshw_df.slc_self]

    slctrk_df = flatten_df(slctrk_df)

    # drop common columns in the slctrk dataframe, keep the trk specific columns 
    test_slctrk_df = slctrk_df.copy()
    for i in np.arange(4,len(slctrk_df.columns)):
        if ("slc" in slctrk_df.columns[i]) | ("pfp" in slctrk_df.columns[i]):
            if ("trk" in slctrk_df.columns[i]):
                continue
            else:
                test_slctrk_df.drop(columns=slctrk_df.columns[i],inplace=True)

    # create list of remaining common columns so we can merge on them 
    col_common = []
    for i in slcshw_df.columns:
        for j in test_slctrk_df.columns:
            if i==j:
                col_common.append(i)
                break

    # merge the two dataframes on the common columns
    slcpfp_df = slcshw_df.merge(test_slctrk_df,on=col_common[:4],how="inner")
    slcpfp_df["fix_trk_len"] = np.sqrt( (slcpfp_df.pfp_trk_start_x - slcpfp_df.pfp_trk_end_x)**2
                                        +(slcpfp_df.pfp_trk_start_y - slcpfp_df.pfp_trk_end_y)**2
                                        +(slcpfp_df.pfp_trk_start_z - slcpfp_df.pfp_trk_end_z)**2)
    return slcpfp_df
