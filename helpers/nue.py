import statsmodels as sm
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
import matplotlib as mpl

# dictionary mapping signal to ints. Signal == 0 is assumed to be the desired topology. 
signal_dict = {"nueCC":0,"numuCC":1,"nuNC":2,"nuOther":3,"dirt":4,"cosmic":5}
# dictionary mapping particle to pdg code
pdg_dict =    {"elec":11, "muon":13, "gamma":22, "proton":2212, "pi0":111, "pi":211, "neu":2112, "other":0}

nue_scale = 1.0 # 1.0 is the default value

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
    Returns a dataframe containing only slices (duplicate slices are not dropped, pfps are dropped)
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.slc__index"])
    return nodup_df

def get_evt(df: pd.DataFrame):
    """
    Returns a dataframe containing only events (duplicate slices are dropped)
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","nu_index","entry"])
    return nodup_df

def get_signal_evt(df: pd.DataFrame):
    """
    Returns a dataframe containing only signal **events** (duplicate slices are dropped).
    Assumes the input dataframe has a "signal" column, where signal=0. 
    """
    # only run on flattened dataframes 
    nodup_df = df.drop_duplicates(subset=["ntuple","nu_index","entry"])
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
    nodup_df = df.drop_duplicates(subset=["ntuple","nu_index","entry"])
    slices_df = nodup_df[nodup_df["signal"]!=0]
    return slices_df

def get_signal_slc(df: pd.DataFrame):
    """
    Returns a dataframe containing only signal slices (Duplicate slices are not dropped, pfps are dropped)
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.slc__index"])
    return nodup_df[nodup_df.signal==0]

def get_backgr_slc(df: pd.DataFrame):
    """
    Returns a dataframe containing only background slices (Duplicate slices are not dropped, pfps are dropped)
    """
    nodup_df = df.drop_duplicates(subset=["ntuple","entry","rec.slc__index"])
    return nodup_df[nodup_df.signal!=0]

def get_slices(df: pd.DataFrame,int_type):
    """
    Returns a dataframe containing slices of a certain interaction type (Duplicate slices are not dropped.)
    
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

def plot_var(df: pd.DataFrame, var: str, bins: np.ndarray, 
             df_add: pd.DataFrame = None,
             label: str = "", title: str = "", 
             cut_val: list = None, mult_factor: float = None):
    """
    Plots a variable for each interaction type in a stacked histogram.
    
    Parameters
    ----------
    df: input dataframe
    var: str
        variable to be plotted, must be a column in the dataframe
    bins: np.ndarray
        histogram binning
    df_add: pd.DataFrame
        if specified, adds this dataframe to the plot with scaling factor set by nue.nue_scale.
    label: str
        x-axis label. If empty, uses the variable name
    title: str
        title of the plot
    cut_val: list
        if specified, plots vertical lines at the specified values
    mult_factor: float
        if specified, multiplies the nueCC histogram by this factor
        
    Notes
    -----
    `df_add` should be the dataframe you are scaling. `df` remains unscaled.
    """
    histtype="step"
    lw = 2
    # plot signal==0 histogram
    if (df_add is None):
        if mult_factor==None:
            plt.hist(get_slices(df,0)[var],bins=bins,histtype=histtype,lw=lw,label="nue cc",zorder=5)
        else:
            hist_values, bin_edges = np.histogram(get_slices(df,0)[var],bins=bins)
            plt.step(x=bin_edges[:-1], y=np.round(hist_values*mult_factor),where="post",linewidth=lw, label=f"nueCC [x{mult_factor:.2f}]",zorder=5)
        # plot other histograms 
        for i, entry in enumerate(signal_dict):
            if i==0:
                continue
            plt.hist(get_slices(df,i)[var],bins=bins,histtype=histtype,lw=lw,label=entry)
    elif (df_add is not None):
        # combine the two dataframes and scale 
        for i, entry in enumerate(signal_dict):
            hist_values, bin_edges = np.histogram(get_slices(df,i)[var],bins=bins)
            add_hist_values, add_bin_edges = np.histogram(get_slices(df_add,i)[var],bins=bins)
            if (i==0 and mult_factor!=None):
                plt.step(x=bin_edges[:-1], y=np.round((hist_values + nue_scale*add_hist_values)*mult_factor),where="post",linewidth=lw, label=entry + f" [x{mult_factor:.2f}]",zorder=5)
            elif (i==0 and mult_factor==None):
                plt.step(x=bin_edges[:-1], y=np.round((hist_values + nue_scale*add_hist_values)),where="post",linewidth=lw, label=entry,zorder=5)
            else:
                plt.step(x=bin_edges[:-1], y=np.round((hist_values + nue_scale*add_hist_values)),where="post",linewidth=lw, label=entry,zorder=5)

    # plot vertical lines
    ymin, ymax = plt.ylim()
    if cut_val != None:
        for i in range(len(cut_val)):
            plt.vlines(cut_val[i],ymin,ymax,lw=2,color="gray",linestyles="--",zorder=6)
    plt.ylim(ymin,ymax)   
    if label == "":
        plt.xlabel(var)
    if title == "":
        plt.title(var+" distribution (for each interaction type)")
    else:
        plt.title(title)
    plt.legend()
    plt.show()
    
def plot_var_pdg(df: pd.DataFrame, var: str, bins: np.ndarray,
             df_add: pd.DataFrame = None,
             label: str = "", title: str = "", 
             cut_val: list = None):
    """
    Plots a variable for each pdg in a stacked histogram.
    
    Parameters
    ----------
    df: pd.DataFrame
        input dataframe
    var: str
        variable to be plotted, must be a column in the dataframe
    bins: np.ndarray
        histogram binning
    df_add: pd.DataFrame
        if specified, adds this dataframe to the plot with scaling factor set by nue.nue_scale
    label: str
        x-axis label
    title: str
        title of the plot
    cut_val: list
        if specified, plots vertical lines at the specified values
    """
    histtype="step"
    lw = 2
    other_df = df.copy()
    add_other_df = df_add.copy() if df_add is not None else None
    for key, value in pdg_dict.items():
        # ignore neutrons and other (take care of "other" more explicitly at the end)
        if (key == "neu" or key == "other"):
            continue
        this_df = df[df.pfp_shw_truth_p_pdg==value]
        # if we have the additional dataframe, need to scale it 
        if (df_add is not None):
            add_df = df_add[df_add.pfp_shw_truth_p_pdg==value]
            hist_values, bin_edges = np.histogram(this_df[var],bins=bins)
            add_hist_values, add_bin_edges = np.histogram(add_df[var],bins=bins)
            plt.step(x=bin_edges[:-1], y=np.round((hist_values + nue_scale*add_hist_values)),where="post",linewidth=lw, label=key,zorder=5)
        # if we don't have the additional dataframe, just plot the histogram
        elif (df_add is None and len(this_df)>0):
            plt.hist(this_df[var],bins=bins,histtype=histtype,lw=lw,label=key)
        
        # if the pdg is not found within the dict, make sure to save it in "other_df"
        other_df = other_df[other_df.pfp_shw_truth_p_pdg!=value]
        if (df_add is not None):
            add_other_df = add_other_df[add_other_df.pfp_shw_truth_p_pdg!=value]    
    
    # if there are other pdgs
    if (len(other_df)>0):
        if (df_add is not None):
            hist_values, bin_edges = np.histogram(other_df[var],bins=bins)
            add_hist_values, add_bin_edges = np.histogram(add_other_df[var],bins=bins)
            plt.step(x=bin_edges[:-1], y=np.round((hist_values + nue_scale*add_hist_values)),where="post",linewidth=lw, label="other",zorder=5)
        else:
            plt.hist(other_df[var],bins=bins,histtype=histtype,lw=lw,label="other")

    ymin, ymax = plt.ylim()
    if cut_val != None:
        for i in range(len(cut_val)):
            plt.vlines(cut_val[i],ymin,ymax,lw=2,color="gray",linestyles="--",zorder=6)
    plt.ylim(ymin,ymax)
    
    if label == "":
        plt.xlabel(var)
    if title == "":
        plt.title(var+" distribution (separated by particle type)")      
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
    maskShower = (df.groupby(["ntuple","entry","rec.slc__index"])[col].transform(lambda x: (x < 0.6).any()))
    return df[maskShower]

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
    df["pfp_cont"] = True
    if whereTrkCont:
        df["pfp_cont"] = np.where(((df.pfp_trackScore >= 0.5) & 
                                   ((abs(df[trk_col+"x"]) > 195) |
                                     (df[trk_col+"y"] < -195)  | (df[trk_col+"y"] > 195) &
                                     (df[trk_col+"z"] < 5)    | (df[trk_col+"z"] > 495))),
                                     False,df["pfp_cont"])
    if whereShwCont: 
        df["pfp_cont"] = np.where(((df.pfp_trackScore < 0.5) & 
                                   ((abs(df[shw_col+"x"]) > 195) |
                                     (df[shw_col+"y"] < -195)  | (df[shw_col+"y"] > 195) &
                                     (df[shw_col+"z"] < 5)    | (df[shw_col+"z"] > 495))),
                                     False,df["pfp_cont"])
    maskPfpCont= ~(df.groupby(["ntuple","entry","rec.slc__index"])['pfp_cont'].transform(lambda x: (x==False).any()))
    return df[maskPfpCont]

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
    # get starting momentum
    flat_nuprim_df["prim_startp"] = np.sqrt(flat_nuprim_df.prim_startp_x**2 + flat_nuprim_df.prim_startp_y**2 + flat_nuprim_df.prim_startp_z**2)
    flat_nuprim_df["prim_abs_pdg"] = abs(flat_nuprim_df["prim_pdg"])
    flat_nuprim_df["prim_depE"] = flat_nuprim_df.prim_startE - flat_nuprim_df.prim_endE

    # set protons with momentum < 200 MeV as other
    # set charged pions with momentum < 100 as other
    flat_nuprim_df["prim_abs_pdg"] = np.where( ((flat_nuprim_df.prim_abs_pdg == 2212) & (flat_nuprim_df.prim_startp < 0.2) ), -1, flat_nuprim_df.prim_abs_pdg)
    flat_nuprim_df["prim_abs_pdg"] = np.where( ((flat_nuprim_df.prim_abs_pdg == 211) &  (flat_nuprim_df.prim_startp < 0.1) ), -1, flat_nuprim_df.prim_abs_pdg)
    # get whether the primaries are contained 
    flat_nuprim_df["prim_exit"] =  np.where((((abs(flat_nuprim_df["prim_end_x"]) > 195) |
                                              (abs(flat_nuprim_df["prim_end_y"]) > 195) &
                                              (flat_nuprim_df["prim_end_z"] < 5)    | (flat_nuprim_df["prim_end_z"] > 495))),
                                            True,False)
    
    # move columns of unspecified pdg into column "nother" and drop original column
    pdg_counts = (flat_nuprim_df.groupby(["ntuple","entry","rec.mc.nu__index"])["prim_abs_pdg"]).value_counts().unstack(fill_value=0).reset_index()
    pdg_counts["nother"] = 0
    
    # get total deposited energy for each neutrino event
    depE_sum =  flat_nuprim_df.groupby(["ntuple","entry","rec.mc.nu__index"])["prim_depE"].sum().reset_index().rename(columns={"prim_depE":"total_prim_depE"})
    
    # get whether the primaries are contained
    prim_cont = (flat_nuprim_df.groupby(["ntuple","entry","rec.mc.nu__index"])["prim_exit"]).sum().reset_index().rename(columns={"prim_exit":"prim_exit_count"})
    prim_cont["prim_cont"] = np.where(prim_cont["prim_exit_count"] == 0,True,False)
    
    # merge the dataframes
    pdg_counts = depE_sum.merge(pdg_counts,on=["ntuple","entry","rec.mc.nu__index"])
    pdg_counts = pdg_counts.merge(prim_cont,on=["ntuple","entry","rec.mc.nu__index"])
    for column in pdg_counts.copy().columns:
        if type(column)!=int:
            continue
        if (abs(column) in pdg_dict.values()) == False:
            pdg_counts["nother"] += pdg_counts[column]
            pdg_counts.drop(columns=[column],inplace=True)
    
    # rename columns from ints to strings
    for pdg_value in pdg_dict.values():
        for col in pdg_counts.columns:
            if col == pdg_value:
                new_col = list(pdg_dict.keys())[list(pdg_dict.values()).index(pdg_value)]
                pdg_counts.rename(columns={col:"n"+new_col},inplace=True)
                
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
    slcshw_df["nu_index"] = slcshw_df["slc_tmatch_idx"]
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
