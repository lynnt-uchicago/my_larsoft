import numpy as np 
import pandas as pd
from scipy.signal import find_peaks
import mplhep as hep
import matplotlib.pyplot as plt

idx_u0 = 0 
idx_v0 = 1984
idx_w0 = 3968
idx_u1 = 5632
idx_v1 = 7616
idx_w1 = 9600
idx_u2 = 11264


def u_ch(input: np.ndarray) -> np.ndarray:
    """
    Returns a np.array of channels that are in the U plane (both TPCs)
    
    Parameters
    ----------
    input : np.ndarray
        Array of channels, should be **per event**, have total shape (3400,11264)
        
    Notes:
    Indexed from 0 to 1984, and 5632 to 7616.
    """
    tpc0_arr = input[0     :idx_v0]
    tpc1_arr = input[idx_u1:idx_v1]
    return np.concatenate([tpc0_arr,tpc1_arr])

def v_ch(input: np.ndarray) -> np.ndarray: 
    """
    Returns a np.array of channels that are in the V plane (both TPCs)
    
    Notes:
    Indexed from 1984 to 3968, and 7616 to 9600.
    """
    tpc0_arr = input[idx_v0:idx_w0]
    tpc1_arr = input[idx_v1:idx_w1]
    return np.concatenate([tpc0_arr,tpc1_arr])

def w_ch(input:np.ndarray) -> np.ndarray: 
    """
    Returns a np.array of channels that are in the W plane (both TPCs)
    
    Notes
    ------
    Indexed from 3968 to 5632, and 9600 to 11264.
    """
    tpc0_arr = input[idx_w0:idx_u1]
    tpc1_arr = input[idx_w1:-1]
    return np.concatenate([tpc0_arr,tpc1_arr])

def u0_ch(input:np.ndarray) -> np.ndarray: 
    """
    Returns a np.array of channels that are in the U plane (TPC0)
    
    Notes:
    Indexed from 0 to 1984.
    """
    tpc0_arr = input[0     :idx_v0]
    return tpc0_arr

def v0_ch(input:np.ndarray) -> np.ndarray: 
    """
    Returns a np.array of channels that are in the V plane (TPC0)
    
    Notes:
    Indexed from 1984 to 3968.
    """
    tpc0_arr = input[idx_v0:idx_w0]
    return tpc0_arr

def w0_ch(input:np.ndarray) -> np.ndarray:
    """
    Returns a np.array of channels that are in the W plane (TPC0)
    
    Notes:
    Indexed from 3968 to 5632.
    """
    tpc0_arr = input[idx_w0:idx_u1]
    return tpc0_arr

def u1_ch(input:np.ndarray) -> np.ndarray:
    """
    Returns a np.array of channels that are in the U plane (TPC1)
    
    Notes:
    Indexed from 5632 to 7616.
    """
    tpc1_arr = input[idx_u1:idx_v1]
    return tpc1_arr

def v1_ch(input:np.ndarray) -> np.ndarray:
    """
    Returns a np.array of channels that are in the V plane (TPC1)
    
    Notes:
    Indexed from 7616 to 9600.
    """
    tpc1_arr = input[idx_v1:idx_w1]
    return tpc1_arr

def w1_ch(input:np.ndarray) -> np.ndarray:
    """
    Returns a np.array of channels that are in the W plane (TPC1)
    
    Notes:
    Indexed from 9600 to 11264.
    """
    tpc1_arr = input[idx_w1:-1]
    return tpc1_arr

def pass_sim(sim_arr:np.ndarray, height_thresh:float = 800, avg_thresh:float = 400) -> np.ndarray:
    """
    Returns a mask of simulated channels with 800 e- peak height **and** average charge above 400 e-
    
    Input
    -----
    sim_arr : np.ndarray
        Array of simulated waveforms, should be **per event**
    height_thresh : float
        Threshold for peak height, defaults to 800 e-
    avg_thresh : float
        Threshold for average charge, defaults to 400 e-/tick
    
    Returns
    -------
    ch_bool : np.array
        Boolean array of channels that pass the cut, aka a mask
    """
    sim_ch = np.where(np.sum(sim_arr,axis=1)>=800)[0]
    pass_sim = []
    for ch in sim_ch:
        sim_peak, sim_prop = find_peaks(sim_arr[ch],height=5,width=0)
        peak_idx = np.argmax(sim_prop["peak_heights"])
        peak_height = sim_prop["peak_heights"][peak_idx]
        peak_charge = sim_arr[ch][sim_prop["left_bases"][peak_idx]:sim_prop["right_bases"][peak_idx]].sum()
        peak_width =  sim_prop["right_bases"][peak_idx]-sim_prop["left_bases"][peak_idx]
        peak_avg = peak_charge / peak_width
        if (peak_height >= height_thresh) and (peak_avg >= avg_thresh): pass_sim.append(ch)
    pass_sim = np.array(pass_sim,dtype=int)
    ch_bool = np.zeros(len(sim_arr),dtype=bool)    
    ch_bool[sim_ch] = True
    return ch_bool

def quartile_reso(data) -> np.ndarray:
    """
    Returns the quartile resolution of a dataset. Rounds the input array to three decimal points to find the median. 
    """
    data = np.round(np.sort(data),3)
    median = np.ma.median(np.ma.array(data, mask=np.isnan(data)))
    median_idx = int(np.round(len(data)/2))
    quar_size =  int(np.round(len(data)*0.5*0.683))
    Q1 = data[median_idx - quar_size]
    Q3 = data[median_idx + quar_size]
    quartile_reso = np.sqrt(0.5*(median-Q1)**2 + 0.5*(median-Q3)**2)
    return quartile_reso


def plot_wvfms(chnum:int, wvfm, sim_names, raw_names, dec_names, evtnum, 
               range: list = [0,3400], textxcoord: float = 0.5):
    """
    Plots waveform of the specified channel number, as well as the neighborhoring two channels.
    
    Needs global variables sim_names, raw_names, dec_names, wvfm.
    """
    fig, axes = plt.subplots(3,1,figsize=(10,10))
    for i, shift in enumerate([-1,0,1]):
        hep.histplot(wvfm[sim_names[evtnum]].values()[chnum+shift]   ,ax=axes[i], label="sim")
        hep.histplot(wvfm[dec_names[evtnum]].values()[chnum+shift]*50,ax=axes[i], label="dec")
        hep.histplot(wvfm[raw_names[evtnum]].values()[chnum+shift]*20,ax=axes[i], label="raw [x20]")
        axes[i].set_title("Channel "+str(chnum+shift))
        axes[i].set_xlim(range)
        axes[i].legend()
        sim_ch_sum = np.sum(wvfm[sim_names[evtnum]].values()[chnum+shift][range[0]:range[1]])
        dec_ch_sum = np.sum(wvfm[dec_names[evtnum]].values()[chnum+shift][range[0]:range[1]]*50)
        diff_ch = (dec_ch_sum - sim_ch_sum)/sim_ch_sum
        axes[i].annotate(f"sim integral: {sim_ch_sum:.0f}",     (textxcoord,0.8) ,xycoords='axes fraction',fontsize=12)
        axes[i].annotate(f"dec integral: {dec_ch_sum:.0f}",     (textxcoord,0.7),xycoords='axes fraction',fontsize=12)
        axes[i].annotate(f"fractional diff: {diff_ch*100:.1f}%",(textxcoord,0.6),xycoords='axes fraction',fontsize=12)
    plt.show()
    
def plot_evt(ch_min, ch_max, sim, raw, dec, 
             title:str=None,ymin:float=0,ymax:float=3400):
    mask = np.where((abs(np.ceil(sim[1]))[:-1]<ch_max) & (abs(np.ceil(sim[1]))[:-1] >= ch_min), True,False)
    xbins = sim[1][(sim[1]<ch_max) & (sim[1] >= ch_min)]
    ybins = sim[2]

    fig, axes = plt.subplots(1,3, figsize=(18,4),sharey=True,sharex=True)
    hep.hist2dplot(H=sim[0][:][mask],xbins=xbins,ybins=ybins,shading="auto",ax=axes[0],cmap="seismic",vmin=-2e4,vmax=2e4)
    hep.hist2dplot(H=raw[0][:][mask],xbins=xbins,ybins=ybins,shading="auto",ax=axes[1],cmap="seismic",vmin=-100,vmax=100,)
    hep.hist2dplot(H=dec[0][:][mask],xbins=xbins,ybins=ybins,shading="auto",ax=axes[2],cmap="seismic",vmin=-5e2,vmax=5e2)
    
    axes[0].set_ylabel("ticks")
    axes[0].set_xlabel("Channel Number");   axes[1].set_xlabel("Channel Number"); axes[2].set_xlabel("Channel Number")
    axes[0].set_title("Simulated (Truth)"); axes[1].set_title("Raw Waveform");    axes[2].set_title("2D Dec+SP")

    plt.xlim(ch_min,ch_max)
    plt.ylim(ymin,ymax)
    plt.suptitle(title,y=1.05,fontsize=15)
    plt.show()
    
def prime_angle(angle):
    """
    Returns the angle (in degrees) in the induction plane reference frame.
    
    Input
    ----------
    angle : float or np.ndarray
        Theta_xz angle in degrees, measured from the z-axis. 
        
    Returns
    -------
    angle_prime : float or np.ndarray
        Induction plane theta_xz angle in degrees, measured from the z-axis.
        Rounded to two decimal points.
        
    """
    return np.round(np.arctan(np.tan(angle*np.pi/180)/np.cos(60*np.pi/180))*180/np.pi,2)
