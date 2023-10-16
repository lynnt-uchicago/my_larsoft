import uproot
import matplotlib.pyplot as plt 
import mplhep as hep
import numpy as np
import fnmatch
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from tqdm import tqdm
import pickle as pkl
import pandas as pd
import glob as glob

idx_u0 = 0 
idx_v0 = 1984
idx_w0 = 3968
idx_u1 = 5632
idx_v1 = 7616
idx_w1 = 9000

def find_u0_ch(ch_num):
    if (((ch_num >= 0 ) & (ch_num < idx_v0))):
        return True
    else:
        return False

def find_v0_ch(ch_num):
    if (((ch_num >= idx_v0 ) & (ch_num < idx_w0))):
        return True
    else:
        return False
    
def find_w0_ch(ch_num):
    if (((ch_num >= idx_w0 ) & (ch_num < idx_u1))):
        return True
    else:
        return False
    
def find_u1_ch(ch_num):
    if (((ch_num >=  idx_u1 ) & (ch_num < idx_v1))):
        return True
    else:
        return False

def find_v1_ch(ch_num):
    if (((ch_num >= idx_v1 ) & (ch_num < idx_w1))):
        return True
    else:
        return False
    
def find_w1_ch(ch_num):
    if (((ch_num >= idx_w1 ))):
        return True
    else:
        return False
    
def u0_ch(input): 
    tpc0_arr = input[0     :idx_v0]
    return tpc0_arr

def v0_ch(input): 
    tpc0_arr = input[idx_v0:idx_w0]
    return tpc0_arr

def w0_ch(input): 
    tpc0_arr = input[idx_w0:idx_u1]
    return tpc0_arr

def u1_ch(input): 
    tpc1_arr = input[idx_u1:idx_v1]
    return tpc1_arr

def v1_ch(input): 
    tpc1_arr = input[idx_v1:idx_w1]
    return tpc1_arr

def w1_ch(input): 
    tpc1_arr = input[idx_w1:]
    return tpc1_arr

def caf_split_tpc(df):
    tpc0 = df[df["rec.mc.nu.prim.start.x"] < 0]
    tpc1 = df[df["rec.mc.nu.prim.start.x"] > 0]
    return tpc0, tpc1

def find_nonzero_ch(wvfm):
    # peak threshold = 800
    mask = (wvfm > 800).sum(axis = 1)
    return np.where(mask>0)[0]
    
with open("proton_caf.list") as f:
    caf_list = f.read().split('\n')

with open("proton_wvfm.list") as f:
    wvfm_list = f.read().split('\n')
    
u_sim_arr = np.array([]); u_dec_arr = np.array([]); u_diff_arr = np.array([]); u_theta_arr = np.array([]); u_depE_arr = np.array([])
v_sim_arr = np.array([]); v_dec_arr = np.array([]); v_diff_arr = np.array([]); v_theta_arr = np.array([]); v_depE_arr = np.array([])
w_sim_arr = np.array([]); w_dec_arr = np.array([]); w_diff_arr = np.array([]); w_theta_arr = np.array([]); w_depE_arr = np.array([])

nfiles = 50

for i in tqdm(range(nfiles)):
    wvfm = uproot.open(wvfm_list[i])
    dec_names = fnmatch.filter(wvfm.keys(), '*run_*_sub_*_evt_*_decon*')
    sim_names = fnmatch.filter(wvfm.keys(), '*run_*_sub_*_evt_*_sim*')

    # caf = uproot.open("/pnfs/sbnd/scratch/users/lynnt/v09_75_03_02/electron_gun/electron_caf/3414904_0/prodsingle_sbnd_SinglesGen-20230812T042829_G4-20230812T043429_WCLS-20230812T050513_405fd978-8f82-4ce9-9dc6-70c89057648f.flat.caf.root:recTree")
    caf = uproot.open(caf_list[i]+":recTree")
    tree = caf.arrays(["rec.hdr.run",
                            "rec.hdr.subrun",
                            "rec.hdr.evt",    
                            "rec.mc.nu.prim.startp.x" ,
                            "rec.mc.nu.prim.startp.y",
                            "rec.mc.nu.prim.startp.z",
                            "rec.mc.nu.prim.startE",
                            "rec.mc.nu.prim.endE"],library='pd')

    tree["theta_xz"] = np.arctan(abs(tree["rec.mc.nu.prim.startp.x"]/tree["rec.mc.nu.prim.startp.z"]))*(180/np.pi)
    tree["depE"] =  tree["rec.mc.nu.prim.startE"] - tree["rec.mc.nu.prim.endE"]
    subrun = tree.iloc[0]["rec.hdr.subrun"]
    
    for event in range(len(dec_names)):
        evt_sim = wvfm[sim_names[event]].values()
        evt_dec = wvfm[dec_names[event]].values()*50
        nonZeroCh = find_nonzero_ch(evt_sim)
        
        theta0= tree.iloc[2*event+1]["theta_xz"]
        theta1= tree.iloc[2*event]["theta_xz"]
        depE0= tree.iloc[2*event+1]["depE"]
        depE1= tree.iloc[2*event]["depE"]
        
        for idx, j in enumerate(nonZeroCh):
            proton_wvfm_dec = wvfm[dec_names[event]].values()[j]*50 # have to re-add the scaling factor
            proton_wvfm_sim = wvfm[sim_names[event]].values()[j]
            # # sim_out = find_peaks(proton_wvfm_sim,height=15,threshold=0,prominence=0,width=5)
            # sim_out = find_peaks(proton_wvfm_sim,height=10,threshold=0,prominence=0,width=0)
            # right_side = 0
            # left_side = 3400
            # for peak in range(len(sim_out[0])):
            #         if sim_out[1]["right_bases"][peak] > right_side:
            #             right_side = sim_out[1]["right_bases"][peak]
            #         if sim_out[1]["left_bases"][peak] < left_side:
            #             left_side = sim_out[1]["left_bases"][peak]
            if ((find_u0_ch(j)) | (find_u1_ch(j))):
                u_sim_arr = np.append(u_sim_arr,np.sum(proton_wvfm_sim))
                u_dec_arr = np.append(u_dec_arr,np.sum(proton_wvfm_dec))
                if (find_u0_ch(j)):
                    u_theta_arr = np.append(u_theta_arr,theta0)
                    u_depE_arr = np.append(u_depE_arr,depE0)
                elif (find_u1_ch(j)):
                    u_theta_arr = np.append(u_theta_arr,theta1)
                    u_depE_arr = np.append(u_depE_arr,depE1)
            if ((find_v0_ch(j)) | (find_v1_ch(j))):
                v_sim_arr = np.append(v_sim_arr,np.sum(proton_wvfm_sim))
                v_dec_arr = np.append(v_dec_arr,np.sum(proton_wvfm_dec))
                if (find_v0_ch(j)):
                    v_theta_arr = np.append(v_theta_arr,theta0)
                    v_depE_arr = np.append(v_depE_arr,depE0)
                elif (find_v1_ch(j)):
                    v_theta_arr = np.append(v_theta_arr,theta1)
                    v_depE_arr = np.append(v_depE_arr,depE1)
            if ((find_w0_ch(j)) | (find_w1_ch(j))):
                w_sim_arr = np.append(w_sim_arr,np.sum(proton_wvfm_sim))
                w_dec_arr = np.append(w_dec_arr,np.sum(proton_wvfm_dec))
                if (find_w0_ch(j)):
                    w_theta_arr = np.append(w_theta_arr,theta0)
                    w_depE_arr = np.append(w_depE_arr,depE0)
                elif (find_w1_ch(j)):
                    w_theta_arr = np.append(w_theta_arr,theta1)
                    w_depE_arr = np.append(w_depE_arr,depE1)
    u_diff_arr = (u_dec_arr - u_sim_arr)/u_sim_arr
    v_diff_arr = (v_dec_arr - v_sim_arr)/v_sim_arr
    w_diff_arr = (w_dec_arr - w_sim_arr)/w_sim_arr
    
    np.savez("npz_files/u_plane/u_subrun_"+str(int(subrun))+".npz",u_sim_arr=u_sim_arr,u_dec_arr=u_dec_arr,u_diff_arr=u_diff_arr,u_theta_arr=u_theta_arr,u_depE_arr=u_depE_arr,)
    np.savez("npz_files/v_plane/v_subrun_"+str(int(subrun))+".npz",v_sim_arr=v_sim_arr,v_dec_arr=v_dec_arr,v_diff_arr=v_diff_arr,v_theta_arr=v_theta_arr,v_depE_arr=v_depE_arr,)
    np.savez("npz_files/w_plane/w_subrun_"+str(int(subrun))+".npz",w_sim_arr=w_sim_arr,w_dec_arr=w_dec_arr,w_diff_arr=w_diff_arr,w_theta_arr=w_theta_arr,w_depE_arr=w_depE_arr,)
    # clear arrays
    u_sim_arr = np.array([]); u_dec_arr = np.array([]); u_diff_arr = np.array([]); u_theta_arr = np.array([]); u_depE_arr = np.array([])
    v_sim_arr = np.array([]); v_dec_arr = np.array([]); v_diff_arr = np.array([]); v_theta_arr = np.array([]); v_depE_arr = np.array([])
    w_sim_arr = np.array([]); w_dec_arr = np.array([]); w_diff_arr = np.array([]); w_theta_arr = np.array([]); w_depE_arr = np.array([])