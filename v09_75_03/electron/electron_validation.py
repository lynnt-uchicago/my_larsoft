import uproot
import mplhep as hep
import numpy as np
import fnmatch
from tqdm import tqdm
import pickle as pkl
import pandas as pd

idx_u0 = 0 
idx_v0 = 1984
idx_w0 = 3968
idx_u1 = 5632
idx_v1 = 7616
idx_w1 = 9000

def u_ch(input): 
    tpc0_arr = input[0     :idx_v0]
    tpc1_arr = input[idx_u1:idx_v1]
    return np.concatenate([tpc0_arr,tpc1_arr])

def v_ch(input): 
    tpc0_arr = input[idx_v0:idx_w0]
    tpc1_arr = input[idx_v1:idx_w1]
    return np.concatenate([tpc0_arr,tpc1_arr])

def w_ch(input): 
    tpc0_arr = input[idx_w0:idx_u1]
    tpc1_arr = input[idx_w1:]
    return np.concatenate([tpc0_arr,tpc1_arr])

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

with open("electron_caf.list") as f:
    caf_list = f.read().split('\n')

with open("electron_wvfm.list") as f:
    wvfm_list = f.read().split('\n')
    
nfiles = 250

N = 20*nfiles*2 # number of events * number of files * 2 particles per event
u_diff_arr = np.zeros(N); v_diff_arr = np.zeros(N); w_diff_arr = np.zeros(N)
u_sim_arr = np.zeros(N); v_sim_arr = np.zeros(N); w_sim_arr = np.zeros(N)
u_dec_arr = np.zeros(N); v_dec_arr = np.zeros(N); w_dec_arr = np.zeros(N)

theta_xz_arr = np.zeros(N)
depE_arr = np.zeros(N)
subrun_arr = np.zeros(N)

for i in tqdm(range(50,nfiles)):
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
    

    for event in range(len(dec_names)):
        evt_sim = wvfm[sim_names[event]].values()
        u0_sim = np.sum(u0_ch(evt_sim)); v0_sim = np.sum(v0_ch(evt_sim)); w0_sim = np.sum(w0_ch(evt_sim))
        u1_sim = np.sum(u1_ch(evt_sim)); v1_sim = np.sum(v1_ch(evt_sim)); w1_sim = np.sum(w1_ch(evt_sim))

        evt_dec = wvfm[dec_names[event]].values()*50
        u0_dec = np.sum(u0_ch(evt_dec)); v0_dec = np.sum(v0_ch(evt_dec)); w0_dec = np.sum(w0_ch(evt_dec))
        u1_dec = np.sum(u1_ch(evt_dec)); v1_dec = np.sum(v1_ch(evt_dec)); w1_dec = np.sum(w1_ch(evt_dec))

        # this_df = tree["rec.hdr.evt" == event]
        # tpc0_df, tpc1_df = caf_split_tpc(this_df)
        
        # tpc 0
        u_diff_arr[i*40 + 2*event+1]   = (u0_dec - u0_sim)/u0_sim
        v_diff_arr[i*40 + 2*event+1]   = (v0_dec - v0_sim)/v0_sim
        w_diff_arr[i*40 + 2*event+1]   = (w0_dec - w0_sim)/w0_sim
        theta_xz_arr[i*40 + 2*event+1] = tree.iloc[2*event+1]["theta_xz"]
        depE_arr[i*40 + 2*event+1] = tree.iloc[2*event+1]["depE"]
        u_sim_arr[i*40 + 2*event+1] = u0_sim
        v_sim_arr[i*40 + 2*event+1] = v0_sim
        w_sim_arr[i*40 + 2*event+1] = w0_sim
        
        u_dec_arr[i*40 + 2*event+1] = u0_dec
        v_dec_arr[i*40 + 2*event+1] = v0_dec
        w_dec_arr[i*40 + 2*event+1] = w0_dec

        # tpc 1
        u_diff_arr[i*40 + 2*event]   = (u1_dec - u1_sim)/u1_sim
        v_diff_arr[i*40 + 2*event]   = (v1_dec - v1_sim)/v1_sim
        w_diff_arr[i*40 + 2*event]   = (w1_dec - w1_sim)/w1_sim
        theta_xz_arr[i*40 + 2*event] = tree.iloc[2*event]["theta_xz"]
        depE_arr[i*40 + 2*event] = tree.iloc[2*event]["depE"]
        u_sim_arr[i*40 + 2*event] = u1_sim
        v_sim_arr[i*40 + 2*event] = v1_sim
        w_sim_arr[i*40 + 2*event] = w1_sim
        
        u_dec_arr[i*40 + 2*event] = u1_dec
        v_dec_arr[i*40 + 2*event] = v1_dec
        w_dec_arr[i*40 + 2*event] = w1_dec
        
        subrun = tree.iloc[0]["rec.hdr.subrun"]
        subrun_arr[i*40 + 2*event] = tree.iloc[0]["rec.hdr.subrun"]
        subrun_arr[i*40 + 2*event+1] = tree.iloc[0]["rec.hdr.subrun"]
        
    df = pd.DataFrame({'sub' : subrun_arr,'u': u_diff_arr, 'v': v_diff_arr, 'w': w_diff_arr,
                       'u_sim': u_sim_arr, 'v_sim': v_sim_arr, 'w_sim': w_sim_arr,
                       'u_dec': u_dec_arr, 'v_dec': v_dec_arr, 'w_dec': w_dec_arr,
                       "theta_xz" : theta_xz_arr, 'depE' : depE_arr})
    df.to_pickle("pkl_files/subrun_"+str(int(subrun))+".pkl")    