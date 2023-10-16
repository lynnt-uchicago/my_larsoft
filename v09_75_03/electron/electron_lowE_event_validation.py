import uproot
import numpy as np
import fnmatch
from tqdm import tqdm
import pandas as pd
from wc_helper import *

start_idx = 165

# with open("electron_caf.list") as f:
with open("electron_caf_edited.list") as f:
    caf_list = f.read().split('\n')

with open("electron_wvfm.list") as f:
    wvfm_list = f.read().split('\n')

print("Checking file access...")    
# check that we can open all the files 
for i in tqdm(range(start_idx,len(wvfm_list)-1)):
    if i == 8: continue
    uproot.open(wvfm_list[i])
    uproot.open(caf_list[i])
print("All files accessible!")


for file_idx in tqdm(range(start_idx,len(wvfm_list)-1)):
    N = 20 # number of events
    if file_idx == 8: continue
        
    file_arr =  np.zeros(shape=(N,2)); subrun_arr = np.zeros(shape=(N,2)); tpc_arr =    np.zeros(shape=(N,2))
    u_diff_arr =np.zeros(shape=(N,2)); v_diff_arr = np.zeros(shape=(N,2)); w_diff_arr = np.zeros(shape=(N,2))
    u_sim_arr = np.zeros(shape=(N,2)); v_sim_arr =  np.zeros(shape=(N,2)); w_sim_arr =  np.zeros(shape=(N,2))
    u_dec_arr = np.zeros(shape=(N,2)); v_dec_arr =  np.zeros(shape=(N,2)); w_dec_arr =  np.zeros(shape=(N,2)) 
    theta_arr = np.zeros(shape=(N,2)); depE_arr =   np.zeros(shape=(N,2))
    
    wvfm = uproot.open(wvfm_list[file_idx])
    dec_names = fnmatch.filter(wvfm.keys(), '*run_*_sub_*_evt_*_decon*')
    sim_names = fnmatch.filter(wvfm.keys(), '*run_*_sub_*_evt_*_sim*')
    caf = uproot.open(caf_list[file_idx]+":recTree")
    tree = caf.arrays(["rec.hdr.run","rec.hdr.subrun","rec.hdr.evt",
                    "rec.mc.nu.prim.startp.x" ,"rec.mc.nu.prim.startp.y","rec.mc.nu.prim.startp.z",
                    "rec.mc.nu.prim.startE","rec.mc.nu.prim.endE"],library='pd')

    tree["theta_xz"] = np.arctan(abs(tree["rec.mc.nu.prim.startp.x"]/tree["rec.mc.nu.prim.startp.z"]))*(180/np.pi)
    tree["depE"] =  tree["rec.mc.nu.prim.startE"] - tree["rec.mc.nu.prim.endE"]
    tree["tpc"] = np.where(tree["rec.mc.nu.prim.startp.x"]>0,1,0)
    subrun = tree.iloc[0]["rec.hdr.subrun"]
    for event in range(len(dec_names)):
        evt_sim = wvfm[sim_names[event]].values()
        evt_dec = wvfm[dec_names[event]].values()*50
        
        for tpc in range(2):
            this_tree = tree[tree.tpc==0] if tpc == 0 else tree[tree.tpc==1]
            u_sim = u0_ch(evt_sim) if tpc == 0 else u1_ch(evt_sim); u_dec = u0_ch(evt_dec) if tpc == 0 else u1_ch(evt_dec)
            v_sim = v0_ch(evt_sim) if tpc == 0 else v1_ch(evt_sim); v_dec = v0_ch(evt_dec) if tpc == 0 else v1_ch(evt_dec)
            w_sim = w0_ch(evt_sim) if tpc == 0 else w1_ch(evt_sim); w_dec = w0_ch(evt_dec) if tpc == 0 else w1_ch(evt_dec)
            
            u_sim_sum = np.sum(u_sim, axis=1); v_sim_sum = np.sum(v_sim, axis=1); w_sim_sum = np.sum(w_sim, axis=1)
            u_dec_sum = np.sum(u_dec, axis=1); v_dec_sum = np.sum(v_dec, axis=1); w_dec_sum = np.sum(w_dec, axis=1)
            
            u_sim_bool = pass_sim(u_sim); v_sim_bool = pass_sim(v_sim); w_sim_bool = pass_sim(w_sim)
            u_sim_pass = np.where( u_sim_bool)[0]; v_sim_pass = np.where( v_sim_bool)[0]; w_sim_pass = np.where( w_sim_bool)[0]

            u_ch_width = np.max(u_sim_pass) - np.min(u_sim_pass); v_ch_width = np.max(v_sim_pass) - np.min(v_sim_pass); w_ch_width = np.max(w_sim_pass) - np.min(w_sim_pass)            
            u_ch_sigma = round(u_ch_width*0.1); v_ch_sigma = round(v_ch_width*0.1); w_ch_sigma = round(w_ch_width*0.1)
            
            u_ch_min = np.min(u_sim_pass) - u_ch_sigma; u_ch_max = np.max(u_sim_pass) + u_ch_sigma
            v_ch_min = np.min(v_sim_pass) - v_ch_sigma; v_ch_max = np.max(v_sim_pass) + v_ch_sigma
            w_ch_min = np.min(w_sim_pass) - w_ch_sigma; w_ch_max = np.max(w_sim_pass) + w_ch_sigma
            
            # the sum we use for dec over the entire event is corrected by only using channels that wouldn't have ghosts in them
            u_dec_evt = np.sum(u_dec_sum[u_ch_min:u_ch_max]); v_dec_evt = np.sum(v_dec_sum[v_ch_min:v_ch_max]); w_dec_evt = np.sum(w_dec_sum[w_ch_min:w_ch_max])
            
            # the sum we use for sim over the entire event should not include sim channels that do not pass 
            u_sim_evt = np.sum(u_sim_sum[u_sim_bool]); v_sim_evt = np.sum(v_sim_sum[v_sim_bool]); w_sim_evt = np.sum(w_sim_sum[w_sim_bool])
            
            u_diff_arr[event][tpc] = (u_dec_evt - u_sim_evt)/u_dec_evt
            v_diff_arr[event][tpc] = (v_dec_evt - v_sim_evt)/v_dec_evt
            w_diff_arr[event][tpc] = (w_dec_evt - w_sim_evt)/w_dec_evt
            u_sim_arr[event][tpc] = u_sim_evt; v_sim_arr[event][tpc] = v_sim_evt; w_sim_arr[event][tpc] = w_sim_evt
            u_dec_arr[event][tpc] = u_dec_evt; v_dec_arr[event][tpc] = v_dec_evt; w_dec_arr[event][tpc] = w_dec_evt
            theta_arr[event][tpc]  = this_tree.iloc[event]["theta_xz"]
            depE_arr[event][tpc]   = this_tree.iloc[event]["depE"]
            tpc_arr[event][tpc]    = tpc
            file_arr[event][tpc]   = file_idx
            subrun_arr[event][tpc] = subrun

    df0 = pd.DataFrame({'file' : file_arr[:,0], 'sub' : subrun_arr[:,0], 'tpc' : tpc_arr[:,0],'u': u_diff_arr[:,0], 'v': v_diff_arr[:,0], 'w': w_diff_arr[:,0],
                        'u_sim' : u_sim_arr[:,0], 'v_sim' : v_sim_arr[:,0], 'w_sim' : w_sim_arr[:,0],
                        'u_dec' : u_dec_arr[:,0], 'v_dec' : v_dec_arr[:,0], 'w_dec' : w_dec_arr[:,0],
                        "theta_xz" : theta_arr[:,0], 'depE' : depE_arr[:,0]})
    df1 = pd.DataFrame({'file' : file_arr[:,1], 'sub' : subrun_arr[:,1], 'tpc' : tpc_arr[:,1],'u': u_diff_arr[:,1], 'v': v_diff_arr[:,1], 'w': w_diff_arr[:,1],
                        'u_sim' : u_sim_arr[:,1], 'v_sim' : v_sim_arr[:,1], 'w_sim' : w_sim_arr[:,1],
                        'u_dec' : u_dec_arr[:,1], 'v_dec' : v_dec_arr[:,1], 'w_dec' : w_dec_arr[:,1],
                        "theta_xz" : theta_arr[:,1], 'depE' : depE_arr[:,1]})
    df = pd.concat([df0,df1])
    df.to_pickle("pkl_lowE_event_files/file_"+str(file_idx)+".pkl")  