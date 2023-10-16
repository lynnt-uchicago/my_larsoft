
from tqdm import tqdm
import fnmatch
import numpy as np 
from wc_helper import * 
import uproot

# with open("electron_caf.list") as f:
with open("proton_caf.list") as f:
    caf_list = f.read().split('\n')

with open("proton_wvfm.list") as f:
    wvfm_list = f.read().split('\n')

# all of these lists contain np.arrays() of a single event, where the plane_sim/dec array has the distributions of the sum of channels
# the theta and depE list have a np.array() of len 1 per event, whereas the plane_sim/dec arr have len=# of channels
# for idx in tqdm(10,range(len(wvfm_list)-1)):
for idx in tqdm(range(len(wvfm_list))):
    N = 20 # number of events
    file_idx = idx
    if file_idx==8: continue
        
    u_sim_sum_arr = np.zeros((N*2,idx_v0)); v_sim_sum_arr = np.zeros((N*2,idx_w0-idx_v0)); w_sim_sum_arr = np.zeros((N*2,idx_u1-idx_w0))
    u_dec_sum_arr = np.zeros((N*2,idx_v0)); v_dec_sum_arr = np.zeros((N*2,idx_w0-idx_v0)); w_dec_sum_arr = np.zeros((N*2,idx_u1-idx_w0))
    u_gho_sum_arr = np.zeros((N*2,idx_v0)); v_gho_sum_arr = np.zeros((N*2,idx_w0-idx_v0)); w_gho_sum_arr = np.zeros((N*2,idx_u1-idx_w0))
    u_sim_rej_arr = np.zeros((N*2,idx_v0)); v_sim_rej_arr = np.zeros((N*2,idx_w0-idx_v0)); w_sim_rej_arr = np.zeros((N*2,idx_u1-idx_w0))
    depE_sum_arr = np.zeros(N*2)
    theta_xz_sum_arr = np.zeros(N*2)
    
    wvfm = uproot.open(wvfm_list[file_idx])
    dec_names = fnmatch.filter(wvfm.keys(), '*run_*_sub_*_evt_*_decon*')
    sim_names = fnmatch.filter(wvfm.keys(), '*run_*_sub_*_evt_*_sim*')
    # raw_names = fnmatch.filter(wvfm.keys(), '*run_*_sub_*_evt_*_raw*')
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
            u_sim_fail = np.where(~u_sim_bool)[0]; v_sim_fail = np.where(~v_sim_bool)[0]; w_sim_fail = np.where(~w_sim_bool)[0]

            if len(u_sim_pass)==0 or len(v_sim_pass)==0 or len(w_sim_pass)==0: continue
            u_ch_width = np.max(u_sim_pass) - np.min(u_sim_pass); v_ch_width = np.max(v_sim_pass) - np.min(v_sim_pass); w_ch_width = np.max(w_sim_pass) - np.min(w_sim_pass)
            u_ch_sigma = round(u_ch_width*0.1); v_ch_sigma = round(v_ch_width*0.1); w_ch_sigma = round(w_ch_width*0.1)
            
            u_ch_min = np.min(u_sim_pass) - u_ch_sigma; u_ch_max = np.max(u_sim_pass) + u_ch_sigma
            v_ch_min = np.min(v_sim_pass) - v_ch_sigma; v_ch_max = np.max(v_sim_pass) + v_ch_sigma
            w_ch_min = np.min(w_sim_pass) - w_ch_sigma; w_ch_max = np.max(w_sim_pass) + w_ch_sigma

            # set everything but ghosts equal to zero in the gho arr
            u_gho_sum = u_dec_sum.copy(); v_gho_sum = v_dec_sum.copy(); w_gho_sum = w_dec_sum.copy()
            u_gho_sum[u_ch_min:u_ch_max] = 0; v_gho_sum[v_ch_min:v_ch_max] = 0; w_gho_sum[w_ch_min:w_ch_max] = 0
                        
            # set ghost to zero in the dec arr 
            u_dec_sum[:u_ch_min] = 0; u_dec_sum[u_ch_max:] = 0
            v_dec_sum[:v_ch_min] = 0; v_dec_sum[v_ch_max:] = 0
            w_dec_sum[:w_ch_min] = 0; w_dec_sum[w_ch_max:] = 0

            # create separate arr for sim channels that didn't pass 
            u_sim_rej = np.where(u_sim_bool==False, u_sim_sum, 0); 
            v_sim_rej = np.where(v_sim_bool==False, v_sim_sum, 0); 
            w_sim_rej = np.where(w_sim_bool==False, w_sim_sum, 0)
                        
            # only include sim channels that pass
            u_sim_sum[u_sim_fail] = 0
            v_sim_sum[v_sim_fail] = 0
            w_sim_sum[w_sim_fail] = 0
                
            u_sim_sum_arr[event*2 + tpc] = u_sim_sum; v_sim_sum_arr[event*2 + tpc] = v_sim_sum; w_sim_sum_arr[event*2 + tpc] = w_sim_sum
            u_sim_rej_arr[event*2 + tpc] = u_sim_rej; v_sim_rej_arr[event*2 + tpc] = v_sim_rej; w_sim_rej_arr[event*2 + tpc] = w_sim_rej
            u_dec_sum_arr[event*2 + tpc] = u_dec_sum; v_dec_sum_arr[event*2 + tpc] = v_dec_sum; w_dec_sum_arr[event*2 + tpc] = w_dec_sum
            u_gho_sum_arr[event*2 + tpc] = u_gho_sum; v_gho_sum_arr[event*2 + tpc] = v_gho_sum; w_gho_sum_arr[event*2 + tpc] = w_gho_sum
            depE_sum_arr[ event*2 + tpc] = this_tree.iloc[event]["depE"]
            theta_xz_sum_arr[event*2 + tpc] = this_tree.iloc[event]["theta_xz"]
            
    np_file_name = "wvfm_arrs_"+str(file_idx)+".npz"
    np.savez("pro_npz_files/"+np_file_name, u_sim_sum_arr=u_sim_sum_arr, v_sim_sum_arr=v_sim_sum_arr, w_sim_sum_arr=w_sim_sum_arr,
                           u_sim_rej_arr=u_sim_rej_arr, v_sim_rej_arr=v_sim_rej_arr, w_sim_rej_arr=w_sim_rej_arr,
                           u_dec_sum_arr=u_dec_sum_arr, v_dec_sum_arr=v_dec_sum_arr, w_dec_sum_arr=w_dec_sum_arr,
                           u_gho_sum_arr=u_gho_sum_arr, v_gho_sum_arr=v_gho_sum_arr, w_gho_sum_arr=w_gho_sum_arr,
                           theta_xz_sum_arr=theta_xz_sum_arr, depE_sum_arr=depE_sum_arr)