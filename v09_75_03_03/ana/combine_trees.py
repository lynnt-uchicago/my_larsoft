import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import nue as nue
import importlib
____ = importlib.reload(nue)
import sys

nu_idx_set = ["ntuple","entry","nu_index"]
slc_idx_set = ["ntuple","entry","rec.slc__index"]

def main(argv):
    hist_path = argv[0]
    file = argv[1]
    
    sample = "sample"
    if "rockbox" in file: 
        sample = "rockbox"
    elif "intrnue" in file:
        sample = "intrnue"
    else: 
        print("Sample type unknown, output files will have name ""sample""")
    print("loading arrays...")

    razzled_col = ["electronScore","muonScore","photonScore","pionScore","protonScore","otherScore","bestScore","bestPDG","truePDG",
               "dazzleMuonScore","dazzlePionScore","dazzleProtonScore","dazzleOtherScore","dazzlePDG","trk_length",
               "razzleElectronScore","razzlePhotonScore","razzleOtherScore","razzlePDG","showerEnergy",
               "pfp_trackScore","unambiguousSlice"]
    
    razzled_df = uproot.open(hist_path+":razzled/pfpTree").arrays(razzled_col, library="pd")
    opt0_df    = uproot.open(hist_path+":opt0finder/flash_match_tree").arrays(['run','subrun','event','tpc','pfpid','score','hypo_pe','flash_pe'],library="pd")
    
    hdr_df_0 = pd.read_hdf(file,key="hdr")
    nuprim_df_0 = pd.read_hdf(file,key="mcnuprim")
    slctrk_df_0 = pd.read_hdf(file,key="slctrk")
    slcshw_df_0 = pd.read_hdf(file,key="slcshw")

    hdr_df = nue.flatten_df(hdr_df_0)[["ntuple","entry","rec_hdr_run","rec_hdr_subrun","rec_hdr_evt","rec_hdr_pot"]].rename(columns={"rec_hdr_subrun":"subrun","rec_hdr_run":"run","rec_hdr_evt":"event"})
    sub_opt0_df = opt0_df[['run','subrun','event']].drop_duplicates()
    sub_df = hdr_df.merge(sub_opt0_df,on=["run","subrun",'event'],how="inner")
    
    hdr_opt0_df = sub_df[["ntuple","entry","run","subrun","event"]].drop_duplicates().merge(opt0_df,how="left",on=["run","subrun","event"])
    hdr_opt0_df = hdr_opt0_df.rename(columns={"pfpid":"slc_self"})
    
    # make dataframes
    nuprim_df = nuprim_df_0.copy()
    whereFV = nue.maskTrueVtxFv(nuprim_df)
    whereSig = ((nuprim_df.iscc==1) & (abs(nuprim_df.pdg)==12) & (abs(nuprim_df.prim.pdg)==11) & (nuprim_df.prim.startE-nuprim_df.prim.endE > 0.2) )
    nuprim_df = nue.defineBackground(nuprim_df)
    nuprim_df["signal"] = np.where(whereFV & whereSig,0,nuprim_df["signal"])
    
    nuprim_df = nuprim_df_0.copy()

    whereFV = nue.maskTrueVtxFv(nuprim_df)
    whereSig = ((nuprim_df.iscc==1) & (abs(nuprim_df.pdg)==12) & (abs(nuprim_df.prim.pdg)==11) & (nuprim_df.prim.startE-nuprim_df.prim.endE > 0.2) )
    nuprim_df = nue.defineBackground(nuprim_df)
    nuprim_df["signal"] = np.where(whereFV & whereSig,0,nuprim_df["signal"])

    nu_df = nuprim_df.loc[:,:,:,0]
    nu_df = nue.flatten_df(nu_df)
    nu_df["nu_index"] = nu_df["rec.mc.nu__index"]

    # get the dataframe that contains the counts of each PDG per event 
    pdg_counts = nue.getPDGCounts(nuprim_df)

    # merge the pdg counts into the full nu_df 
    nu_df = nu_df.merge(pdg_counts,how="left",on=["ntuple","entry","rec.mc.nu__index"])
    pd.to_pickle(nu_df,              "/sbnd/data/users/lynnt/v09_75_03_03/df/"+sample+"_nu.pkl")
    
    if sample == "rockbox":
        rm_nue = nu_df[nu_df.signal==0][["entry","rec.mc.nu__index"]]
        # Number of signal nueCC in bnb sample:  291
        nu_df = (nu_df.merge(rm_nue,how="left",on=["entry","rec.mc.nu__index"], indicator=True)
            .query('_merge == "left_only"')
            .drop('_merge', 1))
    elif sample == "intrnue":
        nu_df = nu_df.query("signal==0")
    
    # get slcpfp dataframe 
    slcpfp_df = nue.getPFP(slcshw_df_0,slctrk_df_0)

    # merge slcpfp with neutrino events 
    slcpfp_nu_df = slcpfp_df.merge(nu_df,on=nu_idx_set,how="left")
    slcpfp_nu_df["signal"] = np.where(slcpfp_df.slc_tmatch_idx==-999,5,slcpfp_nu_df['signal'])
    if sample == "intrnue":
        slcpfp_nu_df = slcpfp_nu_df.query("signal==0")
    
    # fixes for working with updated bdt trees 
    slcpfp_nu_df = slcpfp_nu_df.drop(columns=["pfp_shw_razzle_electronScore","pfp_trk_dazzle_muonScore"])
    pd.to_pickle(slcpfp_nu_df,       "/sbnd/data/users/lynnt/v09_75_03_03/df/"+sample+"_slcpfp_nu.pkl")

    cosmic_df = slcpfp_nu_df.query("slc_is_clear_cosmic==1")
    noncos_df = slcpfp_nu_df.query("slc_is_clear_cosmic==0")
    
    print("performing join for opt0...")
    # perform join for opt0 
    noncos_self_idx = noncos_df.set_index(["ntuple","entry","slc_self"]).sort_index()
    opt_self_idx = hdr_opt0_df.set_index(["ntuple","entry","slc_self"]).sort_index()
    noncos_df = noncos_self_idx.join(opt_self_idx,how="left").reset_index().sort_values(["ntuple","entry","rec.slc__index"])
    
    print("performing join for razzled...")
    # perform join for razzled 
    razzled_df = razzled_df.query("unambiguousSlice==False")
    noncos_df["pfp_idx"] = noncos_df["pfp_trackScore"].apply(str) + "_" + noncos_df["pfp_shw_truth_p_pdg"].apply(str) +"_" + noncos_df["pfp_shw_bestplane_energy"].apply(np.round,args=(4,)).apply(str)
    razzled_df["pfp_idx"] = razzled_df["pfp_trackScore"].apply(str) + "_" + razzled_df["truePDG"].apply(str) +"_" + (1e-3*razzled_df["showerEnergy"]).apply(np.round,args=(4,)).apply(str)
    noncos_pfp_idx = noncos_df.set_index("pfp_idx").sort_index()
    razzled_pfp_idx = razzled_df.set_index("pfp_idx").sort_index()
    
    noncos_opt_raz_df = noncos_pfp_idx.join(razzled_pfp_idx,how="left",lsuffix='_slc',rsuffix='_raz').reset_index().drop(columns={"pfp_idx"})
    slcpfp_opt_raz_df = pd.concat([noncos_opt_raz_df,cosmic_df],ignore_index=True).sort_values(["ntuple","entry","rec.slc__index","rec.slc.reco.pfp__index"])
    
    pd.to_pickle(slcpfp_opt_raz_df,  "/sbnd/data/users/lynnt/v09_75_03_03/df/"+sample+"_slcpfp_opt_raz.pkl")
    print("finished!!")
    
if __name__ == "__main__":
   main(sys.argv[1:])