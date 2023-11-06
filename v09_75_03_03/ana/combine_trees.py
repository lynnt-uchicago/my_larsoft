# %%
import uproot
import numpy as np
import pandas as pd
import nue as nue
from tqdm import tqdm
import importlib
____ = importlib.reload(nue)
import sys

tqdm.pandas()

nu_idx_set = ["ntuple","entry","nu_index"]
slc_idx_set = ["ntuple","entry","rec.slc__index"]
evt_idx_set = ["run","subrun","event"]

def main(argv):
    hist_path = argv[0]
    file = argv[1]
    output_dir = argv[2]
    
    sample = "sample"
    # nu_sample = True
    if "rockbox" in file: 
        sample = "rockbox"
    elif "intrnue" in file:
        sample = "intrnue"
    elif "intime" in file:
        sample = "intime"
    
    print("Sample label is: ", sample)
    
    # if (sample=="intime"):
    #     nu_sample = False
        
    print("loading arrays...")
    
    razz_path = hist_path+":ncpizeroana/events"
    opt0_path = hist_path+":opt0finder/flash_match_tree"

    opt0_col = ['run','subrun','event','tpc','pfpid','score','hypo_pe','flash_pe']

    razz_pfp_col = ['run','subrun','event','slc_pfp_id',
                    'slc_pfp_razzled_electron_score','slc_pfp_razzled_muon_score','slc_pfp_razzled_pdg','slc_pfp_razzled_photon_score','slc_pfp_razzled_pion_score','slc_pfp_razzled_proton_score']

    razz_slc_col = ['run','subrun','event','slc_primary_pfp_id','slc_n_pfps',
                    'slc_n_razzled_electrons','slc_n_razzled_muons','slc_n_razzled_photons','slc_n_razzled_pions','slc_n_razzled_protons',]

    print("opening opt0 file...")
    opt0_df  = uproot.open(opt0_path).arrays(opt0_col,library="pd")
    print("opening razz pfp file...")
    razz_pfp_0 = uproot.open(razz_path,num_workers=64).arrays(razz_pfp_col,library="pd")
    print("opening razz slc file...")
    razz_slc_0 = uproot.open(razz_path,num_workers=64).arrays(razz_slc_col,library="pd")

    # first explode :)
    razz_pfp = razz_pfp_0.copy()
    tqdm.pandas(desc="1st explode")
    razz_pfp = razz_pfp.set_index(['run',"subrun","event"]).progress_apply(pd.Series.explode).reset_index()

    # change the stlvector to np.array
    tqdm.pandas(desc="Apply np array")
    razz_pfp = razz_pfp.set_index(['run',"subrun","event"]).progress_applymap(lambda x: np.array(x))

    # create n_razzled column for sanity check as well as the slc idx for joining
    tqdm.pandas(desc="Create columns")
    razz_pfp["n_razzled"] = razz_pfp["slc_pfp_razzled_pdg"].progress_apply(lambda x: len(x[x>0]))
    razz_pfp = razz_pfp.reset_index()
    razz_pfp["slc_idx"] = razz_pfp.groupby(["run","subrun","event"]).transform("cumcount").add(1) - 1

    # # second explode :)) 
    tqdm.pandas(desc="2nd explode")
    razz_pfp = razz_pfp.set_index(['run',"subrun","event","slc_idx","n_razzled"]).progress_apply(pd.Series.explode).reset_index()

    razz_slc = razz_slc_0.copy()
    # add column for the slc idx for joining
    razz_slc["slc_idx"] = razz_slc.groupby(["run","subrun","event"]).transform("cumcount").add(1) - 1
    # this is for sanity check later 
    razz_slc["n_razzled"] = (razz_slc["slc_n_razzled_electrons"] 
                            + razz_slc["slc_n_razzled_muons"]
                            + razz_slc["slc_n_razzled_photons"]
                            + razz_slc["slc_n_razzled_pions"]
                            + razz_slc["slc_n_razzled_protons"])
    razz_slc = razz_slc.drop(columns=[ 'slc_n_razzled_electrons', 'slc_n_razzled_muons', 'slc_n_razzled_photons', 'slc_n_razzled_pions', 'slc_n_razzled_protons',])

    razz_slc_idx = razz_slc.set_index(["run","subrun","event","slc_idx"]).sort_index()
    razz_pfp_idx = razz_pfp.set_index(['run',"subrun","event","slc_idx"]) .sort_index()
    razz_df      = razz_slc_idx.join(razz_pfp_idx,how="left",lsuffix="_slc",rsuffix="_pfp").reset_index()

    print("reading hdf5 files...")
    hdr_df_0 = pd.read_hdf(file,key="hdr")
    # nuprim_df_0 = pd.read_hdf(file,key="mcnuprim")
    slctrk_df_0 = pd.read_hdf(file,key="slctrk")
    slcshw_df_0 = pd.read_hdf(file,key="slcshw")

    # get run/subrun/event info correlated with entry/ntuple
    hdr_df = nue.flatten_df(hdr_df_0)[["ntuple","entry","rec_hdr_run","rec_hdr_subrun","rec_hdr_evt","rec_hdr_pot"]].rename(columns={"rec_hdr_subrun":"subrun","rec_hdr_run":"run","rec_hdr_evt":"event"})
    sub_opt0_df = opt0_df[evt_idx_set].drop_duplicates()
    sub_opt0_df = hdr_df.merge(sub_opt0_df,on=evt_idx_set,how="inner")
    
    # add frac pe column and keep only one opt0 match per slice 
    opt0_df["frac_pe"] = (opt0_df.hypo_pe - opt0_df.flash_pe)/opt0_df.flash_pe
    opt0_df = opt0_df[opt0_df.pfpid != -1]
    opt0_df = opt0_df[(opt0_df.score == opt0_df.groupby(["run","subrun","event","pfpid"])["score"].transform(max))]
    
    # get opt0 information correlated with entry/ntuple
    opt0_hdr_df = sub_opt0_df[["ntuple","entry","run","subrun","event"]].drop_duplicates().merge(opt0_df,how="left",on= evt_idx_set)
    opt0_hdr_df = opt0_hdr_df.drop(columns=evt_idx_set)
    opt0_hdr_df = opt0_hdr_df.rename(columns={"pfpid":"slc_self"})
        
    # make dataframes
    # if (nu_sample):
    #     nuprim_df = nuprim_df_0.copy()

    #     whereFV = nue.maskTrueVtxFv(nuprim_df)
    #     whereSig = ((nuprim_df.iscc==1) & (abs(nuprim_df.pdg)==12) & (abs(nuprim_df.prim.pdg)==11) & (nuprim_df.prim.startE-nuprim_df.prim.endE > 0.2) )
    #     nuprim_df = nue.defineBackground(nuprim_df)
    #     nuprim_df["signal"] = np.where(whereFV & whereSig,0,nuprim_df["signal"])

    #     nu_df = nuprim_df.loc[:,:,:,0]
    #     nu_df = nue.flatten_df(nu_df)
    #     nu_df["nu_index"] = nu_df["rec.mc.nu__index"]

    #     # get the dataframe that contains the counts of each PDG per event 
    #     pdg_counts = nue.getPDGCounts(nuprim_df)
    #     # nu_df = nu_df[nu_df.signal==0]
    #     # merge the pdg counts into the full nu_df 
    #     nu_df = nu_df.merge(pdg_counts,how="left",on=["ntuple","entry","rec.mc.nu__index"])
    #     pd.to_pickle(nu_df,              "/sbnd/data/users/lynnt/v09_75_03_03/df/ana_"+sample+"_nu.pkl")
    
    # get slcpfp dataframe 
    slcpfp_df = nue.getPFP(slcshw_df_0,slctrk_df_0)

    # merge slcpfp with neutrino events 
    # if (nu_sample):
    #     slcpfp_nu_df = slcpfp_df.merge(nu_df,on=nu_idx_set,how="left")
    #     slcpfp_nu_df["signal"] = np.where(slcpfp_df.slc_tmatch_idx==-999,5,slcpfp_nu_df['signal'])
    # else: 
        # in-time cosmics
    # slcpfp_nu_df = slcpfp_df
    # slcpfp_nu_df["signal"] = 6

    # fixes for working with updated bdt trees 
    slcpfp_df = slcpfp_df.drop(columns=["pfp_shw_razzle_electronScore","pfp_trk_dazzle_muonScore"])
    pd.to_pickle(slcpfp_df,output_dir+"/ana_"+sample+"_slcpfp.pkl")

    cosmic_df = slcpfp_df.query("slc_is_clear_cosmic==1")
    noncos_df = slcpfp_df.query("slc_is_clear_cosmic==0")

    # move run/subrun/event number into the df 
    noncos_idx = noncos_df.set_index(["ntuple","entry"]).sort_index()
    hdr_idx   = hdr_df.set_index(["ntuple","entry"]).sort_index()
    noncos_df = noncos_idx.join(hdr_idx).reset_index()

    # move opt0 information into the df
    noncos_self_idx = noncos_df.set_index(["ntuple","entry","slc_self"]).sort_index()
    opt_self_idx = opt0_hdr_df.set_index(["ntuple","entry","slc_self"]).sort_index()
    noncos_df = noncos_self_idx.join(opt_self_idx,how="left").reset_index()

    # move razz information into the df 
    noncos_evt_idx = noncos_df.set_index(["run","subrun","event","slc_self","pfp_id"]).sort_index() 
    razz_evt_idx   = (razz_df.rename(columns={"slc_primary_pfp_id":"slc_self","slc_pfp_id":"pfp_id"})
                            .set_index(["run","subrun","event","slc_self","pfp_id"])).sort_index()
    noncos_df = noncos_evt_idx.join(razz_evt_idx,how="left").reset_index()

    slcpfp_df = pd.concat([cosmic_df,noncos_df]).sort_values(["ntuple","entry","rec.slc__index"])
    pd.to_pickle(slcpfp_df,output_dir+"/ana_"+sample+"_slcpfp_opt_raz.pkl")
    print("finished!!")

if __name__ == "__main__":
   main(sys.argv[1:])