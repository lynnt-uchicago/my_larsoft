# %%
import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import nue as nue
from tqdm import tqdm
import importlib
____ = importlib.reload(nue)

import dask.dataframe as dd
from tqdm.dask import TqdmCallback
from dask.diagnostics import ProgressBar
ProgressBar().register()

# %%
hist_path  = "/sbnd/data/users/lynnt/v09_75_03_03/rockbox_hists.root" 
razz_path  = hist_path+":razzled/pfpTree"
opt0_path  = hist_path+":opt0finder/flash_match_tree"
# **loading the arrays takes about 2-3 min on build server** 
print("loading arrays...")
razzled_col = ["electronScore","muonScore","photonScore","pionScore","protonScore","otherScore","bestScore","bestPDG","truePDG",
               "dazzleMuonScore","dazzlePionScore","dazzleProtonScore","dazzleOtherScore","dazzlePDG","trk_length",
               "razzleElectronScore","razzlePhotonScore","razzleOtherScore","razzlePDG","showerEnergy",
               "pfp_trackScore","unambiguousSlice"]
razzled_df = uproot.open(razz_path).arrays(razzled_col, library="pd")
opt0_df    = uproot.open(opt0_path).arrays(['run','subrun','event','tpc','pfpid','score','hypo_pe','flash_pe'],library="pd")

file = "/sbnd/data/users/lynnt/v09_75_03_03/rockbox.df"
hdr_df_0 = pd.read_hdf(file,key="hdr")
nuu_df_0 = pd.read_hdf(file,key="mcnu")
nuprim_df_0 = pd.read_hdf(file,key="mcnuprim")
slctrk_df_0 = pd.read_hdf(file,key="slctrk")
slcshw_df_0 = pd.read_hdf(file,key="slcshw")

nu_idx_set = ["ntuple","entry","nu_index"]
slc_idx_set = ["ntuple","entry","rec.slc__index"]

hdr_df = nue.flatten_df(hdr_df_0)[["ntuple","entry","rec_hdr_run","rec_hdr_subrun","rec_hdr_evt","rec_hdr_pot"]].rename(columns={"rec_hdr_subrun":"subrun","rec_hdr_run":"run","rec_hdr_evt":"event"})
sub_opt0_df = opt0_df[['run','subrun','event']].drop_duplicates()
sub_df = hdr_df.merge(sub_opt0_df,on=["run","subrun",'event'],how="inner")

# %%
# make bnb dataframes
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
pd.to_pickle(nu_df,              "/sbnd/data/users/lynnt/v09_75_03_03/df/rockbox_nu.pkl")

# remove nue events from the bnb sample
# ! because we're using a rockbox + intrnue (non rockbox), 
# ! let's try just scaling the number of signal nueCC from 
# ! intrnue by the number of signal nueCC from bnb rockbox.
# ! In this case, we want to know how many signal nueCC are 
# ! in the bnb rockbox sample before removing them. Do not remove
# ! dirt or "other" nue events.
rm_nue = nu_df[nu_df.signal==0][["entry","rec.mc.nu__index"]]
nsignal = len(rm_nue)
print("Number of signal nueCC in bnb sample: ",nsignal)
# Number of signal nueCC in bnb sample:  291

nu_df = (nu_df.merge(rm_nue,how="left",on=["entry","rec.mc.nu__index"], indicator=True)
     .query('_merge == "left_only"')
     .drop('_merge', 1))

# get slcpfp dataframe 
slcpfp_df = nue.getPFP(slcshw_df_0,slctrk_df_0)

# merge slcpfp with neutrino events 
slcpfp_nu_df = slcpfp_df.merge(nu_df,on=nu_idx_set,how="left")
slcpfp_nu_df["signal"] = np.where(slcpfp_df.slc_tmatch_idx==-999,5,slcpfp_nu_df['signal'])

# fixes for working with updated bdt trees 
slcpfp_nu_df = slcpfp_nu_df.drop(columns=["pfp_shw_razzle_electronScore","pfp_trk_dazzle_muonScore"])

# %%
# combine the razzled and opt0 dataframes 
# do this using dask dataframes 
# print("performing preselection...")
# pre_slcpfp_nu_df = nue.cutPreselection(slcpfp_nu_df)
# pre_dd = dd.from_pandas(pre_slcpfp_nu_df,npartitions=30)

slc_dd = dd.from_pandas(slcpfp_nu_df,npartitions=30)

raz_df = razzled_df.query("unambiguousSlice == False").rename(columns={"trk_length":"pfp_trk_len"})
raz_dd = dd.from_pandas(raz_df,npartitions=30)

opt0_df["frac_pe"] = (opt0_df.hypo_pe - opt0_df.flash_pe)/opt0_df.flash_pe
opt0_df = opt0_df[opt0_df.pfpid != -1]
opt0_df = opt0_df[(opt0_df.score == opt0_df.groupby(["run","subrun","event","pfpid"])["score"].transform(max))]
hdr_opt0_df = sub_df[["ntuple","entry","run","subrun","event"]].drop_duplicates().merge(opt0_df,how="left",on=["run","subrun","event"])
hdr_opt0_df = hdr_opt0_df.rename(columns={"pfpid":"slc_self"})
hdr_opt0_df = hdr_opt0_df.astype({'slc_self':'int'})
opt_dd = dd.from_pandas(hdr_opt0_df,npartitions=30)

print("merging slc df and razzled df...")
slc_dd = slc_dd.merge(raz_dd,on=["pfp_trackScore","pfp_trk_len"],how="left")
print("merging slc df and opt0 df...")
slc_dd = slc_dd.merge(opt_dd,on=["ntuple","entry","slc_self"],how="left")

print("converting dask dataframe to pandas dataframe...")
slc_df = slc_dd.compute()
# save everything to data 
pd.to_pickle(nu_df,           "/sbnd/data/users/lynnt/v09_75_03_03/df/rockbox_nu.pkl")
pd.to_pickle(slcpfp_nu_df,    "/sbnd/data/users/lynnt/v09_75_03_03/df/rockbox_slcpfp.pkl")
# pd.to_pickle(pre_slcpfp_nu_df,"/sbnd/data/users/lynnt/v09_75_03_03/df/intrnue_pre_slcpfp.pkl")
pd.to_pickle(slc_df,          "/sbnd/data/users/lynnt/v09_75_03_03/df/rockbox_slc_raz_opt.pkl")
