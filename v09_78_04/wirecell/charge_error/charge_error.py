import uproot
import numpy as np 
import wc_helper
import mplhep as hep
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift, ifft
from tqdm import tqdm

nroi = 500
niter = 1e3

file = uproot.open("sbnd-data-check.root")
tk = file["hu_decon_2D_init0"].to_numpy()[2][:-1] # ticks 

gaus_sigma = 0.1 
def gaussian_filter(f,a):
    return np.exp(-0.5*(abs(f)/a)**2)

roi_array = np.logspace(0.95,2.95,nroi,dtype=int)
roi_array = np.unique(roi_array)

roi_reso = np.zeros((3, len(roi_array)))
slc_reso = np.zeros((3, len(roi_array)))
for iplane in range(3):
    if iplane==0: plane = "u"
    if iplane==1: plane = "v"
    if iplane==2: plane = "w"
    for i,ROI_length in enumerate(tqdm(roi_array)):
        roi_charge_arr = []
        slc_charge_arr = []
        for iteration in range(int(niter)):
            if iplane==0: ch_num = np.random.randint(0,wc_helper.idx_v0)
            if iplane==1: ch_num = np.random.randint(wc_helper.idx_v0,wc_helper.idx_w0)
            if iplane==2: ch_num = np.random.randint(wc_helper.idx_w0,wc_helper.idx_u1)
            dec_mask = np.where(np.ceil(file["h"+plane+"_decon_2D_init0;1"].to_numpy()[1])==ch_num)[0][0]
            adc = file["h"+plane+"_decon_2D_init0"].to_numpy()[0][dec_mask]
            
            # FFT into frequency domain 
            xf = fftfreq(adc.size,0.5)
            yf = fft(adc)
            wvfm = ifft(gaussian_filter(xf,gaus_sigma)*yf).real
            # subtract the pedestal of the wvfm
            wvfm = wvfm - np.median(wvfm) 
            
            # randomly obtain an ROI window 
            beg_tk = np.random.randint(0,3400-ROI_length)
            end_tk = beg_tk + ROI_length
                
            # get the wvfm snippets using the given ROI 
            wvfm_smp = wvfm[beg_tk:end_tk]
            tk_smp = tk[beg_tk:end_tk]
            # !! define the "baseline" of the sample as the average of the two endpoints, induction plane only !! 
            if plane!="w":
                baseline_smp = (wvfm[beg_tk] + wvfm[end_tk])/2
                wvfm_smp= wvfm_smp - baseline_smp
            roi_charge = np.sum(wvfm_smp[wvfm_smp>0])
            roi_charge_arr.append(roi_charge)
            
            # obtain the "central" time slice, time slices are 4 ticks width
            slc = tk_smp[::4] # start:stop:step, sample the wvfm every 4 ticks 
            slc_central = slc[int(len(slc)/2)] 
            if plane != "w":
                slc_smp = (wvfm - baseline_smp)[int(slc_central):int(slc_central+4)]
            else: 
                slc_smp = wvfm[int(slc_central):int(slc_central+4)]
            slc_charge = np.sum(slc_smp[slc_smp>0])
            slc_charge_arr.append(slc_charge)
        # recast arrays into numpy arrays
        roi_charge_arr = np.array(roi_charge_arr)
        slc_charge_arr = np.array(slc_charge_arr)
        
        # get the MPV of the resolution distribution
        roi_reso[iplane,i] = np.median(roi_charge_arr[roi_charge_arr>0])
        slc_reso[iplane,i] = np.median(slc_charge_arr[slc_charge_arr>0])

np.save("roi_reso.npy",roi_reso)
np.save("slc_reso.npy",slc_reso)

