# Directory Contents

This directory contains analysis scripts for nue selection using a rockbox sample.

- `nue_old_selection.ipynb` replicates the selection cuts used for the nue selection on MCP2022A on the existing simulation (as of `sbndcode v09_75_03_02`).
- `nue_selection.ipynb` performs a updated selection given the status of the simulation (`sbndcode v09_75_03_02`). Changes from the last production include:  Pandora track and shower variables being saved for every PFP, updated optical simulation, and more.
- `nue_selection_QL.ipynb` studies the impact of using charge-light calorimetry on the nue selection.

All of these selections use a combination of a BNB rockbox sample and intrnue rockbox sample. In order to properly compare the signal/backgrounds, a scaling factor must be applied to events from the intrnue rockbox sample to account for the difference in POT.
