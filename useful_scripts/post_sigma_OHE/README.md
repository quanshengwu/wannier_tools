Usage: python3   
Copy this file to your work folder where contains wt.in and sigma files
Command: `python post_sigma_OHE.py`  

This python script contains three function, readwtin(), Function readsigma_writerho(band_list, tau_list) and plot_cbar(component). 

1.1 Function readwtin()

    Description: Reads namelist PARAMETERS from wt.in
    Returns: Tmin, Tmax, NumT, OmegaMin, OmegaMax, OmegaNum, BTauMax, BTauNum

1.2 Function readsigma_writerho(band_list, tau_list)

    Description: Reads sigma/tau tensor data from sigma_bands_mu_***eV.dat and calculates the total resistivity.
    The resistivity is stored in rhotau_total_mu_***eV.dat and rho_all.npy.
    Parameters:
    - band_list: A list specifying the bands to be included, in accordance with the SELECTEDBANDS in wt.in
    - tau_list: A list specifying the ratio of relaxation time for the corresponding bands.
    Note: The values in tau_list cannot be zero. To exclude the contribution of certain bands,
    set the corresponding value in tau_list to a very small value (e.g., 0.00001).

1.3 Function plot_cbar(component)

    Description: Plots rhotau(Btau).
    Parameters:
    - component: A list specifying the tensor component to be plotted, such as ['xx']


To make it work, use python3 (from anoconda3) and matplotlib (pip install --upgrade matplotlib).


Author: Hanqi Pi, hqpi1999@gmail.com
2023/06/25


