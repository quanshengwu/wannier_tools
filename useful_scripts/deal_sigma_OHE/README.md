Usage: python3   
Command: `python deal_sigma_OHE.py`  

This python script contains two function, readsigma_writerho() and plot_cbar(component). 
1.function readsigma_writerho() 
  This function read sigma/tau tensor data from every single band, adds them together and calculate the total resistivity. 
  You can modify the script if you need electrons in different bands to have different relaxation time.
  This function reads paramenters from `wt.in` and sigma tensor from `sigma_bands_mu_***eV.dat`.
  The resistivity is stored in `rhotau_total_mu_***eV.dat` and `rho_all.npy`
2.function plot_cbar(component)
  This function is designed to plot rhotau(Btau).
  You should provide a list called `component` that includes the tensor component that you want to plot, such as ['xx']


To make it work, use python3 (from anoconda3) and matplotlib (pip install --upgrade matplotlib).


Author: Hanqi Pi, hqpi1999@gmail.com
2023/05/04


