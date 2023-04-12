import stagpy
from stagpy import stagyydata
sdat = stagyydata.StagyyData('/lfs/jiching/w0201/') # model path
sdat.__dir__() # 
#for ww in 
print(sdat.scales) # acceleration, density, dyn_visc,
# heat_flux, heat_production,length, power,
# sp_heat, stress, temperature, th_cond, th_diff,
# time, velocity 
sdat.refstate #  adiabats [two matrix], systems [one matrix]
sdat.tseries # at_step(setp)[29 個參數], isteps[], 

sdat.par['refstate'] # find input parameters
