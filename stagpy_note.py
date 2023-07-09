import stagpy
from stagpy import field
from stagpy import stagyydata
sdat = stagyydata.StagyyData('/lfs/jiching/thermo_chemical/TC_2D-SPH_Basalt_HR002/') # model path
sdat.__dir__() # 

print(sdat.scales) # acceleration, density, dyn_visc,
# heat_flux, heat_production,length, power,
# sp_heat, stress, temperature, th_cond, th_diff,
# time, velocity 
sdat.refstate #  adiabats [two matrix], systems [one matrix]
sdat.tseries # at_step(setp)[29 個參數], isteps[], 

sdat.par['refstate'] # find input parameters
print(field.get_meshes_fld(sdat.snaps[15],'T'))
