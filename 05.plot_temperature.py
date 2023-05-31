import stagpy
from stagpy import stagyydata
from stagpy import field
path = '/lfs/jiching/'
model = 'w0203'
data = stagyydata.StagyyData(path+ model)
shot = 120
field.get_meshes_fld(data.snaps[shot],'T')
print(data.snaps[shot].fields['T'])
