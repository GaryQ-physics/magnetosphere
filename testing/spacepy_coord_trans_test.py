import spacepy.coordinates as sc
from spacepy.time import Ticktock
import spacepy
print('\nimporting from '+ spacepy.__path__[0])
print('\n')
print('spacepy version ' + spacepy.__version__)

#time = (2019, 1,1,1,1,1)
time = (2019,9,2,6,30,0)
lat = 18.907
lon = 72.815

v_in = [1., lat, lon]

cvals = sc.Coords(v_in, 'GEO', 'sph')
t_str = '%04d-%02d-%02dT%02d:%02d:%02d' %time
cvals.ticks = Ticktock(t_str, 'ISO')
newcoord = cvals.convert('MAG', 'sph')
data = newcoord.data

v_out = data[0, :].tolist()

print(v_out)

sunspot = [0.9999999999999999, 11.016639272301447, 147.323425223922]
local = [1.0, 11.058617393065475, 146.89653169520292]

if sunspot == v_out:
    print('probably correct (was getting on sunspot)')

if local == v_out:
    print('probably incorrect (was getting on local)')
