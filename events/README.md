# Running this code

Download sample data file

```
wget http://mag.gmu.edu/git-data/GaryQ-physics/magnetosphere/3d__var_3_e20031120-070000-000.out.cdf
```

Install Kameleon

```
https://ccmc.gsfc.nasa.gov/Kameleon/Quick_start.html#kameleon-installers
```

Install SciPy compatable with Python 2.7 (Kameleon is not Python 3 compatable)

```
~/kameleon/bin/pip install scipy==0.16
```

Run Kameleon using Python 2.7

```
~/kameleon/bin/python2 cut_plane.py
```
