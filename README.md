# Running this code

## Convert BATSRUS .out files to Kameleon .cdf files using Kameleon-converter

Install (kameleon-converter-5.2)[https://github.com/ccmc/ccmc-software/tree/master/kameleon_converter/tags/kameleon-converter-v5.2.0]

```
curl -O http://mag.gmu.edu/git-data/GaryQ-Physics/magnetosphere/programs/kameleon-converter-v5.2.0.tgz
tar zxvf kameleon-converter-v5.2.0.tgz
cd kameleon-converter-v5.2.0;
```
See README-KAMELEON-CONVERSION-INSTALL-HELP.txt for compiling. After compiling,

```
curl -O http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070000-000.out
./kameleon -v -f cdf -m batsrus -o ~/ ./3d__var_3_e20031120-070000-000.out
```

## Install Kameleon+

Install Kameleon+

```
https://ccmc.gsfc.nasa.gov/Kameleon/Quick_start.html#kameleon-installers
```

Install SciPy compatable with Python 2.7 (Kameleon is not Python 3 compatable)

```
~/kameleon/bin/pip install scipy==0.16
```

## Run code

Download sample data file

```
cd events; wget http://mag.gmu.edu/git-data/GaryQ-Physics/magnetosphere/3d__var_3_e20031120-070000-000.out.cdf
# or
cd events; curl -O http://mag.gmu.edu/git-data/GaryQ-Physics/magnetosphere/3d__var_3_e20031120-070000-000.out.cdf
```

Run Kameleon using Python 2.7

```
cd events; ~/kameleon/bin/python2 cut_plane.py
```

With Anaconda,

```
conda create -n python2.7 python=2.7
conda activate python2.7
conda install scipy=0.16
spyder
```
