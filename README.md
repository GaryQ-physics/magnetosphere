# Running this code

Download sample data file (see next section for information on how this file was created)

```
wget http://mag.gmu.edu/git-data/GaryQ-Physics/magnetosphere/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070000-000.out.cdf
# or
curl -O http://mag.gmu.edu/git-data/GaryQ-Physics/magnetosphere/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070000-000.out.cdf
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

# Converting .out to .cdf

Download kameleon-converter-5.2.0.tgz located at

```
http://mag.gmu.edu/git-data/GaryQ-Physics/magnetosphere/kameleon-converter/
```

into `/tmp`. Follow the compilation instructions in the package. This process is complicated and requires first compiling `libcdf`. Libargtable2 packages may be needed and packages for Ubuntu are at the above URL.

Download the .out files

```
wget -m -nH -R index*,*.gif --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/
```

Convert the files

```
cd SCARR5_GM_IO2/IO2
ls -1 3d_*.out | xargs -i /tmp/kameleon-converter/kameleon-converter-v5.2.0/kameleon -v -f cdf -m batsrus -o . {}
```


