conda --version
#conda 4.9.2
conda create --name tmp37 python=3.7
conda activate tmp37
pip install spacepy
python -c "import sys; print(sys.version); import spacepy; print(spacepy.__version__); import spacepy.pybats.bats"
#3.7.10 (default, Feb 26 2021, 18:47:35) 
#[GCC 7.3.0]
#0.2.2
#/home/gary/miniconda3/envs/tmp37/lib/python3.7/site-packages/spacepy/time.py:2294: UserWarning: Leapseconds may be out of date. Use spacepy.toolbox.update(leapsecs=True)
#  warnings.warn('Leapseconds may be out of date.'
#/home/gary/miniconda3/envs/tmp37/lib/python3.7/site-packages/spacepy/omni.py:415: UserWarning: Qin-Denton/OMNI2 data not found in current format. This module has limited functionality. Run spacepy.toolbox.update(QDomni=True) to download data.
#  "Qin-Denton/OMNI2 data not found in current format."
#Traceback (most recent call last):
#  File "<string>", line 1, in <module>
#  File "/home/gary/miniconda3/envs/tmp37/lib/python3.7/site-packages/spacepy/pybats/__init__.py", line 119, in <module>
#    import spacepy.plot.apionly
#  File "/home/gary/miniconda3/envs/tmp37/lib/python3.7/site-packages/spacepy/plot/__init__.py", line 98, in <module>
#    plt.register_cmap(name='plasma', cmap=_plasma)
#  File "/home/gary/miniconda3/envs/tmp37/lib/python3.7/site-packages/matplotlib/cm.py", line 149, in register_cmap
#    raise ValueError(msg)
#ValueError: Trying to re-register the builtin cmap 'plasma'.
