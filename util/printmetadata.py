import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

def printmetadata(f):
    """Print metadata in a Kameleon .cdf or SWMF .out file
    
    Example
    -------
    import os, sys
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
    from config import conf
    from printmetadata import printmetadata
    
    printmetadata(conf['run_path'] + '3d__var_3_e20031120-070000-000.out.cdf')
    printmetadata(conf['run_path'] + '3d__var_3_e20031120-070000-000.out')
    """

    if not os.path.exists(f):
        raise Exception("File not found: " + f)


    if f[-3:] == "cdf":
        import _CCMC as ccmc

        print('----------------------------------------------------------------------')
        
        kameleon = ccmc.Kameleon()
        print("Opening " + f)
        kameleon.open(f)
        
        for i in range(kameleon.getNumberOfGlobalAttributes()):
            gname = kameleon.getGlobalAttributeName(i)
            gattr = kameleon.getGlobalAttribute(gname)
            if gname != 'README':
                print gname, gattr.toString()
                
        for i in range(kameleon.getNumberOfVariables()):
            varname  = kameleon.getVariableName(i)
            min_attr = kameleon.getVariableAttribute(varname, 'actual_min').getAttributeFloat()
            max_attr = kameleon.getVariableAttribute(varname, 'actual_max').getAttributeFloat()
            units = kameleon.getVisUnit(varname)
            units2 = kameleon.getNativeUnit(varname)
            print varname, '\t', min_attr,'\t', max_attr, units, units2
        
        kameleon.close()    
        print('----------------------------------------------------------------------')

    elif f[-3:] == "out":
        import spacepy.pybats.bats as bats

        print('----------------------------------------------------------------------')
        
        print("Opening " + f)
        
        data3d = bats.Bats2d(f)
        
        # look at keys:
        print(data3d.keys())
        
        print(data3d.attrs)
        print(data3d.meta)
        
        print('----------------------------------------------------------------------')
    else:
        raise ValueError("File name does not end in cdf or out. " + f)
