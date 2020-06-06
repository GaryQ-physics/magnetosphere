import os
import re


def file_list(rootdir, **kwargs):
    """Recursive file list constrained by regular expression
    
    file_list(rootdir)
    file_list(rootdir, regex=regex)

    Example:
    -------
    
    from file_list import file_list
    fl = file_list('.', regex='\.py$') # Find files ending in .py
    print(fl)
    
    """

    file_keep = []
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if 'regex' in kwargs:
                if re.search(kwargs['regex'], file):
                    file_keep.append(file)
            else:
                file_keep.append(file)
    return file_keep
            