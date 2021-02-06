#todo put in seperate repository
import os
import sys

class Fortran_Subroutine:
    """python class allowing for the definition of jit compiled fortran subroutines using f2py

    passed a string that is the regex for a fortran subroutine, along with the 
    subroutines name string. The name passed much match the name in the script 
    in order to function properly

    instance = Fortran_Subroutine(script, subroutine_name)

    the subroutine can be called then via instance.execute(...)
    instance.execute(*args)

    the compiling of the script will happen automatically the first time instance.execute is
    called, but not for sucessive calls so long as the class instance remains in memory. 
    A compilation (or recompilation) can be forced by calling:
    instance.compile()
    """
    def __init__(self, script, subroutineName):#todo, determine name from top of script
        self.script = script
        self.subroutineName = subroutineName
        self.activeSubroutine = None

    def compile(self):
        with open('/tmp/jitFORTRAN_script.f90', 'w') as f:
            f.write(self.script)

        os.system('f2py -c /tmp/jitFORTRAN_script.f90 -m jitFORTRAN_exe -DF2PY_REPORT_ON_ARRAY_COPY=1')
        os.system('mv jitFORTRAN_exe.so /tmp/jitFORTRAN_exe.so')

        sys.path.append('/tmp/')
        import jitFORTRAN_exe
        exec('sub = jitFORTRAN_exe.'+self.subroutineName.lower())
        sys.path.remove('/tmp/')

        self.activeSubroutine = sub

    def execute(self, *args):
        if self.activeSubroutine is None:
            self.compile()

        return self.activeSubroutine(*args)

