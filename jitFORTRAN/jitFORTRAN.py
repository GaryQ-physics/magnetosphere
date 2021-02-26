#todo put in seperate repository
import os
import sys
import shutil

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

    the subroutine must have any array dimensions passed as well, but if you place them
    as the last arguments, then you can ignore them when calling calling instance.execute(*args),
    since the last arguments will be auto filled in by the f2py created module
    """
    def __init__(self, script, subroutineName, include=None):#todo, determine name from top of script
        self.script = script
        self.subroutineName = subroutineName
        self.activeSubroutine = None
        self.include = include
    def compile(self):
        backup_sys_argv = sys.argv[:]
        tempdir = os.tmpnam()

        scriptfile = os.path.join(tempdir,'jitFORTRAN_script.f90')
        os.mkdir(tempdir)
        with open(scriptfile, 'w') as f:
            f.write(self.script)

        ###### set f2py_args , a list such that ' '.join(f2py_args) would be executed in terminal on linux to compile to *.so file
        if isinstance(self.include, str):
            if os.path.exists(self.include+'.o'):
                pass ## TODO: check if existing .o was properly compiled with -fPIC
            elif os.path.exists(self.include+'.f'):
                os.system('gfortran -c -fPIC %s.f -o %s.o'%(self.include,self.include))
            elif os.path.exists(self.include+'.f90'):
                os.system('gfortran -c -fPIC %s.f90 -o %s.o'%(self.include,self.include))
            else:
                raise FileNotFoundError ('no file for %s{.f,.f90,.o}'%(self.include))

            f2py_args = ['f2py','-c',scriptfile, '-I', '%s.o'%(self.include), '-m', 'jitFORTRAN_exe', '-DF2PY_REPORT_ON_ARRAY_COPY=1']
        else:
            f2py_args = ['f2py','-c',scriptfile, '-m', 'jitFORTRAN_exe', '-DF2PY_REPORT_ON_ARRAY_COPY=1']

        ###### try to use numpy.f2py directly to get *.so file, if fails use os.system() to mimick terminal
        try:
            from numpy.f2py.f2py2e import main
            sys.argv = f2py_args[:]
            main()
            sys.argv = backup_sys_argv[:]
        except:
            raise RuntimeWarning ('EXCEPTED')
            print(' '.join(f2py_args))
            sys.argv = backup_sys_argv[:]
            os.system(' '.join(f2py_args))

        ### move *.so to temporary directory, import the library now in the temporary directory, and store callable subroutine
        shutil.move('jitFORTRAN_exe.so', os.path.join(tempdir, 'jitFORTRAN_exe.so'))
        sys.path.append(tempdir)
        import jitFORTRAN_exe
        reload(jitFORTRAN_exe)
        exec('sub = jitFORTRAN_exe.'+self.subroutineName.lower())
        sys.path.remove(tempdir)
        os.remove(os.path.join(tempdir, 'jitFORTRAN_exe.so'))
        del jitFORTRAN_exe
        self.activeSubroutine = sub
        del sub

    def execute(self, *args):
        if self.activeSubroutine is None:
            self.compile()

        return self.activeSubroutine(*args)

