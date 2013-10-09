from distutils.core import setup, Extension

INC_DIR = ['/mnt/big/store/kcorrect/include',
           '/mnt/sda6/Py330/lib/python3.3/site-packages/numpy/core/include']
LIB_DIR = ['/mnt/big/store/kcorrect/lib']
modulekcorrect = Extension('_kcorrect',
                           include_dirs = INC_DIR,
                           libraries = ['kcorrect', 'm'],
                           library_dirs = LIB_DIR,
                           sources = ['src/_kcorrectmodule.c'])

setup (name = 'kcorrect_python',
      version = '20131008',
      description = 'This is kcorrect for Python',
      ext_modules = [modulekcorrect,
                     ])
