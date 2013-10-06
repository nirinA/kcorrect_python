from distutils.core import setup, Extension

INC_DIR = ['/mnt/big/store/kcorrect/include']
LIB_DIR = ['/mnt/big/store/kcorrect/lib']
modulekcorrect = Extension('kcorrect',
                           include_dirs = INC_DIR,
                           libraries = ['kcorrect', 'm'],
                           library_dirs = LIB_DIR,
                           sources = ['kcorrectmodule.c'])

setup (name = 'kcorrect20131006',
      version = '0',
      description = 'This is kcorrect for Python',
      ext_modules = [modulekcorrect,
                     ])
