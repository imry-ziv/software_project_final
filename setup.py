from setuptools import Extension, setup
module = Extension("symnmf_c_api", sources=['ccode/symnmfmodule.c'])
setup(name='symnmf_c_api',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module]
)