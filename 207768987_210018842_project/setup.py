from setuptools import Extension, setup
module = Extension("symnmf_c_api", sources=['symnmfmodule.c', 'symnmf.c'],
                   extra_compile_args=['-g']
)
setup(name='symnmf_c_api',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module],
)