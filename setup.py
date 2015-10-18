from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(
    name='PrimeTest',
    ext_modules=[Extension('C_Shotgun', ['Shotgun/C_Shotgun.pyx']),
                 Extension('C_DeBrujin', ['SBH/C_DeBrujin.pyx'])],
    cmdclass={'build_ext': build_ext}
)