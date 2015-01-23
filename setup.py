from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as _build_ext


# Numpy trick, see:
# http://stackoverflow.com/questions/19919905/how-to-bootstrap-numpy-installation-in-setup-py
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy


setup(
    name='seisflows',
    version='0.1rc1',  # FIXME: dummy version number
    description='a seismic inversion package',
    long_description='',
    url='https://github.com/PrincetonUniversity/seisflows',
    # author='Tromp group, Princeton University',
    # author_email='rmodrak@princeton.edu',
    license='bsd',

    # Classifiers, see list at:
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Physics'
    ],

    keywords='seismic modeling migration inversion fwi',

    packages=find_packages(exclude=['docs', 'tests']),

    # Trick numpy distutils
    cmdclass={'build_ext':build_ext},
    setup_requires=['numpy'],

    install_requires=[
        'numpy',
        'scipy',
        'matplotlib'
    ],

    scripts=[
        'scripts/plotgrid',
        'scripts/plotmesh',
        'scripts/plotsegy',
        'scripts/plotseis',
        'scripts/sfclean',
        'scripts/sfexamples',
        'scripts/sfexamples-research',
        'scripts/sfrun'
    ]

)
