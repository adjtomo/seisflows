from setuptools import setup, find_packages

setup(
    name='seisflows',
    version='0.1rc1',  # FIXME: dummy version number
    description='a seismic inversion package',
    long_description='',
    url='https://github.com/PrincetonUniversity/seisflows',
    author='Tromp group, Princeton University',
    license='bsd-2',

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

    install_requires=[
        'numpy'
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
