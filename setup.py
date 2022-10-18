from setuptools import setup, find_packages

setup(name="seisflows",
      version="2.1.0",
      description="SeisFlows: A seismic inversion package",
      url="https://github.com/adjtomo/seisflows",
      author="Seisflows Development Team",
      packages=find_packages(),
      entry_points={
          "console_scripts": ["seisflows=seisflows.seisflows:main",]},
      license="GPL",
      install_requires=[
          "obspy>=1.3.0",
          "pyyaml>=5.3.1",
          "IPython>=7.31.1",
          "dill>=0.3.5.1",
          ],
      zip_save=False
      )
