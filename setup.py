from setuptools import setup, find_packages

setup(name="seisflows",
      version="1.0.0",
      description="SeisFlows3: A seismic inversion package",
      url="https://github.com/seisflows/seisflows",
      author="Seisflows Development Team",
      packages=find_packages(),
      entry_points={
          "console_scripts": ["seisflows=seisflows.seisflows:main",]},
      license="GPL",
      install_requires=[
          "obspy>=1.2.2",
          "pyyaml>=5.3.1",
          "IPython>=7.31.1"
          ],
      zip_save=False
      )
