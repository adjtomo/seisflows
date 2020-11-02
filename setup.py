from setuptools import setup, find_packages

setup(name="seisflows3",
      version="1.0.0",
      description="SeisFlows3: A seismic inversion package",
      url="https://github.com/seisflows/seisflows",
      author="Seisflows Development Team",
      packages=find_packages(),
      entry_points={
          "console_scripts": ["seisflows=seisflows.scripts.seisflows:main",]},
      license="GPL",
      zip_save=False
      )
