from setuptools import setup, find_packages

setup(name="seisflows3",
      version="1.0.0",
      description="SeisFlows3: A seismic inversion package",
      url="https://github.com/seisflows/seisflows",
      author="Seisflows Development Team",
      packages=find_packages(),
      entry_points={
        "console_scripts": ["seisflows=seisflows.scripts.seisflows:main",
                            "sfsubmit=seisflows.scripts.seisflows:submit",
                            "sfresume=seisflows.scripts.seisflows:resume",
                            "sfclean=seisflows.scripts.seisflows:clean",
                            "sfdebug=seisflows.scripts.seisflows:debug",
                            "sfrestart=seisflows.scripts.seisflows:restart",
                            "sfcheck=seisflows.scripts.sfcheck:main"]},
      license="GPL",
      zip_save=False
      )
