from setuptools import setup

setup(name='seapipy',
      version='0.1.0',
      author='Augusto Borges',
      author_email='borges.augustoar@gmail.com',
      packages=['seapipy'],
      url='https://github.com/borgesaugusto/seapi',
      license='LICENSE',
      description='A Surface evolver API for python',
      long_description=open('README.md').read(),
      install_requires=[
          "pytest",
      ],
      platforms=["nt", "posix"]
      )
