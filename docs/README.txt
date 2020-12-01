Docs are generated using sphinx.

Three files are critical to generating the docs:
1.) source/conf.py
2.) source/index.rst
3.) source/build_index.py

conf.py speficies configuration variables for Sphinx. It lists necessary extensions and sets default variables.

index.rst specifies the sections in the sphinx output.

build_index.py generates the index.rst file. build_index iterates recursively across all modules and extracts the name
of each class or function. Each class is then listed as its own section in the docs.