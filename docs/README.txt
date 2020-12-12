Docs are generated using sphinx.

Three files are critical to generating the docs:
1.) source/conf.py
2.) source/index.rst
3.) source/build_index.py

conf.py speficies configuration variables for Sphinx. It lists necessary extensions and sets default variables.

index.rst specifies the sections in the sphinx output.

build_index.py generates the index.rst file. build_index iterates recursively across all modules and extracts the name
of each class or function. Each class is then listed as its own section in the docs.

build_index.py was written in order to include classes in the table of contents. While the autosummary extension can
summarize all of the classes in a module at the beginning of each module section, it will not generate separate section
headers for each class that can be included in the table of contents on the first page.

In order to generate the index.rst file, navigate to docs/source and execute the build_index.py file.

In order to build the sphinx docs in latex, navigate to the docs/build directory and execute this command:
    sphinx-build -b latex ../source .

In order to generate a pdf of the docs, remain in the docs/build directory and execute this command:
    pdflatex feast.tex
