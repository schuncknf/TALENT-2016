/*******************************************/
/                                           /
/             Damn Ducks Towels'            /
/                                           /
/       Spherical Hartree-Fock project      /
/                                           /
/*******************************************/


1_ INSTALL

Before using our solver, please install the following softwares and
libraries using your distribution's package manager:
    Armadillo (>= 6.500)
    Boost
    Doxygen
    Python 2.7
    Graphviz
    Dia
    gsl
    swig3.0 (>= 3.0.8)
    astyle
    Python distutils (python-setuptools with apt)
    gcc and g++ (>= 5.4.0)
    


2_ COMPILE

You can compile the project by using 
    make
in the main directory. Use
    make mrproper
if you want to clean everything.


3_ RUN

To run the solver, please modify the /python/config.py file according
to the physics system you want to study. Once done, execute
    python ./python/iter.py


4_ DOCUMENTATION

To produce a complete documentation on the Hartree-Fock solver and of
the different functions implemented, you can simply execute in the
main directory :
    make doc
This produces a html documentation in /doc/html/index.html and a LaTeX
file in /doc/latex/refman.tex, which you can further compile using
latex, or better pdflatex.
An UML diagram is available in /doc using Dia.
