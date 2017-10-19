#Complilation instructions

For python use
`python Q2_Perturbation.py`
The plot will be generated and displayed.

For FORTRAN95 use
`gfortran Q2_Pert.f95`
Provide the the range and step-size of omega in REAL format. The data is
generated in 'Q2_Pert.dat'. Plot using in gnuplot
`plot 'Q2_Pert.dat' u 1:2 w l`
