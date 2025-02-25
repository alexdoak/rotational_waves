# rotational_waves
Codes to compute free-surface waves with arbitrary vorticity. This is the code used in the paper ***DOI HERE***. 

There are a few notation changes from the paper, in particular 
  (alpha,beta) --> (s,t)
  gamma --> vort_fun


The code is designed to solve for the independent variables psi and y, both two-dimensional fields in the rectangular domain (s,t). They solve the coupled PDE systems given by equations (2.3) and (2.4). The code makes heavy use of MATLAB struct's, a form of data storage. The labelling of these is roughly

pa - Parameters of the flow.
un - Unknowns of the flow 
me - Variables related to the computational mesh
fe - Choices for the solver. (see below).
ja - Variables related to the Jacobian
de - Derivatives of variables
li - Lists, data storage along bifurcation curves

Note that in fe, there are a number of features that do not concern the paper. In this sense, please set:
fe.Bathymetry=0, fe.Stratified=0, fe.Embed=0, fe.Embedloc=1, fe.Scont = 0, fe.Freesurface=1, fe.Mapping=0.

Meanwhile, the other fe values are
fe.Fixdepth = 0,1: Choose whether you fix the value of Y at the M*N meshpoint as pa.H
fe.Fixwavelength = 0,1: Choose whether you fix the wavelength of the solution to pa.wavelength
fe.Fix = 'Q', 'amplitude': Choose whether you are fixing the flux to pa.fixQ, or are fixing the amplitude to Amp=**.
fe.Ghostpoints = 0,1: Choose ...




