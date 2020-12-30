.TH RECIPES LOCAL
.SH NAME
.I Numerical Recipes
.SH SYNOPSIS
librecipes_f.a - Numerical recipes in Fortran
.br
librecipes_c.a - Numerical recipes in C (ANSI version)
.br
librecipes_ckr.a - Numerical recipes in C (K&R version)
.SH DESCRIPTION
The
.I Numerical Recipes
is a collection of numerical routines which are useful for many
scientific and engineering applications. Versions exist for both the C
and Fortran languages. The book
.I Numerical Recipes
will tell you everything you want to know.
.PP
The routines are called from the user's program.
At compilation, the
"-lrecipes_f"
or
"-lrecipes_c -lm"
or
"-lrecipes_ckr -lm"
option is added. Thus a program
.I myprog.f,
which contains calls to the recipes, is compiled with:
.I
f77 myprog.f -lrecipes_f
.br
Note that the C compiler does not automatically include the math library
the way the Fortran compiler does, so -lm is often required.
.br
.SH FILES
/usr/local/lib/librecipes_f.a
.br
/usr/local/lib/librecipes_c.a
.br
/usr/local/lib/librecipes_ckr.a
.br
/usr/local/src/recipes
.br
/usr/local/src/recipes/recipes_f
.br
/usr/local/src/recipes/recipes_c-ansi
.br
/usr/local/src/recipes/recipes_c-kr
.PP
The directories /usr/local/src and /usr/local/lib are generic examples.
The actual locations may vary on your system.
.SH RECIPES
Numerical Recipes in the recipes directory (all .f)
.br
(main book section numbers are in parentheses)
.br
.TP 8
.I ADDINT
interpolate and add, used by \fIMGLIN\fR (19.6)
.br
.TP 8
.I AIRY
Airy functions (6.7)
.br
.TP 8
.I AMEBSA
simulated annealing in continuous spaces (10.9)
.br
.TP 8
.I AMOEBA
minimize in N\-dimensions by downhill simplex method (10.4)
.br
.TP 8
.I AMOTRY
evaluate a trial point, used by \fIAMOEBA\fR (10.4)
.br
.TP 8
.I AMOTSA
evaluate a trial point, used by \fIAMEBSA\fR (10.9)
.br
.TP 8
.I ANNEAL
traveling salesman problem by simulated annealing (10.9)
.br
.TP 8
.I ANORM2
utility used by \fIMGFAS\fR (19.6)
.br
.TP 8
.I ARCMAK
construct an arithmetic code (20.5)
.br
.TP 8
.I ARCODE
encode or decode a character using arithmetic coding (20.5)
.br
.TP 8
.I ARCSUM
add integer to byte string, used by \fIARCODE\fR (20.5)
.br
.TP 8
.I ASOLVE
used by \fILINBCG\fR for preconditioner (2.7)
.br
.TP 8
.I ATIMES
used by \fILINBCG\fR for sparse multiplication (2.7)
.br
.TP 8
.I AVEVAR
calculate mean and variance of a data set (14.2)
.br
.TP 8
.I BADLUK
Friday the 13th when the moon is full (1.1)
.br
.TP 8
.I BALANC
balance a nonsymmetric matrix (11.5)
.br
.TP 8
.I BANBKS
band diagonal systems, backsubstitution (2.4)
.br
.TP 8
.I BANDEC
band diagonal systems, decomposition (2.4)
.br
.TP 8
.I BANMUL
multiply vector by band diagonal matrix (2.4)
.br
.TP 8
.I BCUCOF
construct two\-dimensional bicubic (3.6)
.br
.TP 8
.I BCUINT
two\-dimensional bicubic interpolation (3.6)
.br
.TP 8
.I BESCHB
Chebyshev expansion used by \fIBESSJY\fR (6.7)
.br
.TP 8
.I BESSI
modified Bessel function I of integer order (6.6)
.br
.TP 8
.I BESSI0
modified Bessel function I_0 (6.6)
.br
.TP 8
.I BESSI1
modified Bessel function I_1 (6.6)
.br
.TP 8
.I BESSIK
modified Bessel functions of fractional order (6.7)
.br
.TP 8
.I BESSJ
Bessel function J of general integer order (6.5)
.br
.TP 8
.I BESSJ0
Bessel function J_0 (6.5)
.br
.TP 8
.I BESSJ1
Bessel function J_1 (6.5)
.br
.TP 8
.I BESSJY
Bessel functions of fractional order (6.7)
.br
.TP 8
.I BESSK
modified Bessel function K of integer order (6.6)
.br
.TP 8
.I BESSK0
modified Bessel function K_0 (6.6)
.br
.TP 8
.I BESSK1
modified Bessel function K_1 (6.6)
.br
.TP 8
.I BESSY
Bessel function Y of general integer order (6.5)
.br
.TP 8
.I BESSY0
Bessel function Y_0 (6.5)
.br
.TP 8
.I BESSY1
Bessel function Y_1 (6.5)
.br
.TP 8
.I BETA
beta function (6.1)
.br
.TP 8
.I BETACF
continued fraction used by \fIBETAI\fR (6.4)
.br
.TP 8
.I BETAI
incomplete beta function (6.4)
.br
.TP 8
.I BICO
binomial coefficients function (6.1)
.br
.TP 8
.I BKSUB
backsubstitution, used by \fISOLVDE\fR (17.3)
.br
.TP 8
.I BNLDEV
binomial distributed random deviates (7.3)
.br
.TP 8
.I BRENT
find minimum of a function by Brent's method (10.2)
.br
.TP 8
.I BROYDN
secant method for systems of equations (9.7)
.br
.TP 8
.I BSSTEP
integrate ODEs, Bulirsch\-Stoer step (16.4)
.br
.TP 8
.I CALDAT
calendar date from Julian day number (1.1)
.br
.TP 8
.I CHDER
derivative of a function already Chebyshev fitted (5.9)
.br
.TP 8
.I CHEBEV
Chebyshev polynomial evaluation (5.8)
.br
.TP 8
.I CHEBFT
fit a Chebyshev polynomial to a function (5.8)
.br
.TP 8
.I CHEBPC
polynomial coefficients from a Chebyshev fit (5.10)
.br
.TP 8
.I CHINT
integrate a function already Chebyshev fitted (5.9)
.br
.TP 8
.I CHIXY
used by \fIFITEXY\fR to calculate a chi^2 (15.3)
.br
.TP 8
.I CHOLDC
Cholesky decomposition (2.9)
.br
.TP 8
.I CHOLSL
Cholesky backsubstitution (2.9)
.br
.TP 8
.I CHSONE
chi\-square test for difference between data and model (14.3)
.br
.TP 8
.I CHSTWO
chi\-square test for difference between two data sets (14.3)
.br
.TP 8
.I CISI
cosine and sine integrals Ci and Si (6.9)
.br
.TP 8
.I CNTAB1
contingency table analysis using chi\-square (14.4)
.br
.TP 8
.I CNTAB2
contingency table analysis using entropy measure (14.4)
.br
.TP 8
.I CONVLV
convolution or deconvolution of data using FFT (13.1)
.br
.TP 8
.I COPY
utility used by \fIMGLIN\fR, \fIMGFAS\fR (19.6)
.br
.TP 8
.I CORREL
correlation or autocorrelation of data using FFT (13.2)
.br
.TP 8
.I COSFT1
fast cosine transform with endpoints (12.3)
.br
.TP 8
.I COSFT2
``staggered'' fast cosine transform (12.3)
.br
.TP 8
.I COVSRT
rearrange covariance matrix, used by \fILFIT\fR (15.4)
.br
.TP 8
.I CRANK
replaces array elements by their rank (14.6)
.br
.TP 8
.I CYCLIC
solution of cyclic tridiagonal systems (2.7)
.br
.TP 8
.I DAUB4
Daubechies 4\-coefficient wavelet filter (13.10)
.br
.TP 8
.I DAWSON
Dawson's integral (6.10)
.br
.TP 8
.I DBRENT
find minimum of a function using derivative information (10.3)
.br
.TP 8
.I DDPOLY
evaluate a polynomial and its derivatives (5.3)
.br
.TP 8
.I DECCHK
decimal check digit calculation or verification (20.3)
.br
.TP 8
.I DERIVS
sample derivatives routine for \fISTIFF\fR (16.6)
.br
.TP 8
.I DF1DIM
alternative function used by \fILINMIN\fR (10.6)
.br
.TP 8
.I DFPMIN
minimize in N\-dimensions by variable metric method (10.7)
.br
.TP 8
.I DFRIDR
numerical derivative by Ridders' method (5.7)
.br
.TP 8
.I DFTCOR
compute endpoint corrections for Fourier integrals (13.9)
.br
.TP 8
.I DFTINT
high\-accuracy Fourier integrals (13.9)
.br
.TP 8
.I DIFEQ
spheroidal matrix coefficients, used by \fISFROID\fR (17.4)
.br
.TP 8
.I DPYTHAG, DREALFT, DSPRSAX, DSPRSTX, DSVBKSB, DSVDCMP
double precision versions of routines minus initial D
.br
.TP 8
.I ECLASS
determine equivalence classes from list (8.6)
.br
.TP 8
.I ECLAZZ
determine equivalence classes from procedure (8.6)
.br
.TP 8
.I EI
exponential integral Ei (6.3)
.br
.TP 8
.I EIGSRT
eigenvectors, sorts into order by eigenvalue (11.1)
.br
.TP 8
.I ELLE
Legendre elliptic integral of the second kind (6.11)
.br
.TP 8
.I ELLF
Legendre elliptic integral of the first kind (6.11)
.br
.TP 8
.I ELLPI
Legendre elliptic integral of the third kind (6.11)
.br
.TP 8
.I ELMHES
reduce a general matrix to Hessenberg form (11.5)
.br
.TP 8
.I ERF
error function (6.2)
.br
.TP 8
.I ERFC
complementary error function (6.2)
.br
.TP 8
.I ERFCC
complementary error function, concise routine (6.2)
.br
.TP 8
.I EULSUM
sum a series by Euler\-van Wijngaarden algorithm (5.1)
.br
.TP 8
.I EVLMEM
power spectral estimation from MEM coefficients (13.7)
.br
.TP 8
.I EXPDEV
exponential random deviates (7.2)
.br
.TP 8
.I EXPINT
exponential integral E_n (6.3)
.br
.TP 8
.I F1DIM
function used by \fILINMIN\fR (10.5)
.br
.TP 8
.I FACTLN
logarithm of factorial function (6.1)
.br
.TP 8
.I FACTRL
factorial function (6.1)
.br
.TP 8
.I FASPER
power spectrum of unevenly sampled larger data sets (13.8)
.br
.TP 8
.I FDJAC
finite\-difference Jacobian, used by \fINEWT\fR (9.7)
.br
.TP 8
.I FGAUSS
fit a sum of Gaussians using \fIMRQMIN\fR (15.5)
.br
.TP 8
.I FILL0
utility used by \fIMGLIN\fR (19.6)
.br
.TP 8
.I FIT
least\-squares fit data to a straight line (15.2)
.br
.TP 8
.I FITEXY
fit data to a straight line, errors in both x and y (15.3)
.br
.TP 8
.I FIXRTS
reflect roots of a polynomial into unit circle (13.6)
.br
.TP 8
.I FLEG
fit a Legendre polynomial using \fILFIT\fR or \fISVDFIT\fR (15.4)
.br
.TP 8
.I FLMOON
calculate phases of the moon by date (1.0)
.br
.TP 8
.I FMIN
norm of a vector function, used by \fINEWT\fR (9.7)
.br
.TP 8
.I FOUR1
fast Fourier transform (FFT) in one dimension (12.2)
.br
.TP 8
.I FOUREW
rewind and permute files, used by \fIFOURFS\fR (12.6)
.br
.TP 8
.I FOURFS
FFT for huge data sets on external media (12.6)
.br
.TP 8
.I FOURN
fast Fourier transform in multidimensions (12.4)
.br
.TP 8
.I FPOLY
fit a polynomial using \fILFIT\fR or \fISVDFIT\fR (15.4)
.br
.TP 8
.I FRED2
solve linear Fredholm equations of the second kind (18.1)
.br
.TP 8
.I FREDEX
example of solving a singular Fredholm equation (18.3)
.br
.TP 8
.I FREDIN
interpolate solutions obtained with \fIfred2\fR (18.1)
.br
.TP 8
.I FRENEL
Fresnel integrals S(x) and C(x) (6.9)
.br
.TP 8
.I FRPRMN
minimize in N\-dimensions by conjugate gradient (10.6)
.br
.TP 8
.I FTEST
F\-test for difference of variances (14.2)
.br
.TP 8
.I GAMDEV
gamma\-law distribution random deviates (7.3)
.br
.TP 8
.I GAMMLN
logarithm of gamma function (6.1)
.br
.TP 8
.I GAMMP
incomplete gamma function (6.2)
.br
.TP 8
.I GAMMQ
complement of incomplete gamma function (6.2)
.br
.TP 8
.I GASDEV
normally distributed random deviates (7.2)
.br
.TP 8
.I GAUCOF
quadrature weights from orthogonal polynomials (4.5)
.br
.TP 8
.I GAUHER
Gauss\-Hermite weights and abscissas (4.5)
.br
.TP 8
.I GAUJAC
Gauss\-Jacobi weights and abscissas (4.5)
.br
.TP 8
.I GAULAG
Gauss\-Laguerre weights and abscissas (4.5)
.br
.TP 8
.I GAULEG
Gauss\-Legendre weights and abscissas (4.5)
.br
.TP 8
.I GAUSSJ
Gauss\-Jordan matrix inversion and linear equation solution (2.1)
.br
.TP 8
.I GCF
continued fraction used by \fIGAMMP\fR and \fIGAMMQ\fR (6.2)
.br
.TP 8
.I GOLDEN
find minimum of a function by golden section search (10.1)
.br
.TP 8
.I GSER
series used by \fIGAMMP\fR and \fIGAMMQ\fR (6.2)
.br
.TP 8
.I HPSEL
find M largest values, without altering an array (8.5)
.br
.TP 8
.I HPSORT
sort an array by heapsort method (8.3)
.br
.TP 8
.I HQR
eigenvalues of a Hessenberg matrix (11.6)
.br
.TP 8
.I HUFAPP
append bits to a Huffman code, used by \fIHUFMAK\fR (20.4)
.br
.TP 8
.I HUFDEC
use Huffman code to decode and decompress a character (20.4)
.br
.TP 8
.I HUFENC
use Huffman code to encode and compress a character (20.4)
.br
.TP 8
.I HUFMAK
construct a Huffman code (20.4)
.br
.TP 8
.I HUNT
search a table when calls are correlated (3.4)
.br
.TP 8
.I HYPDRV
complex hypergeometric function, derivative of (6.12)
.br
.TP 8
.I HYPGEO
complex hypergeometric function (6.12)
.br
.TP 8
.I HYPSER
complex hypergeometric function, series evaluation (6.12)
.br
.TP 8
.I ICRC
cyclic redundancy checksum (20.3)
.br
.TP 8
.I ICRC1
cyclic redundancy checksum, used by \fIICRC\fR (20.3)
.br
.TP 8
.I IGRAY
Gray code and its inverse (20.2)
.br
.TP 8
.I IINDEXX
integer array version of \fIINDEXX\fR (8.4)
.br
.TP 8
.I INDEXX
construct an index for an array (8.4)
.br
.TP 8
.I INTERP
bilinear prolongation, used by \fIMGLIN\fR, \fIMGFAS\fR (19.6)
.br
.TP 8
.I IRBIT1
random bit sequence (7.4)
.br
.TP 8
.I IRBIT2
random bit sequence (7.4)
.br
.TP 8
.I JACOBI
eigenvalues and eigenvectors of a symmetric matrix (11.1)
.br
.TP 8
.I JACOBN
sample Jacobian routine for \fISTIFF\fR (16.6)
.br
.TP 8
.I JULDAY
Julian Day number from calendar date (1.1)
.br
.TP 8
.I KENDL1
correlation between two data sets, Kendall's tau (14.6)
.br
.TP 8
.I KENDL2
contingency table analysis using Kendall's tau (14.6)
.br
.TP 8
.I KERMOM
sample routine for moments of a singular kernel (18.3)
.br
.TP 8
.I KS2D1S
K\-S test in two dimensions, data vs. model (14.7)
.br
.TP 8
.I KS2D2S
K\-S test in two dimensions, data vs. data (14.7)
.br
.TP 8
.I KSONE
Kolmogorov\-Smirnov test of data against model (14.3)
.br
.TP 8
.I KSTWO
Kolmogorov\-Smirnov test between two data sets (14.3)
.br
.TP 8
.I LAGUER
find a root of a polynomial by Laguerre's method (9.5)
.br
.TP 8
.I LFIT
general linear least\-squares fit by normal equations (15.4)
.br
.TP 8
.I LINBCG
biconjugate gradient solution of sparse systems (2.7)
.br
.TP 8
.I LINMIN
minimum of a function along a ray in N\-dimensions (10.5)
.br
.TP 8
.I LNSRCH
search along a line, used by \fINEWT\fR (9.7)
.br
.TP 8
.I LOCATE
search an ordered table by bisection (3.4)
.br
.TP 8
.I LOP
applies nonlinear operator, used by \fIMGFAS\fR (19.6)
.br
.TP 8
.I LUBKSB
linear equation solution, backsubstitution (2.3)
.br
.TP 8
.I LUDCMP
linear equation solution, LU decomposition (2.3)
.br
.TP 8
.I MACHAR
diagnose computer's floating arithmetic (20.1)
.br
.TP 8
.I MALOC
memory allocation utility used by \fIMGLIN\fR, \fIMGFAS\fR (19.6)
.br
.TP 8
.I MATADD
utility used by \fIMGFAS\fR (19.6)
.br
.TP 8
.I MATSUB
utility used by \fIMGFAS\fR (19.6)
.br
.TP 8
.I MEDFIT
fit data to a straight line robustly, least absolute deviation (15.7)
.br
.TP 8
.I MEMCOF
evaluate maximum entropy (MEM) coefficients (13.6)
.br
.TP 8
.I METROP
Metropolis algorithm, used by \fIANNEAL\fR (10.9)
.br
.TP 8
.I MGFAS
nonlinear elliptic PDE solved by multigrid method (19.6)
.br
.TP 8
.I MGLIN
linear elliptic PDE solved by multigrid method (19.6)
.br
.TP 8
.I MIDEXP
integrate a function that decreases exponentially (4.4)
.br
.TP 8
.I MIDINF
integrate a function on a semi\-infinite interval (4.4)
.br
.TP 8
.I MIDPNT
extended midpoint rule (4.4)
.br
.TP 8
.I MIDSQL
integrate a function with lower square\-root singularity (4.4)
.br
.TP 8
.I MIDSQU
integrate a function with upper square\-root singularity (4.4)
.br
.TP 8
.I MISER
recursive multidimensional Monte Carlo integration (7.8)
.br
.TP 8
.I MMID
integrate ODEs by modified midpoint method (16.3)
.br
.TP 8
.I MNBRAK
bracket the minimum of a function (10.1)
.br
.TP 8
.I MNEWT
Newton's method for systems of equations (9.6)
.br
.TP 8
.I MOMENT
calculate moments of a data set (14.1)
.br
.TP 8
.I MP2DFR
multiple precision conversion to decimal base (20.6)
.br
.TP 8
.I MPDIV
multiple precision divide and remainder (20.6)
.br
.TP 8
.I MPINV
multiple precision reciprocal (20.6)
.br
.TP 8
.I MPMUL
multiple precision multiply, using FFT methods (20.6)
.br
.TP 8
.I MPOPS
multiple precision arithmetic, simpler operations (20.6)
.br
.TP 8
.I MPPI
multiple precision example, compute many digits of pi (20.6)
.br
.TP 8
.I MPROVE
linear equation solution, iterative improvement (2.5)
.br
.TP 8
.I MPSQRT
multiple precision square root (20.6)
.br
.TP 8
.I MRQCOF
used by \fIMRQMIN\fR to evaluate coefficients (15.5)
.br
.TP 8
.I MRQMIN
nonlinear least\-squares fit, Marquardt's method (15.5)
.br
.TP 8
.I NEWT
globally convergent multi\-dimensional Newton's method (9.7)
.br
.TP 8
.I ODEINT
integrate ODEs with accuracy monitoring (16.2)
.br
.TP 8
.I ORTHOG
construct nonclassical orthogonal polynomials (4.5)
.br
.TP 8
.I PADE
Pade approximant from power series coefficients (5.12)
.br
.TP 8
.I PCCHEB
inverse of \fICHEBPC\fR; use to economize power series (5.11)
.br
.TP 8
.I PCSHFT
polynomial coefficients of a shifted polynomial (5.10)
.br
.TP 8
.I PEARSN
Pearson's correlation between two data sets (14.5)
.br
.TP 8
.I PERIOD
power spectrum of unevenly sampled data (13.8)
.br
.TP 8
.I PIKSR2
sort two arrays by straight insertion (8.1)
.br
.TP 8
.I PIKSRT
sort an array by straight insertion (8.1)
.br
.TP 8
.I PINVS
diagonalize a sub\-block, used by \fISOLVDE\fR (17.3)
.br
.TP 8
.I PLGNDR
Legendre polynomials, associated (spherical harmonics) (6.8)
.br
.TP 8
.I POIDEV
Poisson distributed random deviates (7.3)
.br
.TP 8
.I POLCOE
polynomial coefficients from table of values (3.5)
.br
.TP 8
.I POLCOF
polynomial coefficients from table of values (3.5)
.br
.TP 8
.I POLDIV
divide one polynomial by another (5.3)
.br
.TP 8
.I POLIN2
two\-dimensional polynomial interpolation (3.6)
.br
.TP 8
.I POLINT
polynomial interpolation (3.1)
.br
.TP 8
.I POWELL
minimize in N\-dimensions by Powell's method (10.5)
.br
.TP 8
.I PREDIC
linear prediction using MEM coefficients (13.6)
.br
.TP 8
.I PROBKS
Kolmogorov\-Smirnov probability function (14.3)
.br
.TP 8
.I PSDES
pseudo\-DES hashing of 64 bits (7.5)
.br
.TP 8
.I PWT
partial wavelet transform (13.10)
.br
.TP 8
.I PWTSET
initialize coefficients for \fIPWT\fR (13.10)
.br
.TP 8
.I PYTHAG
calculate (a^2+b^2)^1/2 without overflow (2.6)
.br
.TP 8
.I PZEXTR
polynomial extrapolation, used by \fIBSSTEP\fR (16.4)
.br
.TP 8
.I QGAUS
integrate a function by Gaussian quadratures (4.5)
.br
.TP 8
.I QRDCMP
QR decomposition (2.10)
.br
.TP 8
.I QROMB
integrate using Romberg adaptive method (4.3)
.br
.TP 8
.I QROMO
integrate using open Romberg adaptive method (4.4)
.br
.TP 8
.I QROOT
complex or double root of a polynomial, Bairstow (9.5)
.br
.TP 8
.I QRSOLV
QR backsubstitution (2.10)
.br
.TP 8
.I QRUPDT
update a QR decomposition (2.10)
.br
.TP 8
.I QSIMP
integrate using Simpson's rule (4.2)
.br
.TP 8
.I QTRAP
integrate using trapezoidal rule (4.2)
.br
.TP 8
.I QUAD3D
integrate a function over a three\-dimensional space (4.6)
.br
.TP 8
.I QUADCT
count points by quadrants, used by \fIks2d1s\fR (14.7)
.br
.TP 8
.I QUADMX
sample routine for a quadrature matrix (18.3)
.br
.TP 8
.I QUADVL
quadrant probabilities, used by \fIks2d1s\fR (14.7)
.br
.TP 8
.I RAN0
random deviate by Park and Miller minimal standard (7.1)
.br
.TP 8
.I RAN1
random deviate, minimal standard plus shuffle (7.1)
.br
.TP 8
.I RAN2
random deviate by L'Ecuyer long period plus shuffle (7.1)
.br
.TP 8
.I RAN3
random deviate by Knuth subtractive method (7.1)
.br
.TP 8
.I RAN4
random deviates from DES\-like hashing (7.5)
.br
.TP 8
.I RANK
construct a rank table for an array (8.4)
.br
.TP 8
.I RANPT
get random point, used by \fIMISER\fR (7.8)
.br
.TP 8
.I RATINT
rational function interpolation (3.2)
.br
.TP 8
.I RATLSQ
rational fit by least\-squares method (5.13)
.br
.TP 8
.I RATVAL
evaluate a rational function (5.3)
.br
.TP 8
.I RC
Carlson's degenerate elliptic integral (6.11)
.br
.TP 8
.I RD
Carlson's elliptic integral of the second kind (6.11)
.br
.TP 8
.I REALFT
fast Fourier transform of a single real function (12.3)
.br
.TP 8
.I REBIN
sample rebinning used by \fIVEGAS\fR (7.8)
.br
.TP 8
.I RED
reduce columns of a matrix, used by \fISOLVDE\fR (17.3)
.br
.TP 8
.I RELAX
Gauss\-Seidel relaxation, used by \fIMGLIN\fR (19.6)
.br
.TP 8
.I RELAX2
Gauss\-Seidel relaxation, used by \fIMGFAS\fR (19.6)
.br
.TP 8
.I RESID
calculate residual, used by \fIMGLIN\fR (19.6)
.br
.TP 8
.I REVCST
cost of a reversal, used by \fIANNEAL\fR (10.9)
.br
.TP 8
.I REVERS
do a reversal, used by \fIANNEAL\fR (10.9)
.br
.TP 8
.I RF
Carlson's elliptic integral of the first kind (6.11)
.br
.TP 8
.I RJ
Carlson's elliptic integral of the third kind (6.11)
.br
.TP 8
.I RK4
integrate one step of ODEs, fourth\-order Runge\-Kutta (16.1)
.br
.TP 8
.I RKCK
Cash\-Karp\-Runge\-Kutta step used by \fIRKQS\fR (16.2)
.br
.TP 8
.I RKDUMB
integrate ODEs by fourth\-order Runge\-Kutta (16.1)
.br
.TP 8
.I RKQS
integrate one step of ODEs with accuracy monitoring (16.2)
.br
.TP 8
.I RLFT3
FFT of real data in two or three dimensions (12.5)
.br
.TP 8
.I ROFUNC
fit data robustly, used by \fIMEDFIT\fR (15.7)
.br
.TP 8
.I ROTATE
Jacobi rotation used by \fIQRUPDT\fR (2.10)
.br
.TP 8
.I RSOLV
right triangular backsubstitution (2.10)
.br
.TP 8
.I RSTRCT
half\-weighting restriction, used by \fIMGLIN\fR, \fIMGFAS\fR (19.6)
.br
.TP 8
.I RTBIS
find root of a function by bisection (9.1)
.br
.TP 8
.I RTFLSP
find root of a function by false\-position (9.2)
.br
.TP 8
.I RTNEWT
find root of a function by Newton\-Raphson (9.4)
.br
.TP 8
.I RTSAFE
find root of a function by Newton\-Raphson and bisection (9.4)
.br
.TP 8
.I RTSEC
find root of a function by secant method (9.2)
.br
.TP 8
.I RZEXTR
rational function extrapolation, used by \fIBSSTEP\fR (16.4)
.br
.TP 8
.I SAVGOL
Savitzky\-Golay smoothing coefficients (14.8)
.br
.TP 8
.I SCRSHO
graph a function to search for roots (9.0)
.br
.TP 8
.I SELECT
find the Nth largest in an array (8.5)
.br
.TP 8
.I SELIP
find the Nth largest, without altering an array (8.5)
.br
.TP 8
.I SFROID
spheroidal functions by method of \fISOLVDE\fR (17.4)
.br
.TP 8
.I SHELL
sort an array by Shell's method (8.1)
.br
.TP 8
.I SHOOT
solve two point boundary value problem by shooting (17.1)
.br
.TP 8
.I SHOOTF
ditto, by shooting to a fitting point (17.2)
.br
.TP 8
.I SIMP1
linear programming, used by \fISIMPLX\fR (10.8)
.br
.TP 8
.I SIMP2
linear programming, used by \fISIMPLX\fR (10.8)
.br
.TP 8
.I SIMP3
linear programming, used by \fISIMPLX\fR (10.8)
.br
.TP 8
.I SIMPLX
linear programming maximization of a linear function (10.8)
.br
.TP 8
.I SIMPR
integrate stiff ODEs by semi\-implicit midpoint rule (16.6)
.br
.TP 8
.I SINFT
fast sine transform (12.3)
.br
.TP 8
.I SLVSM2
solve on coarsest grid, used by \fIMGFAS\fR (19.6)
.br
.TP 8
.I SLVSML
solve on coarsest grid, used by \fIMGLIN\fR (19.6)
.br
.TP 8
.I SNCNDN
Jacobian elliptic functions (6.11)
.br
.TP 8
.I SNRM
used by \fILINBCG\fR for vector norm (2.7)
.br
.TP 8
.I SOBSEQ
Sobol's quasi\-random sequence (7.7)
.br
.TP 8
.I SOLVDE
two point boundary value problem, solve by relaxation (17.3)
.br
.TP 8
.I SOR
elliptic PDE solved by successive overrelaxation method (19.5)
.br
.TP 8
.I SORT
sort an array by quicksort method (8.2)
.br
.TP 8
.I SORT2
sort two arrays by quicksort method (8.2)
.br
.TP 8
.I SORT3
sort, use an index to sort 3 or more arrays (8.4)
.br
.TP 8
.I SPCTRM
power spectrum estimation using FFT (13.4)
.br
.TP 8
.I SPEAR
Spearman's rank correlation between two data sets (14.6)
.br
.TP 8
.I SPHBES
spherical Bessel functions j_n and y_n (6.7)
.br
.TP 8
.I SPHFPT
spheroidal functions by method of \fISHOOTF\fR (17.4)
.br
.TP 8
.I SPHOOT
spheroidal functions by method of \fISHOOT\fR (17.4)
.br
.TP 8
.I SPLIE2
construct two\-dimensional spline (3.6)
.br
.TP 8
.I SPLIN2
two\-dimensional spline interpolation (3.6)
.br
.TP 8
.I SPLINE
construct a cubic spline (3.3)
.br
.TP 8
.I SPLINT
cubic spline interpolation (3.3)
.br
.TP 8
.I SPREAD
extirpolate value into array, used by \fIFASPER\fR (13.8)
.br
.TP 8
.I SPRSAX
product of sparse matrix and vector (2.7)
.br
.TP 8
.I SPRSIN
convert matrix to sparse format (2.7)
.br
.TP 8
.I SPRSPM
pattern multiply two sparse matrices (2.7)
.br
.TP 8
.I SPRSTM
threshold multiply two sparse matrices (2.7)
.br
.TP 8
.I SPRSTP
transpose of sparse matix (2.7)
.br
.TP 8
.I SPRSTX
product of transpose sparse matrix and vector (2.7)
.br
.TP 8
.I STIFBS
integrate stiff ODEs, Bulirsch\-Stoer step (16.6)
.br
.TP 8
.I STIFF
integrate stiff ODEs by fourth\-order Rosenbrock (16.6)
.br
.TP 8
.I STOERM
integrate conservative second\-order ODEs (16.5)
.br
.TP 8
.I SVBKSB
singular value backsubstitution (2.6)
.br
.TP 8
.I SVDCMP
singular value decomposition of a matrix (2.6)
.br
.TP 8
.I SVDFIT
linear least\-squares fit by singular value decomposition (15.4)
.br
.TP 8
.I SVDVAR
variances from singular value decomposition (15.4)
.br
.TP 8
.I TOEPLZ
solve Toeplitz systems (2.8)
.br
.TP 8
.I TPTEST
Student's t\-test for means, case of paired data (14.2)
.br
.TP 8
.I TQLI
eigensolution of a symmetric tridiagonal matrix (11.3)
.br
.TP 8
.I TRAPZD
trapezoidal rule (4.2)
.br
.TP 8
.I TRED2
Householder reduction of a real, symmetric matrix (11.2)
.br
.TP 8
.I TRIDAG
solution of tridiagonal systems (2.4)
.br
.TP 8
.I TRNCST
cost of a transposition, used by \fIANNEAL\fR (10.9)
.br
.TP 8
.I TRNSPT
do a transposition, used by \fIANNEAL\fR (10.9)
.br
.TP 8
.I TTEST
Student's t\-test for difference of means (14.2)
.br
.TP 8
.I TUTEST
Student's t\-test for means, case of unequal variances (14.2)
.br
.TP 8
.I TWOFFT
fast Fourier transform of two real functions (12.3)
.br
.TP 8
.I VANDER
solve Vandermonde systems (2.8)
.br
.TP 8
.I VEGAS
adaptive multidimensional Monte Carlo integration (7.8)
.br
.TP 8
.I VOLTRA
linear Volterra equations of the second kind (18.2)
.br
.TP 8
.I WT1
one\-dimensional discrete wavelet transform (13.10)
.br
.TP 8
.I WTN
multidimensional discrete wavelet transform (13.10)
.br
.TP 8
.I WWGHTS
quadrature weights for an arbitrarily singular kernel (18.3)
.br
.TP 8
.I ZBRAC
outward search for brackets on roots (9.1)
.br
.TP 8
.I ZBRAK
inward search for brackets on roots (9.1)
.br
.TP 8
.I ZBRENT
find root of a function by Brent's method (9.3)
.br
.TP 8
.I ZRHQR
roots of a polynomial by eigenvalue methods (9.5)
.br
.TP 8
.I ZRIDDR
find root of a function by Ridders' method (9.2)
.br
.TP 8
.I ZROOTS
roots of a polynomial by Laguerre's method with deflation (9.5)
.br
.SH "SEE ALSO"
.I
Numerical Recipes in Fortran, The Art of Scientific Computing (Second Edition),
W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling,
Cambridge University Press, 1992
.br
.I
Numerical Recipes in C, The Art of Scientific Computing (Second Edition),
W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling,
Cambridge University Press, 1992
.SH BUGS
See associated documentation.
