#!/usr/bin/env python3

# MCFIT    NDD     21/09/2019

# Python alternative to the Fortran bootstrap Monte Carlo fitting
# program.  Much shorter, possibly easier to use, and much slower.  To
# insert a new model, parameters or set of derived properties, look
# for the text "YOU NEED TO CHANGE THIS IF YOU HAVE A NEW MODEL".

import sys
import textwrap as tw
import numpy as np
from scipy import optimize
from scipy import special

usage="Usage: mcfit.py [datafile]"

# YOU NEED TO CHANGE THIS IF YOU HAVE A NEW MODEL.
# Names of the parameters and an array of their (first-guess) initial values.
paramnames=("a","b","alpha")
initparams=np.array((0.0,2.0,1.0))
# Names of derived properties.  Avoid slashes, etc.: it should be
# possible for these names to be part of filenames.  If there are no
# derived properties, then just set propnames equal to ().
propnames=("a+b",)

# Number of samples for Monte Carlo resampling of data.
nsamples=2000
# Create a histogram of each derived property.
make_histogram=True


def model_y(x,a,b,alpha):
    """The model function.  This should be an elemental function of x, i.e., if
    x is an array then the function should return an array of model values.
    YOU NEED TO CHANGE THIS IF YOU HAVE A NEW MODEL."""
    return a+b*np.exp(-alpha*x)


def derived_props(a,b,alpha):
    """Return an array of parameter-dependent properties to average.
    YOU NEED TO CHANGE THIS IF YOU HAVE A NEW MODEL."""
    # Just return np.array(()) if there are no derived properties.
    return np.array((a+b,))


def errstop(message):
    """Report an error and halt."""
    print(f"\n{usage}\n")
    sys.exit(tw.fill(message))


def chisqfn(model_y,xx,yy,erryy,params):
    """The (non-reduced) chi^2 value."""
    return np.sum(((yy-model_y(xx,*params))/erryy)**2)


def displayparams(paramnames,params,errparams=None):
    """Write out a parameter set (with error bars)."""
    namewidth=max(len(p) for p in paramnames)
    if errparams is not None:
        for pn,p,errp in zip(paramnames,params,errparams):
            print(f"  {pn:>{namewidth}s}: {p} +/- {errp}")
    else:
        for pn,p in zip(paramnames,params):
            print(f"  {pn:>{namewidth}s}: {p}")
    print()


def displaychisq(chisq,ndof):
    """Report the reduced chi^2 function."""
    print(f"            Reduced chi^2: {chisq/ndof}")
    print(f"  Prob. chi^2 is this bad: {special.gammaincc(0.5*ndof,0.5*chisq)}"
          "\n")


def displayprops(propnames,props,errprops=None):
    """Write out a set of derived properties (with error bars)."""
    if len(propnames)>0:
        namewidth=max(len(p) for p in propnames)
        if errprops is not None:
            for pn,p,errp in zip(propnames,props,errprops):
                print(f"  {pn:>{namewidth}s}: {p} +/- {errp}")
        else:
            for pn,p in zip(propnames,props):
                print(f"  {pn:>{namewidth}s}: {p}")
        print()


def bootstrapfit(model_y,xx,yy,erryy,initparams,nprops):
    """Bootstrap Monte Carlo fitting."""
    mcavparams=np.zeros(len(initparams))
    mcavparamssq=np.zeros(len(initparams))
    mcavprops=np.zeros(nprops) ; mcavpropssq=np.zeros(nprops)
    if make_histogram:
        histfiles=[open(f"histogram_{p}.dat","w") for p in propnames]
    for _i in range(nsamples):
        yyp=np.random.normal(yy,erryy)
        popt,_pcov=optimize.curve_fit(model_y,xx,yyp,sigma=erryy,p0=initparams)
        props=derived_props(*popt)
        mcavparams+=popt ; mcavparamssq+=popt**2
        mcavprops+=props ; mcavpropssq+=props**2
        if make_histogram:
            for f,prop in zip(histfiles,props):
                f.write(str(prop)+"\n")
    if make_histogram:
        for f in histfiles:
            f.close()
    rec_nsamples=1.0/nsamples
    mcavparams*=rec_nsamples ; mcavparamssq*=rec_nsamples
    mcavprops*=rec_nsamples ; mcavpropssq*=rec_nsamples
    bessel_corr=float(nsamples)/float(nsamples-1)
    errmcavparams=np.sqrt((mcavparamssq-mcavparams**2)*bessel_corr)
    errmcavprops=np.sqrt((mcavpropssq-mcavprops**2)*bessel_corr)
    return mcavparams,errmcavparams,mcavprops,errmcavprops


def sortdata(xx,yy,erryy):
    """Sort data into ascending order.  Warn if the data are not already in
    ascending order (including repeated data points)."""
    tol=1.E-14
    if any(xx[i+1]-xx[i]<tol for i in range(len(xx)-1)):
        print("Warning: input data are not in ascending order.")
    ira=np.argsort(xx)
    xx=xx[ira] ; yy=yy[ira] ; erryy=erryy[ira]
    if any(abs(xx[i+1]-xx[i])<tol for i in range(len(xx)-1)):
        print("Warning: there are repeated data points.")


# Main program starts here.
print("\nMCFIT")
print("=====\n")

# Get filename and read data.
if len(sys.argv)>1:
    if len(sys.argv)>2:
        errstop("Only one or zero command-line arguments should be present.")
    fname=sys.argv[1]
else:
    print("Please enter a filename.")
    fname=sys.stdin.readline()
print(tw.fill(f"Reading data from {fname}."))
try:
    xx,yy,erryy=np.loadtxt(fname,unpack=True)
except (IOError,ValueError):
    errstop(f"Unable to read {fname}.")
if any(erryy<=0.0): errstop("Error bars should all be strictly positive.")

# Ensure data are in ascending order; warn if not.  Warn about repeated data.
sortdata(xx,yy,erryy)

# Sort out number of data points and number of parameters.
if len(initparams)!=len(paramnames):
    errstop("Initial parameter vector and name vector have different lengths.")
print(f"Number of data lines: {len(xx)}")
print(f"Number of parameters: {len(initparams)}\n")
ndof=len(xx)-len(initparams)
if ndof<0: errstop("Have more parameters than data points.")

print("Initial parameters:")
displayparams(paramnames,initparams)
displaychisq(chisqfn(model_y,xx,yy,erryy,initparams),ndof)

# Perform initial fit.
print("Optimised parameters:")
optparams,pcov=optimize.curve_fit(model_y,xx,yy,sigma=erryy,p0=initparams)
optprops=derived_props(*optparams)
displayparams(paramnames,optparams)
displaychisq(chisqfn(model_y,xx,yy,erryy,optparams),ndof)
displayprops(propnames,optprops)

# Perform bootstrap Monte Carlo fit.
print(f"Monte Carlo-averaged parameters (using {nsamples} samples):")
mcavparams,errmcavparams,mcavprop,errmcavprop\
    =bootstrapfit(model_y,xx,yy,erryy,optparams,len(propnames))
displayparams(paramnames,mcavparams,errmcavparams)
displaychisq(chisqfn(model_y,xx,yy,erryy,mcavparams),ndof)
displayprops(propnames,mcavprop,errmcavprop)
if make_histogram:
    for propname in propnames:
        print(f"  Histogram of {propname} written to histogram_{propname}.dat.")
    print()

print("Program finished.\n")
