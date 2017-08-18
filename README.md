# Bspline.py

Python/Numpy implementation of Bspline basis functions via Cox - de Boor algorithm.

Also provided are higher-order differentiation, collocation matrix generation, and a minimal procedural API (mainly for dealing with knot vectors) which may help in converting MATLAB codes.

# Usage

```python
import numpy
import bspline
import splinelab

## Spline setup and evaluation

p = 3              # order of spline (as-is; 3 = cubic)
nknots = 11        # number of knots to generate (here endpoints count only once)
tau = [0.1, 0.33]  # collocation sites (i.e. where to evaluate)

knots = numpy.linspace(0,1,nknots)  # create a knot vector without endpoint repeats
k     = splinelab.augknt(knots, p)  # add endpoint repeats as appropriate for spline order p
B     = bspline.Bspline(k, p)       # create spline basis of order p on knots k

A0 = B.collmat(tau)                 # collocation matrix for function value at sites tau
A2 = B.collmat(tau, deriv_order=2)  # collocation matrix for second derivative at sites tau

print( A0 )
print( A2 )

D3 = B.diff(order=3)  # third derivative of B as lambda x: ...
print( D3(0.4) )

D  = numpy.array( [D3(t) for t in tau], dtype=numpy.float64 )  # third derivative of B at sites tau


## Spline setup by defining collocation sites

ncolloc = 7
tau = numpy.linspace(0,1,ncolloc)  # These are the sites to which we would like to interpolate
k   = splinelab.aptknt(tau, p)     # Given the collocation sites, generate a knot vector
                                   # (incl. endpoint repeats). To get meaningful results,
                                   # here one must choose ncolloc such that  ncolloc >= p+1.
B   = bspline.Bspline(k, p)

A0  = B.collmat(tau)

print( A0 )


## Evaluate a function expressed in the spline basis:

# set up coefficients (in a real use case, fill this with something sensible,
#                      e.g. with an L2 projection onto the spline basis)
#
nbasis = A0.shape[1]  # A0.shape = (num_collocation_sites, num_basis_functions)
c = numpy.ones( (nbasis,), dtype=numpy.float64 )

# evaluate f(0.4)
y1 = numpy.sum( B(0.4) * c )

# evaluate at each tau[k]
y2 = numpy.array( [numpy.sum( B(t) * c ) for t in tau], dtype=numpy.float64 )

# equivalent, using the collocation matrix
#
# NOTE: the sites tau are built into the matrix when collmat() is called.
#
y3 = numpy.sum( A0 * c, axis=-1 )
```

# Installation

## From PyPI

Install as user:

```bash
pip install bspline --user
```

Install as admin:

```bash
sudo pip install bspline
```

## From GitHub

As user:

```bash
git clone https://github.com/johntfoster/bspline.git
cd bspline
python setup.py install --user
```

As admin, change the last command to

```bash
sudo python setup.py install
```

## Old method

Copy `bspline.py` and `splinelab.py` files from the `bspline` subdirectory next to your source code,
or leave them there and call it as a module.


# Tested on

 - Python 2.7 and 3.4.
 - Linux Mint.


# Dependencies

* [NumPy](http://www.numpy.org)
* [Matplotlib](http://matplotlib.org/) (for demo script)

# License

[MIT](LICENSE.md)
