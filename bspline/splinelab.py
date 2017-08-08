# -*- coding: utf-8 -*-
"""Utility functions for B-splines, potentially useful for converting MATLAB codes to Python.

CAVEATS:
    - Only a very minimal set of functionality is implemented here.
    - Some technical details differ from the MATLAB equivalents.

Particularly, we use spline order `p` **as-is** instead of MATLAB's `k` parameter (`k` = `p` + 1)
in the function parameters.

Created on Fri Mar 24 13:52:37 2017

@author: Juha Jeronen, juha.jeronen@tut.fi
"""

from __future__ import division, print_function, absolute_import

import numpy as np

import bspline.bspline


def augknt(knots, order):
    """Augment a knot vector.

Parameters:
    knots:
        Python list or rank-1 array, the original knot vector (without endpoint repeats)
    order:
        int, >= 0, order of spline

Returns:
    list_of_knots:
        rank-1 array that has (`order` + 1) copies of ``knots[0]``, then ``knots[1:-1]``, and finally (`order` + 1) copies of ``knots[-1]``.

Caveats:
    `order` is the spline order `p`, not `p` + 1, and existing knots are never deleted.
    The knot vector always becomes longer by calling this function.
"""
    if isinstance(knots, np.ndarray)  and  knots.ndim > 1:
        raise ValueError("knots must be a list or a rank-1 array")
    knots = list(knots)  # ensure Python list

    # One copy of knots[0] and knots[-1] will come from "knots" itself,
    # so we only need to prepend/append "order" copies.
    #
    return np.array( [knots[0]] * order  +  knots  +  [knots[-1]] * order )


def aveknt(t, k):
    """Compute the running average of `k` successive elements of `t`. Return the averaged array.

Parameters:
    t:
        Python list or rank-1 array
    k:
        int, >= 2, how many successive elements to average

Returns:
    rank-1 array, averaged data. If k > len(t), returns a zero-length array.

Caveat:
    This is slightly different from MATLAB's aveknt, which returns the running average
    of `k`-1 successive elements of ``t[1:-1]`` (and the empty vector if  ``len(t) - 2 < k - 1``).

"""
    t = np.atleast_1d(t)
    if t.ndim > 1:
        raise ValueError("t must be a list or a rank-1 array")

    n = t.shape[0]
    u = max(0, n - (k-1))  # number of elements in the output array
    out = np.empty( (u,), dtype=t.dtype )

    for j in range(u):
        out[j] = sum( t[j:(j+k)] ) / k

    return out


def aptknt(tau, order):
    """Create an acceptable knot vector.

Minimal emulation of MATLAB's ``aptknt``.

The returned knot vector can be used to generate splines of desired `order`
that are suitable for interpolation to the collocation sites `tau`.

Note that this is only possible when ``len(tau)`` >= `order` + 1.

When this condition does not hold, a valid knot vector is returned,
but using it to generate a spline basis will not have the desired effect
(the spline will return a length-zero array upon evaluation).

Parameters:
    tau:
        Python list or rank-1 array, collocation sites

    order:
        int, >= 0, order of spline

Returns:
    rank-1 array, `k` copies of ``tau[0]``, then ``aveknt(tau[1:-1], k-1)``,
    and finally `k` copies of ``tau[-1]``, where ``k = min(order+1, len(tau))``.
"""
    tau = np.atleast_1d(tau)
    k   = order + 1

    if tau.ndim > 1:
        raise ValueError("tau must be a list or a rank-1 array")

    # emulate MATLAB behavior for the "k" parameter
    #
    # See
    #   https://se.mathworks.com/help/curvefit/aptknt.html
    #
    if len(tau) < k:
        k = len(tau)

    if not (tau == sorted(tau)).all():
        raise ValueError("tau must be nondecreasing")

    # last processed element needs to be:
    #     i + k - 1 = len(tau)- 1
    # =>  i + k = len(tau)
    # =>  i = len(tau) - k
    #
    u = len(tau) - k
    for i in range(u):
        if tau[i+k-1] == tau[i]:
            raise ValueError("k-fold (or higher) repeated sites not allowed, but tau[i+k-1] == tau[i] for i = %d, k = %d" % (i,k))

    # form the output sequence
    #
    prefix = [ tau[0]  ] * k
    suffix = [ tau[-1] ] * k

    # https://se.mathworks.com/help/curvefit/aveknt.html
    # MATLAB's aveknt():
    #  - averages successive k-1 entries, but ours averages k
    #  - seems to ignore the endpoints
    #
    tmp    = aveknt(tau[1:-1], k-1)
    middle = tmp.tolist()
    return np.array( prefix + middle + suffix, dtype=tmp.dtype )


def knt2mlt(t):
    """Count multiplicities of elements in a sorted list or rank-1 array.

Minimal emulation of MATLAB's ``knt2mlt``.

Parameters:
    t:
        Python list or rank-1 array. Must be sorted!

Returns:
    out
        rank-1 array such that
        out[k] = #{ t[i] == t[k] for i < k }

Example:
    If ``t = [1, 1, 2, 3, 3, 3]``, then ``out = [0, 1, 0, 0, 1, 2]``.

Caveat:
    Requires input to be already sorted (this is not checked).
"""
    t = np.atleast_1d(t)
    if t.ndim > 1:
        raise ValueError("t must be a list or a rank-1 array")

    out   = []
    e     = None
    for k in range(t.shape[0]):
        if t[k] != e:
            e     = t[k]
            count = 0
        else:
            count += 1
        out.append(count)

    return np.array( out )


def spcol(knots, order, tau):
    """Return collocation matrix.

Minimal emulation of MATLAB's ``spcol``.

Parameters:
    knots:
        rank-1 array, knot vector (with appropriately repeated endpoints; see `augknt`, `aptknt`)
    order:
        int, >= 0, order of spline
    tau:
        rank-1 array, collocation sites

Returns:
    rank-2 array A such that

        A[i,j] = D**{m(i)} B_j(tau[i])

    where
        m(i) = multiplicity of site tau[i]

        D**k  = kth derivative (0 for function value itself)
"""
    m = knt2mlt(tau)
    B = bspline.Bspline(knots, order)

    dummy = B(0.)
    nbasis = len(dummy)  # perform dummy evaluation to get number of basis functions

    A = np.empty( (tau.shape[0], nbasis), dtype=dummy.dtype )
    for i,item in enumerate(zip(tau,m)):
        taui,mi = item
        f       = B.diff(order=mi)
        A[i,:]  = f(taui)

    return A
