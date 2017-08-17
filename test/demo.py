# -*- coding: utf-8 -*-
"""Demonstration / usage example for the new features.

Created on Fri Mar 24 13:58:36 2017

@author: Juha Jeronen, juha.jeronen@tut.fi
"""

from __future__ import division, print_function, absolute_import

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors

import bspline
import bspline.splinelab as splinelab


def test():
    ########################
    # config
    ########################

    p      = 3  # order of spline basis (as-is! 3 = cubic)

    nknots = 5  # for testing: number of knots to generate (here endpoints count only once)

    tau = [0.1, 0.33]  # collocation sites (i.e. where to evaluate)

    ########################
    # usage example
    ########################

    knots = np.linspace(0,1,nknots)     # create a knot vector without endpoint repeats
    k     = splinelab.augknt(knots, p)  # add endpoint repeats as appropriate for spline order p
    B     = bspline.Bspline(k, p)       # create spline basis of order p on knots k

    # build some collocation matrices:
    #
    A0 = B.collmat(tau)                 # function value at sites tau
    A2 = B.collmat(tau, deriv_order=2)  # second derivative at sites tau

    print( A0 )
    print( A2 )

    ########################
    # tests
    ########################

    # number of basis functions
    n_interior_knots = len(knots) - 2
    n_basis_functions_expected = n_interior_knots + (p + 1)
    n_basis_functions_actual = len(B(0.))  # perform dummy evaluation to get number of basis functions
    assert n_basis_functions_actual == n_basis_functions_expected, "something went wrong, number of basis functions is incorrect"

    # partition-of-unity property of the spline basis
    assert np.allclose( np.sum(A0, axis=1), 1.0 ), "something went wrong, the basis functions do not form a partition of unity"


def main():
    """Demonstration: plot a B-spline basis and its first three derivatives."""

    #########################################################################################
    # config
    #########################################################################################

    # Order of spline basis.
    #
    p = 3

    # Knot vector, including the endpoints.
    #
    # For convenience, endpoints are specified only once, regardless of the value of p.
    #
    # Duplicate knots *in the interior* are allowed, with the standard meaning for B-splines.
    #
    knots = [0., 0.25, 0.5, 0.75, 1.]

    # How many plotting points to use on each subinterval [knots[i], knots[i+1]).
    #
    # Only intervals with length > 0 are actually plotted.
    #
    nt_per_interval = 101

    #########################################################################################
    # the demo itself
    #########################################################################################

    # The evaluation algorithm used in bspline.py uses half-open intervals  t_i <= x < t_{i+1}.
    #
    # This causes the right endpoint of each interval to actually be the start point of the next interval.
    #
    # Especially, the right endpoint of the last interval is the start point of the next (nonexistent) interval,
    # so the basis will return a value of zero there.
    #
    # We work around this by using a small epsilon to avoid evaluation exactly at t_{i+1} (for each interval).
    #
    epsrel = 1e-10
    epsabs = epsrel * (knots[-1] - knots[0])

    original_knots = knots
    knots = splinelab.augknt( knots, p )  # add repeated endpoint knots for splines of order p

    # treat each interval separately to preserve discontinuities
    #
    # (useful especially when plotting the highest-order nonzero derivative)
    #
    B   = bspline.Bspline(knots, p)
    xxs = []
    for I in zip( knots[:-1], knots[1:] ):
        t_i   = I[0]
        t_ip1 = I[1] - epsabs
        if t_ip1 - t_i > 0.:  # accept only intervals of length > 0 (to skip higher-multiplicity knots in the interior)
            xxs.append( np.linspace(t_i, t_ip1, nt_per_interval) )

    # common settings for all plotted lines
    settings = { "linestyle" : 'solid',
                 "linewidth" : 1.0 }

    # create a list of unique colors for plotting
    #
    # http://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib
    #
    NUM_COLORS = nbasis = len( B(0.) )  # perform dummy evaluation to get number of basis functions
    cm         = plt.get_cmap('gist_rainbow')
    cNorm      = matplotlib.colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap  = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    colors     = [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]

    labels = [ r"$B$",
               r"$\mathrm{d}B\,/\,\mathrm{d}x$",
               r"$\mathrm{d}^2B\,/\,\mathrm{d}x^2$",
               r"$\mathrm{d}^3B\,/\,\mathrm{d}x^3$" ]

    # for plotting the knot positions:
    unique_knots_xx = np.unique(original_knots)
    unique_knots_yy = np.zeros_like(unique_knots_xx)

    # plot the basis functions B(x) and their first three derivatives
    plt.figure(1)
    plt.clf()
    for k in range(4):
        ax = plt.subplot(2,2, k+1)

        # place the axis label where it fits
        if k % 2 == 0:
            ax.yaxis.set_label_position("left")
        else:
            ax.yaxis.set_label_position("right")

        # plot the kth derivative; each basis function gets a unique color
        f = B.diff(order=k)  # order=0 is a passthrough
        for xx in xxs:
            yy = np.array( [f(x) for x in xx] )  # f(scalar) -> rank-1 array, one element per basis function
            for i in range(nbasis):
                settings["color"] = colors[i]
                plt.plot( xx, yy[:,i], **settings )
            plt.ylabel( labels[k] )

        # show knot positions
        plt.plot( unique_knots_xx, unique_knots_yy, "kx" )

    plt.suptitle(r"$B$-spline basis functions, $p=%d$" % (p))


if __name__ == '__main__':
    test()
    main()
    plt.show()
