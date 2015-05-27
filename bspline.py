from functools import partial
import numpy as np
import matplotlib.pyplot as plt
np.seterr(divide='ignore', invalid='ignore')

class memoize(object):
    """cache the return value of a method
    
    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments passed to a method decorated with memoize must
    be hashable.
    
    If a memoized method is invoked directly on its class the result will not
    be cached. Instead the method will be invoked like a static method:
    class Obj(object):
        @memoize
        def add_to(self, arg):
            return self + arg
    Obj.add_to(1) # not enough arguments
    Obj.add_to(1, 2) # returns 3, result is not cached
    
    Script borrowed from here:
    MIT Licensed, attributed to Daniel Miller, Wed, 3 Nov 2010
    http://code.activestate.com/recipes/577452-a-memoize-decorator-for-instance-methods/
    """
    def __init__(self, func):
        self.func = func
    def __get__(self, obj, objtype=None):
        if obj is None:
            return self.func
        return partial(self, obj)
    def __call__(self, *args, **kw):
        obj = args[0]
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        key = (self.func, args[1:], frozenset(kw.items()))
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kw)
        return res



class Bspline():
    
    def __init__(self, knot_vector, order):
        
        self.knot_vector = np.array(knot_vector)
        self.p = order
        
        
    def __basis0(self, xi):
        
        return np.where(np.all([self.knot_vector[:-1] <=  xi, 
                                xi < self.knot_vector[1:]],axis=0), 1.0, 0.0)
    
    def __basis(self, xi, p, compute_derivatives=False):
        
        if p == 0:
            return self.__basis0(xi)
        else:
            basis_p_minus_1 = self.__basis(xi, p - 1)
        
        first_term_numerator = xi - self.knot_vector[:-p] 
        first_term_denominator = self.knot_vector[p:] - self.knot_vector[:-p]
        
        second_term_numerator = self.knot_vector[(p + 1):] - xi
        second_term_denominator = (self.knot_vector[(p + 1):] - 
                                   self.knot_vector[1:-p])
                
        
        first_term = np.where(first_term_denominator != 0.0, 
                              (first_term_numerator / 
                               first_term_denominator), 0.0)
        second_term = np.where(second_term_denominator != 0.0,
                               (second_term_numerator / 
                                second_term_denominator), 0.0)
        
        if compute_derivatives and p == self.p:
            
            first_term = np.where(first_term_denominator != 0.0, 
                                  p / first_term_denominator, 0.0)
            second_term = np.where(second_term_denominator != 0.0,
                                   -p / second_term_denominator, 0.0)
        
        return  (first_term[:-1] * basis_p_minus_1[:-1] + 
                 second_term * basis_p_minus_1[1:])
            
    
    @memoize
    def __call__(self, xi):
        
        return self.__basis(xi, self.p, compute_derivatives=False)
    
    @memoize
    def d(self, xi):
        
        return self.__basis(xi, self.p, compute_derivatives=True)
    
    def plot(self):
        
        x_min = np.min(self.knot_vector)
        x_max = np.max(self.knot_vector)
        
        x = np.linspace(x_min, x_max, num=1000)
        
        N = np.array([self(i) for i in x]).T;
        
        for n in N:
            
            plt.plot(x,n)
            
        return plt.show()
    
    def dplot(self):
        
        x_min = np.min(self.knot_vector)
        x_max = np.max(self.knot_vector)
        
        x = np.linspace(x_min, x_max, num=1000)
        
        N = np.array([self.d(i) for i in x]).T;
        
        for n in N:
            
            plt.plot(x,n)
            
        return plt.show()
