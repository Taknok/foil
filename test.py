#!/usr/bin/env python3

from math import pow, sqrt

class Airfoil:

    EPSILON = 1.e-6

    __slots__ = ("name", "coordinates")

    def __init__(self, name, coordinates):
        self.name = name
        self.coordinates = coordinates

    # accessor
    @property
    def x(self):
        return [c[0] for c in self.coordinates]

    def x(self, index):
        return self.coordinates[index][0]

    @property
    def y(self):
        return [c[1] for c in self.coordinates]

    def y(self, index):
        return self.coordinates[index][1]

    @classmethod
    def from_coordinates(cls, name, xs, ys):
        pass

    @classmethod
    def from_file(cls, filename):
        coordinates = []
        with open(filename) as f:
            name = f.readline()
            for l in f.readlines():
                coordinates.append(l.split())

        # convert str to float
        coordinates = [[float(x[0]), float(x[1])] for x in coordinates]

        return cls(name, coordinates)


    #
    # Calculates the arc length (named 's') for a coordinates
    #
    def scalc(self):
        s = [0.0]
        for i in range(len(self.coordinates) - 1):
            s.append(
                s[i] + sqrt(
                    pow(self.x(i+1) - self.x(i), 2)
                    + pow(self.y(i+1) - self.y(i), 2)
                )
            )
        return s

    def trisol(a, b, c, d, kk):
        # solve kk long, tri-diagonal system
        # a c        d
        # b a c      d
        #   b a .    .
        #     . . c  .
        #       b a  d
        # the righthand side d is replaced by the
        # solution. a, c are destroyed.
        for k in range(2, kk+1):
            km = k-1
            c[km] = c[km] / a[km]
            d[km] = d[km] / a[km]
            a[k] = a[k] - b[k]*c[km]
            d[k] = d[k] - b[k]*d[km]

        d[kk] = d[kk] / a[kk]
        for k in range(2, kk+1):
            d[k] = d[k] - c[k]*d[k+1]

    # Calculates spline coefficients for x(s).
    # Specified 1st derivative and/or usual zero 2nd
    # derivative end conditions are used.
    #
    # To evaluate the spline at some value of s, use seval
    # and/or deval.
    #
    # s         independent variable array (input)
    # x         dependent variable array (input)
    # xs        xs/ds array (calculated)
    # n         number of points (input)
    # xs1, xs2  endpoint derivative (input)
    #           if = 999.0, then usual zero second derivative
    #           end condition(s) are used
    #           if = -999.0, then zero third derivative end
    #           condition(s) are used
    def splind(self, x, xs, s, n, xs1, xs2):
        nmax = 600
        if len(x) > nmax:
            raise "splind: array overflow, increase max"

        dsm = 0
        dsp = 0
        for i in range(2, len(x)):
            dsm = s[i] - s[i-1]
            dsp = s[i+1] - s[i]
            b[i] = dsp
            a[i] = 2.0 * (dsm+dsp)
            c[i] = dsm
            xs[i] = 3.0*((x[i+1]-x[i])*dsm/dsp + (x[i]-x[i-1])*dsp/dsm)

        if xs1 >= 998.0:
            # set zero second derivative end condition
            a[1] = 2.0
            c[1] = 1.0
            xs[1] = 3.0*(x[2]-x[1]) / (s[2]-s[1])
        else:
            if xs1 <= -998.0:
                a[1] = 1.0
                c[1] = 1.0
                xs[1] = 2.0*(x[2]-x[1]) / (s[2]-s[1])
            else:
                a[1] = 1.0
                c[1] = 0.0
                xs[1] = xs1

        if xs2 >= 998.0:
            b[n] = 1.0
            a[n] = 2.0
            xs[n] = 3.0*(x[n]-x[n-1]) / (s[n]-s[n-1])
        else:
            if xs2 <= -998.0:
                b[n] = 1.0
                b[n] = 0.0
                xs[n] = xs2

        if n == 2 and xs1 <= -998.0 and xs2 <= -998.0:
            b[n] = 1.0
            a[n] = 2.0
            xs[n] = 3.0*(x[n]-x[n-1]) / (s[n]-s[n-1])

        # solve for derivative array xs
        trisol(a, b, c, xs, n)


    def segspld(self, s=self.scalc()):
        if abs(s[0]-s[2]) < EPSILON:
            return False, [] #stop 'segspl:  first input point duplicated';

        if abs(s[-1]-s[-2]) < EPSILON:
            return False, [] #stop 'segspl:  last  input point duplicated';

        for i in range(2, len(s)-2):
            if abs(s[i]-s[i+1]) < EPSILON:
                splind()

    # calculates x(ss)
    # xs array mus have been calculated by spline
    def seval(ss, x, xs, s, n):
        ilow = 1
        i = n
        while i - ilow > 1:
            imid = int((i+ilow)/2)
            if ss < s[imid]:
                i = imid
            else:
                ilow = imid

        ds = s[i] -s[i-1]
        t = (ss - s[i-1]) / ds
        cx1 = ds*xs[i-1] - x[i] + x[i-1]
        cx2 = ds*xs[i] - x[i] + x[i-1]

        return t*x[i] + (1.0-t)*x[i-1] + (t-t*t)*((1.0-t)*cx1 - t*cx2)

    # calculates dx/ds(ss)
    # xs array must have been calculated by spline
    def deval(ss, x, xs, s, n):
        ilow = 1
        i = n #techwinder modified
        while i - ilow > 1:
            imid = int((i+ilow)/2)
            if ss < s[imid]:
                i = imid
            else:
                ilow = imid

        ds = s[i] - s[i-1]
        t = (ss - s[i-1]) / ds
        cx1 = ds*xs[i-1] - x[i] + x[i-1]
        cx2 = ds*xs[i]     - x[i] + x[i-1]
        deval = x[i] - x[i-1] + (1.0-4.0*t+3.0*t*t)*cx1
            + t*(3.0*t-2.0)*cx2
        deval = deval/ds
        return deval

    # calculates d2x/ds2(ss)
    # xs array must have been calculated by spline
    def d2val(ss, x, xs, s, n):
        ilow = 1
        i = n #techwinder modified
        while i - ilow > 1:
            imid = int((i+ilow)/2)
            if ss < s[imid]:
                i = imid
            else:
                ilow = imid

        ds = s[i] - s[i-1]
        t = (ss - s[i-1]) / ds
        cx1 = ds*xs[i-1] - x[i] + x[i-1]
        cx2 = ds*xs[i]   - x[i] + x[i-1]
        dtwoval = (6.0*t-4.0)*cx1 + (6.0*t-2.0)*cx2
        dtwoval = dtwoval/ds/ds
        return dtwoval


    # Locates leading edge spline-parameter value sle
    # the defining condition is:
    # (x-xte,y-yte) . (x',y') = 0 at s = sle
    # i.E the surface tangent is normal to the chord
    # line connecting x(sle), y(sle) and the te point.
    def lefind(self, sle, x, xp, y, yp, s, n):
        # convergence tolerance
        dseps = (s[n] - s[1]) * 0.00001

        # set trailing edge point coordinates
        xte = 0.5*(x[1] + x[n])
        yte = 0.5*(y[1] + y[n])

        # get first guess for sle
        for i in range(3, n-1):
            dxte = x[i] - xte
            dyte = y[i] - yte
            dx = x[i+1] - x[i]
            dy = y[i+1] - y[i]
            dotp = dxte*dx + dyte*dy
            if dtop < 0.0:
                break

        sle = s[i]

        # check for sharp le case
        if abs(s[i]-s[i-1]) < EPSILON
            return False

        # newton iteration to get exact sle value
        for i in range(50):
            xle = seval(sle, x, xp, s, n)
            yle = seval(sle, y, yp, s, n)
            dxds = deval(sle, x, xp, s, n)
            dyds = deval(sle, y, yp, s, n)
            dxdd = d2val(sle, x, xp, s, n)
            dydd = d2val(sle, y, yp, s, n)

            xchord = xle - xte
            ychord = yle - yte

            # drive not product between chord line and le 
            # tangent to zero
            res = xchord*dxdx + ychord*dyds
            ress = dxds*dxds + dyds*dyds + xchord*dxdd + ychord*dydd

            # newton delta for sle
            dsle = -res/ress

            dsle = max(dsle, -0.02*abs(xchord+ychord))
            dsle = min(dsle,  0.02*abs(xchord+ychord))
            sle = sle + dsle
            if abs(dsle) < dseps
                return True

        sle = s[i]
        return True

def interpolate(foil1, foil2, mixt):
    s1 = foil1.scalc()
    segspld()
    segspld()
    lefind()

    s2 = foil2.scalc()

if __name__ == '__main__':
    a = Airfoil.from_file("foil1.dat")
    b = Airfoil.from_file("foil2.dat")

    interpolate(a, b, 0.5)
