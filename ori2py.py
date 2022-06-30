#!/usr/bin/env python3

IBX = 302  #< 300 = number of surface panel nodes + 6 
EPSILON = 1.e-6

# -----------------------------------------
#      calculates the arc length array s  |
#      for a 2-d array of points (x,y).   |
# -----------------------------------------
def scalc(x: list, y: list , s: list, n: int) -> list:
    s[1] = 0.0
    for i in range(2, n+1):
        s[i] = s[i-1] + sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]))

    return s

#-----------------------------------------
#     solves kk long, tri-diagonal system |
#                                         |
#             a c          d              |
#             b a c        d              |
#               b a .      .              |
#                 . . c    .              |
#                   b a    d              |
#                                         |
#     the righthand side d is replaced by |
#     the solution.  a, c are destroyed.  |
#-----------------------------------------
def trisol(a: list, b: list, c: list, d: list, kk: int):
    for k in range(2, kk+1):
        km = k - 1
        c[km] = c[km] / a[km]
        d[km] = d[km] / a[km]
        a[k] = a[k] - b[k]*c[km]
        d[k] = d[k] - b[k]*d[km]

    d[kk] = d[kk]/a[kk]

    for k in range(kk-1, 0):
        d[k] = d[k] - c[k]*d[k+1]

    return d


# -------------------------------------------------------
#      Calculates spline coefficients for x(s).          |
#      Specified 1st derivative and/or usual zero 2nd    |
#      derivative end conditions are used.               |
#                                                        |
#      To evaluate the spline at some value of s,        |
#      use seval and/or deval.                           |
#                                                        |
#      s        independent variable array (input)       |
#      x        dependent variable array   (input)       |
#      xs       dx/ds array                (calculated)  |
#      n        number of points           (input)       |
#      xs1,xs2  endpoint derivatives       (input)       |
#               if = 999.0, then usual zero second       |
#               derivative end condition(s) are used     |
#               if = -999.0, then zero third             |
#               derivative end condition(s) are used     |
#                                                        |
# -------------------------------------------------------
def splind(x: list, xs: list, s: list, n: int, xs1: int, xs2: int):
    nmax = 600
    a = [None for i in range(601)]
    b = [None for i in range(601)]
    c = [None for i in range(601)]
    dsm = 0
    dsp = 0

    if n > nmax
        print("splind: array overflow, increase nmax")
        return False

    for i in range(2, n):
        dsm = s[i] - s[i-1]
        dsp = s[i+1] - s[i]
        b[i] = dsp
        a[i] = 2.0*(dsm+dsp)
        c[i] = dsm
        xs[i] = 3.0*((x[i+1]-x[i])*dsm/dsp + (x[i]-x[i-1])*dsp/dsm)

    if xs1 => 998.0:
        # set zero second derivative end condition
        a[1] = 2.0
        c[1] = 1.0
        xs[1] = 3.0*(x[2]-x[1]) / (s[2]-s[1])
    else:
        if xs1 <= -998.0:
            # set third derivative end condition
            a[1] = 1.0
            c[1] = 1.0
            xs[1] = 2.0*(x[2]-x[1]) / (s[2]-s[1])
        else:
            # set specified first derivative end condition
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
            a[n] = 1.0
            xs[n] = 2.0*(x[n]-x[n-1]) / (s[n]-s[n-1])
        else:
            a[n] = 1.0
            b[n] = 0.0
            xs[n] = xs2

    if n == 2 && xs1 <= -998.0 && xs2 <= -998.0:
        b[n] = 1.0
        a[n] = 2.0
        xs[n] = 3.0*(x[n]-x[n-1]) / (s[n]-s[n-1])

    # solve for derivative array xs
    xs = trisol(a, b, c, xs, n)
    return xs


# -----------------------------------------------
#     splines x(s) array just like splind,      |
#     but allows derivative discontinuities     |
#     at segment joints.  segment joints are    |
#     defined by identical successive s values. |
#------------------------------------------------
def segspld(x: list, xs: list, s: int, n: int, xs1: int, xs2: int):
    nseg = 0
    iseg = 0
    iseg0 = 0

    if abs(s[1]-s[2]) < EPSILON:
        return False
    if abs(s[n]-s[n-1]) < EPSILON:
        return False

    iseg0 = 1
    for iseg in range(2, n+1-2):
        if abs(s[iseg]-s[iseg+1]) < EPSILON:
            nseg = iseg - iseg0 + 1
            # CONTINUE: 1
            tmp = splind(x[iseg0 - 1:], xs[iseg0 -1:], s[iseg0 - 1:], nseg, xs1, xs2)
            xs = xs[..iseg0 - 1] + tmp
            iseg0 = iseg + 1

    nseg = n - iseg0 + 1
    tmp = splind(x + iseg0 - 1, xs + iseg0 - 1, s + iseg0 - 1, nseg, xs1, xs2)
    xs = xs[iseg0 - 1] + tmp
    return xs

def interpolate(xf1, yf1, n1,
                xf2, yf2, n2,
                mixt):
    i = 0
    x1 = [None for i in range(IBX)]
    y1 = [None for i in range(IBX)]
    x2 = [None for i in range(IBX)]
    y2 = [None for i in range(IBX)]
    xp1 = [None for i in range(IBX)]
    yp1 = [None for i in range(IBX)]
    xp2 = [None for i in range(IBX)]
    yp2 = [None for i in range(IBX)]
    s1 = [None for i in range(IBX)]
    s2 = [None for i in range(IBX)]
    sleint1=0
    sleint2=0
    for i in range(n1):
        x1[i+1] = xf1[i]
        y1[i+1] = yf1[i]

    for i in range(n2):
        x2[i+1] = xf2[i]
        y2[i+1] = yf2[i]

    s1 = scalc(x1,y1,s1,n1)
    xp1 = segspld(x1,xp1,s1,n1, -999.0, -999.0)
    yp1 = segspld(y1,yp1,s1,n1, -999.0, -999.0)
    lefind(sleint1, x1, xp1, y1, yp1, s1, n1)

    scalc(x2,y2,s2,n2)
    segspld(x2,xp2,s2,n2, -999.0, -999.0)
    segspld(y2,yp2,s2,n2, -999.0, -999.0)
    lefind(sleint2, x2, xp2, y2, yp2, s2, n2)

    inter(x1, xp1, y1, yp1, s1, n1, sleint1,
          x2, xp2, y2, yp2, s2, n2, sleint2,
          xb,yb,nb,mixt);

    scalc(xb,yb,sb,nb)
    segspl(xb,xbp,sb,nb)
    segspl(yb,ybp,sb,nb)

    geopar(xb,xbp,yb,ybp,sb,nb,w1,sble,chordb,areab,radble,angbte,
           ei11ba,ei22ba,apx1ba,apx2ba,ei11bt,ei22bt,apx1bt,apx2bt)
