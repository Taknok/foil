#!/usr/bin/env python3

import math

IBX = 302  #< 300 = number of surface panel nodes + 6 
EPSILON = 1.e-6

# -----------------------------------------
#      calculates the arc length array s  |
#      for a 2-d array of points (x,y).   |
# -----------------------------------------
def scalc(x: list, y: list , s: list, n: int) -> list:
    s[1] = 0.0
    for i in range(2, n+1):
        s[i] = s[i-1] + math.sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]))

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

    if n > nmax:
        print("splind: array overflow, increase nmax")
        return False

    for i in range(2, n):
        dsm = s[i] - s[i-1]
        dsp = s[i+1] - s[i]
        b[i] = dsp
        a[i] = 2.0*(dsm+dsp)
        c[i] = dsm
        xs[i] = 3.0*((x[i+1]-x[i])*dsm/dsp + (x[i]-x[i-1])*dsp/dsm)

    if xs1 >= 998.0:
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

    if n == 2 and xs1 <= -998.0 and xs2 <= -998.0:
        b[n] = 1.0
        a[n] = 2.0
        xs[n] = 3.0*(x[n]-x[n-1]) / (s[n]-s[n-1])

    # solve for derivative array xs
    xs = trisol(a, b, c, xs, n)
    return xs




#-----------------------------------------------
#      Splines x(s) array just like spline,      |
#      but allows derivative discontinuities     |
#      at segment joints.  Segment joints are    |
#      defined by identical successive s values. |
#----------------------------------------------- */
def segspl(x: list, xs: list, s: list, n: int):
    nseg=0
    iseg=0
    iseg0=0

    if abs(s[1]-s[2]) < EPSILON:
        return false #stop 'segspl:  first input point duplicated'
    if abs(s[n]-s[n-1]) < EPSILON:
        return false #stop 'segspl:  last  input point duplicated'

    iseg0 = 1
    for iseg in range(2, n+1-2):
        if abs(s[iseg]-s[iseg+1]) < EPSILON:
            nseg = iseg - iseg0 + 1
            tmp = splind(x[iseg0 - 1:], xs[iseg0 -1:], s[iseg0 - 1:], nseg, -999, -999)
            xs = xs[:iseg0-1] + tmp
            iseg0 = iseg + 1

    nseg = n - iseg0 + 1
    tmp = splind(x[iseg0 - 1:], xs[iseg0 - 1:], s[iseg0 - 1:], nseg, -999, -999)
    xs = xs[:iseg0-1] + tmp
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
            tmp = splind(x[iseg0 - 1:], xs[iseg0 -1:], s[iseg0 - 1:], nseg, xs1, xs2)
            xs = xs[:iseg0-1] + tmp
            iseg0 = iseg + 1

    nseg = n - iseg0 + 1
    tmp = splind(x[iseg0 - 1:], xs[iseg0 - 1:], s[iseg0 - 1:], nseg, xs1, xs2)
    xs = xs[:iseg0-1] + tmp
    return xs


#--------------------------------------------------
 #       calculates dx/ds(ss)                         |
 #       xs array must have been calculated by spline |
#-------------------------------------------------- */
def deval(ss: int, x: list, xs: list, s: list, n: int):
    ilow = 1;
    #    i = nc;
    i = n; # techwinder modified


    while i-ilow > 1:
        # chelou ici
        imid = int((i+ilow)/2)
        if ss < s[imid]:
            i = imid
        else:
            ilow = imid

    ds = s[i] - s[i-1]
    t = (ss - s[i-1]) / ds
    cx1 = ds*xs[i-1] - x[i] + x[i-1]
    cx2 = ds*xs[i] - x[i] + x[i-1]
    deval = x[i] - x[i-1] + (1.0-4.0*t+3.0*t*t)*cx1 + t*(3.0*t-2.0)*cx2
    deval = deval/ds
    return deval

# --------------------------------------------------
 #      calculates d2x/ds2(ss)                       |
 #      xs array must have been calculated by spline |
 # --------------------------------------------------- */
def d2val(ss: float, x: list, xs: list, s: list, n: int):

    ilow = 1;
    i = n;
    while i-ilow > 1:
        # chelou ici aussi
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
    dtwoval = dtwoval/ds/ds;
    return dtwoval


#------------------------------------------------------
 #     locates leading edge spline-parameter value sle
 #
 #     the defining condition is
 #
 #      (x-xte,y-yte) . (x',y') = 0     at  s = sle
 #
 #     i.e. the surface tangent is normal to the chord
 #     line connecting x(sle),y(sle) and the te point.
#------------------------------------------------------ */
def lefind(sle: int, x: list, xp: list, y: list, yp: list, s: list, n: int):
    # convergence tolerance
    dseps = (s[n]-s[1]) * 0.00001

    # set trailing edge point coordinates
    xte = 0.5*(x[1] + x[n])
    yte = 0.5*(y[1] + y[n])

    # get first guess for sle
    for i in range(3, n - 1):
        dxte = x[i] - xte
        dyte = y[i] - yte
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        dotp = dxte*dx + dyte*dy
        if dotp < 0.0:
          break

    sle = s[i]

    # check for sharp le case
    if abs(s[i]-s[i-1]) < EPSILON:
      return False

    # newton iteration to get exact sle value
    for iter in range(1, 51):
        xle  = seval(sle,x,xp,s,n)
        yle  = seval(sle,y,yp,s,n)
        dxds = deval(sle,x,xp,s,n)
        dyds = deval(sle,y,yp,s,n)
        dxdd = d2val(sle,x,xp,s,n)
        dydd = d2val(sle,y,yp,s,n)

        xchord = xle - xte
        ychord = yle - yte

        # drive dot product between chord line and le tangent to zero
        res  = xchord*dxds + ychord*dyds
        ress = dxds*dxds + dyds*dyds + xchord*dxdd + ychord*dydd;

        # newton delta for sle
        dsle = -res/ress

        dsle = max(dsle, -0.02*abs(xchord+ychord))
        dsle = min(dsle,  0.02*abs(xchord+ychord))
        sle = sle + dsle
        if abs(dsle) < dseps:
          return sle

    sle = s[i]
    return sle

#      Calculates x(ss)
#       xs array must have been calculated by spline
def seval(ss: int, x: list, xs: list, s: list, n: int):
    ilow=0
    i=0
    imid=0
    ds=0
    t=0
    cx1=0
    cx2=0

    ilow = 1
    i = n

    while i-ilow > 1:
        imid = int((i+ilow)/2)
        if ss < s[imid]:
            i = imid;
        else:
            ilow = imid;

    ds = s[i] - s[i-1]
    t = (ss - s[i-1]) / ds
    if xs[i-1] == None:
        import pdb; pdb.set_trace()
    cx1 = ds*xs[i-1] - x[i] + x[i-1]
    cx2 = ds*xs[i]   - x[i] + x[i-1]
    return  t*x[i] + (1.0-t)*x[i-1] + (t-t*t)*((1.0-t)*cx1 - t*cx2)



#.....................................................................
#
#     interpolates two source airfoil shapes into an "intermediate" shape.
#
#     procedure:
#        The interpolated x coordinate at a given normalized spline
#        parameter value is a weighted average of the two source
#        x coordinates at the same normalized spline parameter value.
#        ditto for the y coordinates. The normalized spline parameter
#        runs from 0 at the leading edge to 1 at the trailing edge on
#        each surface.
#     .....................................................................
def inter(x0: list, xp0: list, y0: list, yp0: list, s0: list, n0: int, sle0: int, x1: list, xp1: list, y1: list, yp1: list, s1: list, n1: int, sle1: int, x: list, y: list, n: int, frac: int):
    f0=0
    f1=0
    tops0=0
    tops1=0
    bots0=0
    bots1=0
    sn=0
    st0=0
    st1=0
    # number of points in interpolated airfoil is the same as in airfoil 0
    n = n0

    # interpolation weighting fractions
    f0 = 1.0 - frac
    f1 = frac

    # top side spline parameter increments
    tops0 = s0[1] - sle0
    tops1 = s1[1] - sle1

    # bottom side spline parameter increments
    bots0 = s0[n0] - sle0
    bots1 = s1[n1] - sle1

    for i in range(1, n+1):

        # normalized spline parameter is taken from airfoil 0 value
        if s0[i]< sle0:
            sn = (s0[i] - sle0) / tops0  # top side
        else:
            sn = (s0[i] - sle0) / bots0  # bottom side

        # set actual spline parameters
        st0 = s0[i]
        if st0 < sle0:
            st1 = sle1 + tops1 * sn
        #        if(st0>=sle0) st1 = sle1 + bots1 * sn;
        else:
            st1 = sle1 + bots1 * sn

        # set interpolated x,y coordinates
        x[i] = f0*seval(st0,x0,xp0,s0,n0) + f1*seval(st1,x1,xp1,s1,n1)
        y[i] = f0*seval(st0,y0,yp0,s0,n0) + f1*seval(st1,y1,yp1,s1,n1)
        print(st0)
        
    return x, y, n

#------------------------------------------------------
#      sets geometric parameters for airfoil shape
 # ------------------------------------------------------ */
#def geopar(x: list, xp: list, y: list, yp: list, s: list, n: int, t: list, sle: int, chord: )



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
    xb = [None for i in range(IBX)]
    yb = [None for i in range(IBX)]
    sb = [None for i in range(IBX)]
    xbp = [None for i in range(IBX)]
    ybp = [None for i in range(IBX)]
    nb = IBX
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
    sleint1 = lefind(sleint1, x1, xp1, y1, yp1, s1, n1)

    s2 = scalc(x2,y2,s2,n2)
    xp2 = segspld(x2,xp2,s2,n2, -999.0, -999.0)
    yp2 = segspld(y2,yp2,s2,n2, -999.0, -999.0)
    sleint2 = lefind(sleint2, x2, xp2, y2, yp2, s2, n2)

    # where is defined xb ?
    xb, yb, nb = inter(x1, xp1, y1, yp1, s1, n1, sleint1, x2, xp2, y2, yp2, s2, n2, sleint2, xb, yb, nb, mixt);

    #sb = scalc(xb,yb,sb,nb)
    #xbp = segspl(xb,xbp,sb,nb)
    #ybp = segspl(yb,ybp,sb,nb)

    return xb, yb, nb
    # geopar(xb,xbp,yb,ybp,sb,nb,w1,sble,chordb,areab,radble,angbte, ei11ba,ei22ba,apx1ba,apx2ba,ei11bt,ei22bt,apx1bt,apx2bt)


if __name__ == "__main__":
    import sys
    mixt = float(sys.argv[1]) if len(sys.argv) > 1 else 20
    mixt = mixt / 100.0
    
    with open("foil1.dat") as f1, open("foil2.dat") as f2:
          lines = f1.readlines()
          name0 = lines[0].strip()
          x1 = []
          y1 = []
          for i in lines[1:]:
              x, y = i.split()
              x1.append(float(x))
              y1.append(float(y))
          n1 = len(x1)
          
          lines = f2.readlines()
          name1 = lines[0].strip()
          x2 = []
          y2 = []
          for i in lines[1:]:
              x, y = i.split()
              x2.append(float(x))
              y2.append(float(y))
          n2 = len(x2)
          
          xb, yb, nb = interpolate(x1, y1, n1, x2, y2, n2, mixt)
          
          for i in range(nb):
              print(f"{xb[i+1]:.5f} {yb[i+1]:.5f} ", end="")
          name = f"{name0}-{100-mixt* 100}_{name1}-{mixt*100}"
          with open(f"{name}.dat", "w") as f:
            f.write(f"{name}\n")
            for i in range(nb):
                f.write(f" {xb[i+1]:.5f}    {yb[i+1]:.5f}\n")

            f.write("\n")