/* Copyright (c) 2012,2013 Genome Research Ltd.
 *
 * Author: Stephan Schiffels <stephan.schiffels@sanger.ac.uk>
 *
 * This file is part of msmc.
 * msmc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
import std.algorithm;
import std.math;
import std.stdio;

class Amoeba {
  double ftol;
  int nfunc;
  int mpts;
  int ndim;
  double fmin;
  double[] y;
  double[][] p;

  this(double ftol=3.0e-8) {
    this.ftol = ftol;
  }
  
  double[] minimize(T)(in double[] point, double del, T func, bool verbose=false) {
    double[] dels = new double[point.length];
    dels[] = del;
    return minimize(point, dels, func, verbose);
  }      
    
  double[] minimize(T)(in double[] point, in double[] dels, T func, bool verbose=false) {
    int ndim = cast(int)point.length;
    double[][] pp = new double[][](ndim + 1, ndim);
    for (int i=0;i<ndim+1;i++) {
      for (int j=0;j<ndim;j++)
        pp[i][j]=point[j];
      if (i !=0 ) pp[i][i-1] += dels[i-1];
    }
    return minimize(pp,func, verbose);
  }

  double[] minimize(T)(double[][] pp, T func, bool verbose=false) {
    immutable int NMAX = 5000;
    immutable double TINY = 1.0e-10;
    int ihi, ilo, inhi;
    mpts = cast(int)pp.length;
    ndim = cast(int)pp[0].length;
    double[] psum = new double[ndim];
    double[] pmin = new double[ndim];
    double[] x = new double[ndim];
    this.p = pp.dup;
    y.length = mpts;
    for(int i = 0; i < mpts; ++i) {
      for(int j = 0; j < ndim; ++j)
        x[j] = p[i][j];
      y[i] = func(x);
    }
    nfunc = 0;
    get_psum(p, psum);
    for (;;) {
      if(nfunc % 100 == 0 && verbose)
          stderr.writefln("amoeba point #%s: %s, val=%s", nfunc,
              p[ilo], y[ilo]);
      ilo = 0;
      ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
      for(int i = 0; i < mpts; ++i) {
        if(y[i] <= y[ilo]) ilo = i;
        if(y[i] > y[ihi]) {
          inhi = ihi;
          ihi = i;
        }
        else if (y[i] > y[inhi] && i != ihi) inhi = i;
      }
      double rtol = 2.0 * abs(y[ihi] - y[ilo]) / (abs(y[ihi]) + abs(y[ilo]) +
          TINY);
      if(rtol < ftol || nfunc >= NMAX) {
        swap(y[0], y[ilo]);
        for(int i = 0;i < ndim; ++i) {
          swap(p[0][i], p[ilo][i]);
          pmin[i] = p[0][i];
        }
        fmin = y[0];
        if(nfunc >= NMAX) {
          // throw new Exception("NMAX exceeded");
          stderr.writeln("amoeba: WARNING, NMAX exceeded");
        }
        return pmin;
      }
      nfunc += 2;
      double ytry=amotry(p,y,psum,ihi,-1.0,func);
      if (ytry <= y[ilo])
        ytry = amotry(p,y,psum,ihi,2.0,func);
      else if (ytry >= y[inhi]) {
        double ysave=y[ihi];
        ytry = amotry(p,y,psum,ihi,0.5,func);
        if (ytry >= ysave) {
          for (int i=0;i<mpts;i++) {
            if (i != ilo) {
              for (int j=0;j<ndim;j++)
                p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
              y[i]=func(psum);
            }
          }
          nfunc += ndim;
          get_psum(p,psum);
        }
      } else --nfunc;
    }
  }
    
  void get_psum(in double[][] p, double[] psum) {
    for (int j=0;j<ndim;j++) {
      double sum=0.0;
      for (int i=0;i<mpts;i++)
        sum += p[i][j];
      psum[j]=sum;
    }
  }
    
  double amotry(T)(double[][] p, double[] y, double[] psum, int ihi,
      double fac, T func) {
    double[] ptry = new double[ndim];
    double fac1=(1.0-fac)/ndim;
    double fac2=fac1-fac;
    for (int j=0;j<ndim;j++)
      ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    double ytry=func(ptry);
    if (ytry < y[ihi]) {
      y[ihi]=ytry;
      for (int j=0;j<ndim;j++) {
        psum[j] += ptry[j]-p[ihi][j];
        p[ihi][j]=ptry[j];
      }
    }
    return ytry;
  }
}

unittest {
  writeln("test amoeba.minimize");
  double myfunc(double[] x) {
    return pow(x[0] - 3.0, 2.0) + pow(x[1] + 6.0, 2.0);
  }

  auto am = new Amoeba(1.0E-5);
  
  auto pmin = am.minimize([0, 0], 1, &myfunc);
  assert(approxEqual(pmin[0], 3.0) && approxEqual(pmin[1], -6.0));
}
