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
 
import std.math;
import std.stdio;
import std.string;
import brent;
import logger;

class F1dim(T) {
  double[] p, xi;
  int n;
  T func;
  double[] xt;

  this(double[] p, double[] xi, T func) {
    this.p = p;
    this.xi = xi;
    n = cast(int)p.length;
    this.func = func;
    xt.length = n;
  }
  
  double opCall(double x) {
    for(int j=0;j<n;j++)
      xt[j]=p[j]+x*xi[j];
    return func(xt);
  }
};

class Linemethod(T) {
  double[] p, xi;
  T func;
  int n;

  this(T func) {
    this.func = func;
  }
  
  double linmin() {
    double ax,xx,xmin;
    n=cast(int)p.length;
    auto f1dim = new F1dim!T(p,xi,func);
    ax=0.0;
    xx=1.0;
    auto brent = new Brent();
    brent.bracket(ax,xx,f1dim);
    xmin=brent.minimize(f1dim);
    for(int j=0;j<n;j++) {
      xi[j] *= xmin;
      p[j] += xi[j];
    }
    return brent.fmin;
  }
};

class Powell(T) : Linemethod!T {
  int iter;
  double fret;
  double ftol;
  this(T func, double ftol=3.0e-8) {
    super(func);
    this.ftol = ftol;
  }

  double[] minimize(double[] pp) {
    int n=cast(int)pp.length;
    auto ximat = new double[][](n, n);
    foreach(ref row; ximat) row[] = 0.0;
    for(int i=0;i<n;i++) ximat[i][i]=1.0;
    return minimize(pp,ximat);
  }
  
  double[] minimize(double[] pp, double[][] ximat) {
    immutable int ITMAX=200;
    immutable double TINY=1.0e-25;
    double fptt;
    int n=cast(int)pp.length;
    p=pp.dup;
    auto pt = new double[n];
    auto ptt = new double[n];
    xi.length = n;
    fret=func(p);
    for(int j=0;j<n;j++) pt[j]=p[j];
    for(iter=0;;++iter) {
      double fp=fret;
      logInfo(format("\r  * [%s/200(max)] Maximization Step", iter));
      int ibig=0;
      double del=0.0;
      for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) xi[j]=ximat[j][i];
        fptt=fret;
        fret=linmin();
        if(fptt-fret > del) {
          del=fptt-fret;
          ibig=i+1;
        }
      }
      if(2.0*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) {
        return p;
      }
      if(iter == ITMAX)
        throw new Exception("powell exceeding maximum iterations.");
      for(int j=0;j<n;j++) {
        ptt[j]=2.0*p[j]-pt[j];
        xi[j]=p[j]-pt[j];
        pt[j]=p[j];
      }
      fptt=func(ptt);
      if(fptt < fp) {
        double t=2.0*(fp-2.0*fret+fptt)*(fp-fret-del)^^2-del*(fp-fptt)^^2;
        if(t < 0.0) {
          fret=linmin();
          for(int j=0;j<n;j++) {
            ximat[j][ibig-1]=ximat[j][n-1];
            ximat[j][n-1]=xi[j];
          }
        }
      }
    }
  }
};

// class Powell(T) : Linemethod!T {
//   int iter;
//   double fret;
//   double ftol;
//   
//   // these are needed during the iteration
//   double fptt;
//   double[] pt, ptt;
//   double[][] ximat;
//   double fp;
//   
//   immutable int ITMAX=200;
//   immutable double TINY=1.0e-25;
//   
//   this(T func, double ftol=3.0e-8) {
//     super(func);
//     this.ftol = ftol;
//   }
// 
//   double[] step() {
//     fp=fret;
//     int ibig=0;
//     double del=0.0;
//     for (int i=0;i<n;i++) {
//       for (int j=0;j<n;j++) xi[j]=ximat[j][i];
//       fptt=fret;
//       fret=linmin();
//       if (fptt-fret > del) {
//         del=fptt-fret;
//         ibig=i+1;
//       }
//     }
//     for (int j=0;j<n;j++) {
//       ptt[j]=2.0*p[j]-pt[j];
//       xi[j]=p[j]-pt[j];
//       pt[j]=p[j];
//     }
//     fptt=func(ptt);
//     if (fptt < fp) {
//       double t=2.0*(fp-2.0*fret+fptt)*(fp-fret-del)^^2-del*(fp-fptt)^^2;
//       if (t < 0.0) {
//         fret=linmin();
//         for (int j=0;j<n;j++) {
//           ximat[j][ibig-1]=ximat[j][n-1];
//           ximat[j][n-1]=xi[j];
//         }
//       }
//     }
//     ++iter;
//     return p;
//   }
//   
//   bool finished() {
//     return 2.0*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY;
//   }
//   
//   void init(double[] pp, double[][] ximat) {
//     //copy ximat into this.ximat
//     this.ximat.length = 0;
//     foreach(row; ximat) {
//       this.ximat ~= row.dup;
//     }
// 
//     n=cast(int)pp.length;
//     p=pp.dup;
//     pt = new double[n];
//     ptt = new double[n];
//     xi.length = n;
//     fret=func(p);
//     for (int j=0;j<n;j++) pt[j]=p[j];
//   }
//   
//   void init(double[] pp) {
//     n=cast(int)pp.length;
//     auto ximat = new double[][](n, n);
//     foreach(ref row; ximat) row[] = 0.0;
//     for(int i=0;i<n;i++) ximat[i][i]=1.0;
//     init(pp, ximat);
//   }
//   
//   double[] minimize() {
//     for (;;) {
//       step();
//       if (finished()) {
//         return p;
//       }
//       if(iter == ITMAX) {
//         stderr.writeln("powell exceeding maximum iterations.");
//         return p;
//       }
//     }
//   }
//   
//   double[] minimize(double[] pp) {
//     init(pp);
//     return minimize();
//   }
//   
//   double[] minimize(double[] pp, double[][] ximat) {
//     init(pp, ximat);
//     return minimize();
//   }
// };

unittest {
  double myfunc(double[] x) {
    return pow(x[0] - 3.0, 2.0) + pow(x[1] + 6.0, 2.0);
  }

  auto pw = new Powell!(typeof(&myfunc))(&myfunc);
  
  auto pmin = pw.minimize([0, 0]);
  assert(approxEqual(pmin[0], 3.0) && approxEqual(pmin[1], -6.0));
  writeln(pmin);
}

