//
//  mds.cpp
//  MR-MEGA
//
//  Created by reedik on 03/08/2016.
//  Copyright © 2016 Estonian Genome Center. All rights reserved.
//

#include <vector>
#include <math.h>
#include <iostream>
#include <map>
#include "mds.h"
#include "structures.h"




double
pythag(const double a, const double b)
{
    double absa,absb;
    
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+pow(absb/absa,2));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+pow(absa/absb,2)));
}



matrixD
calculateMDS(matrixD D, int pcCount){
    int nc =D.getCols();
    double mean = 0;
    arrayD M(nc);

    for (int c1 = 0 ; c1 < nc ; c1++)
    {
        for (int c2 = 0 ; c2 < nc ; c2++)
        {
            M.put(c1, M.get(c1)+ D.get(c1,c2));
        }
        M.put(c1,M.get(c1) /(double)nc);
        mean += M.get(c1);
    }
    mean /= (double)nc;
   // cout << "mean: " << mean << endl;
    // For each element for D, double center

    //cout << "D:"<<endl;
    //D.print();
    for (int c1 = 0 ; c1 < nc ; c1++)
    for (int c2 = c1 ; c2 < nc ; c2++)
    {
//test change here!
        double putter = - 0.5 * ( D.get(c1,c2) - M.get(c1) - M.get(c2) + mean );
        D.put(c1,c2, putter);
        D.put(c2,c1, putter);
//        D.put(c1,c2, - 0.5 * ( D.get(c1,c2) - M.get(c1) - M.get(c2) + mean ));
//        D.put(c2,c1, - 0.5 * ( D.get(c1,c2) - M.get(c1) - M.get(c2) + mean ));
    }
  //  cout << "D double centered:"<<endl;
  //  D.print();

    
    arrayD eigenvalue(nc);
    matrixD eigenvector(nc,nc);

    bool flag = svd(D,eigenvalue,eigenvector);
    
//    cout << "Eigenvalue:"<<endl;
//    eigenvalue.print();
//    cout << "Eigenvector:"<<endl;
//    eigenvector.print();
    
    map<double,int> emap;
    for (int i=0; i<nc; i++)
        emap.insert(std::make_pair( eigenvalue.get(i) , i ) );
    
    map<double,int>::reverse_iterator e = emap.rbegin();
    int inc = pcCount;
    
    vector<int> elist;
    while ( e != emap.rend() && inc > 0 )
    {
        elist.push_back(e->second);
        inc--;
        e++;
    }
    
    if (pcCount < 1)
        pcCount = 1;
    if (pcCount > nc)
        pcCount = nc;
    
    if ( elist.size() != pcCount )
    {
        cerr << "Internal problem extracting MDS solution\n";
        elist.resize(pcCount);
    }
    
    
    // Sqrt(D)
    arrayD sqrt_eigenvalue(nc);
    
    for (int i=0; i<nc; i++)
    {
        if (eigenvalue.get(i)>=0)sqrt_eigenvalue.put(i,sqrt(eigenvalue.get(i)));
        else sqrt_eigenvalue.put(i,0.0);
    }
    
    // Make solution
    // EVEC * sqrt(EVAL)  but filter on rows that are in solution
    // with EVAL as diagonal matrix
    
    matrixD mds(nc,pcCount);
    
    for (int c1=0; c1<nc; c1++)
        for (int c2=0; c2<pcCount; c2++)
            // ERROR: *** in 1.02 and below was ***   for (int c3=0; c3<nc; c3++)
            for (int c3=0; c3<pcCount; c3++)
            {
                int i2 = elist[c2];
                int i3 = elist[c3];
                if ( i3 == i2 )
                  //  mds.put(c1,c2,mds.get(c1,c2)+eigenvector.get(c1,i3) * sqrt_eigenvalue.get(i2));
                    mds.put(c1,c2,mds.get(c1,c2)+eigenvector.get(i3,c1) * sqrt_eigenvalue.get(i2));
            }
    

    return mds;
}

vector<double> eigenvalues(vector<vector<double> > & a)
{
    
    // 'a' should be a square, symmetric matrix
    int n=a.size();
    vector<double> e(n);
    vector<double> d(n);
    tred2(a,d,e);
    vector<vector<double> > z; // dummy
    tqli(d,e,z);
    return d;
}

// Householder method to reduce real, symmetric matrix
// to tridiagonal form
// Modified to return only eigenvalues.
void tred2(vector<vector<double> > & a,
           vector<double> & d,
           vector<double> &e)
{
    int l,k,j,i;
    double scale,hh,h,g,f;
    
    int n=d.size();
    for (i=n-1;i>0;i--) {
        l=i-1;
        h=scale=0.0;
        if (l > 0) {
            for (k=0;k<l+1;k++)
                scale += fabs(a[i][k]);
            if (scale == 0.0)
                e[i]=a[i][l];
            else {
                for (k=0;k<l+1;k++) {
                    a[i][k] /= scale;
                    h += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i]=scale*g;
                h -= f*g;
                a[i][l]=f-g;
                f=0.0;
                for (j=0;j<l+1;j++) {
                    // Next statement can be omitted if eigenvectors not wanted
                    // 	  a[j][i]=a[i][j]/h;
                    g=0.0;
                    for (k=0;k<j+1;k++)
                        g += a[j][k]*a[i][k];
                    for (k=j+1;k<l+1;k++)
                        g += a[k][j]*a[i][k];
                    e[j]=g/h;
                    f += e[j]*a[i][j];
                }
                hh=f/(h+h);
                for (j=0;j<l+1;j++) {
                    f=a[i][j];
                    e[j]=g=e[j]-hh*f;
                    for (k=0;k<j+1;k++)
                        a[j][k] -= (f*e[k]+g*a[i][k]);
                }
            }
        } else
            e[i]=a[i][l];
        d[i]=h;
    }
    // Next statement can be omitted if eigenvectors not wanted
    //   d[0]=0.0;
    e[0]=0.0;
    // Contents of this loop can be omitted if eigenvectors not
    //	wanted except for statement d[i]=a[i][i];
    for (i=0;i<n;i++) {
        //     l=i;
        //     if (d[i] != 0.0) {
        //       for (j=0;j<l;j++) {
        // 	g=0.0;
        // 	for (k=0;k<l;k++)
        // 	  g += a[i][k]*a[k][j];
        // 	for (k=0;k<l;k++)
        // 	  a[k][j] -= g*a[k][i];
        //       }
        //     }
        d[i]=a[i][i];
        //     a[i][i]=1.0;
        //     for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
    }
}

// Modified to return only eigenvalues.
void tqli(vector<double> &d, vector<double>&e,
          vector<vector<double> > &z)
{
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    double volatile temp;
    int n=d.size();
    for (i=1;i<n;i++) e[i-1]=e[i];
    e[n-1]=0.0;
    for (l=0;l<n;l++) {
        iter=0;
        do {
            for (m=l;m<n-1;m++) {
                dd=fabs(d[m])+fabs(d[m+1]);
                temp=fabs(e[m])+dd;
                if (temp == dd) break;
            }
            if (m != l) {
                if (iter++ == 30) cerr << "Internal problem in tqli routine" << endl;
                g=(d[l+1]-d[l])/(2.0*e[l]);
                r=pythag(g,1.0);
                g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=l;i--) {
                    f=s*e[i];
                    b=c*e[i];
                    e[i+1]=(r=pythag(f,g));
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m]=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2.0*c*b;
                    d[i+1]=g+(p=s*r);
                    g=c*r-b;
                    // Next loop can be omitted if eigenvectors not wanted
                    /* for (k=0;k<n;k++) {
                     f=z[k][i+1];
                     z[k][i+1]=s*z[k][i]+c*f;
                     z[k][i]=c*z[k][i]-s*f;
                     } */
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l]=g;
                e[m]=0.0;
            }
        } while (m != l);
    }
}

bool svd(matrixD & u, arrayD &w, matrixD &v)
{
    
    // #ifdef WITH_LAPACK
    //   matrix_t u2;
    //   svd_lapack(u,w,u2,v);
    //   u = u2;
    //#else
    
    const double eps = 1e-12;
    
    if (u.size() == 0)
        cerr << "Internal problem: matrix with no rows in svd()" << endl;
    
 /* SHOULD BE DONE BEFORE
    int r = u.size();
    int c = u[0].size();
    w.resize(c);
    sizeMatrix(v,c,c);
 */

    vector <vector <double> > u1 = u.toVector();
    vector <double> w1 = w.toVector();
    vector <vector <double> > v1 = v.toVector();
    
    
    bool flag = svdcmp(u1,w1,v1);
    
    u.putFromVector(u1);
    w.putFromVector(w1);
    v.putFromVector(v1);
    
//    v.print();
    
    return flag;
    
    // Look for singular values
    //   double wmax = 0;
    //   for (int i=0; i<n; i++)
    //     wmax = w[i] > wmax ? w[i] : wmax;
    //   double wmin = wmax * eps;
    //   for (int i=0; i<n; i++)
    //     {
    //       w[i] = w[i] < wmin ? 0 : 1/w[i];
    //     }  
    
    
    // #endif
}


void svdvar(vector<vector<double> > & v,
            vector<double> & w,
            vector<vector<double> > & cvm)
{
    int i,j,k;
    double sum;
    
    int ma=w.size();
    vector<double> wti(ma);
    for (i=0;i<ma;i++) {
        wti[i]=0.0;
        if (w[i] != 0.0) wti[i]=1.0/(w[i]*w[i]);
    }
    for (i=0;i<ma;i++) {
        for (j=0;j<i+1;j++) {
            sum=0.0;
            for (k=0;k<ma;k++)
                sum += v[i][k]*v[j][k]*wti[k];
            cvm[j][i]=cvm[i][j]=sum;
        }
    }
}

bool svdcmp(vector<vector<double> > & a,
            vector<double> & w,
            vector<vector<double> > &v)
{
    bool flag;
    int i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;
    double volatile temp;
    
    int m=a.size();
    if (m==0) cerr<<"Internal problem in SVD function (no observations left?)"<<endl;
    int n=a[0].size();
    
    vector<double> rv1(n);
    g=scale=anorm=0.0;
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += fabs(a[k][i]);
            if (scale != 0.0) {
                for (k=i;k<m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l-1;j<n;j++) {
                    for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i+1 <= m && i+1 != n) {
            for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
            if (scale != 0.0) {
                for (k=l-1;k<n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l-1];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l-1]=f-g;
                for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
                for (j=l-1;j<m;j++) {
                    for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
                    for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l-1;k<n;k++) a[i][k] *= scale;
            }
        }
        anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n-1;i>=0;i--) {
        if (i < n-1) {
            if (g != 0.0) {
                for (j=l;j<n;j++)
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=MIN(m,n)-1;i>=0;i--) {
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) a[i][j]=0.0;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {
                for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<m;j++) a[j][i] *= g;
        } else for (j=i;j<m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
    for (k=n-1;k>=0;k--) {
        for (its=0;its<30;its++) {
            flag=true;
            for (l=k;l>=0;l--) {
                nm=l-1;
                temp=fabs(rv1[l])+anorm;
                if (temp == anorm) {
                    flag=false;
                    break;
                }
                temp=fabs(w[nm])+anorm;
                if (temp == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    temp = fabs(f)+anorm;
                    if (temp == anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=0;j<m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j=0;j<n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 29)
                return false; // cannot converge: multi-collinearity?
            x=w[l];
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=0;jj<n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z;
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    return true;
}





/*
bool svdcmp(matrixD & a,
            arrayD & w,
            matrixD &v)
{
    
    cout << "A,W,V" << endl;
    a.print();
    w.print();
    v.print();
    
    
    bool flag;
    int i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;
    double volatile temp;
    
    int m=a.size();
    if (m==0) cerr << "Internal problem in SVD function (no observations left?)" << endl;
    int n=a.getCols();  //TO CHECK IF COLS, NOT ROWS
    
    vector<double> rv1(n);
    g=scale=anorm=0.0;
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=(scale*g);
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += fabs(a.get(k,i));
            if (scale != 0.0) {
                for (k=i;k<m;k++) {
                    a.put(k,i,a.get(k,i) / scale);
                    s += pow(a.get(k,i),2);
                }
                f=a.get(i,i);
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a.put(i,i,f-g);
                for (j=l-1;j<n;j++) {
                    for (s=0.0,k=i;k<m;k++) s += a.get(k,i)*a.get(k,j);
                    f=s/h;
                    for (k=i;k<m;k++) a.put(k,j, a.get(k,j)+f*a.get(k,i));
                }
                for (k=i;k<m;k++) a.put(k,i,a.get(k,i)*scale);
            }
        }
        w.put(i,(scale *g));
        g=s=scale=0.0;
        if (i+1 <= m && i+1 != n) {
            for (k=l-1;k<n;k++) scale += fabs(a.get(i,k));
            if (scale != 0.0) {
                for (k=l-1;k<n;k++) {
                    a.put(i,k,a.get(i,k) / scale);
                    s += pow(a.get(i,k),2);
                }
                f=a.get(i,l-1);
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a.put(i,l-1,f-g);
                for (k=l-1;k<n;k++) rv1[k]=a.get(i,k)/h;
                for (j=l-1;j<m;j++) {
                    for (s=0.0,k=l-1;k<n;k++) s += a.get(j,k)*a.get(i,k);
                    for (k=l-1;k<n;k++) a.put(j,k,a.get(j,k)+s*rv1[k]);
                }
                for (k=l-1;k<n;k++) a.put(i,k,a.get(i,k)* scale);
            }
        }
        anorm=MAX(anorm,(fabs(w.get(i))+fabs(rv1[i])));
    }
    for (i=n-1;i>=0;i--) {
        if (i < n-1) {
            if (g != 0.0) {
                for (j=l;j<n;j++)
                    v.put(j,i,(a.get(i,j)/a.get(i,l))/g);
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s += a.get(i,k)*v.get(k,j);
                    for (k=l;k<n;k++) v.put(k,j,v.get(k,j) + s*v.get(k,i));
                }
            }
            for (j=l;j<n;j++){ v.put(i,j,0.0);v.put(j,i,0.0);}
        }
        v.put(i,i,1.0);
        g=rv1[i];
        l=i;
    }
    for (i=MIN(m,n)-1;i>=0;i--) {
        l=i+1;
        g=w.get(i);
        for (j=l;j<n;j++) a.put(i,j,0.0);;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {
                for (s=0.0,k=l;k<m;k++) s += a.get(k,i)*a.get(k,j);
                f=(s/a.get(i,i))*g;
                for (k=i;k<m;k++) a.put(k,j,a.get(k,j)+f*a.get(k,i));
            }
            for (j=i;j<m;j++) a.put(j,i,a.get(j,i)*g);
        } else for (j=i;j<m;j++) a.put(j,i,0.0);
        a.put(i,i,a.get(i,i)+1);
    }
    for (k=n-1;k>=0;k--) {
        for (its=0;its<30;its++) {
            flag=true;
            for (l=k;l>=0;l--) {
                nm=l-1;
                temp=fabs(rv1[l])+anorm;
                if (temp == anorm) {
                    flag=false;
                    break;
                }
                temp=fabs(w.get(nm))+anorm;
                if (temp == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=l;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    temp = fabs(f)+anorm;
                    if (temp == anorm) break;
                    g=w.get(i);
                    h=pythag(f,g);
                    w.put(i,h);
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=0;j<m;j++) {
                        y=a.get(j,nm);
                        z=a.get(j,i);
                        a.put(j,nm,y*c+z*s);
                        a.put(j,i,z*c-y*s);
                    }
                }
            }
            z=w.get(k);
            if (l == k) {
                if (z < 0.0) {
                    w.put(k, -z);
                    for (j=0;j<n;j++) v.put(j,k,  -v.get(j,k));
                }
                break;
            }
            if (its == 29) 
                return false; // cannot converge: multi-collinearity?
            x=w.get(l);
            nm=k-1;
            y=w.get(nm);
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w.get(i);
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=0;jj<n;jj++) {
                    x=v.get(jj,j);
                    z=v.get(jj,i);
                    v.put(jj,j,x*c+z*s);
                    v.put(jj,i,z*c-x*s);
                }
                z=pythag(f,h);
                w.put(j,z);
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=a.get(jj,j);
                    z=a.get(jj,i);
                    a.put(jj,j,y*c+z*s);
                    a.put(jj,i,z*c-y*s);
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w.put(k,x);
        }
    }
    return true;
}
*/
