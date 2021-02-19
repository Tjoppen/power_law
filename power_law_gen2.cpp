//v=300; q=10; w=30; d=30; o=30; A=sparse(v+w+o,v+o); for k=1:v; for l=1:q; A(1+floor(rand()*k),1+floor(rand()*k)) = -rand()/q; endfor; A(k,k)=1; endfor; for ww=1:w; for dd=1:d; A(v+ww,1+floor(rand()*v)) = rand(); endfor; endfor;  A((v+w+1):(v+w+o),1:v) = -1; A((v+w+1):(v+w+o),(v+1):(v+o)) = speye(o);
//d1 = sum(A!=0,1);
//d2 = sum(A!=0,2);
//scatter(d1(1:n),d2(1:n))

//
/*
x = A\b;
x /= min(A(1:(v+w),:)*x./b(1:(v+w),:));
x *= 1.1;
x((v+1):n,:) = -A((v+w+1):m,1:v)*x(1:v,:) * 1.1;

if min(x) <= 0 || min(A*x - b) <= 0
  error('bad x');
endif

l = [1./b(1:(v+w)); 0.9*ones(o,1)];

if min(l) <= 0 || min(c'-l'*A) <= 0
  error('bad l');
endif

l'*b
c'*x
*/
#include <map>
#include <set>
#include <vector>
#include <deque>
#include <stdio.h>
#include <stdlib.h>
#include <matio.h>

using namespace std;

#define SCALE 1000000

typedef float real;
#define MAT_REAL MAT_T_SINGLE

static void swap(mat_int32_t *Ai, mat_int32_t *Aj, real *Ad, int a, int b) {
  if (a != b) {
    mat_int32_t x = Ai[a];
    Ai[a] = Ai[b];
    Ai[b] = x;

    x = Aj[a];
    Aj[a] = Aj[b];
    Aj[b] = x;

    real f = Ad[a];
    Ad[a] = Ad[b];
    Ad[b] = f;
  }
}

static int compare(mat_int32_t ai, mat_int32_t aj, mat_int32_t bi, mat_int32_t bj) {
  if (ai < bi) return -1;
  if (ai > bi) return 1;
  if (aj < bj) return -1;
  if (aj > bj) return 1;
  return 0;
}

static int part(mat_int32_t *Ai, mat_int32_t *Aj, real *Ad, int lo, int hi) {
  int m = (lo+hi)/2;
  mat_int32_t pi = Ai[m];
  mat_int32_t pj = Aj[m];

  int i = lo - 1;
  int j = hi + 1;

  for (;;) {
    do {
      i++;
    } while (compare(Ai[i], Aj[i], pi, pj) < 0);
    do {
      j--;
    } while (compare(Ai[j], Aj[j], pi, pj) > 0);
    if (i >= j) {
      return j;
    }
    swap(Ai, Aj, Ad, i, j);
  }
}

static void qsort2(mat_int32_t *Ai, mat_int32_t *Aj, real *Ad, int lo, int hi, int level) {
  if (level < 6) {
    fprintf(stderr, ".");
    fflush(stderr);
  }
  if (lo < hi) {
    int p = part(Ai, Aj, Ad, lo, hi);
    
    qsort2(Ai, Aj, Ad, lo, p, level+1);
    qsort2(Ai, Aj, Ad, p+1, hi, level+1);
  }
}

static int deduplicate(mat_int32_t *Ai, mat_int32_t *Aj, real *Ad, int n) {
  int j = 0;
  for (int i = 1; i < n;) {
    if (Ai[j] == Ai[i] && Aj[j] == Aj[i]) {
      if (Ad[i] > Ad[j]) {
        Ad[j] = Ad[i];
      }
      i++;
    } else {
      j++;
      swap(Ai, Aj, Ad, i, j);
      i++;
    }
  }
  return j+1;
}

int main(int argc, char **argv) {
  if (argc < 6) {
    fprintf(stderr, "USAGE: %s num_variables inputs_per_industry num_baskets inputs_per_basket num_opt [Hawkins-Simon_constant]\n", argv[0]);
    return 1;
  }
  uint32_t v = atoi(argv[1]); // number of variables
  uint32_t q = atoi(argv[2]); // how many inputs to each industry, on average
  uint32_t w = atoi(argv[3]); // number of baskets
  uint32_t d = atoi(argv[4]); // number of inputs to each basket
  uint32_t o = atoi(argv[5]); // number of balance equations we're optimizing over

  // sometimes it is necessary to scale the coefficients down somewhat to satisfy the Hawkins-Simon condition
  // or, at least scale them down enough that bicgstab(A(1:v,1:v), b(1:v)) gives a strictly positive x
  int hk = SCALE/q;
  if (argc >= 7) {
    hk *= atof(argv[6]);
  } else if (v > 10000) {
    hk /= 2;  // determined impirically
  }

  uint32_t m = v+w+o;
  uint32_t n = v+o;

  mat_t *matfp = Mat_CreateVer("program.mat", NULL, MAT_FT_MAT5);
  enum matio_compression comp = MAT_COMPRESSION_NONE; //or MAT_COMPRESSION_ZLIB
  matvar_t *var;

  size_t dims1[1] = {1};
  size_t dimsc[1] = {n};
  size_t dimsb[1] = {m};

#define WRITE_UINT32(name) do {\
    var = Mat_VarCreate(#name, MAT_C_UINT32, MAT_T_UINT32, 1, dims1, &name, MAT_F_DONT_COPY_DATA); \
    Mat_VarWrite(matfp, var, comp); \
    Mat_VarFree(var); \
  } while (0)
  
  WRITE_UINT32(m);
  WRITE_UINT32(n);
  WRITE_UINT32(v);
  WRITE_UINT32(q);
  WRITE_UINT32(w);
  WRITE_UINT32(d);
  WRITE_UINT32(o);

  {
    vector<real> b;
    for (uint32_t k = 0; k < m; k++) {
      if (k < v) {
        b.push_back(SCALE);
      } else if (k < v+w) {
        b.push_back(SCALE);
      } else {
        b.push_back(0);
      }
    }
    var = Mat_VarCreate("b", MAT_C_DOUBLE, MAT_REAL, 1, dimsb, b.data(), MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
  }


  {
    vector<real> c;
    for (uint32_t k = 0; k < v+o; k++) {
      c.push_back(k < v ? 0 : 1);
    }
    var = Mat_VarCreate("c", MAT_C_DOUBLE, MAT_REAL, 1, dimsc, c.data(), MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
  }

  {
    int cap = v*(1+q) + w*d + o*(v+1);
    vector<mat_int32_t> Ai;
    vector<mat_int32_t> Aj;
    vector<real> Ad;
    Ai.reserve(cap);
    Aj.reserve(cap);
    Ad.reserve(cap);

    //2221060
    //2242120

    // fill in technical coefficients
    fprintf(stderr, "tech\n");
    for (uint32_t k = 0; k < v; k++) {
      if (k % 10000 == 0) {
        fprintf(stderr, "%i/%i\n", k/10000, v/10000);
      }
      Ai.push_back(k);
      Aj.push_back(k);
      Ad.push_back(SCALE);
      for (uint32_t l = 0; l < q; l++) {
        Ai.push_back(rand() % (k+1));
        Aj.push_back(rand() % (k+1));
        Ad.push_back(-1 - (real)(rand() % hk));
      }
    }
    // fill in baskets
    fprintf(stderr, "baskets\n");
    for (uint32_t ww = 0; ww < w; ww++) {
      if (ww % 10 == 0 && ww > 0) {
        fprintf(stderr, "%3i/%3i\n", ww, w);
      }
      for (uint32_t dd = 0; dd < d; dd++) {
        Ai.push_back(rand() % v);
        Aj.push_back(v+ww);
        Ad.push_back(1 + (real)(rand() % SCALE));
      }
    }
    // fill in balance equations
    fprintf(stderr, "balance\n");
    for (uint32_t k = 0; k < o; k++) {
      if (k % 10 == 0 && k > 0) {
        fprintf(stderr, "%3i/%3i\n", k, o);
      }
      for (uint32_t i = 0; i < v; i++) {
        Ai.push_back(i);
        Aj.push_back(v+w+k);
        Ad.push_back(-SCALE);
      }
      Ai.push_back(v+k);
      Aj.push_back(v+w+k);
      Ad.push_back(SCALE);
    }

    /*int nnz = 0;
    for (const auto& colit : Ai) {
      nnz += colit.size();
    }*/
    int nnz = Ai.size();
    fprintf(stderr, "        %i\n", nnz*(2*sizeof(uint32_t)+sizeof(real))/1024);


    fprintf(stderr, "sorting\n");
    qsort2(Ai.data(), Aj.data(), Ad.data(), 0, nnz-1, 0);

    /*for (int x = 0; x < nnz; x++) {
      fprintf(stderr, "%i %i\n", Ai[x], Aj[x]);
    }*/

    // verify sort
    fprintf(stderr, "\nverify\n");
    for (int x = 1; x < nnz; x++) {
      if (compare(Ai[x-1], Aj[x-1], Ai[x], Aj[x]) > 0) {
        fprintf(stderr, "bad @ %i: %i,%i %i,%i\n", x, Ai[x-1], Aj[x-1], Ai[x], Aj[x]);
        return 1;
      }
    }

    fprintf(stderr, "deduplicate\n");
    int nnz2 = deduplicate(Ai.data(), Aj.data(), Ad.data(), nnz);
    fprintf(stderr, "%i -> %i\n", nnz, nnz2);
    nnz = nnz2;
    /*for (int x = 0; x < nnz; x++) {
      fprintf(stderr, "%i %i %f\n", Ai[x], Aj[x], Ad[x]);
    }*/

    fprintf(stderr, "verify2\n");
    for (int x = 1; x < nnz; x++) {
      if (compare(Ai[x-1], Aj[x-1], Ai[x], Aj[x]) >= 0) {
        fprintf(stderr, "bad @ %i: %i,%i %i,%i\n", x, Ai[x-1], Aj[x-1], Ai[x], Aj[x]);
        return 1;
      }
    }

    vector<mat_int32_t> jc;
    int ofs = 0;
    for (int col = 0; col < n; col++) {
      //fprintf(stderr, "%i @ %i\n", col, ofs);
      jc.push_back(ofs);
      while (Ai[ofs] == col && ofs < nnz) {
        ofs++;
      }
    }
    if (ofs != nnz) {
      fprintf(stderr, "nnz mismatch %i %i\n", ofs, nnz);
    }
    jc.push_back(ofs);

    mat_sparse_t s;
    s.ndata = s.nir = s.nzmax = nnz;
    s.ir = Aj.data();
    s.jc = jc.data();
    s.data = Ad.data();
    s.njc = n+1;

    size_t dims2[2] = {m, n};
    var = Mat_VarCreate("A", MAT_C_SPARSE, MAT_REAL, 2, dims2, &s, MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
  }

  Mat_Close(matfp);

  return 0;
}
