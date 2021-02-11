//v=300; q=10; w=30; d=30; o=30; A=sparse(v+w+o,v+o); for k=1:v; for l=1:q; A(1+floor(rand()*k),1+floor(rand()*k)) = -rand()/q; endfor; A(k,k)=1; endfor; for ww=1:w; for dd=1:d; A(v+ww,1+floor(rand()*v)) = rand(); endfor; endfor;  A((v+w+1):(v+w+o),1:v) = -1; A((v+w+1):(v+w+o),(v+1):(v+o)) = speye(o);
//d1 = sum(A!=0,1);
//d2 = sum(A!=0,2);
//scatter(d1(1:n),d2(1:n))
#include <map>
#include <set>
#include <vector>
#include <deque>
#include <stdio.h>
#include <stdlib.h>
#include <matio.h>

using namespace std;

#define SCALE 1000000

int main(int argc, char **argv) {
  if (argc < 6) {
    fprintf(stderr, "USAGE: %s n m b d o\n", argv[0]);
    return 1;
  }
  uint32_t v = atoi(argv[1]); // number of variables
  uint32_t q = atoi(argv[2]); // how many inputs to each industry, on average
  uint32_t w = atoi(argv[3]); // number of baskets
  uint32_t d = atoi(argv[4]); // number of inputs to each basket
  uint32_t o = atoi(argv[5]); // number of balance equations we're optimizing over

  uint32_t m = v+w+o;
  uint32_t N = v+o;

  mat_t *matfp = Mat_CreateVer("program.mat", NULL, MAT_FT_MAT5);
  enum matio_compression comp = MAT_COMPRESSION_NONE; //or MAT_COMPRESSION_ZLIB
  matvar_t *var;

  size_t dims1[1] = {1};
  size_t dimsc[1] = {N};
  size_t dimsb[1] = {m};

  var = Mat_VarCreate("m", MAT_C_UINT32, MAT_T_UINT32, 1, dims1, &m, MAT_F_DONT_COPY_DATA);
  Mat_VarWrite(matfp, var, comp);
  Mat_VarFree(var);
  var = Mat_VarCreate("n", MAT_C_UINT32, MAT_T_UINT32, 1, dims1, &N, MAT_F_DONT_COPY_DATA);
  Mat_VarWrite(matfp, var, comp);
  Mat_VarFree(var);
  /*var = Mat_VarCreate("b", MAT_C_DOUBLE, MAT_T_DOUBLE, 1, dimsb, bd.data(), MAT_F_DONT_COPY_DATA);
  Mat_VarWrite(matfp, var, comp);
  Mat_VarFree(var);*/


  {
    vector<float> c;
    for (uint32_t k = 0; k < v+o; k++) {
      c.push_back(k < v ? 0 : 1);
    }  var = Mat_VarCreate("c", MAT_C_DOUBLE, MAT_T_DOUBLE, 1, dimsc, c.data(), MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
  }

  {
    // A is indexed like A[col][row]
    deque<map<uint32_t, float> > A(v+o, map<uint32_t, float>());

    // fill in technical coefficients
    fprintf(stderr, "tech\n");
    for (uint32_t k = 0; k < v; k++) {
      if (k % 10000 == 0) {
        fprintf(stderr, "%i/%i\n", k/10000, v/10000);
      }
      for (uint32_t l = 0; l < q; l++) {
        A[rand() % (k+1)][rand() % (k+1)] = -1 - (rand() % (SCALE/q));
      }
      A[k][k] = SCALE;
    }
    // fill in baskets
    fprintf(stderr, "baskets\n");
    for (uint32_t ww = 0; ww < w; ww++) {
      for (uint32_t dd = 0; dd < d; dd++) {
        A[rand() % (v+1)][v+ww] = 1 + (rand() % SCALE);
      }
    }
    // fill in balance equations
    fprintf(stderr, "balance\n");
    for (uint32_t k = 0; k < o; k++) {
      for (uint32_t i = 0; i < v; i++) {
        A[i][v+w+k] = -1;
      }
      A[v+k][v+w+k] = 1;
    }

    int nnz = 0;
    for (const auto& colit : A) {
      nnz += colit.size();
    }
    fprintf(stderr, "        %i\n", nnz*(sizeof(_Rb_tree_node_base)+sizeof(uint32_t)+sizeof(float))/1024);

    // destructively convert A to Matlab format
    // this should minimize memory consumption
    vector<mat_int32_t> ir;
    vector<mat_int32_t> jc;
    vector<float> data;

    mat_int32_t j = 0;
    while (A.size() > 0) {
      auto& col = A[0];
      jc.push_back(j);
      j += col.size();
      for (auto it = col.begin(); it != col.end(); it++) {
        ir.push_back(it->first);
        data.push_back(it->second);
      }
      col.clear();
      A.pop_front();
    }
    jc.push_back(j);

    mat_sparse_t s;
    s.ndata = s.nir = s.nzmax = j;
    s.ir = ir.data();
    s.jc = jc.data();
    s.data = data.data();
    s.njc = N+1;

    size_t dims2[2] = {m, N};
    var = Mat_VarCreate("A", MAT_C_SPARSE, MAT_T_SINGLE, 2, dims2, &s, MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
  }

  Mat_Close(matfp);

  return 0;
}
