#include <map>
#include <set>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <matio.h>

using namespace std;

int main(int argc, char **argv) {
  if (argc < 4) return 1;
  uint32_t m = atoi(argv[1]);
  uint32_t n = atoi(argv[2]);
  int k = atoi(argv[3]);
  int fmt = 0;
  if (argc >= 5 && argv[4][0] == 'm') fmt = 1;
  if (argc >= 5 && argv[4][0] == 'M') fmt = 2;
  map<int, map<int,int> > A;
  vector<int32_t> b, c;
  vector<double> bd, cd;
  long nnz = 0;

  for (int i = 0; i < n; i++) {
    c.push_back(i+1);
    cd.push_back(c[i]);
  }

  for (int i = 0; i < m; i++) {
    if (i < n) {
      A[i][i] = 1 + rand() / m;
    }
    b.push_back(1 + rand());
    bd.push_back(b[i]);

    for (int j = 0; j < k; j++) {
      int f = i < n ? i : n;
      int l = rand() % (f+1);
      if (l >= n) l = n-1;
      A[i][l] = 1 + rand() / m;
    }

    nnz += A[i].size();
  }
  
  if (fmt == 0) {
    printf("max: ");
    for (int i = 0; i < n; i++) {
      if (i > 0) printf(" + ");
      printf("%i x%i", c[i], i);
    }
    printf(";\n");

    printf("// nnz=%li\n", nnz);
    for (int i = 0; i < m; i++) {
      auto& row = A[i];
      for (auto it2 = row.begin(); it2 != row.end(); it2++) {
        if (it2 != row.begin()) printf(" + ");
        printf("%i x%i", it2->second, it2->first);
      }
      printf(" <= %i;\n", b[i]);
    }
  } else if (fmt == 1) {
    printf("c = [");
    for (int i = 0; i < n; i++) {
      printf("%i;", c[i]);
    }
    printf("];\n");
    printf("b = [");
    for (int i = 0; i < m; i++) {
      printf("%i;", b[i]);
    }
    printf("];\n");
    printf("m = %i;\n", m);
    printf("n = %i;\n", n);
    printf("A = sparse(m,n);\n");
    for (int i = 0; i < m; i++) {
      auto& row = A[i];
      if (i % 1000 == 999) {
        printf("printf(\"loaded %i/%i rows so far\\n\");\n", i+1, m);
      }
      for (auto it2 = row.begin(); it2 != row.end(); it2++) {
        //printf("%i;", it2->second);
        printf("A(%i,%i)=%i;\n", i+1, it2->first+1, it2->second);
      }
    }
  } else if (fmt == 2) {
    mat_t *matfp = Mat_CreateVer("program.mat", NULL, MAT_FT_MAT5);

    enum matio_compression comp = MAT_COMPRESSION_NONE; //or MAT_COMPRESSION_ZLIB
    size_t dims1[1] = {1};
    size_t dimsc[1] = {n};
    size_t dimsb[1] = {m};

    matvar_t *var = Mat_VarCreate("m", MAT_C_UINT32, MAT_T_UINT32, 1, dims1, &m, MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
    var = Mat_VarCreate("n", MAT_C_UINT32, MAT_T_UINT32, 1, dims1, &n, MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
    var = Mat_VarCreate("b", MAT_C_DOUBLE, MAT_T_DOUBLE, 1, dimsb, bd.data(), MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
    var = Mat_VarCreate("c", MAT_C_DOUBLE, MAT_T_DOUBLE, 1, dimsc, cd.data(), MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);

    mat_sparse_t s;
    s.ndata = s.nir = s.nzmax = nnz;
    s.ir = (mat_int32_t*)malloc(sizeof(mat_int32_t)*(nnz+1));
    s.jc = (mat_int32_t*)malloc(sizeof(mat_int32_t)*(m+1));
    double *d = (double*)malloc(sizeof(double)*(nnz+1));
    s.data = d;
    s.njc = m+1;


    int j = 0;
    for (int i = 0; i < m; i++) {
      auto& row = A[i];
      //int k = j;
      s.jc[i] = j;
      for (auto it2 = row.begin(); it2 != row.end(); it2++) {
        //printf("%i;", it2->second);
        //printf("A(%i,%i)=%i;\n", i+1, it2->first+1, it2->second);
        s.ir[j] = it2->first;
        d[j] = it2->second;
        j++;
      }
    }
    s.jc[m] = nnz;
    s.ir[nnz] = n;

    size_t dims2[2] = {n, m};
    // transposed because matlab is CSC
    // gets un-transposed in code
    var = Mat_VarCreate("A", MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims2, &s, MAT_F_DONT_COPY_DATA);
    Mat_VarWrite(matfp, var, comp);
    Mat_VarFree(var);
    Mat_Close(matfp);
  }

  return 0;
}
