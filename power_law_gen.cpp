#include <map>
#include <set>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
int main(int argc, char **argv) {
  if (argc < 4) return 1;
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int k = atoi(argv[3]);
  int fmt = 0;
  if (argc >= 5 && argv[4][0] == 'm') fmt = 1;
  map<int, map<int,int> > A;
  map<int,int> b, c;
  long nnz = 0;

  for (int i = 0; i < n; i++) {
    c[i] = i+1;
  }

  for (int i = 0; i < m; i++) {
    if (i < n) {
      A[i][i] = 1 + rand() / m;
    }
    b[i] = 1 + rand();

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
  } else {
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
    printf("Ai = [");
    for (int i = 0; i < m; i++) {
      auto& row = A[i];
      for (size_t j = 0; j < row.size(); j++) {
        printf("%i;", i+1);
      }
    }
    printf("];\n");
    printf("Aj = [");
    for (int i = 0; i < m; i++) {
      auto& row = A[i];
      for (auto it2 = row.begin(); it2 != row.end(); it2++) {
        printf("%i;", it2->first+1);
      }
    }
    printf("];\n");
    printf("Av = [");
    for (int i = 0; i < m; i++) {
      auto& row = A[i];
      for (auto it2 = row.begin(); it2 != row.end(); it2++) {
        printf("%i;", it2->second);
      }
    }
    printf("];\n");
    printf("m = %i;\n", m);
    printf("n = %i;\n", n);
    printf("A = sparse(Ai,Aj,Av,m,n);\n");
    
  }

  return 0;
}
