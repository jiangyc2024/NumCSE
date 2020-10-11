#include "ColumnMajorMatrix.hpp"

// Constructor:
ColumnMajorMatrix::ColumnMajorMatrix(int _n, int _m) : n(_n), m(_m) {
  data = new Real[n * m];
  for (int i = 0; i < n * m; ++i) {
    data[i] = 0;
  }
}

// Copy Constructor:
ColumnMajorMatrix::ColumnMajorMatrix(const ColumnMajorMatrix &B)
    : n(B.n), m(B.m) {
  data = new Real[n * m];
  for (int i = 0; i < n * m; ++i) {
    data[i] = B.data[i];
  }
}

// Destructor:
ColumnMajorMatrix::~ColumnMajorMatrix() { delete[] data; }

// Assignment operator:
ColumnMajorMatrix &ColumnMajorMatrix::operator=(const ColumnMajorMatrix &B) {
  assert(this->n == B.n && this->m == B.m);
  for (int i = 0; i < n * m; ++i) {
    data[i] = B.data[i];
  }
  return *this;
}

// Function to initialize:
void ColumnMajorMatrix::initGrow() {
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      operator()(i, j) = (Real)(n * j + i);
    }
  }
}

// Function to initialize:
void ColumnMajorMatrix::initRand() {
  srand(42);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      operator()(i, j) = (Real)(arc4random());
    }
  }
}
Real ColumnMajorMatrix::CalcErr(const ColumnMajorMatrix &B) {
  Real Err(0);
  for (int i = 0; i < n * m; ++i) {
    Err += fabs(data[i] - B.data[i]);
  }
  return Err;
}

void ColumnMajorMatrix::print() {
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      printf("%g ", operator()(i, j));
    }
    printf("\n");
  }
}
