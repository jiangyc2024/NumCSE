#include "ColumnMajorMatrix.hpp"
#include <random>

// Constructor:
ColumnMajorMatrix::ColumnMajorMatrix(int _n, int _m) : n(_n), m(_m) {
  data.resize(static_cast<size_t>(_n) * static_cast<size_t>(_m), 0);
}

// Copy Constructor:
ColumnMajorMatrix::ColumnMajorMatrix(const ColumnMajorMatrix &B) = default;

//Move Constructor
ColumnMajorMatrix::ColumnMajorMatrix(ColumnMajorMatrix &&B) noexcept
    : data(std::move(B.data)), n(B.n), m(B.m) {
}

// Copy assignment operator:
ColumnMajorMatrix &ColumnMajorMatrix::operator=(const ColumnMajorMatrix &B) {
  if(this == &B) {
    return *this;
  }
  assert(this->n == B.n && this->m == B.m);
  for (int i = 0; i < n * m; ++i) {
    data[i] = B.data[i];
  }
  return *this;
}

// Move assignment operator:
ColumnMajorMatrix & ColumnMajorMatrix::operator=(ColumnMajorMatrix &&B) noexcept {
  assert(this->n == B.n && this->m == B.m);
  data = std::move(B.data);
  return *this;
}

// Function to initialize:
void ColumnMajorMatrix::initGrow() {
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      operator()(i, j) = static_cast<Real>(n * j + i);
    }
  }
}

// Function to initialize:
void ColumnMajorMatrix::initRand() {
  
  std::default_random_engine engine;
  engine.seed(42);
  std::uniform_real_distribution<Real> dist(0.0,1.0);

  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      operator()(i, j) = dist(engine);
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
      std::cout << operator()(i, j);
    }
    std::cout << "\n";
  }
}
