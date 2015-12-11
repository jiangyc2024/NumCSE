# include <mylibrary.hpp>

std::vector<double> sample()
{
  const std::size_t N = 10;
  std::vector<double> v(N);
  for (std::size_t i = 0; i < N; ++i)
    v[i] = i*i;
  return v;
}
