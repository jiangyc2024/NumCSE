#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <stdint.h>
#include <string.h>

//using gsl::owner without actual Guideline Support Library
namespace gsl { 

  template<class T>
  using owner = T; 
} //namespace gsl

namespace imread {


inline void fread_or_exit(void * const ptr, size_t const size, size_t const nmemb, FILE * const stream, const std::string & filename){

  const auto bytes_read = fread(ptr, size, nmemb, stream);
  if (bytes_read < 54) {
    std::cerr << "Error reading file " << filename << std::endl;
    std::quick_exit(1);
  }
}

// reads a bmp file and returns a Eigen::MatrixXi formated like:
// 00000000bbbbbbbbggggggggrrrrrrrr modfied version of
// http://stackoverflow.com/questions/9296059/read-pixel-value-in-bmp-file
// See also https://solarianprogrammer.com/2018/11/19/cpp-reading-writing-bmp-images/
inline Eigen::MatrixXi readBMP(const std::string & filename) {
  gsl::owner<FILE *> f = fopen(filename.c_str(), "rbe");

  if (f == nullptr) {
    std::cerr << "Error reading file " << filename << std::endl;
    std::quick_exit(1);
  }

  //unsigned char info[54];
  std::array<uint8_t, 54> info{};
  fread_or_exit(info.data(), sizeof(uint8_t), 54, f, filename); // read the 54-byte header

  // extract image height and width from header
  //NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
  const auto width = static_cast<Eigen::Index>(*reinterpret_cast<int *>(&info[18]));
  const auto height = static_cast<Eigen::Index>(*reinterpret_cast<int *>(&info[22]));
  //NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

  Eigen::MatrixXi mat(height, width);
  
  //round up row size to multiples of 4, 24 bits per pixel
  const size_t alignment_mask = 0b11;
  const size_t color_depth = 24/8;
  const size_t row_size = width * color_depth;
  const size_t row_size_padded = (row_size + alignment_mask) & ~alignment_mask;

  std::vector<uint8_t> data(row_size_padded);

  for (Eigen::Index i = 0; i < height; i++) {
    fread_or_exit(data.data(), sizeof(uint8_t), row_size_padded, f, filename);
    for (Eigen::Index j = 0; j < width * 3; j += 3) {
      const size_t blue = data[j] << 16U;
      const size_t green = data[j + 1] << 8U;
      const size_t red = data[j + 2];
      mat(height - i - 1, j / 3) = static_cast<int>(blue | green | red);
    }
  }

  static_cast<void>(fclose(f));
  return mat;
}

// extracts a color channel from an image
inline Eigen::MatrixXd getcolor(Eigen::MatrixXi img, uint8_t channel) {
  assert(channel >= 0 && channel <= 2);
  Eigen::MatrixXd ret(img.rows(), img.cols());
  for (int i = 0; i < img.rows(); ++i) {
    for (int j = 0; j < img.cols(); ++j) {
      const size_t p = (static_cast<size_t>(img(i, j)) >> (channel * 8U)) & 0xffU;
      ret(i, j) = static_cast<double>(p);
    }
  }

  return ret;
}

// converts rgb to greyscale
inline Eigen::MatrixXd greyscale(Eigen::MatrixXi img) {
  Eigen::MatrixXd grey(img.rows(), img.cols());
  for (int i = 0; i < img.rows(); ++i) {
    for (int j = 0; j < img.cols(); ++j) {
      const size_t p = img(i, j);
      grey(i, j) =
          static_cast<double>(((p >> 16U) & 0xffU) + ((p >> 8U) & 0xffU) + (p & 0xffU)) / 3.;
    }
  }

  return grey;
}


} //namespace imread