#ifndef _FFT_
#define _FFT_

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

int binary_inversion(int N, int digit) {
  int a = 0;
  for (int i = 0; i < digit; ++i) {
    a = a << 1;
    a += (N >> i) & 1;
  }
  return a;
}

std::complex<double> Exp_i(double theta) {
  double Re = std::cos(theta);
  double Im = std::sin(theta);
  return {Re, Im};
}

template <class T> T power(T const &x, long long N) {
  T result = 1;
  for (int i = 0; i < N; ++i) {
    result *= x;
  }
  return result;
}

template <class T>
std::vector<std::complex<double>> FFT_impl(int n_sample, T create_input) {
  // ------ size investigation ------ //
  assert((n_sample & (n_sample - 1)) == 0);
  auto log2 = [](int N) {
    int count = 0;
    while (N != 1) {
      N /= 2;
      ++count;
    }
    return count;
  };

  const int log2_of_n_sample = log2(n_sample);

  // ------ Rotation factor ---- //
  const double pi = 3.14159265358979; // pi
  std::vector<std::complex<double>> W(n_sample);
  for (int i = 0; i < n_sample; i++) {
    W[i] = Exp_i(-i * 2 * pi / n_sample);
  }

  // ---- Reserve array for output of butterfly operation ---- //
  std::array<std::vector<std::complex<double>>, 2> aft_btfly = {
      create_input(n_sample, log2_of_n_sample),
      std::vector<std::complex<double>>(n_sample),
  }; // include input data itself at[0];

  // ------ FFT ------ //
  for (int i = 1; i <= log2_of_n_sample; ++i) {
    int num_btfly = power(2, i - 1);
    int rot_factor = power(2, log2_of_n_sample - i);
    int num_group = rot_factor;
    for (int j = 0; j < num_btfly; ++j) {
      for (int k = 0; k < num_group; ++k) {
        int top_index = 2 * num_btfly * k + j;
        int bottom_index = 2 * num_btfly * k + j + num_btfly;
        aft_btfly[1][top_index] =
            aft_btfly[0][top_index] +
            W[rot_factor * j] * aft_btfly[0][bottom_index];
        aft_btfly[1][bottom_index] =
            aft_btfly[0][top_index] -
            W[rot_factor * j] * aft_btfly[0][bottom_index];
      } // k
    }   // j
    aft_btfly[0].swap(aft_btfly[1]);
  } // i : 1 -> log2_of_n_sample

  return std::move(aft_btfly[0]);
}

template <class InputIterator, class OutputIterator>
void FFT(int n_sample, InputIterator data_first, OutputIterator real_first,
         OutputIterator imag_first) {
  const auto result =
      FFT_impl(n_sample, [data_first](int n_sample, int log2_of_n_sanmple) {
        std::vector<std::complex<double>> input(n_sample);
        for (int i = 0; i < n_sample; ++i) {
          int bit_inv_i = binary_inversion(i, log2_of_n_sanmple);
          input[i] = {static_cast<double>(*(data_first + bit_inv_i)), 0};
        }
        return input;
      });

  // ------ COPY ------ //
  for (auto &&r : result) {
    *real_first = r.real();
    *imag_first = r.imag();
    ++real_first;
    ++imag_first;
  }

  return;
}

template <class InputIterator, class OutputIterator>
void InverseFFT(int n_sample, InputIterator in_first_Re,
                InputIterator in_first_Im, OutputIterator out_first_Re,
                OutputIterator out_first_Im) {
  const auto result = FFT_impl(
      n_sample, [in_first_Re, in_first_Im](int n_sample, int log2_of_n_sample) {
        std::vector<std::complex<double>> input(n_sample);
        for (int i = 0; i < n_sample; i++) {
          int bit_inv_i = binary_inversion(i, log2_of_n_sample);
          input[i] = {static_cast<double>(*(in_first_Re + bit_inv_i)),
                      static_cast<double>(*(in_first_Im + bit_inv_i)) * (-1)};
        }
        return input;
      });

  // ------ COPY ------ //
  for (auto &&r : result) {
    *out_first_Re = r.real() / n_sample;
    *out_first_Im = r.imag() / n_sample;
    ++out_first_Re;
    ++out_first_Im;
  }

  return;
}
#endif
