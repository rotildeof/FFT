#ifndef _FFT_
#define _FFT_

#include <iostream>
#include <cassert>
#include <complex>
#include <cmath>
#include <vector>

int binary_inversion(int N, int digit){
  int a = 0;
  for(int i = 0 ; i < digit ; i++){
    a = a << 1;
    a += (N >> i) & 1;
  }
  return a;
}

std::complex<double> Exp_i(double theta){
  double Re = std::cos(theta);
  double Im = std::sin(theta);
  return {Re, Im};
}

template <class T>
T power(T const &x, long long N){
  T result = 1;
  for(int i = 0 ; i < N ; i++){
    result *= x;
  }
  return result;
}

template <class InputIterator, class OutputIterator>
  void FFT(int n_sample, InputIterator  data_first,
	   OutputIterator real_first, OutputIterator imag_first){
  // ------ size investigation ------ //
  assert( (n_sample & (n_sample - 1) ) == 0);
  auto log2 = [](int N)
    {
     int count = 0;
     while(N != 1){
       N /= 2;
       ++count;
     }
     return count;
    };
  
  int n = log2(n_sample);
  
  // ------ Rotation factor ---- //
  const double pi = 3.14159265358979; // pi
  std::vector<std::complex<double> > W(n_sample);
  for(int i = 0; i < n_sample ; i++){
    W[i] = Exp_i( - i * 2 * pi / n_sample);
  }
  
  // ---- Reserve array for output of butterfly operation ---- //
  // --> The results of butterfly operations are input to
  // --> "aft_btfly[0][]" and "aft_btfly[1][]" alternately for saving memory.
  
  std::vector<std::vector<std::complex<double> > > aft_btfly(2, std::vector<std::complex<double> >(n_sample) );
  
  // ---- Data as complex numbers -------//
  // --> Input data to aft_btfly[0][i] 
  for(int i = 0 ; i < n_sample ; i++){
    int bit_inv_i = binary_inversion(i, n);
    aft_btfly[0][i] = { static_cast<double>(*(data_first + bit_inv_i)), 0 };
  }
  
  // ------ FFT ------ //
  for(int i = 1 ; i <= n ; i++){
    int i_next =  i      % 2;
    int i_prev = (i - 1) % 2;
    int num_btfly  = power(2, i - 1);
    int rot_factor = power(2, n - i);
    int num_group  = rot_factor;
    for(int j = 0 ; j < num_btfly ; ++j){
      for(int k = 0 ; k < num_group ; ++k){
	int top_index      = 2 * num_btfly * k + j;
	int bottom_index   = 2 * num_btfly * k + j + num_btfly;
	aft_btfly[i_next][top_index]    = aft_btfly[i_prev][top_index] + W[rot_factor * j] * aft_btfly[i_prev][bottom_index];
	aft_btfly[i_next][bottom_index] = aft_btfly[i_prev][top_index] - W[rot_factor * j] * aft_btfly[i_prev][bottom_index];
      }// k
    }// j
  }// i : 1 -> n
  
  // ------ COPY ------ //
  for(auto &&r : aft_btfly[n % 2]){
    *real_first = r.real();
    *imag_first = r.imag();
    ++real_first;
    ++imag_first;
  }
  
  return;
}

template <class InputIterator, class OutputIterator>
void InverseFFT(int n_sample, InputIterator in_first_Re, InputIterator in_first_Im,
		OutputIterator out_first_Re, OutputIterator out_first_Im){
  // ------ size investigation ------ //
  assert( ( n_sample & (n_sample - 1) ) == 0);
  auto log2 = [](int N)
	      {
		int count = 0;
		while(N != 1){
		  N /= 2;
		  ++count;
		}
		return count;
	      };
  
  int n = log2(n_sample);
  
  // ------ Rotation factor ---- //
  const double pi = 3.14159265358979; // pi
  std::vector<std::complex<double> > W(n_sample);
  for(int i = 0; i < n_sample ; i++){
    W[i] = Exp_i( - i * 2 * pi / n_sample);
  }
  
  // ---- Reserve array for output of butterfly operation ---- //
  // --> include input data itself at[0];
  std::vector<std::vector< std::complex<double> > > aft_btfly(2, std::vector<std::complex<double> >(n_sample) );
  
  // ---- Data as complex numbers -------//
  /////// Input data to aft_btfly[0][i] ///////
  for(int i = 0 ; i < n_sample ; i++){
    int bit_inv_i = binary_inversion(i, n);
    aft_btfly[0][i] = { static_cast<double>(*(in_first_Re + bit_inv_i) ), static_cast<double>( *(in_first_Im + bit_inv_i) ) * (-1)};
  }
  
  // ------ Inverse FFT ------ //
  for(int i = 1 ; i <= n ; i++){
    int i_next =  i      % 2;
    int i_prev = (i - 1) % 2;
    int num_btfly  = power(2, i - 1);
    int rot_factor = power(2, n - i);
    int num_group  = rot_factor;
    for(int j = 0 ; j < num_btfly ; ++j){
      for(int k = 0 ; k < num_group ; ++k){
	int top_index      = 2 * num_btfly * k + j;
	int bottom_index   = 2 * num_btfly * k + j + num_btfly;
	aft_btfly[i_next][top_index]    = aft_btfly[i_prev][top_index] + W[rot_factor * j] * aft_btfly[i_prev][bottom_index];
	aft_btfly[i_next][bottom_index] = aft_btfly[i_prev][top_index] - W[rot_factor * j] * aft_btfly[i_prev][bottom_index];
      }// k
    }// j
  }// i : 1 -> n
  
  // ------ COPY ------ //
  for(auto &&r : aft_btfly[n % 2]){
    *out_first_Re = r.real() / n_sample;
    *out_first_Im = r.imag() / n_sample;
    ++out_first_Re;
    ++out_first_Im;
  }
  
  return;
}
#endif
