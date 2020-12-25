# FFT
Fast Fourier Transform

How to Use
--

高速フーリエ変換のライブラリ。
サンプリング数が2の累乗じゃないと使えないです。

FFT 定義

```c++
template <class InputIterator, class OutputIterator>  
void FFT(int n_sample, 
         InputIterator  data_first, 
         OutputIterator real_first, 
         OutputIterator imag_first);
```

フーリエ変換は入力関数を複素数関数に変換するので、出力には実部と虚部それぞれ関する配列を用意する必要があります。
入力に関する引数が1つしかないのは、実際に得られるデータ形式が大体の場合において実数の配列だからです。(というか、データが複素数として手に入ることってあるのでしょうか？？)
したがって、引数の意味は順番に 「サンプリング数(need to be power of 2)」, 「入力データ配列のイテレータ(begin)」, 「出力データ配列(実部)のイテレータ(begin)」, 「出力データ配列(虚部)のイテレータ(begin)」
です。出力用の配列は事前に n_sample 分だけ確保しておく必要があります。


InverseFFT 定義

```c++
template <class InputIterator, class OutputIterator>  
void InverseFFT(int n_sample, 
                InputIterator in_first_Re, 
                InputIterator in_first_Im, 
                OutputIterator out_first_Re, 
                OutputIterator out_first_Im);
```

逆フーリエ変換用の関数です。
引数の意味は順番に 「サンプリング数(need to be power of 2)」, 「入力データ配列(実部)のイテレータ(begin)」,「入力データ配列(虚部)のイテレータ(begin)」, 「出力データ配列(実部)のイテレータ(begin)」, 「出力データ配列(虚部)のイテレータ(begin)」。
FFTの時とは違い、入力には実部と虚部それぞれに関する配列が必要となります。これは、フーリエ変換を行なった後のデータを解析などで変更を加えたあと、また逆変換で元に戻すことを想定しているからです。例えば、パワースペクトルで高周波成分のノイズを除去した後、逆変換で除去後の波形を得る、など。
また、FFTの時と同様出力用の配列は事前に n_sample 分だけ確保しておく必要があります。

多分色々と高速にできるところはあると思います。
付け加えるとしたら、配列を複素数型で管理している時用に別途オーバーロードしたものを作っておけば便利かも。
あとはクラスで管理した方が使いやすいかもしれない。。
