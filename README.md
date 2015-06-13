競技プログラミングで使えるアルゴリズム C++編
====

## 数学？

### 高速フーリエ変換

* 参考1: http://atc001.contest.atcoder.jp/tasks/fft_c
* 参考2: http://www.prefield.com/algorithm/math/fft.html

```cpp
// 高速フーリエ変換
// inverse = 1 でフーリエ変換、inverse = -1 で逆フーリエ変換
const double PI = 4.0*atan(1.0);
const complex<double> I(0,1);
void fft(complex<double> a[], int n, int inverse) {
    double theta = 2 * inverse * PI / n;
    
    for (int m = n; m >= 2; m >>= 1) {
        int mh = m >> 1;
        for (int i = 0; i < mh; i++) {
            complex<double> w = exp(i*theta*I);
            for (int j = i; j < n; j += m) {
                int k = j + mh;
                complex<double> x = a[j] - a[k];
                a[j] += a[k];
                a[k] = w * x;
            }
        }
        theta *= 2;
    }
    int i = 0;
    for (int j = 1; j < n - 1; j++) {
        for (int k = n >> 1; k > (i ^= k); k >>= 1);
        if (j < i) swap(a[i], a[j]);
    }
    
    if(inverse == -1){
        complex<double> d(n,0);
        REP(i,n){
            a[i] = a[i] / d;
        }
    }
}
```

```cpp
int main(int argc, const char * argv[]){
    int n;
    cin >> n;
    int g[100000], h[100000];
    REP(i,n){
        cin >> g[i] >> h[i];
    }
    
    int nn = 1;
    while(nn <= n + n - 1){
        nn *= 2;
    }
    
    complex<double> *gg = new complex<double>[nn];
    complex<double> *hh = new complex<double>[nn];
    REP(i,nn){
        if(i < n){
            gg[i] = g[i];
            hh[i] = h[i];
        }else{
            gg[i] = 0;
            hh[i] = 0;
        }
    }
    fft(gg, nn, 1);
    fft(hh, nn, 1);
    
    complex<double> *ff = new complex<double>[nn];
    REP(i,nn){
        ff[i] = gg[i] * hh[i];
    }
    fft(ff, nn, -1);
    cout << 0 << endl;
    REP(i, 2*n-1){
        cout << (ff[i].real()) << endl;
    }
}
```	
