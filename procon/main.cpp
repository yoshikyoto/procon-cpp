#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include <stack>
#include <queue>
#include <set>
#include <cstdio>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <map>
#include <complex>
using namespace std;

#define REP(i, n) for(int i = 0; i < n; i++)
#define RREP(i, n) for(int i=(n)-1;i>=0;i--)
#define FOR(i, b, e) for(int i = b; i < e; i++)
#define to_bit(i) static_cast< bitset<8> >(i)
#define INF (1<<28)
#define EPS 1e-9





typedef long long ll;
typedef unsigned long long ull;
typedef vector<int> VI;
typedef vector<string> VS;
typedef pair<int, int> PII;
typedef pair<long long, long long>PLL;
typedef queue<int> QI;
typedef priority_queue<int> maxpq;
typedef priority_queue<int, vector<int>, greater<int> > minpq;

struct edge{int to, cost;};

int gcd(int a, int b){if(a%b==0){return(b);}else{return(gcd(b,a%b));}};
int lcm(int m, int n){if((0 == m)||(0 == n)){return 0;} return ((m / gcd(m, n)) * n);};
unsigned long long comb(long n, long m){unsigned long long c = 1; m = (n - m < m ? n - m : m);
    for(long ns = n - m + 1, ms = 1; ms <= m; ns ++, ms ++){c *= ns; c /= ms;} return c;};

void cp(int a1[], int a2[], int l){REP(i, l) a2[i] = a1[i];};
void cp(string a1[], string a2[], int l){REP(i, l) a2[i] = a1[i];};
double sq(double d){return d*d;};
int sq(int i){return i*i;};
double sqdist(int x1, int y1, int x2, int y2){
    double dx = x1 - x2, dy = y1 - y2; return dx*dx + dy*dy;
};
bool inside(int y, int x, int h, int w){return 0 <= y && y <= (h-1) && 0 <= x && x <= (w-1);};


// 線分の交差判定
bool isCross(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4){
    // 並行な場合
    int m = (x2-x1)*(y4-y3)-(y2-y1)*(x4-x3);
    if(m == 0) return false;
    // 媒介変数の値が0より大きく1以下なら交差する、これは問題の状況によって変わるかも。
    double r =(double)((y4-y3)*(x3-x1)-(x4-x3)*(y3-y1))/m;
    double s =(double)((y2-y1)*(x3-x1)-(x2-x1)*(y3-y1))/m;
    return (0 < r && r <= 1 && 0 < s && s <= 1);
};

// 外積の計算 AB CD の内積を求める
int crossProdct(int ax, int ay, int bx, int by, int cx, int cy, int dx, int dy){
    int abx = bx - ax; int aby = by - ay;
    int cdx = dx - cx; int cdy = dy - cy;
    return abx * cdy - cdx * aby;
};
double crossProdct(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy){
    double abx = bx - ax; double aby = by - ay;
    double cdx = dx - cx; double cdy = dy - cy;
    return abx * cdy - cdx * aby;
};
double innerProduct(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy){
    double abx = bx - ax; double aby = by - ay;
    double cdx = dx - cx; double cdy = dy - cy;
    return abx * cdx + aby * cdy;
};

// 三角形の内部判定 ABCの中にPがあるか判定
bool inTriangle(int ax, int ay, int bx, int by, int cx, int cy, int px, int py){
    int c1 = crossProdct(ax, ay, bx, by, bx, by, px, py);
    int c2 = crossProdct(bx, by, cx, cy, cx, cy, px, py);
    int c3 = crossProdct(cx, cy, ax, ay, ax, ay, px, py);
    return (c1 > 0 && c2 > 0 && c3 > 0) || (c1 < 0 && c2 < 0 && c3 < 0);
};
bool inTriangle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py){
    double c1 = crossProdct(ax, ay, bx, by, bx, by, px, py);
    double c2 = crossProdct(bx, by, cx, cy, cx, cy, px, py);
    double c3 = crossProdct(cx, cy, ax, ay, ax, ay, px, py);
    return (c1 > 0 && c2 > 0 && c3 > 0) || (c1 < 0 && c2 < 0 && c3 < 0);
};

// 三角形の外心
void circumcenter(double ax, double ay, double bx, double by, double cx, double cy, double res[3]){
    // AB AC を求める
    double abx = bx - ax; double aby = by - ay;
    double acx = cx - ax; double acy = cy - ay;
    double m = abx * acy - acx * aby; // 分母 m = 0 となるのは3点が一直線になるとき
    double s = (abx * acx + aby * acy - sq(acx) - sq(acy)) / m; // 媒介変数s
    res[0] = abx / 2 + s * aby / 2; res[1] = aby / 2 - s * abx / 2;
    // cout << res[0] << " " << res[1] << endl;
    res[2] = sqrt(sqdist(0, 0, res[0], res[1]));
    res[0] += ax; res[1] += ay;
};



class BinaryIndexedTree{
public:
    int n; vector<int> bit;
    BinaryIndexedTree(int _n){
        n = _n;
        bit.resize(n+1,0);
    }
    int sum(int i){
        int sum = 0;
        while(i > 0){
            sum += bit[i];
            i -= i & -i;
        }
        return sum;
    }
    void add(int i, int v){
        while(i <= n){
            bit[i] += v;
            i += i & -i;
        }
    }
};

class UnionFindTree{
public:
    vector<int> parent, rank;
    int n;
    std::set<int> set;
    // 初期化
    UnionFindTree(int _n){
        n = _n;
        for(int i = 0; i < n; i++){
            parent.resize(n);
            rank.resize(n);
            parent[i] = i;
            rank[i] = 0;
        }
    }
    // 根を求める
    int find(int x){
        if(parent[x] == x){
            return x;
        }else{
            return parent[x] = find(parent[x]);
        }
    }
    // xとyの集合を結合
    void unite(int x, int y){
        x = find(x);
        y = find(y);
        if(x == y){
            set.insert(x);
            return;
        }
        if(rank[x] < rank[y]){
            parent[x] = y;
            set.erase(x);
            set.insert(y);
        }else{
            parent[y] = x;
            set.erase(y);
            set.insert(x);
            if(rank[x] == rank[y]) rank[x]++;
        }
    }
    // xとyが同じ集合か
    bool same(int x, int y){
        return find(x) == find(y);
    }
    // 集合の数を数える
    int count(){
        return (int)set.size();
    }
};

// ワーシャルフロイド法
void warshallFloyd(int graph[100][100], int graph_size){
    for(int mid_node = 0; mid_node < graph_size; mid_node++)
        for(int s_node = 0; s_node < graph_size; s_node++)
            for(int g_node = 0; g_node < graph_size; g_node++)
                if(s_node == g_node) graph[s_node][g_node] = 0;
                else graph[s_node][g_node] = min(graph[s_node][g_node], graph[s_node][mid_node] + graph[mid_node][g_node]);
};


// d:隣接行列  n:グラフのサイズ  s:始点  dist:始点からの距離が入る配列
void dijkstra(int graph[1000][1000], int node_count, int start_node, int distances[1000]){
    REP(i, node_count) distances[i] = -1;
    distances[start_node] = 0;
    priority_queue<PII, vector<PII>, greater<PII> > dijkstra_pq;
    dijkstra_pq.push(PII(0, start_node));
    while (!dijkstra_pq.empty()) {
        PII p = dijkstra_pq.top(); dijkstra_pq.pop();
        int i = p.second;
        if(distances[i] < p.first) continue;
        for(int j = 0; j < node_count; j++){
            if(graph[i][j] == -1) continue;
            if(distances[j] == -1){
                distances[j] = distances[i] + graph[i][j];
                dijkstra_pq.push(PII(distances[j], j));
            }else if(distances[j] > distances[i] + graph[i][j]){
                distances[j] = distances[i] + graph[i][j];
                dijkstra_pq.push(PII(distances[j], j));
            }
        }
    }
};


// return とかの位置によってうまい具合にする
int vi[4] = {-1, 1, 0, 0}, vj[4] = {0, 0, -1, 1};
void dfs2d(int i, int j, int r, int c){
    // 終了条件とかをここに書く
    REP(k, 4){
        int ni = i + vi[k];
        int nj = j + vj[k];
        if(inside(ni, nj, r, c)){
            dfs2d(ni, nj, r, c);
        }
    }
};


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