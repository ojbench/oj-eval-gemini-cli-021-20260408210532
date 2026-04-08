#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int INF = 1e9 + 7;

struct Point {
    int x, y;
};

int N;
vector<Point> pts;

// Fenwick tree for prefix max
int tree_pmax[200005];
int version_pmax[200005];
int cur_ver_pmax = 0;

void update_pmax(int idx, int val) {
    for (; idx <= N; idx += idx & -idx) {
        if (version_pmax[idx] != cur_ver_pmax) {
            tree_pmax[idx] = -INF;
            version_pmax[idx] = cur_ver_pmax;
        }
        tree_pmax[idx] = max(tree_pmax[idx], val);
    }
}

int query_pmax(int idx) {
    int res = -INF;
    for (; idx > 0; idx -= idx & -idx) {
        if (version_pmax[idx] == cur_ver_pmax) {
            res = max(res, tree_pmax[idx]);
        }
    }
    return res;
}

// Fenwick tree for suffix min (implemented as prefix min on N - idx + 1)
int tree_smin[200005];
int version_smin[200005];
int cur_ver_smin = 0;

void update_smin(int idx, int val) {
    idx = N - idx + 1;
    for (; idx <= N; idx += idx & -idx) {
        if (version_smin[idx] != cur_ver_smin) {
            tree_smin[idx] = INF;
            version_smin[idx] = cur_ver_smin;
        }
        tree_smin[idx] = min(tree_smin[idx], val);
    }
}

int query_smin(int idx) {
    idx = N - idx + 1;
    int res = INF;
    for (; idx > 0; idx -= idx & -idx) {
        if (version_smin[idx] == cur_ver_smin) {
            res = min(res, tree_smin[idx]);
        }
    }
    return res;
}

// Fenwick tree for sum
int tree_sum[200005];
int version_sum[200005];
int cur_ver_sum = 0;

void update_sum(int idx, int val) {
    for (; idx <= N; idx += idx & -idx) {
        if (version_sum[idx] != cur_ver_sum) {
            tree_sum[idx] = 0;
            version_sum[idx] = cur_ver_sum;
        }
        tree_sum[idx] += val;
    }
}

int query_sum(int idx) {
    int res = 0;
    for (; idx > 0; idx -= idx & -idx) {
        if (version_sum[idx] == cur_ver_sum) {
            res += tree_sum[idx];
        }
    }
    return res;
}

struct L_Point {
    int y, y_max;
};

struct R_Point {
    int y, y_min;
};

long long solve(int left, int right) {
    if (left >= right) return 0;
    int mid = left + (right - left) / 2;
    long long ans = solve(left, mid) + solve(mid + 1, right);

    cur_ver_smin++;
    vector<L_Point> L_arr;
    L_arr.reserve(mid - left + 1);
    for (int i = mid; i >= left; --i) {
        int y = pts[i].y;
        int ym = query_smin(y + 1);
        L_arr.push_back({y, ym});
        update_smin(y, y);
    }

    cur_ver_pmax++;
    vector<R_Point> R_arr;
    R_arr.reserve(right - mid);
    for (int i = mid + 1; i <= right; ++i) {
        int y = pts[i].y;
        int ym = query_pmax(y - 1);
        R_arr.push_back({y, ym});
        update_pmax(y, y);
    }

    sort(L_arr.begin(), L_arr.end(), [](const L_Point& a, const L_Point& b) {
        return a.y_max > b.y_max;
    });
    sort(R_arr.begin(), R_arr.end(), [](const R_Point& a, const R_Point& b) {
        return a.y > b.y;
    });

    cur_ver_sum++;
    int ptr = 0;
    for (const auto& B : R_arr) {
        while (ptr < L_arr.size() && L_arr[ptr].y_max > B.y) {
            update_sum(L_arr[ptr].y, 1);
            ptr++;
        }
        int q_high = B.y - 1;
        int q_low = (B.y_min == -INF) ? 0 : B.y_min;
        if (q_high >= q_low) {
            ans += query_sum(q_high) - query_sum(q_low);
        }
    }

    return ans;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    if (!(cin >> N)) return 0;
    pts.resize(N);
    vector<int> ys(N);
    for (int i = 0; i < N; ++i) {
        cin >> pts[i].x >> pts[i].y;
        ys[i] = pts[i].y;
    }

    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());

    for (int i = 0; i < N; ++i) {
        pts[i].y = lower_bound(ys.begin(), ys.end(), pts[i].y) - ys.begin() + 1;
    }

    sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
        return a.x < b.x;
    });

    cout << solve(0, N - 1) << "\n";

    return 0;
}
