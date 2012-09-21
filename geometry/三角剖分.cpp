/*
给出平面上的一个点集，求该点集的三角剖分。

此程序是利用平面三角剖分求平面欧几里德最小生成树。
可以证明，平面欧几里德最小生成树的树边一定在平面
三角剖分的边中。经过 LiZhiXu 的精心简化，已经达
到能在现场赛中输入的程度了

——Felicia
*/

#include <iostream>
#include <cmath>
using namespace std;

#define Oi(e) ((e)->oi)
#define Dt(e) ((e)->dt)
#define On(e) ((e)->on)
#define Op(e) ((e)->op)
#define Dn(e) ((e)->dn)
#define Dp(e) ((e)->dp)
#define Other(e, p) ((e)->oi == p ? (e)->dt : (e)->oi)
#define Next(e, p) ((e)->oi == p ? (e)->on : (e)->dn)
#define Prev(e, p) ((e)->oi == p ? (e)->op : (e)->dp)
#define V(p1, p2, u, v) (u = p2->x - p1->x, v = p2->y - p1->y)
#define C2(u1, v1, u2, v2) (u1 * v2 - v1 * u2)
#define C3(p1, p2, p3) ((p2->x - p1->x) * (p3->y - p1->y) - (p2->y - p1->y) * (p3->x - p1->x))
#define Dot(u1, v1, u2, v2) (u1 * u2 + v1 * v2)

#define MAXN 100001

struct point {
    int x, y;
    struct edge *in;
    bool operator < (const point &p1) const {
        return x < p1.x || (x == p1.x && y < p1.y);
    }
};

struct edge {
    point *oi, *dt;
    edge *on, *op, *dn, *dp;
};

struct gEdge {
    int u, v;
    double w;
    bool operator < (const gEdge &e1) const {return w < e1.w;}
}E[3 * MAXN], MST[MAXN];

int N, M;
int f[MAXN], d[MAXN];

void Init() {
    for(int i = 0; i < N; ++i) f[i] = i, d[i] = 0;
}

int Find(int x) {
    int i = x, t = x;
    while(f[i] != i) i = f[i];
    while(f[t] != i) {
        int temp = f[t]; f[t] = i; t = temp;
    }
    return i;
}

void Make(int x, int y) {
    x = Find(x); y = Find(y);
    if(d[x] > d[y]) f[y] = x;
    else {
        f[x] = y;
        if(d[x] == d[y]) d[y]++;
    }
}

void Kruskal() {
    double length = 0.0;
    Init();
    sort(E, E + M);
    for(int i = 0, k = 0; i < M && k < N - 1; ++i) {
        if(Find(E[i].u) != Find(E[i].v)) {
            Make(E[i].u, E[i].v); MST[k++] = E[i]; length += E[i].w;
        }
    }
    printf("%.4lf\n", length);
}

point p[MAXN], *Q[MAXN];
edge mem[3 * MAXN], *elist[3 * MAXN];
int nfree;

void Alloc_memory() {
    nfree = 3 * N;
    edge *e = mem;
    for(int i = 0; i < nfree; ++i) elist[i] = e++;
}

void Splice(edge *a, edge *b, point *v) {
    edge *next;
    if(Oi(a) == v) next = On(a), On(a) = b;
    else next = Dn(a), Dn(a) = b;
    if(Oi(next) == v) Op(next) = b;
    else Dp(next) = b;
    if(Oi(b) == v) On(b) = next, Op(b) = a;
    else Dn(b) = next, Dp(b) = a;
}

edge *Make_edge(point *u, point *v) {
    edge *e = elist[--nfree];
    e->on = e->op = e->dn = e->dp = e; e->oi = u; e->dt = v;
    if(u->in == NULL) u->in = e; if(v->in == NULL) v->in = e;
    return e;
}

edge *Join(edge *a, point *u, edge *b, point *v, int side) {
    edge *e = Make_edge(u, v);
    if(side == 1) {
        if(Oi(a) == u) Splice(Op(a), e, u);
        else Splice(Dp(a), e, u);
        Splice(b, e, v);
    }
    else {
        Splice(a, e, u);
        if(Oi(b) == v) Splice(Op(b), e, v);
        else Splice(Dp(b), e, v);
    }
    return e;
}

void Remove(edge *e) {
    point *u = Oi(e), *v = Dt(e);
    if(u->in == e) u->in = e->on; if(v->in == e) v->in = e->dn;
    if(Oi(e->on) == u) e->on->op = e->op;
    else e->on->dp = e->op;
    if(Oi(e->op) == u) e->op->on = e->on;
    else e->op->dn = e->on;
    if(Oi(e->dn) == v) e->dn->op = e->dp;
    else e->dn->dp = e->dp;
    if(Oi(e->dp) == v) e->dp->on = e->dn;
    else e->dp->dn = e->dn;
    elist[nfree++] = e;
}

void Make_Graph() {
    for(int i = 0; i < N; ++i) {
        point *u = &p[i];
        edge *start = u->in, *e = u->in;
        do {
            point *v = Other(e, u);
            if(u < v) {
                E[M].u = u - p, E[M].v = v - p;
                E[M++].w = hypot(u->x - v->x, u->y - v->y);
            }
            e = Next(e, u);
        } while(e != start);
    }
}

void Low_tangent(edge *e_l, point *o_l, edge *e_r, point *o_r, edge **l_low, point **OL, edge **r_low, point **OR) {
    point *d_l = Other(e_l, o_l), *d_r = Other(e_r, o_r);
    while(true) {
        if(C3(o_l, o_r, d_l) < 0.0) {
            e_l = Prev(e_l, d_l);
            o_l = d_l; d_l = Other(e_l, o_l);
        }
        else if(C3(o_l, o_r, d_r) < 0.0) {
            e_r = Next(e_r, d_r);
            o_r = d_r; d_r = Other(e_r, o_r);
        }
        else break;
    }
    *OL = o_l, *OR = o_r;
    *l_low = e_l, *r_low = e_r;
}

void Merge(edge *lr, point *s, edge *rl, point *u, edge **tangent) {
    double l1, l2, l3, l4, r1, r2, r3, r4, cot_L, cot_R, u1, v1, u2, v2, N1, cot_N, P1, cot_P;
    point *O, *D, *OR, *OL;
    edge *B, *L, *R;
    Low_tangent(lr, s, rl, u, &L, &OL, &R, &OR);
    *tangent = B = Join(L, OL, R, OR, 0);
    O = OL, D = OR;
    do {
        edge *El = Next(B, O), *Er = Prev(B, D), *next, *prev;
        point *l = Other(El, O), *r = Other(Er, D);
        V(l, O, l1, l2); V(l, D, l3, l4); V(r, O, r1, r2); V(r, D, r3, r4);
        double cl = C2(l1, l2, l3, l4), cr = C2(r1, r2, r3, r4);
        bool BL = cl > 0.0, BR = cr > 0.0;
        if(!BL && !BR) break;
        if(BL) {
            double dl = Dot(l1, l2, l3, l4);
            cot_L = dl / cl;
            do {
                next = Next(El, O);
                V(Other(next, O), O, u1, v1); V(Other(next, O), D, u2, v2);
                N1 = C2(u1, v1, u2, v2);
                if(!(N1 > 0.0)) break;
                cot_N = Dot(u1, v1, u2, v2) / N1;
                if(cot_N > cot_L) break;
                Remove(El);
                El = next;
                cot_L = cot_N;
            } while(true);
        }
        if(BR) {
            double dr = Dot(r1, r2, r3, r4);
            cot_R = dr / cr;
            do {
                prev = Prev(Er, D);
                V(Other(prev, D), O, u1, v1); V(Other(prev, D), D, u2, v2);
                P1 = C2(u1, v1, u2, v2);
                if(!(P1 > 0.0)) break;
                cot_P = Dot(u1, v1, u2, v2) / P1;
                if(cot_P > cot_R) break;
                Remove(Er);
                Er = prev;
                cot_R = cot_P;
            } while(true);
        }
        l = Other(El, O); r = Other(Er, D);
        if(!BL || (BL && BR && cot_R < cot_L)) { B = Join(B, O, Er, r, 0); D = r; }
        else { B = Join(El, l, B, D, 0); O = l; }
    } while(true);
}

void Divide(int s, int t, edge **L, edge **R) {
    edge *a, *b, *c, *ll, *lr, *rl, *rr, *tangent;
    int n = t - s + 1;
    if(n == 2) *L = *R = Make_edge(Q[s], Q[t]);
    else if(n == 3) {
        a = Make_edge(Q[s], Q[s + 1]), b = Make_edge(Q[s + 1], Q[t]);
        Splice(a, b, Q[s + 1]);
        double v = C3(Q[s], Q[s + 1], Q[t]);
        if(v > 0.0) {
            c = Join(a, Q[s], b, Q[t], 0);
            *L = a; *R = b;
        }
        else if(v < 0.0) {
            c = Join(a, Q[s], b, Q[t], 1);
            *L = c; *R = c;
        }
        else { *L = a; *R = b; }
    }
    else if(n > 3) {
        int split = (s + t) / 2;
        Divide(s, split, &ll, &lr); Divide(split + 1, t, &rl, &rr);
        Merge(lr, Q[split], rl, Q[split + 1], &tangent);
        if(Oi(tangent) == Q[s]) ll = tangent;
        if(Dt(tangent) == Q[t]) rr = tangent;
        *L = ll; *R = rr;
    }
}

int main() {
    while(scanf("%d", &N) != EOF) {
        if(N == 1) {
            printf("0.0000\n");
            continue;
        }
        Alloc_memory();
        for(int i = 0; i < N; ++i) {
            scanf("%d %d", &p[i].x, &p[i].y);
            p[i].in = NULL;
        }
        sort(p, p + N);
        for(int i = 0; i < N; i++) Q[i] = p + i;
        edge *L, *R;
        Divide(0, N - 1, &L, &R);
        M = 0;
        Make_Graph();
        Kruskal();
    }
    return 0;
}
