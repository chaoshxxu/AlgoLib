#include<stdio.h>
#include<math.h>
#include<algorithm>
#define PI 3.14159265358979323846
#define eps 1e-8
using namespace std;

inline int SGN(double x){return x < -eps ? -1 : x < eps ? 0 : 1;}

struct pt
{
	double x, y;
	pt(){}
	pt(double _x, double _y):x(_x), y(_y){}
	pt operator - (const pt p1){return pt(x - p1.x, y - p1.y);}
	pt operator + (const pt p1){return pt(x + p1.x, y + p1.y);}
	pt operator * (double r){return pt(x * r, y * r);}
	pt operator / (double r){return pt(x / r, y / r);}
	void read(){scanf("%lf %lf", &x, &y);}
};

inline double cpr(const pt &a,const pt &b){return a.x*b.y-a.y*b.x;}
inline double dis(const pt &a, const pt &b){return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));}

///////////////////////////////////////////////////

int n;
pt o[1010];
double r[1010];

pair<double, int> e[2010];
int cnt;

double ans[1010];

//0:a包含b	1:b包含a	2:相交	3:相离
inline int rlt(int a, int b)
{
	double d = dis(o[a], o[b]), d1 = SGN(d - r[a] + r[b]), d2 = SGN(d - r[b] + r[a]);
	if (d1 < 0 || !d1 && (d > eps || a > b))return 0;
	if (d2 < 0 || !d2 && (d > eps || a < b))return 1;
	return d < r[a] + r[b] - eps ? 2 : 3;
}

inline double arcArea(pt &o, double r, double ang1, double ang2)
{
	pt a(o.x + r * cos(ang1), o.y + r * sin(ang1));
	pt b(o.x + r * cos(ang2), o.y + r * sin(ang2));
	double dif = ang2 - ang1;
	return (cpr(a, b) + (dif - sin(dif)) * r * r) * 0.5;
}

void solve()
{
	double last, center, d2, ang, angX, angY;
	pt X, Y;
	
	for (int i = 0; i < n; i++) if (r[i] > eps)
	{
		int acc = 0;
		cnt = 0;
		e[cnt++] = make_pair(-PI, 1);
		e[cnt++] = make_pair(PI, -1);
		for (int j = 0; j < n; j++) if (i != j && r[j] > eps)
		{
			int rel = rlt(i, j);
			if (rel == 1)
			{
				e[cnt++] = make_pair(-PI, 1);
				e[cnt++] = make_pair(PI, -1);
			} else if (rel == 2)
			{
				center = atan2(o[j].y - o[i].y, o[j].x - o[i].x);
				d2 = (o[i].x - o[j].x) * (o[i].x - o[j].x) + (o[i].y - o[j].y) * (o[i].y - o[j].y);
				ang = acos((r[i] * r[i] + d2 - r[j] * r[j]) / (2 * r[i] * sqrt(d2)));
				angX = center + ang;
				angY = center - ang;
				if (angX > PI)angX -= 2*PI;
				if (angY < -PI)angY += 2*PI;
				if (angX < angY) acc++;
				e[cnt++] = make_pair(angY, 1);
				e[cnt++] = make_pair(angX, -1);
			}
		}
		sort(e, e + cnt);
		last = -PI;
		for (int j = 0; j < cnt; j++)
		{
			double tmp = arcArea(o[i], r[i], last, e[j].first);
			ans[acc] += tmp;
			ans[acc - 1] -= tmp;
			acc += e[j].second;
			last = e[j].first;
		}
	}
}

int main()
{
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
	{
		o[i].read();
		scanf("%lf", &r[i]);
	}
	solve();
	for (int i = 1; i <= n; i++)
	{
		printf("[%d] = %.3f\n", i, ans[i]);
	}
	return 0;
}
