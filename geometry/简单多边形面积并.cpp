#include<stdio.h>
#include<math.h>
#include<algorithm>
#include<vector>
#define PI 3.14159265358979323846
#define eps 1e-8
#define SGN(x) ((x)>eps?1:((x)>-eps?0:-1))
using namespace std;

struct pt
{
	double x, y;
	pt(){}
	pt(double _x, double _y):x(_x), y(_y){}
	pt operator - (const pt p1){return pt(x - p1.x, y - p1.y);}
	pt operator + (const pt p1){return pt(x + p1.x, y + p1.y);}
	pt operator * (double s){return pt(x * s, y * s);}
	pt operator / (double s){return pt(x / s, y / s);}
	bool operator < (const pt p1)const{return y < p1.y-eps || y < p1.y+eps && x < p1.x;}
};

inline double cpr(const pt &a,const pt &b,const pt &c){return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);}
inline double cpr(const pt &a,const pt &b){return a.x*b.y-a.y*b.x;}
inline double dpr(const pt &a,const pt &b,const pt &c){return (b.x-a.x)*(c.x-a.x)+(b.y-a.y)*(c.y-a.y);}
inline double dpr(const pt &a,const pt &b){return a.x*b.x+a.y*b.y;}

pt its(const pt &a, const pt &b, const pt &c, const pt &d)
{
	pt ret = a;
	double t =  ((c.x - a.x)*(d.y - c.y) - (c.y - a.y)*(d.x - c.x))/
				((b.x - a.x)*(d.y - c.y) - (b.y - a.y)*(d.x - c.x));
	ret.x += (b.x - a.x) * t;
	ret.y += (b.y - a.y) * t;
	return ret;
}

//////////////////////////////////////////////////////////

pair<double, int> e[510];
int cnt;

inline void insert(pt &s, pt &t, pt X, int inc)
{
	double ratio = SGN(t.x - s.x) ? (X.x - s.x) / (t.x - s.x) : (X.y - s.y) / (t.y - s.y);
	if (ratio > 1.0)ratio = 1.0;
	if (ratio < 0.0)ratio = 0.0;
	e[cnt++] = make_pair(ratio, inc);
}

double poly_union(vector<vector <pt> > &p)
{
	double ans = 0.0;
	int cp0, cp1, cp2, cp3;
	
	for (int i = 0; i < p.size(); i++)
	{
		for (int k = 0; k < p[i].size(); k++)
		{
			pt &s = p[i][k], &t = p[i][(k + 1) % p[i].size()];
			if (fabs(cpr(s, t)) < eps)continue;
			cnt = 0;
			e[cnt++] = make_pair(0.0, 1);
			e[cnt++] = make_pair(1.0, -1);
			for (int j = 0; j < p.size(); j++) if (i != j)
			{
				for (int l = 0; l < p[j].size(); l++)
				{
					pt &a = p[j][l], &b = p[j][(l + 1) % p[j].size()];
					cp0 = SGN(cpr(s, t, p[j][(l + p[j].size() - 1) % p[j].size()]));
					cp1 = SGN(cpr(s, t, a));
					cp2 = SGN(cpr(s, t, b));
					if (cp1 * cp2 < 0)
						insert(s, t, its(s, t, a, b), -cp2);
					else if (!cp1 && cp0 * cp2 < 0)
						insert(s, t, a, -cp2);
					else if (!cp1 && !cp2)
					{
						cp3 = SGN(cpr(s, t, p[j][(l + 2) % p[j].size()]));
						int dp = SGN(dpr(t - s, b - a));
						if (dp && cp0)insert(s, t, a, dp > 0 ? cp0 * (j > i ^ cp0 < 0) : -(cp0 < 0));
						if (dp && cp3)insert(s, t, b, dp > 0 ? -cp3 * (j > i ^ cp3 < 0) : cp3 < 0);
					}
				}
			}
			sort(e, e + cnt);
			int acc = 0;
			double total = 0.0, last;
			for (int j = 0; j < cnt; j++)
			{
				if (acc == 1)
					total += e[j].first - last;
				acc += e[j].second;
				last = e[j].first;
			}
			ans += cpr(s, t) * total;
		}
	}
	return fabs(ans) * 0.5;
}

int main()
{
	return 0;
}
