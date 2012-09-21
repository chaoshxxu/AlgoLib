#include<stdio.h>
#include<math.h>
#include<algorithm>
#define eps 1e-8
#define SGN(x) ((x)>eps?1:((x)>-eps?0:-1))
using namespace std;

struct pt
{
	double x, y;
	pt(){}
	pt(double _x, double _y):x(_x), y(_y){}
	pt operator - (const pt p1){return pt(x - p1.x, y - p1.y);}
	void read(){scanf("%lf %lf", &x, &y);}
};

double cpr(const pt &a,const pt &b,const pt &c){return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);}
double cpr(const pt &a,const pt &b){return a.x*b.y-a.y*b.x;}
double dpr(const pt &a,const pt &b){return a.x*b.x+a.y*b.y;}

pt its(const pt &a, const pt &b, const pt &c, const pt &d)
{
	pt ret = a;
	double t = ((c.x - a.x)*(d.y - c.y) - (c.y - a.y)*(d.x - c.x)) / ((b.x - a.x)*(d.y - c.y) - (b.y - a.y)*(d.x - c.x));
	ret.x += (b.x - a.x) * t;
	ret.y += (b.y - a.y) * t;
	return ret;
}

////////////////////////////////////////////////

int n;
pt p[110][3];

pair<double, int> e[210];
int cnt;

inline void insert(pt &s, pt &t, pt X, int inc)
{
	double ratio = SGN(t.x - s.x) ? (X.x - s.x) / (t.x - s.x) : (X.y - s.y) / (t.y - s.y);
	e[cnt++] = make_pair(ratio, inc);
}

double solve()
{
	double ans = 0.0;
	int cp0, cp1, cp2;
	
	for (int i = 0; i < n; i++) if (SGN(cpr(p[i][0], p[i][1], p[i][2])))
	{
		for (int k = 0; k < 3; k++)
		{
			pt &s = p[i][k], &t = p[i][k==2?0:k+1]; 
			cnt = 0;
			e[cnt++] = make_pair(0.0, 1);
			e[cnt++] = make_pair(1.0, -1);
			for (int j = 0; j < n; j++) if (i != j && SGN(cpr(p[j][0], p[j][1], p[j][2])))
			{
				for (int l = 0; l < 3; l++)
				{
					cp0 = SGN(cpr(s, t, p[j][l==0?2:l-1]));
					cp1 = SGN(cpr(s, t, p[j][l]));
					cp2 = SGN(cpr(s, t, p[j][l==2?0:l+1]));
					if (cp1 * cp2 < 0)
						insert(s, t, its(s, t, p[j][l], p[j][l==2?0:l+1]), cp2);
					else if (!cp1 && cp0 * cp2 < 0)
						insert(s, t, p[j][l], cp2);
					else if (!cp1 && !cp2 && j > i && dpr(t - s, p[j][l==2?0:l+1] - p[j][l]) > -eps)
					{
						insert(s, t, p[j][l], -1);
						insert(s, t, p[j][l==2?0:l+1], 1);
					}
				}
			}
			sort(e, e + cnt);
			int acc = 0;
			double total = 0.0, last;
			for (int j = 0; j < cnt; j++)
			{
				acc += e[j].second;
				if (acc == 0 && e[j].second < 0)
					total += e[j].first - last;
				last = e[j].first;
			}
			ans += cpr(s, t) * total;
		}
	}
	return ans * 0.5;
}

int main()
{
	int tc;
	scanf("%d", &tc);
	for (int t = 1; t <= tc; t++)
	{
		scanf("%d", &n);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 3; j++)
				p[i][j].read();
			if (cpr(p[i][0], p[i][1], p[i][2]) < 0)
				swap(p[i][1], p[i][2]);
		}
		printf("Case %d: %.3f\n", t, solve());
		
	}
	return 0;
}
