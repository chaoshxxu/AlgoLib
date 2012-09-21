#include<stdio.h>
#include<math.h>
#include<algorithm>
#define MAXV 2000
#define PI 3.14159265358979323846
#define eps 1e-8
#define zero(x) (fabs(x)<eps)
#define SGN(x) ((x)>eps?1:((x)>-eps?0:-1))
#define _sign(x) ((x)>eps?1:((x)<-eps?2:0))
using namespace std;

//��ά��
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
	bool operator == (const pt p1)const{return !SGN(y - p1.y) && !SGN(x - p1.x);}
	bool operator != (const pt p1)const{return SGN(y - p1.y) || SGN(x - p1.x);}
	void read(){scanf("%lf %lf", &x, &y);}
};


//��� ���(�����) 
double cpr(const pt &a,const pt &b,const pt &c){return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);}
double dpr(const pt &a,const pt &b,const pt &c){return (b.x-a.x)*(c.x-a.x)+(b.y-a.y)*(c.y-a.y);}

//��� ���(��������) 
double cpr(const pt &a,const pt &b){return a.x*b.y-a.y*b.x;}
double dpr(const pt &a,const pt &b){return a.x*b.x+a.y*b.y;}

//������� 
double vlen(const pt &a){return sqrt(a.x*a.x+a.y*a.y);}
double dis(const pt &a, const pt &b){return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));}

//�ж�ֱ��ab���߶�cd�ϸ��ཻ
bool over(const pt &a, const pt &b, const pt &c, const pt &d)
{
	double p, q;
	p = cpr(a, b, c);
	q = cpr(a, b, d);
	return p>eps && q<-eps || p<-eps && q>eps;
}

//�ж��߶�ab���߶�cd�ϸ��ཻ
bool cross(const pt &a, const pt &b, const pt &c, const pt &d)
{
	double p, q;
	p = cpr(a, b, c);
	q = cpr(a, b, d);
	if (!(p>eps && q<-eps || p<-eps && q>eps))
		return 0;
	p = cpr(c, d, a);
	q = cpr(c, d, b);
	return p>eps && q<-eps || p<-eps && q>eps;
}

//��ֱ��ab��ֱ��cd�Ľ���
pt its(pt a, pt b, pt c, pt d)
{
	double v1 = cpr(a, b, c), v2 = cpr(a, b, d);
	return (c * v2 - d * v1) / (v2 - v1);
}

//����������Ǻ�,��ʱ��Ϊ�� 
double angsum(pt p[], int n)
{
	double ret = 0, tmp;
	for (int i = 0; i < n; i++)
	{
		pt &A = p[i], &B = p[(i+1)%n], &C = p[(i+2)%n];
		tmp = PI - acos(dpr(B, A, C) / dis(A, B) / dis(B, C));
		if (cpr(A, B, C) < 0)
			tmp = -tmp;
		ret += tmp;
	}
	return ret;
}

//�ж�͹�����,���㰴˳ʱ�����ʱ�����,�������ڱ߹���(n>=3)
int is_convex(pt p[], int n)
{
	int s[3] = {1,1,1};
	for (int i = 0; i<n && s[1]|s[2]; i++)
		s[_sign(cpr(p[i], p[(i+1)%n], p[(i+2)%n]))] = 0;
	return s[1]|s[2];
}

//�ж�͹�����,���㰴˳ʱ�����ʱ�����,���������ڱ߹���
int is_convex_v2(pt p[], int n)
{
	int s[3] = {1,1,1};
	for (int i = 0; i<n && s[0] && s[1]|s[2]; i++)
		s[_sign(cpr(p[i], p[(i+1)%n], p[(i+2)%n]))] = 0;
	return s[0] && s[1]|s[2];
}

//�е���͹������ڻ����α���,���㰴˳ʱ�����ʱ�����
int inside_convex(pt q, pt p[], int n)
{
	int s[3] = {1,1,1};
	for (int i = 0; i<n && s[1]|s[2]; i++)
		s[_sign(cpr(p[i], p[(i+1)%n],q))] = 0;
	return s[1]|s[2];
}

//�е���͹�������,���㰴˳ʱ�����ʱ�����,�ڶ���α��Ϸ���0
int inside_convex_v2(pt q, pt p[], int n)
{
	int s[3] = {1,1,1};
	for (int i = 0; i<n && s[0] && s[1]|s[2]; i++)
		s[_sign(cpr(p[i], p[(i+1)%n],q))] = 0;
	return s[0] && s[1]|s[2];
}

//���������
pt barycenter(int n, pt *p)
{
	pt ret(0, 0), t;
	double t1 = 0, t2;
	for (int i = 1; i < n - 1; i++)
	{
		if (fabs(t2 = cpr(p[i+1], p[0], p[i])) > eps)
		{
			t.x = (p[0].x + p[i].x + p[i+1].x) /3.0;
			t.y = (p[0].y + p[i].y + p[i+1].y) /3.0;
			ret.x += t.x*t2;
			ret.y += t.y*t2;
			t1 += t2;
		}
	}
	if (fabs(t1) > eps)
		ret.x /= t1, ret.y /= t1;
	return ret;
}

//��ƽ�潻(�и�ֱ��Ϊab,pol[]Ϊԭʼ�����,������Ϊpolcnt 
//exch[]��ʱ������ŵ㼯,���������Է���pol[])
//����ε���Ϊ��ʱ��,�и������Ϊ��Ч�� 
void half_its(pt &a, pt &b, pt pol[], int &polcnt, pt exch[])
{
	int i, p2 = 0;
	bool now, last = cpr(a, b, pol[polcnt-1]) > -eps;
	for (i = 0; i < polcnt; i++)
	{
		now = cpr(a, b, pol[i]) > -eps;
		if (now ^ last)
			exch[p2++] = its(a, b, pol[i], pol[(i+polcnt-1)%polcnt]);				
		if (now)
			exch[p2++] = pol[i];
		last = now;
	}
	polcnt = p2;
	for (i = 0; i < p2; i++)
		pol[i] = exch[i];
}


//��p��ֱ��ab�ϵ������
pt ptoline(pt p, pt a, pt b)
{
	pt t = p;
	t.x += a.y - b.y;
	t.y += b.x - a.x;
	return its(p, t, a, b);
}

//��p��ֱ��ab���� 
double disptoline(const pt &p, const pt &a, const pt &b)
{
	return fabs(cpr(p, a, b)) / dis(a, b);
}

//��p���߶�ab�ϵ������
pt ptoseg(pt p, pt a, pt b)
{
	pt t = p;
	t.x += a.y - b.y;
	t.y += b.x - a.x;
	if (cpr(p,a,t) * cpr(p,b,t) > eps)
		return dis(p,a) < dis(p,b) ? a : b;
	return its(p, t, a, b);
}

//�㵽�߶ξ���
double disptoseg(pt p, pt a, pt b)
{
	pt t = p;
	t.x += a.y - b.y;
	t.y += b.x - a.x;
	if (cpr(p,a,t) * cpr(p,b,t) > eps)
		return min(dis(p,a), dis(p,b));
	return fabs(cpr(b,p,a)) / dis(a,b);
}

//��v���ŵ�p��ʱ����תangle���Ŵ�scale��
pt rotate(pt v, pt p, double angle, double scale)
{
	pt ret = p;
	v.x -= p.x;
	v.y -= p.y;
	p.x = scale * cos(angle);
	p.y = scale * sin(angle);
	ret.x += v.x * p.x - v.y * p.y;
	ret.y += v.x * p.y + v.y * p.x;
	return ret;
}

//����������
pt circumcenter(pt a, pt b, pt c)
{
	pt u1, u2, v1, v2;
	u1.x = (a.x+b.x)/2;
	u1.y = (a.y+b.y)/2;
	u2.x = u1.x-a.y+b.y;
	u2.y = u1.y+a.x-b.x;
	v1.x = (a.x+c.x)/2;
	v1.y = (a.y+c.y)/2;
	v2.x = v1.x-a.y+c.y;
	v2.y = v1.y+a.x-c.x;
	return its(u1, u2, v1, v2);
}

//����������
pt incenter(pt a, pt b, pt c)
{
	pt u1, u2, v1, v2;
	double m, n;
	u1 = a;
	m = atan2(b.y-a.y, b.x-a.x);
	n = atan2(c.y-a.y, c.x-a.x);
	u2.x = u1.x + cos((m+n)/2);
	u2.y = u1.y + sin((m+n)/2);
	v1 = b;
	m = atan2(a.y-b.y, a.x-b.x);
	n = atan2(c.y-b.y, c.x-b.x);
	v2.x = v1.x + cos((m+n)/2);
	v2.y = v1.y + sin((m+n)/2);
	return its(u1, u2, v1, v2);
}

//����
pt perpencenter(pt a, pt b, pt c)
{
	pt u1, u2, v1, v2;
	u1 = c;
	u2.x = u1.x-a.y+b.y;
	u2.y = u1.y+a.x-b.x;
	v1 = b;
	v2.x = v1.x-a.y+c.y;
	v2.y = v1.y+a.x-c.x;
	return its(u1, u2, v1, v2);
}

//����
//������������������ƽ������С�ĵ�
//�������ڵ����߾���֮�����ĵ�
pt barycenter(pt a, pt b, pt c)
{
	pt ret;
	ret.x = (a.x + b.x + c.x) / 3;
	ret.y = (a.y + b.y + c.y) / 3;
	return ret;
}

//�����
//�����������������֮����С�ĵ�
pt fermentpoint(pt a, pt b, pt c)
{
	pt u, v;
	double step = fabs(a.x)+fabs(a.y)+fabs(b.x)+fabs(b.y)+fabs(c.x)+fabs(c.y);
	int i, j, k;
	u.x = (a.x+b.x+c.x)/3;
	u.y = (a.y+b.y+c.y)/3;
	while (step > 1e-10)
	{
		for (k = 0; k < 10; step/=2, k++)
		for (i = -1; i <= 1; i++)
		for (j = -1; j <= 1; j++)
		{
			v.x = u.x + step*i;
			v.y = u.y + step*j;
			if (dis(u,a) + dis(u,b) + dis(u,c) > dis(v,a) + dis(v,b) + dis(v,c))
				u = v;
		}
	}
	return u;
}

//����ֱ����Բ�Ľ���,��ֱ֤����Բ�н���
//�����߶���Բ�Ľ����������������е��Ƿ����߶���
void intersection_line_circle(pt c, double r, pt l1, pt l2, pt &p1, pt &p2)
{
	pt p = c;
	p.x += l1.y - l2.y;
	p.y += l2.x - l1.x;
	p = its(p, c, l1, l2);
	double d = dis(p,c), t = sqrt(r*r - d*d) / dis(l1,l2);
	p2.x = p.x + (l2.x-l1.x)*t;
	p2.y = p.y + (l2.y-l1.y)*t;
	p1.x = p.x - (l2.x-l1.x)*t;
	p1.y = p.y - (l2.y-l1.y)*t;
}


//����Բ��Բ�Ľ���,��֤Բ��Բ�н���,Բ�Ĳ��غ�
void intersection_circle_circle(pt c1, double r1, pt c2, double r2, pt &p1, pt &p2)
{
	double d2 = (c1.x - c2.x) * (c1.x - c2.x) + (c1.y - c2.y) * (c1.y - c2.y);
	double cos = (r1 * r1 + d2 - r2 * r2) / (2 * r1 * sqrt(d2));
	pt v1 = (c2 - c1) / dis(c1, c2), v2 = pt(-v1.y, v1.x) * (r1 * sqrt(1 - cos * cos));
	pt X = c1 + v1 * (r1 * cos);
	p1 = X + v2;
	p2 = X - v2;
}


//��Բ�������е��Ӧ�ĽǶ�

//�������� 
void find_tp(double a, double b, double c, double &ang1, double &ang2)
{
	double v1, v2;
	v1 = fabs(c) > eps ? atan2(b*c, a*c) : atan2(b, a);
	v2 = acos(fabs(c)/sqrt(a*a+b*b));
	ang1 = v1 - v2;
	ang2 = v1 + v2;
}

//�⹫����(����Ƕ�t1 t2������Բ������)
void tangent1(pt c1, double r1, pt c2, double r2, double &t1, double &t2)
{
	find_tp(c2.x-c1.x, c2.y-c1.y, r1-r2, t1, t2);
}

//�ڹ�����(����Ƕ�t1 t2����Բc1����,��Բc2�����(��)PI)
void tangent2(pt c1, double r1, pt c2, double r2, double &t1, double &t2)
{
	find_tp(c2.x-c1.x, c2.y-c1.y, r1+r2, t1, t2);
} 

//ƽ��͹��
//��ʼ�㣺p[],������n	Ŀ��㣺f[],������top	 
void make_ch(pt *p, pt *f, int n, int &top)
{
	top = 0;
	sort(p, p + n);
	for (int i = 0; i < 2*n-1; i++)
	{
		int j = (i < n) ? i : 2*(n-1)-i;
		while (top > 1 && cpr(f[top-2], f[top-1], p[j]) < -eps)
			top--;
		f[top++] = p[j];
	}
	top--;
}

//ƽ��͹��
//��ʼ�㣺p[],������n	Ŀ���±����У�s[],������top	 
void make_ch(pt *p, int *s, int n, int &top)
{
	top = 0;
	sort(p, p + n);
	for (int i = 0; i < 2*n-1; i++)
	{
		int j = (i < n) ? i : 2*(n-1)-i;
		while (top > 1 && cpr(p[s[top-2]], p[s[top-1]], p[j]) < -eps)
			top--;
		s[top++] = j;
	}
	top--;
}

int main()
{}
