/*
  Name: 3D Computing Geometry
  Copyright: 
  Author: 
  Date: 05-11-09 11:59
  Description: 
*/

#include<stdio.h>
#include<math.h>
#include<algorithm>
#define MAXV 2000
#define PI 3.14159265358979323846
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
using namespace std;

//��ά��
struct pt3
{
	double x, y, z;
	pt3(){}
	pt3(double _x, double _y, double _z): x(_x), y(_y), z(_z){}
	void read()
	{
		scanf("%lf %lf %lf", &x, &y, &z);
	}
	void disp()
	{
		printf("%lf %lf %lf\n", x, y, z);
	}
	pt3 operator - (const pt3 p1){return pt3(x - p1.x, y - p1.y, z - p1.z);}
	pt3 operator + (const pt3 p1){return pt3(x + p1.x, y + p1.y, z + p1.z);}
	pt3 operator * (pt3 p){return pt3(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);}		//��� 
	double operator ^ (pt3 p){return x*p.x+y*p.y+z*p.z;}							//��� 
	pt3 operator / (double r){return pt3(x / r, y / r, z / r);}
	pt3 operator * (double r){return pt3(x * r, y * r, z * r);}
	bool operator == (const pt3 p1)const{return fabs(x-p1.x)<eps && fabs(y-p1.y)<eps && fabs(z-p1.z)<eps;}
	bool operator != (const pt3 p1)const{return fabs(x-p1.x)>eps || fabs(y-p1.y)>eps || fabs(z-p1.z)>eps;}
};

//������� 
double vlen(const pt3 &a){return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);}
double dis(const pt3 &a, const pt3 &b){return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));}

//�㵽ƽ����� 
double ptoplane(pt3 p, pt3 s1, pt3 s2, pt3 s3)
{
	pt3 norm = (s2 - s1) * (s3 - s1);
	return fabs(norm ^ (p - s1)) / vlen(norm);
}

//�㵽ֱ�߾���
double ptoline(pt3 p, pt3 l1, pt3 l2){
	return vlen((p-l1)*(l2-l1)) / dis(l1, l2);
}

//ֱ�ߵ�ֱ�߾���
double linetoline(pt3 u1, pt3 u2, pt3 v1, pt3 v2){
	pt3 n = (u1 - u2) * (v1 - v2);
	return fabs((u1 - v1) ^ n) / vlen(n);
}

//�е��Ƿ��ڿռ���������,�����߽�,���㹲��������
bool dot_intri_in(pt3 p, pt3 s1, pt3 s2, pt3 s3)
{
	return zero(vlen((s1-s2)*(s1-s3))-vlen((p-s1)*(p-s2))-vlen((p-s2)*(p-s3))-vlen((p-s3)*(p-s1)));
}

//�е��Ƿ��ڿռ���������,�������߽�,���㹲��������
bool dot_intri_ex(pt3 p, pt3 s1, pt3 s2, pt3 s3)
{
	return dot_intri_in(p,s1,s2,s3)&&
	vlen((p-s1)*(p-s2))>eps&&vlen((p-s2)*(p-s3))>eps&&vlen((p-s3)*(p-s1))>eps;
}

//����ֱ����ƽ�潻��,ע�������ж��Ƿ�ƽ��,����֤���㲻����!
//�߶κͿռ������ν����������ж�
pt3 intersection(pt3 l1, pt3 l2, pt3 s1, pt3 s2, pt3 s3)
{
	pt3 ret = (s1 - s2) * (s2 - s3);
	double t = (ret ^ (s1 - l1)) / (ret ^ (l2 - l1));
	return l1 + (l2 - l1) * t;
}

//�ж�ֱ��ab�Ƿ񴩹��ռ�������s1s2s3 
bool line_throughtri(pt3 a, pt3 b, pt3 s1, pt3 s2, pt3 s3)
{
	pt3 norm = (s2 - s1) * (s3 - s1);
	if (((a - s1)^norm) * ((b - s1)^norm) > 0 || fabs((a - b)^norm) < eps)
		return 0;
	pt3 X = intersection(a, b, s1, s2, s3);
	return dot_intri_ex(X, s1, s2, s3);
}

//��p�ƹ�ԭ�������v��ת˳ʱ���A��õ��ĵ�
pt3 rotate(pt3 p, pt3 v, double A)
{
	double len = sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
	double x = v.x/len, y = v.y/len, z = v.z/len;
	double M[][3] = 
	{
		cos(A)+(1-cos(A))*x*x, (1-cos(A))*x*y-sin(A)*z, (1-cos(A))*x*z+sin(A)*y,
		(1-cos(A))*y*x+sin(A)*z, cos(A)+(1-cos(A))*y*y, (1-cos(A))*y*z-sin(A)*x,
		(1-cos(A))*z*x-sin(A)*y, (1-cos(A))*z*y+sin(A)*x, cos(A)+(1-cos(A))*z*z
	};
	return pt3 (p.x * M[0][0] + p.y * M[1][0] + p.z * M[2][0],
				p.x * M[0][1] + p.y * M[1][1] + p.z * M[2][1],
				p.x * M[0][2] + p.y * M[1][2] + p.z * M[2][2]);
}

//ȡƽ�淨����
pt3 pvec(pt3 s1, pt3 s2, pt3 s3){
	return (s1 - s2) * (s2 - s3);
}

//�����㹲��
bool dots_inline(pt3 p1, pt3 p2, pt3 p3){
	return vlen((p1 - p2) * (p2 - p3)) < eps;
}

//���ĵ㹲��
bool dots_onplane(pt3 a, pt3 b, pt3 c, pt3 d){
	return zero(pvec(a, b, c) ^ (d - a));
}

//�е��Ƿ����߶���,�����˵�
bool dot_online_in(pt3 p, pt3 l1, pt3 l2){
	return vlen((p - l1) * (p - l2)) < eps && ((p - l1) ^ (p - l2)) < eps;
}

//�е��Ƿ����߶���,�������˵�
bool dot_online_ex(pt3 p, pt3 l1, pt3 l2){
	return vlen((p - l1) * (p - l2)) < eps && ((p - l1) ^ (p - l2)) < -eps;
}

//���������߶�ͬ��,�����߶��Ϸ���0,������������
bool same_side(pt3 p1, pt3 p2, pt3 l1, pt3 l2){
	return ((l1 - l2) * (p1 - l2) ^ (l1 - l2) * (p2 - l2)) > eps;
}

//���������߶����,�����߶��Ϸ���0,������������
bool opposite_side(pt3 p1,pt3 p2,pt3 l1,pt3 l2){
	return ((l1 - l2) * (p1 - l2) ^ (l1 - l2) * (p2 - l2)) < -eps;
}

//��������ƽ��ͬ��,����ƽ���Ϸ���0
bool same_side(pt3 p1, pt3 p2, pt3 s1, pt3 s2, pt3 s3){
	return (pvec(s1,s2,s3) ^ (p1-s1)) * (pvec(s1,s2,s3) ^ (p2 - s1)) > eps;
}

//��������ƽ�����,����ƽ���Ϸ���0
int opposite_side(pt3 p1, pt3 p2, pt3 s1, pt3 s2, pt3 s3){
	return (pvec(s1,s2,s3) ^ (p1-s1)) * (pvec(s1,s2,s3) ^ (p2 - s1)) < -eps;
}

//����ֱ��ƽ��
int parallel(pt3 u1, pt3 u2, pt3 v1, pt3 v2){
	return vlen((u1 - u2) * (v1 - v2)) < eps;
}

//����ƽ��ƽ��
int parallel(pt3 u1, pt3 u2, pt3 u3, pt3 v1, pt3 v2, pt3 v3){
	return vlen(pvec(u1,u2,u3) * pvec(v1,v2,v3)) < eps;
}

//��ֱ����ƽ��ƽ��
int parallel(pt3 l1, pt3 l2, pt3 s1, pt3 s2, pt3 s3){
	return zero((l1 - l2) ^ pvec(s1,s2,s3));
}

//����ֱ�ߴ�ֱ
int perpendicular(pt3 u1, pt3 u2, pt3 v1, pt3 v2){
	return zero((u1 - u2) ^ (v1 - v2));
}

//����ƽ�洹ֱ
int perpendicular(pt3 u1, pt3 u2, pt3 u3, pt3 v1, pt3 v2, pt3 v3){
	return zero(pvec(u1, u2, u3) ^ pvec(v1, v2, v3));
}

//��ֱ����ƽ�洹ֱ
int perpendicular(pt3 l1, pt3 l2, pt3 s1, pt3 s2, pt3 s3){
	return vlen((l1 - l2) * pvec(s1, s2, s3)) < eps;
}

//�����߶��ཻ,�����˵�Ͳ����غ�
int intersect_in(pt3 u1, pt3 u2, pt3 v1, pt3 v2){
	if (!dots_onplane(u1, u2, v1, v2))
		return 0;
	if (!dots_inline(u1, u2, v1) || !dots_inline(u1, u2, v2))
		return !same_side(u1, u2, v1, v2) && !same_side(v1, v2, u1, u2);
	return  dot_online_in(u1, v1, v2) || dot_online_in(u2, v1, v2) || 
			dot_online_in(v1, u1, u2) || dot_online_in(v2,u1,u2);
}

//�����߶��ཻ,�������˵�Ͳ����غ�
int intersect_ex(pt3 u1, pt3 u2, pt3 v1, pt3 v2){
	return  dots_onplane(u1, u2, v1, v2) &&
			opposite_side(u1, u2, v1, v2) && opposite_side(v1, v2, u1, u2);
}

//���߶���ռ��������ཻ,�������ڱ߽��(����)����
int intersect_in(pt3 l1, pt3 l2, pt3 s1, pt3 s2, pt3 s3){
	return  !same_side(l1, l2, s1, s2, s3) &&
			!same_side(s1, s2, l1, l2, s3) &&
			!same_side(s2, s3, l1, l2, s1) &&
			!same_side(s3, s1, l1, l2, s2);
}

//���߶���ռ��������ཻ,���������ڱ߽��(����)����
int intersect_ex(pt3 l1, pt3 l2, pt3 s1, pt3 s2, pt3 s3){
	return  opposite_side(l1, l2, s1, s2, s3) && 
			opposite_side(s1, s2, l1, l2, s3) &&
			opposite_side(s2, s3, l1, l2, s1) &&
			opposite_side(s3, s1, l1, l2, s2);
}
 

//������ֱ�߽���,ע�������ж�ֱ���Ƿ����ƽ��!
//�߶ν������������߶��ཻ(ͬʱ����Ҫ�ж��Ƿ�ƽ��!)
pt3 intersection(pt3 u1, pt3 u2, pt3 v1, pt3 v2){
	pt3 ret = u1;
	double t =  ((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))/
   				((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
	return ret + (u2 - u1) * t;
}

//������ƽ�潻��,ע�������ж��Ƿ�ƽ��,����֤���㲻����!
void intersection(pt3 u1, pt3 u2, pt3 u3, pt3 v1, pt3 v2, pt3 v3, pt3 &a, pt3 &b){
	a = parallel(v1,v2,u1,u2,u3)?intersection(v2,v3,u1,u2,u3):intersection(v1,v2,u1,u2,u3);
	b = parallel(v3,v1,u1,u2,u3)?intersection(v2,v3,u1,u2,u3):intersection(v3,v1,u1,u2,u3);
}


//��ֱ�߼н�cosֵ
double angle_cos(pt3 u1, pt3 u2, pt3 v1, pt3 v2){
	return ((u1-u2)^(v1-v2)) / vlen(u1-u2) / vlen(v1-v2);
}

//��ƽ��н�cosֵ
double angle_cos(pt3 u1, pt3 u2, pt3 u3, pt3 v1, pt3 v2, pt3 v3){
	return (pvec(u1,u2,u3) ^ pvec(v1,v2,v3)) / vlen(pvec(u1,u2,u3)) / vlen(pvec(v1,v2,v3));
}

//ֱ�ߺ�ƽ��н�sinֵ
double angle_sin(pt3 l1, pt3 l2, pt3 s1, pt3 s2, pt3 s3){
	return ((l1-l2)^pvec(s1,s2,s3)) / vlen(l1-l2) / vlen(pvec(s1,s2,s3));
}

////////////////////////////////////////////////////////////////////

int main()
{
	return 0;
}
