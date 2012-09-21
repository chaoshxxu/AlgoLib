#include <cstdio>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <cstdlib>
#include<algorithm>
#include<time.h>

#define eps 1e-7
using namespace std;

FILE *visualizer;
int tc;

double xminv, xmaxv;
double yminv, ymaxv;


struct pt
{
	double x, y;
	pt(){}
	pt(double _x, double _y)
	{
		x = _x;
		y = _y;
	}
	pt operator - (const pt p1){return pt(x - p1.x, y - p1.y);}
	pt operator + (const pt p1){return pt(x + p1.x, y + p1.y);}
	pt operator / (double r){return pt(x / r, y / r);}
	pt operator * (double r){return pt(x * r, y * r);}
	bool operator == (const pt &p1)const{return fabs(y-p1.y)<eps && fabs(x-p1.x)<eps;}
	bool operator != (const pt &p1)const{return fabs(y-p1.y)>eps || fabs(x-p1.x)>eps;}
};

double cpr(const pt &a,const pt &b,const pt &c){return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);}
double dpr(const pt &a,const pt &b,const pt &c){return (b.x-a.x)*(c.x-a.x)+(b.y-a.y)*(c.y-a.y);}
double dis(const pt &a, const pt &b){return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));}

///////////////////////////////////////////////////////////////////////

void mp_begin() {
    char s[30];
    tc++;
	sprintf(s, "Graphic_%d.mp", tc);
    visualizer = fopen(s, "w");
    assert(visualizer != NULL);
    fprintf(visualizer, "beginfig(1);\n");
    fprintf(visualizer, "u := 2cm;\n");
}

void mp_end() {
    fprintf(visualizer, "endfig;\n");
    fprintf(visualizer, "end\n");
    fclose(visualizer);
    char s[30];
	sprintf(s, "mpost Graphic_%d.mp", tc);	system(s);
	sprintf(s, "epstopdf Graphic_%d.1", tc);	system(s);
	sprintf(s, "del Graphic_%d.1", tc);		system(s);
	sprintf(s, "del Graphic_%d.log", tc);	system(s);
	sprintf(s, "del Graphic_%d.mp", tc);		system(s);
}

void mp_print_segment(pt a, pt b, const int& pen_width = 1, const string& color = "black") {
	a.x = -50 + 100.0 / (xmaxv - xminv) * (a.x - xminv);
	a.y = -50 + 100.0 / (ymaxv - yminv) * (a.y - yminv);
	b.x = -50 + 100.0 / (xmaxv - xminv) * (b.x - xminv);
	b.y = -50 + 100.0 / (ymaxv - yminv) * (b.y - yminv);
    fprintf(visualizer, "draw ");
    fprintf(visualizer, "(%.2lfu, %.2lfu) -- (%.2lfu, %.2lfu) withpen pencircle scaled %dpt withcolor %s;\n", a.x, a.y, b.x, b.y, pen_width, color.c_str());
}

//【注意】边太多的话会挂掉!最好一条边一条边自己画。 
void mp_print_polygon(const vector<point>& p, const int& pen_width = 1, const string& color = "black") {
    fprintf(visualizer, "draw ");
    for (int i = 0; i < p.size(); ++i) {
        fprintf(visualizer, "(%.2lfu, %.2lfu) -- ", p[i].x, p[i].y);
    }
    fprintf(visualizer, "cycle withpen pencircle scaled %dpt withcolor %s;\n", pen_width, color.c_str());
}


void mp_print_point(pt a, const int& pen_width = 2, const string& color = "black") {
	a.x = -50 + 100.0 / (xmaxv - xminv) * (a.x - xminv);
	a.y = -50 + 100.0 / (ymaxv - yminv) * (a.y - yminv);

    fprintf(visualizer, "draw ");
    fprintf(visualizer, "(%.2lfu, %.2lfu) withpen pencircle scaled %dpt withcolor %s;\n", a.x, a.y, pen_width, color.c_str());
}

void mp_print_circle(const pt& a, const double& r, const int& pen_width = 1, const string& color = "black") {
    fprintf(visualizer, "draw ");
    fprintf(visualizer, "fullcircle scaled %.2lfu shifted (%.2lfu, %.2lfu) withpen pencircle scaled %dpt withcolor %s;\n", r * 2, a.x, a.y, pen_width, color.c_str());
}

/////////////////////////////////////////////////

int main()
{
	return 0;
}
