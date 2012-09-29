#include<stdio.h>
#include<math.h>
#include<algorithm>
using namespace std;

const int N = 3;
const double eps = 1e-8;

double a[][N] = {
	2, 3, 4,
	3, 4, 5,
	8, 6, 5
};

double b[] = {
	29,
	38,
	54
};

void debug() {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			printf("%g ", a[i][j]);
		}
		printf("- %g\n", b[i]);
	}
	puts("");
}

//返回主元个数
int gause(double a[N][N], double b[N]) {
	int majorNumber = N;
	for (int i = 0; i < N; ++i) {
		int major = i;
		for (int j = i; j < N; ++j) {
			if (fabs(a[j][i]) > fabs(a[major][i])) {
				major = j;
			}
		}
		for (int k = 0; k < N; ++k) {
			swap(a[major][k], a[i][k]);
		}
		swap(b[major], b[i]);
		if (fabs(a[i][i]) < eps) {
			--majorNumber;
			continue;
		}
		b[i] /= a[i][i];
		for (int k = N - 1; k >= i; --k) {
			a[i][k] /= a[i][i];
		}
		for (int j = 0; j < N; ++j) if (j != i && fabs(a[j][i]) > eps){
			b[j] -= b[i] * a[j][i];
			for (int k = N - 1; k >= i; --k) {
				a[j][k] -= a[i][k] * a[j][i];
			}
		}
	}
	return majorNumber;
}


int main() {
	printf("%d\n", gause(a, b));
	debug();
}
