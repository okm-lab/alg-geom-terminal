#define _CRT_SECURE_NO_WARNINGS
#include "Main.h"
#include <stdio.h>
#include <math.h>

double powerD(double a, int n) {
	double t = 1;
	for (int i = 0; i < n; i++) {
		t *= a;
	}
	return t;
}

int fact(int n) {
	int result = 1;
	for (int i = 1; i <= n; i++) {
		result *= i;
	}
	return result;
}

void algToTrig(double* r, double* a, double* b) {
	*r = sqrt(*a * *a + *b * *b);
	*a /= *r;
	*b /= *r;
}

void trigToAlg(double* r, double* a, double* b) {
	*a *= *r;
	*b *= *r;
}

void revTrig(double* r, double* a, double* b) {
	*r = 1 / *r;
	*b = -*b;
}

void revAlg(double* a, double* b) {
	double t = *a;
	*a = *a / (*a * *a + *b * *b);
	*b = -*b / (t * t + *b * *b);
}

void powerAlg(double* a, double* b, int n) {
	int ka = n;
	int kb = 0;
	int suma = 0;
	int sumb = 0;
	while (ka >= 0 && kb <= n) {
		int c_n_k = fact(n) / (fact(kb) * fact(n - kb));

		switch (kb % 4) {
		case 0:
			suma += c_n_k * powerD(*a, ka) * (powerD(*b, (n - ka)));
			break;
		case 1:
			sumb += c_n_k * powerD(*a, ka) * (powerD(*b, (n - ka)));
			break;
		case 2:
			suma -= c_n_k * powerD(*a, ka) * (powerD(*b, (n - ka)));
			break;
		case 3:
			sumb -= c_n_k * powerD(*a, ka) * (powerD(*b, (n - ka)));
			break;
		}
		ka--;
		kb++;
	}
	*a = suma;
	*b = sumb;
}

void powerTrig(double* r, double* a, double* b, int n) {
	*r = pow(*r, n);
	*a = cos(acos(*a) * n);
	*b = sin(asin(*b) * n);
}

void rootTrig(double r, double a, double b, int n) {
	printf("Roots:\n");
	for (int i = 0; i < n; i++) {
		double r0 = pow(r, 1.0 / n);
		double a0 = cos((acos(a) + 2 * acos(-1.0) * i) / n);
		double b0 = sin((asin(b) + 2 * acos(-1.0) * i) / n);
		if (b0 >= 0)
			printf("%d: %.4g(cos(%.4gpi) + isin(%.4gpi))     ", i + 1, r0, acos(a0)/acos(-1.0), acos(b0)/acos(-1.0));
		else
			printf("%d: %.4g(cos(%.4gpi) - isin(%.4gpi)     ", i + 1, r0, acos(a0)/acos(-1.0), -asin(b0)/acos(-1.0));
		trigToAlg(&r0, &a0, &b0);
		if (b0 >= 0)
			printf("(%.4g + %.4gi)\n", a0, b0);
		else
			printf("(%.4g - %.4gi)\n", a0, -b0);
	}
	printf("\n\n");
}

int main() {
	double a, b, r, r2;
	int selection = -1, state = -1, n, sqrt_state;
	printf("Select number format:\n1. Algebraic form\n2. Trigonometric form\n");
	while (selection != 1 && selection != 2)
		scanf("%d", &selection);

	system("cls");

	if (selection == 1) {
		printf("Input a and b (z = a + bi)\n");
		scanf("%lf%lf", &a, &b);
		state = 1;
	}
	else {
		printf("Input r, cosa, sina (z = r(cosa + sina))\n");
		scanf("%lf%lf%lf", &r, &a, &b);
		r2 = r * r;
		state = 2;
	}

	system("cls");

	while (selection != -1) {
		if (state == 1) {
			if (b >= 0)
				printf("Current number: %.4g + %.4gi", a, b);
			else
				printf("Current number: %.4g - %.4gi", a, -b);
			printf("\n\n\nSelect:\n1. Reverse\n2. Transform to trigonometric form\n3. Power\n4. Root\n-1. Exit\n\n\n\n");
			scanf("%d", &selection);

			while ((selection < 1 && selection != -1) || selection > 4)
				scanf("%d", &selection);

			system("cls");

			switch (selection) {
			case 1:
				revAlg(&a, &b);
				break;
			case 2:
				algToTrig(&r, &a, &b);
				state = 2;
				break;
			case 3:
				printf("Enter power:\n");
				scanf("%d", &n);
				system("cls");
				powerAlg(&a, &b, n);
				break;
			case 4:
				printf("Enter root:\n");
				scanf("%d", &n);
				system("cls");
				algToTrig(&r, &a, &b);
				rootTrig(r, a, b, n);
				trigToAlg(&r, &a, &b);
				break;
			}
		}
		else {
			if (b >= 0)
				printf("Current number: %.4g(cos(%.4gpi) + isin(%.4gpi)) ", r, acos(a)/acos(-1.0), asin(b)/acos(-1.0));
			else
				printf("Current number: %.4g(cos(%.4g) - isin(%.4g))", r, acos(a)/acos(-1.0), -asin(b)/acos(-1.0));
			printf("\n\n\nSelect:\n1. Reverse\n2. Transform to algebraic form\n3. Power\n4. Root\n-1. Exit\n\n\n\n");
			scanf("%d", &selection);

			while ((selection < 1 && selection != -1) || selection > 4)
				scanf("%d", &selection);

			system("cls");
			switch (selection) {
			case 1:
				revTrig(&r, &a, &b);
				break;
			case 2:
				trigToAlg(&r, &a, &b);
				state = 1;
				break;
			case 3:
				printf("Enter power:\n");
				scanf("%d", &n);
				system("cls");
				powerTrig(&r, &a, &b, n);
				break;
			case 4:
				printf("Enter root:\n");
				scanf("%d", &n);
				system("cls");
				rootTrig(r, a, b, n);
				break;
			}
		}
	}
}