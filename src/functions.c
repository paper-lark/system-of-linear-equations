#include <stdlib.h>
#include <math.h>
#include <assert.h>
/*
 * Functions for tasks in Appendix 2.
 * Tasks 1-4
 */

extern double **matrix_create(unsigned n, unsigned m);

static unsigned M = 0;
static unsigned n = 0;
double *(*matrix_fill_vector)(double x) = NULL;

static double *f1(double x) {
	double *vector;
	if ((vector = malloc(n * sizeof(vector[0]))) == NULL) {
		return NULL;
	}
	for (unsigned i = 0; i < n; i++) {
		vector[i] = n * exp(x / i) * cos(x);
	}

	return vector;
}

static double *f2(double x) {
	double *vector;
	if ((vector = malloc(n * sizeof(vector[0]))) == NULL) {
		return NULL;
	}
	for (unsigned i = 0; i < n; i++) {
		vector[i] = fabs(x - (double) n / 10) * i * sin(x);
	}

	return vector;
}

static double *f3(double x) {
	double *vector;
	if ((vector = malloc(n * sizeof(vector[0]))) == NULL) {
		return NULL;
	}
	for (unsigned i = 0; i < n; i++) {
		vector[i] = x * exp(x / i) * cos(x / i);
	}

	return vector;
}

static double *f4(double x) {
	double *vector;
	if ((vector = malloc(n * sizeof(vector[0]))) == NULL) {
		return NULL;
	}
	for (unsigned i = 0; i < n; i++) {
		vector[i] = n * exp(x / i) * cos(x);
	}

	return vector;
}

/*
 * Function initializes the module for using formula-determined matrices
 */ 
unsigned functions_init(unsigned task) {
	switch(task) {
		case 1:
			M = 1;
			n = 50;
			matrix_fill_vector = f1;
			break;
		case 2:
			M = 2;
			n = 40;
			matrix_fill_vector = f2;
			break;
		case 3:
			M = 3;
			n = 30;
			matrix_fill_vector = f3;
			break;
		case 4:
			M = 4;
			n = 100;
			matrix_fill_vector = f4;
			break;
		default:
			assert(0);
	}

	return n;
}

double **matrix_fill() {
	const double q = 1.001 - 2 * M * 0.001;

	double **matrix;
	if ((matrix = matrix_create(n, n)) == NULL) {
		return NULL;
	}
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < n; j++) {
			if (i == j) {
				matrix[i][j] = pow(q - 1, i + j);
			} else {
				matrix[i][j] = pow(q, i + j) + 0.1 * (j - i);
			}
		}
	}

	return matrix;
}
