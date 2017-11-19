#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

extern double **matrix_create(unsigned n, unsigned m);
extern void matrix_destroy(double **matrix, unsigned n);

/* Function normalizes lower triangular n * m matrix */
void normalize_reversed(double **matrix, unsigned n, unsigned m) {
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < i; j++) {
            matrix[i][j] /= matrix[i][i];
        }
		matrix[i][m - 1] /= matrix[i][i];
		matrix[i][i] = 1;
    }
}

/* Function performs back substitution */
void back_reversed(double **matrix, unsigned n, unsigned m) {
    for (unsigned i = 0; i < n; i++) {
        for (unsigned next = i + 1; next < n; next++) {
            // subtract current row from previous
            double multiplier = matrix[next][i];
            if (multiplier == 0) {
                continue;
            }
            for (unsigned j = 0; j < m; j++) {
                matrix[next][j] -= matrix[i][j] * multiplier;
            }
        }
    }
}


/*
 * Function calculates next element of the iteration sequence.
 */
void matrix_iteration_next(const double **a, double **b, const double *f, unsigned n, double omega, double *current) {
	// create new right part in the augmented B matrix
	for (unsigned i = 0; i < n; i++) {
		b[i][n] = f[i];
		for (unsigned j = i; j < n; j++) {
			if (i == j) {
				b[i][n] -= (1 - 1 / omega) * a[i][j] * current[i];
			} else {
				b[i][n] -= a[i][j] * current[j];
			}
		}
	}

	// solve matrix equation with the new f vector: Bx = f
	normalize_reversed(b, n, n + 1);
	back_reversed(b, n, n + 1);
	
	// copy the result
	for (unsigned i = 0; i < n; i++) {
		current[i] = b[i][n];
	}
}

/*
 * Function calculates vector norm
 */
double matrix_vector_diff(const double *a, const double *b, unsigned n) {
	double result = 0;
	for (unsigned i = 0; i < n; i++) {
		result += fabs(a[i] - b[i]);
	}

	return result;
}

/*
 * Functions calculates an estimate solution of the matrix equation 
 * Ax = f using succesive over-relaxation method. Function returns a handle to the result. 
 */
double *matrix_iteration_solve(const double **a, const double *f, unsigned n, double omega, const double *start, double precision) {
	// allocate vectors for current and previous elements of iteration sequence
	double *current, *previous;
	if ((current = malloc(sizeof(current[0]) * n)) == NULL) {
		return NULL;
	}
	if ((previous = malloc(sizeof(previous[0]) * n)) == NULL) {
		free(current);
		return NULL;
	}

	// create matrix B
	double **aug;
	if ((aug = matrix_create(n, n + 1)) == NULL) {
		free(current);
		free(previous);
		return NULL;
	}

	// start iteration process
	unsigned iterationCount = 0;
	memcpy(current, start, n * sizeof(start[0]));
	do {
		// construct B matrix
		for (unsigned i = 0; i < n; i++) {
			for (unsigned j = 0; j < n; j++) {
				if (i < j) {
					aug[i][j] = 0;
				} else if (i == j) {
					aug[i][j] = (1 / omega) * a[i][j];
				} else {
					aug[i][j] = a[i][j];
				}
			}
		}

		// move result of the previous iteration and perform iteration
		memcpy(previous, current, n * sizeof(current[0]));
		matrix_iteration_next(a, aug, f, n, omega, current);
		iterationCount++;
	} while(matrix_vector_diff(current, previous, n) > precision);

	printf("> log: %d iterations made\n", iterationCount);

	// release memory
	matrix_destroy(aug, n);
	free(previous);

	return current;
}