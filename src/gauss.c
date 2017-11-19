#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>

/* 
 * Function allocates memory for a matrix of size n * m on the heap and returns a handle to it.
 */
double **matrix_create(unsigned n, unsigned m) {
    double **matrix;
    if ((matrix = malloc(n * sizeof(matrix[0]))) == NULL) {
        // failed to allocate memory
        return NULL;
    }
    for (unsigned i = 0; i < n; i++) {
        if ((matrix[i] = malloc(m * sizeof(matrix[i][0]))) == NULL) {
            // failed to allocate memory
            for (unsigned j = 0; j < i; j++) {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

/*
 * Function releases memory occupied by matrix of size n * n.
 */
void matrix_destroy(double **matrix, unsigned n) {
    for (unsigned i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

/*
 * Function reads a matrix of size n * n and returns a handle to it.
 */

double **matrix_read(unsigned n) {
    double **matrix;
    if ((matrix = matrix_create(n, n)) == NULL) {
        return NULL;
    }
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            if (scanf("%lf", &matrix[i][j]) == EOF) {
                // failed to read
                matrix_destroy(matrix, n);
                errno = EINVAL;
                return NULL;
            }
        }
    }
    return matrix;
}

/*
 * Function reads a matrix of size n * 1 and returns a handle to it.
 */
double *matrix_read_vector(unsigned n) {
    double *vector;
    if ((vector = malloc(n * sizeof(vector[0]))) == NULL) {
        return NULL;
    }
    for (unsigned i = 0; i < n; i++) {
        if (scanf("%lf", &vector[i]) == EOF) {
            // falied to read
            free(vector);
            errno = EINVAL;
            return NULL;
        }
    }
    return vector;
}

/*
 * Function prints a matrix of size n * n.
 */
int matrix_print(double **matrix, unsigned n) {
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j  = 0; j < n; j++) {
            // substitute negative zero with zero
            if (matrix[i][j] == 0) {
                matrix[i][j] = +0.0;
            }
            // print
            if (printf("%15.10g ", matrix[i][j]) == 0) {
                // failed to print
                return -1;
            }
        }
        // put 'end of line'
        if (putchar('\n') == EOF) {
            // failed to print    
            return -1;
        }
    }
    return 0;
}

/*
 * Function prints a matrix of size n * 1.
 */
int matrix_print_vector(double *vector, unsigned n) {
    for (unsigned i = 0; i < n; i++) {
        // substitute negative zero with zero
        if (vector[i] == 0) {
            vector[i] = +0.0;
        }
        // print
        if (printf("%15.10g ", vector[i]) == 0) {
            // failed to print
            return -1;
        }
    }
    // put 'end of line'
    if (putchar('\n') == EOF) {
        return -1;
    }
    return 0;
}

/*
 * Function solves matrix equation Ax = f using Gaussian Elimination method
 */
void matrix_forward(double **matrix, unsigned n, unsigned m);
void matrix_back(double **matrix, unsigned n, unsigned m); 
void matrix_normalize(double **matrix, unsigned n, unsigned m);

double *matrix_gauss_solve(double **matrix, const double *f, unsigned n) {
   // create augmented matrix
   double **aug = matrix_create(n, n + 1);
   if (aug == NULL) {
       return NULL;
   }
   for (unsigned i = 0; i < n; i++) {
       for (unsigned j = 0; j < n; j++) {
           aug[i][j] = matrix[i][j];
       }
       aug[i][n] = f[i];
   }

   // perform forward elimination and normalization
   matrix_forward(aug, n, n + 1);
   matrix_normalize(aug, n, n + 1);
   matrix_back(aug, n, n + 1);

   // calculate result
   double *result;
   if ((result = malloc(n * sizeof(result[0]))) == NULL) {
       matrix_destroy(aug, n);
       return NULL;
   }
   for (unsigned i = 0; i < n; i++) {
       result[i] = aug[i][n];
   }

   // free memory
   matrix_destroy(aug, n);
   return result;
}

/* Function returns the number of the row with the greatest primary element */
unsigned matrix_find_greatest(double **matrix, unsigned current, unsigned n) {
    unsigned result = current;
    for (unsigned j = current + 1; j < n; j++) {
        if (matrix[j][current] != 0 && (!matrix[result][current] || fabs(matrix[result][current]) < fabs(matrix[j][current])) ) {
            result = j;
        }
    }
    return result;
}

/* Swap lines */
void matrix_swap_rows(double **matrix, unsigned i, unsigned j) {
    double *temp = matrix[i];
    matrix[i] = matrix[j];
    matrix[j] = temp;
}

/* Function subtracts current line from the following */
void matrix_subtract(double **matrix, unsigned current, unsigned n, unsigned m) {
    for (unsigned i = current + 1; i < n; i++) {
        double multiplier = matrix[i][current] / matrix[current][current];
        if (multiplier != 0) {
            for (unsigned j = current; j < m; j++) {
                matrix[i][j] -= multiplier * matrix[current][j];
            }
        }
    }
}

/* Function performs Forward Elimination of the augmented matrix n * m (n >= m) */
void matrix_forward(double **matrix, unsigned n, unsigned m) {
    for (unsigned i = 0; i < n; i++) {
        // subtract greatest from others
        unsigned greatest = matrix_find_greatest(matrix, i, n);
        if (i != greatest) {
            matrix_swap_rows(matrix, i, greatest);
        }
        matrix_subtract(matrix, i, n, m);
    }
}

/* Function normalizes upper triangular n * m matrix */
void matrix_normalize(double **matrix, unsigned n, unsigned m) {
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = i + 1; j < m; j++) {
            matrix[i][j] /= matrix[i][i];
        }
        matrix[i][i] = 1;
    }
}

/* Function performs back substitution */
void matrix_back(double **matrix, unsigned n, unsigned m) {
    for (unsigned i = n - 1; i > 0; i--) {
        for (unsigned prev = 0; prev < i; prev++) {
            // subtract current row from previous
            double multiplier = matrix[prev][i];
            if (multiplier == 0) {
                continue;
            }
            for (unsigned j = i; j < m; j++) {
                matrix[prev][j] -= matrix[i][j] * multiplier;
            }
        }
    }
}

/*
 * Function returns determinant of an n * n matrix.
 */
double matrix_determinant(double **matrix, unsigned n) {
    
    // calculate upper triangulal matrix
    double **aug = matrix_create(n, n);
    if (aug == NULL) {
        return 0; //::TODO how to handle?
    }
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            aug[i][j] = matrix[i][j];
        }
    }
    matrix_forward(aug, n, n);

    // calculate determinant
    double result = 1;
    for (unsigned i = 0; i < n; i++) {
        result *= aug[i][i];
    }
    return result;
}

/*
 *  Function returns a handle to te inverse of an n * n matrix.
 */
double **matrix_inverse(double **matrix, unsigned n) {
    // create augmented matrix
    double **aug = matrix_create(n, 2 * n);
    if (aug == NULL) {
        return NULL;
    }
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            aug[i][j] = matrix[i][j];
        }
        for (unsigned j = n; j < n * 2; j++) {
            aug[i][j] = (i == (j - n) ? 1 : 0);
        }
    }

    // perform forward elimination and back substitution
    matrix_forward(aug, n, 2 * n);
    matrix_normalize(aug, n, 2 * n);
    matrix_back(aug, n, 2 * n);

    // copy the result
    double **result = matrix_create(n, n);
    if (result == NULL) {
        return NULL;
    }
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            result[i][j] = aug[i][j + n];
        }
    } 

    // free memory
    matrix_destroy(aug, n);

    return result;
}