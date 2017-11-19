#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

/* Prototypes */
double **matrix_create(unsigned n, unsigned m);
double **matrix_read(unsigned n);
double *matrix_read_vector(unsigned n);
int matrix_print(double **matrix, unsigned n);
int matrix_print_vector(double *vector, unsigned n);
void matrix_destroy(double **matrix, unsigned n);
double matrix_determinant(double **matrix, unsigned n);
double **matrix_inverse(double **matrix, unsigned n);
double *matrix_gauss_solve(double **matrix, double *f, unsigned n);


/*
 * Entry Point
 */
int main(void) {
    // read matrix
    unsigned n;
    scanf("%u", &n);
    double **matrix, *f;
    if ((matrix = matrix_read(n)) == NULL || (f = matrix_read_vector(n)) == NULL) {
        fprintf(stderr, "> error: reading failed with message: %s\n", strerror(errno));
        exit(1);
    }
   
    // print processed matrix
    printf("> matrix:\n");
    matrix_print(matrix, n);
    putchar('\n');

    // calculate inverse matrix
    double **inverse;
    if ((inverse = matrix_inverse(matrix, n)) == NULL) {
        fprintf(stderr, "> error: invertion failed\n");
        exit(1);
    }
    printf("> inverse matrix:\n");
    matrix_print(inverse, n);
    matrix_destroy(inverse, n);
    putchar('\n');

    // calculate condition number
    //::TODO

    // calculate solution using Gaussian elimination method
    double *solution = matrix_gauss_solve(matrix, f, n); //::TODO error handling
    printf("> solution:\n");
    matrix_print_vector(solution, n);
    putchar('\n');
    
    // calculate determinant
    double determinant = matrix_determinant(matrix, n);
    printf("> determinant:\n %5.10g\n", determinant);
    putchar('\n');

    // free memory
    free(solution);
    free(f);
    matrix_destroy(matrix, n);

    // completion message
    printf("Program sucessfully finished!\n");
}

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
 * Function releases memory occupied by matrix of size n * n.
 */
void matrix_destroy(double **matrix, unsigned n) {
    for (unsigned i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}


/*
 * Function solves matrix equation Ax = f using Gaussian Elimination method
 */
void forward(double **matrix, unsigned n, unsigned m);
void back(double **matrix, unsigned n, unsigned m); 
void normalize(double **matrix, unsigned n, unsigned m);

double *matrix_gauss_solve(double **matrix, double *f, unsigned n) {
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
   forward(aug, n, n + 1);
   normalize(aug, n, n + 1);

   // calculate result
   double *result;
   if ((result = malloc(n * sizeof(result[0]))) == NULL) {
       return NULL;
   }
   back(aug, n, n + 1);
   for (unsigned i = 0; i < n; i++) {
       result[i] = aug[i][n];
   }

   // free memory
   matrix_destroy(aug, n);

   return result;
}

/* Function returns the number of the row with the greatest primary element */
unsigned find_greatest(double **matrix, unsigned current, unsigned n) {
    unsigned result = current;
    for (unsigned j = current + 1; j < n; j++) {
        if (matrix[j][current] != 0 && (!matrix[result][current] || fabs(matrix[result][current]) < fabs(matrix[j][current])) ) {
            result = j;
        }
    }
    return result;
}

/* Swap lines */
void swap(double **matrix, unsigned i, unsigned j) {
    double *temp = matrix[i];
    matrix[i] = matrix[j];
    matrix[j] = temp;
}

/* Function subtracts current line from the following */
void subtract(double **matrix, unsigned current, unsigned n, unsigned m) {
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
void forward(double **matrix, unsigned n, unsigned m) {
    for (unsigned i = 0; i < n; i++) {
        // subtract greatest from others
        unsigned greatest = find_greatest(matrix, i, n);
        if (i != greatest) {
            swap(matrix, i, greatest);
        }
        subtract(matrix, i, n, m);
    }
}

/* Function normalizes upper triangular n * m matrix */
void normalize(double **matrix, unsigned n, unsigned m) {
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = i + 1; j < m; j++) {
            matrix[i][j] /= matrix[i][i];
        }
        matrix[i][i] = 1;
    }
}

/* Function performs back substitution */
void back(double **matrix, unsigned n, unsigned m) {
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
    forward(aug, n, n);

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
    forward(aug, n, 2 * n);
    normalize(aug, n, 2 * n);
    back(aug, n, 2 * n);

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
