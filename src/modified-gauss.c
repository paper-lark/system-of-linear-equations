#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

extern double **matrix_create(unsigned n, unsigned m);
extern void matrix_destroy(double **matrix, unsigned n);
extern void matrix_swap_rows(double **matrix, unsigned i, unsigned j);

void matrix_modified_forward(double **matrix, unsigned *transform, unsigned n, unsigned m);
void matrix_modified_back(double **matrix, const unsigned *transform, unsigned n, unsigned m); 
void matrix_modified_normalize(double **matrix, const unsigned *transform, unsigned n, unsigned m);

typedef struct {
	unsigned i;
	unsigned j;
} Point;

/*
 * Function solves matrix equation Ax = f using modified Gaussian Elimination method.
 */
double *matrix_modified_solve(double **matrix, const double *f, unsigned n) {
   // create augmented matrix and a transformation vector
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
   unsigned transform[n];
   for (unsigned i = 0; i < n; i++) {
	   transform[i] = i;
   }

   // perform forward elimination and normalization
   matrix_modified_forward(aug, transform, n, n + 1);

   // calculate result
   double *result;
   if ((result = malloc(n * sizeof(result[0]))) == NULL) {
       matrix_destroy(aug, n);
       return NULL;
   }
   for (unsigned i = 0; i < n; i++) {
       // find effective x
       unsigned col = 0;
       for (unsigned j = 0; j < n; j++) {
           if (transform[j] == i) {
               col = j;
               break;
           }
       }
       // write out the result
       result[col] = aug[i][n];
   }

   // free memory
   matrix_destroy(aug, n);
   return result;
}

/* Function returns the number of the row and column of the greatest element */
Point matrix_modified_find_greatest(double **matrix, const unsigned *transform,
		unsigned current, unsigned n) {
    Point result = {current, 0};
    bool isFound = false;
	for (unsigned i = current; i < n; i++) {
    	for (unsigned j = 0; j < n; j++) {
			if (transform[j] < current) {
				continue;
			}
    	    if (matrix[i][j] != 0 && 
					(!isFound || fabs(matrix[result.i][result.j]) < fabs(matrix[i][j])) ) {
    	        result.i = i;
				result.j = j;
                isFound = true;
    	    }
    	}
	}
    return result;
}

/* Function performs Forward Elimination of the augmented matrix n * m */
void matrix_modified_forward(double **matrix, unsigned *transform, unsigned n, unsigned m) {
    for (unsigned i = 0; i < n; i++) {
        // find greatest from others and put it to front
        Point greatest = matrix_modified_find_greatest(matrix, transform, i, n);
        if (i != greatest.i) {
            matrix_swap_rows(matrix, i, greatest.i);
        }
        if (i != transform[greatest.j]) {
            unsigned col = 0;
            for (unsigned j = 0; j < n; j++) {
                if (transform[j] == i) {
                    col = j;
                    break;
                }
            }
		    unsigned temp = transform[greatest.j];
		    transform[greatest.j] = transform[col];
            transform[col] = temp;
        }

		// divide current row
		for (unsigned j = 0; j < m; j++) {
            if (j == greatest.j) {
                continue;
            }
			matrix[i][j] /= matrix[i][greatest.j];
		}
        matrix[i][greatest.j] = 1;
		
        // reduce other rows
        const unsigned col = greatest.j;
        for (unsigned row = 0; row < n; row++) {
            if (row == i) {
                continue;
            }
		    double multiplier = matrix[row][col];
            if (multiplier != 0) {
                for (unsigned j = 0; j < m; j++) {
                    matrix[row][j] -= multiplier * matrix[i][j];
                }
            }
        }
    }
}