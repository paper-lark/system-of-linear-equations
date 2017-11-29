#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

/* Constants */
typedef enum {
    SOR, GAUSS, MODIFIED
} Method;

#ifndef OMEGA
    #define OMEGA  1.3333333333333333
#endif
#define PRECISION 0.00000001


/* Prototypes */
extern double **matrix_create(unsigned n, unsigned m);
extern double **matrix_read(unsigned n);
extern double *matrix_read_vector(unsigned n);
extern int matrix_print(double **matrix, unsigned n);
extern int matrix_print_vector(double *vector, unsigned n);
extern void matrix_destroy(double **matrix, unsigned n);
extern double matrix_determinant(double **matrix, unsigned n);
extern double **matrix_inverse(double **matrix, unsigned n);

extern double *matrix_gauss_solve(double **matrix, const double *f, unsigned n);
extern double *matrix_modified_solve(double **matrix, const double *f, unsigned n);
extern double *matrix_iteration_solve(double **a, const double *f, unsigned n, double omega, const double *start, double precision);

extern unsigned functions_init(unsigned task);
extern double **matrix_fill();
extern double *(*matrix_fill_vector)(double x);


/*
 * Entry Point.
 */
int main(int argc, char *argv[]) {
    // read matrix
    unsigned n;
    double **matrix, *f;

    // parse options
    if (argc == 2) {
        // use standard input stream
        scanf("%u", &n);
        if ((matrix = matrix_read(n)) == NULL || (f = matrix_read_vector(n)) == NULL) {
            fprintf(stderr, "> error: reading failed with message: %s\n", strerror(errno));
            exit(1);
        }
    } else if (argc == 4) {
        // use formula
        unsigned long formula = strtoul(argv[2], NULL, 10);
        if (errno != 0 || formula > 4) {
            fprintf(stderr, "> error: invalid option\n");
            exit(2);
        }
        double x = atof(argv[3]);
        n = functions_init((unsigned) formula);
        matrix = matrix_fill();
        f = matrix_fill_vector(x);
    } else {
        fprintf(stderr, "> error: wrong number arguments: %s gauss|mod|sor [formula x]\n",
            argv[0]);
        exit(2);
    }

    // choose method
    Method method;
    double start[n];
    if (strcmp(argv[1], "gauss") == 0) {
        // Gaussian elimination method
        method = GAUSS;
    } else if (strcmp(argv[1], "mod") == 0) {
        // Modified Gaussian elimination method
        method = MODIFIED;
    } else if (strcmp(argv[1], "sor") == 0) {
        // SOR method
        method = SOR;
        // create first approximation vector
        for (unsigned i = 0; i < n; i++) {
            start[i] = 0;
        }
    } else {
        fprintf(stderr, "> error: unknown method\n");
        exit(2);
    }
   
    // print processed matrix
    printf("> matrix:\n");
    matrix_print(matrix, n);
    putchar('\n');

    // calculate determinant
    double determinant = matrix_determinant(matrix, n);
    printf("> determinant:\n %.10g\n", determinant);
    putchar('\n');
    
    // process if posssible if possible
    if (determinant == 0) {
        printf("> inverse matrix:\n does not exist\n\n");
        printf("> solution:\n cannot be found using Gaussian elimination method\n\n");
    } else {
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

        // calculate solution using Gaussian elimination method
        double *solution;
        switch (method) {
            case GAUSS:
                solution = matrix_gauss_solve(matrix, f, n);
                break;
            case MODIFIED:
                solution = matrix_modified_solve(matrix, f, n);
                break;
            case SOR:
                solution = matrix_iteration_solve(matrix, f, n, OMEGA, start, PRECISION);
                break;
            default:
                fprintf(stderr, "> error: unknown method option\nterminationg...\n");
                exit(1);
        }

        if (solution == NULL) {
            fprintf(stderr, "> error: failed to calculate solution: %s\n", strerror(errno));
            exit(1);
        }

        printf("> solution:\n");
        matrix_print_vector(solution, n);
        putchar('\n');

        // free memory
        free(solution);
    }

    // free memory
    free(f);
    matrix_destroy(matrix, n);

    // completion message
    printf("> program sucessfully finished!\n");
}