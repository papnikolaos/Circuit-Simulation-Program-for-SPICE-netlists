#include <gsl/gsl_linalg.h> 
#include<stdio.h>
#include <stdbool.h>
#include"csparse.h"
#include<string.h>
#include<stdlib.h>
#include<complex.h>
#include"CXSparse/Sparse/CXSparse/Include/cs.h"

#ifndef __sol_methods_H_
#define __sol_methods_H_ 1
gsl_vector *LU_solution(double **LeftPart, double *RightPart, int nodes, int col2,char **nodes_sorted_buffer);
gsl_vector *Cholesky_solution(double **LeftPart, double *RightPart, int nodes, int col2,char **nodes_sorted_buffer);
gsl_vector *CG(double **LeftPart, double *RightPart, int nodes, int col2,double itol);
gsl_vector *Bi_CG(double **LeftPart, double *RightPart, int nodes, int col2,double itol);
double **matrix_transpose(int m_row, int m_col, double **m);
double *matrix_x_vector(int m_row, int m_col, double **m, double *v);
double *constant_x_vector(int m_row, double *vector, double constant);
double vector_inner_product(int m_row, double *first, double *second);
double norm(int m_row, int m_col, double **m, double *v);
double *inverse_diagonal(int m_row, int m_col, double **matrix);
double *solve_diagonal(int m_row, double *left, double *right);
double *vector_add(int m_row, double *first, double *second);
double *vector_sub(int m_row, double *first, double *second);
double *vector_copy(int m_row, double *initial);
gsl_vector *method(double **LeftPart, double *RightPart, int nodes,char **nodes_sorted_buffer, int col2,double itol,bool iter,bool cholesky);
double sparse_norm(int m_row, int m_col, cs *m, double *v);
double *sparse_CG(cs *LeftPart, double *RightPart, int nodes, int col2,double itol);
double *sparse_Bi_CG(cs *LeftPart, double *RightPart, int nodes, int col2,double itol);
double *method_cs(cs *C,double *RightPart,bool cholesky,bool iter, int n,int col2,double itol);
double **matrix_add(int m_row,int m_col,double **m1,double **m2);
double **matrix_sub(int m_row,int m_col,double **m1,double **m2);
double *extract_gsl_vector(gsl_vector *x,int size);
gsl_vector_complex *AC_LU_solution(gsl_matrix_complex *LeftPart, gsl_vector_complex *RightPart, int nodes, int col2);
double complex *AC_Bi_CG(double complex **LeftPart, double complex *RightPart, int nodes, int col2,double complex itol);
double complex **matrix_transpose_conjugate(int m_row, int m_col, double complex **m);
double complex *solve_diagonal_complex(int m_row, double complex **left, double complex *right);
double complex *solve_conjugate_diagonal_complex(int m_row,double complex **LeftPart,double complex *RightPart);
double complex *inverse_diagonal_conjugate(int m_row, int m_col, double complex **matrix);
double complex vector_inner_product_conjugate(int m_row, double complex *first, double complex *second);
double complex norm_complex(int m_row, int m_col, double complex **m, double complex *v);
double complex *vector_copy_complex(int m_row, double complex *initial);
double complex *vector_add_complex(int m_row, double complex *first, double complex *second);
double complex *vector_sub_complex(int m_row, double complex *first, double complex *second);
double complex **matrix_add_complex(int m_row,int m_col,double complex **m1,double complex **m2);
double complex ** matrix_sub_complex(int m_row,int m_col,double complex **m1,double complex **m2);
double complex *constant_x_vector_complex(int m_row, double complex *vector, double complex constant);
double complex *matrix_x_vector_complex(int m_row, int m_col, double complex **m, double complex *v);
double vector_complex_norm(int size,double complex *v);
double complex *AC_sparse_Bi_CG(cs_ci *LeftPart, double complex*RightPart, int nodes, int col2,double itol);
#endif