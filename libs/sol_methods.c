#include"sol_methods.h"
#include <gsl/gsl_linalg.h> 
#include<stdio.h>
#include"csparse.h"
#include<string.h>
#include <stdbool.h>
#include<stdlib.h>
#include <complex.h>
#include"CXSparse/Sparse/CXSparse/Include/cs.h"

// ********************************************************Standard Implementation******************************************************************************
gsl_vector *LU_solution(double **LeftPart, double *RightPart, int nodes, int col2,char **nodes_sorted_buffer){
    
    int array_size = nodes + col2;
    double left[array_size*array_size];
    int k=0;
    
    for(int i=0;i<array_size;i++){ 
		for(int j=0;j<array_size;j++){ 
			left[k]=LeftPart[i][j]; 
			k++; 
		}
	}
	
    gsl_matrix_view A = gsl_matrix_view_array (left, array_size, array_size);
    gsl_vector_view b = gsl_vector_view_array (RightPart, array_size);
    gsl_vector *x = gsl_vector_alloc (array_size);
    
    int s;
    
    gsl_permutation * p = gsl_permutation_alloc (array_size);
    
    gsl_linalg_LU_decomp (&A.matrix, p, &s);    
    gsl_linalg_LU_solve (&A.matrix, p, &b.vector, x);
    
    gsl_permutation_free (p);

    return x;
}

gsl_vector *Cholesky_solution(double **LeftPart, double *RightPart, int nodes, int col2,char **nodes_sorted_buffer){
    
    int array_size = nodes + col2;
    double left[array_size*array_size];
    int k=0;
    
    gsl_set_error_handler_off(); 

    for(int i=0;i<array_size;i++){ 
		for(int j=0;j<array_size;j++){ 
			left[k]=LeftPart[i][j]; 
			k++; 
		}
	}
    gsl_matrix_view A = gsl_matrix_view_array (left, array_size, array_size);
    gsl_vector_view b = gsl_vector_view_array (RightPart, array_size);
    gsl_vector *x = gsl_vector_alloc (array_size);
    if(gsl_linalg_cholesky_decomp (&A.matrix) != GSL_EDOM){
        gsl_linalg_cholesky_solve (&A.matrix,&b.vector,x);
        
    }
    else{
        printf("[ERROR] Matrix is not positive definite, proceeding with LU decomposition.\n\n");
        x = LU_solution(LeftPart,RightPart,nodes,col2,nodes_sorted_buffer);
    }
    return x;
}

gsl_vector *CG(double **LeftPart, double *RightPart, int nodes, int col2,double itol){
    int array_size = nodes + col2;
    double *Ax, *r, *z, *M, *p, *q;
    double b_norm, r_norm, rho, rho1 = 1, beta, alpha;
    int iter=0;
    
    p = (double *) calloc(array_size, sizeof(double));
    
    double* x = (double *) calloc(array_size, sizeof(double));
   
    
    M = inverse_diagonal(array_size, array_size, LeftPart);

    Ax = (double *) calloc(array_size, sizeof(double));

    r = vector_sub(array_size, RightPart, Ax);
    
    r_norm = b_norm = norm(array_size, array_size, LeftPart, RightPart);
    
    if(isnan(fabs(b_norm)) || b_norm == 0){
      r_norm = b_norm = 1;  
    }
    
    while(((r_norm/b_norm)>itol) && (iter<array_size)){
        iter++;
        
        z = solve_diagonal(array_size, M, r);
        rho = vector_inner_product(array_size, r, z);
        
        if(iter == 1){
           p = vector_copy(array_size, z);  
        }
        else{
           beta = rho/rho1;
           p = vector_add(array_size, z, constant_x_vector(array_size, p, beta));
        }
        
        rho1 = rho;
        q = matrix_x_vector(array_size, array_size, LeftPart, p);
        
        alpha = rho/vector_inner_product(array_size, p, q);
        x = vector_add(array_size, x, constant_x_vector(array_size, p, alpha));
        r = vector_sub(array_size, r, constant_x_vector(array_size, q, alpha));
        r_norm = norm(array_size,array_size,LeftPart,r);
    }
    
    
    gsl_vector *res = gsl_vector_alloc (array_size);
    for(int i=0;i<array_size;i++){
        gsl_vector_set(res, i, x[i]);
    }
    
    return res;
}
gsl_vector *Bi_CG(double **LeftPart, double *RightPart, int nodes, int col2,double itol){
    int array_size = nodes + col2;
    double *Ax, *r, *z, *M, *p, *q, *_r, *_z, *_q, *_p;
    double b_norm, r_norm, rho, rho1 = 1, beta, alpha, omega;
    int iter=0;
    double **transpose_left = matrix_transpose(array_size, array_size, LeftPart);
    
    p = (double *) calloc(array_size, sizeof(double));

    double* x = (double *) calloc(array_size, sizeof(double));
    
    M = inverse_diagonal(array_size, array_size, LeftPart);
    
    Ax = (double *) calloc(array_size, sizeof(double));
    r = vector_sub(array_size, RightPart, Ax);
    _r = vector_sub(array_size, RightPart, Ax);
    
    r_norm = b_norm = norm(array_size, array_size, LeftPart, RightPart);
    
    if(isnan(fabs(b_norm))|| b_norm == 0){
      r_norm = b_norm = 1;  
    }
    
    
    while(((r_norm/b_norm)>itol) && (iter<array_size)){
        iter++;
        
        z = solve_diagonal(array_size, M, r);
        _z = solve_diagonal(array_size, M, _r);
        
        rho = vector_inner_product(array_size, _r, z);
        
        if(isnan(fabs(rho))){
           break; 
        }
        if(iter == 1){
           p = vector_copy(array_size, z); 
           _p = vector_copy(array_size, _z); 
        }
        else{
           beta = rho/rho1;
           p = vector_add(array_size, z, constant_x_vector(array_size, p, beta));
           _p = vector_add(array_size, _z, constant_x_vector(array_size, _p, beta));
        }
        
        rho1 = rho;
        q = matrix_x_vector(array_size, array_size, LeftPart, p);
        _q = matrix_x_vector(array_size, array_size, transpose_left, _p);
        
        omega = vector_inner_product(array_size, _p, q);
        if(isnan(fabs(omega))){
           break; 
        }
        alpha = rho/omega;
        x = vector_add(array_size, x, constant_x_vector(array_size, p, alpha));
        r = vector_sub(array_size, r, constant_x_vector(array_size, q, alpha));
        _r = vector_sub(array_size, _r, constant_x_vector(array_size, _q, alpha));
        r_norm = norm(array_size,array_size,LeftPart,r);
        
    }
    
    gsl_vector *res = gsl_vector_alloc (array_size);
    for(int i=0;i<array_size;i++){
        gsl_vector_set(res, i, x[i]);
    }
    
    return res;
}

double **matrix_add(int m_row,int m_col,double **m1,double **m2){
    double **res = calloc(m_row,sizeof(double **));
    for(int i = 0; i < m_row; i++){
        res[i] = calloc(m_col,sizeof(double));
    }

    for(int i = 0; i < m_row; i++){
        for(int j = 0; j < m_col; j++){
            res[i][j] = m1[i][j] + m2[i][j];
        }
    }

    return res;
}

double **matrix_sub(int m_row,int m_col,double **m1,double **m2){
    double **res = calloc(m_row,sizeof(double **));
    for(int i = 0; i < m_row; i++){
        res[i] = calloc(m_col,sizeof(double));
    }

    for(int i = 0; i < m_row; i++){
        for(int j = 0; j < m_col; j++){
            res[i][j] = m1[i][j] - m2[i][j];
        }
    }

    return res;
}

double **matrix_transpose(int m_row, int m_col, double **m){

    double **transpose = (double **) calloc(m_col, sizeof(double *));
    for(int i = 0; i < m_col; i++)
        transpose[i] = (double *) calloc(m_row, sizeof(double));

    for(int i = 0; i < m_col; i++)
       for(int j = 0; j < m_row; j++){
          transpose[j][i] = m[i][j];
    }
    
    return(transpose);
}

double *matrix_x_vector(int m_row, int m_col, double **m, double *v){
  int i, j; 
  double *res = (double *) calloc(m_row, sizeof(double));

  for(i=0; i<m_row; i++){
     for(j=0; j<m_col; j++){
      res[i] += m[i][j] * v[j];
    }
  }
  
  return(res);

}

double *constant_x_vector(int m_row, double *vector, double constant){
  int i; 
  double *res = (double *) calloc(m_row, sizeof(double));

  for(i=0; i<m_row; i++){
      res[i] = constant * vector[i];
  }
  
  return(res);

}

double vector_inner_product(int m_row, double *first, double *second){

  int i; 
  double res;
  
  for(i=0; i<m_row; i++){
     res += first[i] * second[i];
  }
  
  return(res);
}

double *inverse_diagonal(int m_row, int m_col, double **matrix){
    
    int i, j; 
    double *res = (double *) calloc(m_row, sizeof(double));

    int k =0;
    for(i=0; i<m_row; i++){
        for(j=0; j<m_col; j++){
            if(i == j){
               if(matrix[i][j]==0){
                 res[k]=1;
               }
               else{
                  res[k] = 1/matrix[i][j]; 
               }
            k++;
            }
        }
    }
  
  return(res);
}

double *solve_diagonal(int m_row, double *left, double *right){    
    
    int i; 
    double *res = (double *) calloc(m_row, sizeof(double));
    
    for(i=0; i<m_row; i++){
        res[i] = right[i]*left[i];
    }
    
    return(res);
    
}

double norm(int m_row, int m_col, double **m, double *v){
    
   double* between =  matrix_x_vector(m_row, m_col, m, v);
   double between2 = vector_inner_product(m_row, v, between);
   
   double absolute = fabs(between2);
   double res = sqrt(absolute);
   
   return(res);
}

double *vector_copy(int m_row, double *initial){

  int i; 
  double *res = (double *) calloc(m_row, sizeof(double));
  
  for(i=0; i<m_row; i++){
     res[i] = initial[i];
  }
  
  return(res);
}

double *vector_add(int m_row, double *first, double *second){

  int i; 
  double *res = (double *) calloc(m_row, sizeof(double));
  
  for(i=0; i<m_row; i++){
     res[i] = first[i] + second[i];
  }
  
  return(res);
}

double *vector_sub(int m_row, double *first, double *second){

  int i; 
  double *res = (double *) calloc(m_row, sizeof(double));
  
  for(i=0; i<m_row; i++){
     res[i] = first[i] - second[i];
  }
  
  return(res);
}

gsl_vector *method(double **LeftPart, double *RightPart, int nodes,char **nodes_sorted_buffer, int col2,double itol,bool iter,bool cholesky){
  if(!iter && !cholesky){return LU_solution(LeftPart,RightPart,nodes,col2,nodes_sorted_buffer);}
  else if(!iter && cholesky){return Cholesky_solution(LeftPart,RightPart,nodes,col2,nodes_sorted_buffer);}
  else if(iter && cholesky){return CG(LeftPart,RightPart,nodes,col2,itol);}
  else{return Bi_CG(LeftPart,RightPart,nodes,col2,itol);}
}

double *extract_gsl_vector(gsl_vector *x,int size){
    double *ext = calloc(size,sizeof(double));
    for(int i = 0; i < size; i++){
        ext[i] = gsl_vector_get(x,i);
    }

    return ext;
}

//************************************************************************************************************************************************************


// **************************************************************************** SPARSE ***********************************************************************

double sparse_norm(int m_row, int m_col, cs *m, double *v){
    
    double* between = (double *) calloc(m_row, sizeof(double));
    cs_gaxpy(m,v,between);
    double between2 = vector_inner_product(m_row, v, between);
   
    double absolute = fabs(between2);
    double res = sqrt(absolute);
    free(between);
    return(res);
}

double *sparse_CG(cs *LeftPart, double *RightPart, int nodes, int col2,double itol){ 


    int array_size = nodes + col2;
    double *Ax, *r, *z, *M, *p, *q;
    double b_norm, r_norm, rho, rho1 = 1, beta, alpha;
    int iter=0;
    q = (double *) calloc(array_size, sizeof(double));
    M = (double *) calloc(array_size, sizeof(double));
    
    double* x = (double *) calloc(array_size, sizeof(double));
   
    for(int j = 0 ; j < array_size ; j++){
        for (int k = LeftPart->p[j] ; k < LeftPart->p[j+1] ; k++){
            if (LeftPart->i[k] == j){
                M[j] += LeftPart->x[k];
            }
        }
    }
    
    for(int k=0;k<array_size;k++){
       if(M[k]==0){
         M[k] = 1;  
       }
       else{
         M[k] = 1/M[k];    
       }
    }

    Ax = (double *) calloc(array_size, sizeof(double));

    r = vector_sub(array_size, RightPart, Ax);
    
    r_norm = b_norm = sparse_norm(array_size, array_size, LeftPart, RightPart);
    
    if(isnan(fabs(b_norm)) || b_norm == 0){
      r_norm = b_norm = 1;  
    }
    
    while(((r_norm/b_norm)>itol) && (iter<array_size)){
        iter++;
        
        z = solve_diagonal(array_size, M, r);
        rho = vector_inner_product(array_size, r, z);
        
        if(iter == 1){
           p = vector_copy(array_size, z);  
        }
        else{
           beta = rho/rho1;
           double *temp = constant_x_vector(array_size, p, beta);
           free(p);
           p = vector_add(array_size, z, temp);
           free(temp);
        }
        
        rho1 = rho;
        for(int i=0;i<array_size;i++){
            q[i]=0;
        }
        cs_gaxpy(LeftPart,p,q);
        
        alpha = rho/vector_inner_product(array_size, p, q);
        double *temp1 = constant_x_vector(array_size, p, alpha);
        double *temp2 = constant_x_vector(array_size, q, alpha);
        double *temp3 = x;
        double *temp4 = r;
        x = vector_add(array_size, temp3, temp1);
        r = vector_sub(array_size, temp4, temp2);
        free(temp1);
        free(temp2);
        free(temp3);
        free(temp4);
        r_norm = sparse_norm(array_size,array_size,LeftPart,r);
        free(z);
    }

    free(Ax);
    free(r);
    //free(z);
    free(M);
    free(p);
    free(q);
  
    return x;
}

double *sparse_Bi_CG(cs *LeftPart, double *RightPart, int nodes, int col2,double itol){
    int array_size = nodes + col2;
    double *Ax, *r, *z, *M, *p, *q, *_r, *_z, *_q, *_p;
    double b_norm, r_norm, rho, rho1 = 1, beta, alpha, omega;
    int iter=0;
    
    //p = (double *) calloc(array_size, sizeof(double));
    //_p = (double *) calloc(array_size, sizeof(double));
    q = (double *) calloc(array_size, sizeof(double));
    _q = (double *) calloc(array_size, sizeof(double));
    M = (double *) calloc(array_size, sizeof(double));

    double* x = (double *) calloc(array_size, sizeof(double));
    
    for(int j = 0 ; j < array_size ; j++){
        for (int k = LeftPart->p[j] ; k < LeftPart->p[j+1] ; k++){
            if (LeftPart->i[k] == j){
                M[j] += LeftPart->x[k];
            }
        }
    }
    
    for(int k=0;k<array_size;k++){
       if(M[k]==0){
         M[k] = 1;  
       }
       else{
         M[k] = 1/M[k];    
       }
    }
    
    Ax = (double *) calloc(array_size, sizeof(double));
    r = vector_sub(array_size, RightPart, Ax);
    _r = vector_sub(array_size, RightPart, Ax);
    
    r_norm = b_norm = sparse_norm(array_size, array_size, LeftPart, RightPart);
    
    if(isnan(fabs(b_norm))|| b_norm == 0){
      r_norm = b_norm = 1;  
    }
    
    
    while(((r_norm/b_norm)>itol) && (iter<array_size)){
        iter++;
        
        z = solve_diagonal(array_size, M, r);
        _z = solve_diagonal(array_size, M, _r);
        rho = vector_inner_product(array_size, _r, z);
        
        if(isnan(fabs(rho))){
           break; 
        }
        if(iter == 1){
           p = vector_copy(array_size, z); 
           _p = vector_copy(array_size, _z); 
        }
        else{
           beta = rho/rho1;
           double *temp1 = constant_x_vector(array_size, p, beta);
           double *temp2 = constant_x_vector(array_size, _p, beta);
           free(p);
           free(_p);
           p = vector_add(array_size, z, temp1);
           _p = vector_add(array_size, _z, temp2);
           free(temp1);
           free(temp2);
        }
        
        rho1 = rho;
        for(int i=0;i<array_size;i++){
            q[i]=0;
        }
        cs_gaxpy(LeftPart,p,q);
        for (int i = 0 ; i < array_size ; i++){
            _q[i] = 0.0;
            for (int j = LeftPart->p[i] ; j < LeftPart->p[i+1] ; j++){
                _q[i] = _q[i] + LeftPart->x[j] * _p[LeftPart->i[j]];
            }
        }
        
        omega = vector_inner_product(array_size, _p, q);
        if(isnan(fabs(omega))){
           break; 
        }
        alpha = rho/omega;
        double *temp1 = constant_x_vector(array_size, p, alpha);
        double *temp2 = constant_x_vector(array_size, q, alpha);
        double *temp3 = constant_x_vector(array_size, _q, alpha);
        double *temp4 = x;
        x = vector_add(array_size, temp4, temp1);
        double *temp5 = r;
        double *temp6 = _r;
        r = vector_sub(array_size, temp5, temp2);
        _r = vector_sub(array_size, temp6, temp3);
        r_norm = sparse_norm(array_size,array_size,LeftPart,r);
        free(temp1);
        free(temp2);
        free(temp3);
        free(temp4);
        free(temp5);
        free(temp6);
        free(z);
        free(_z);
    }
    
    free(p);
    free(_p);
    free(r);
    free(_r);
    free(q);
    free(_q);
    free(M);
    free(Ax);
    

    return x;
}

double *method_cs(cs *C,double *RightPart,bool cholesky,bool iter, int n,int col2,double itol){

  double *b;
  if(!cholesky && !iter){
        b = (double *) calloc(n + col2, sizeof(double));
        memcpy(b,RightPart,(n + col2)* sizeof(double));
        int ok = cs_lusol(2,C,b,1);
        if(!ok){
            printf("Cannot solve using LU!\n");
            exit(3);
        }
    }
    else if(cholesky && !iter){
        b = (double *) calloc(n + col2, sizeof(double));
        memcpy(b,RightPart,(n + col2)* sizeof(double));
        int ok = cs_lusol(1,C,b,1);
        if(!ok){
            printf("Cannot solve with Cholesky.\nTrying with LU!\n");
            int ok1 = cs_lusol(2,C,b,1);
            if(!ok1){
                printf("Cannot solve using LU!\n");
                exit(3);
            }
        }
    }

    else if(!cholesky && iter){
        b = sparse_Bi_CG(C,RightPart,n,col2,itol);
    }

    else{
        b = sparse_CG(C, RightPart,n,col2,itol);
    }

    return b;
}
// ********************************************************AC Implementation******************************************************************************

gsl_vector_complex *AC_LU_solution(gsl_matrix_complex *LeftPart, gsl_vector_complex *RightPart, int nodes, int col2){
    
    int array_size = nodes + col2;
    int signum = 0;
    gsl_permutation * p = gsl_permutation_alloc (array_size);
    gsl_vector_complex *x = gsl_vector_complex_alloc (array_size);
    gsl_linalg_complex_LU_decomp(LeftPart, p, &signum);
    gsl_linalg_complex_LU_solve(LeftPart,p,RightPart,x);
    
    
    gsl_permutation_free (p);

    return x;
}

double complex *AC_Bi_CG(double complex **LeftPart, double complex *RightPart, int nodes, int col2,double complex itol){

    int array_size = nodes + col2;
    double complex *Ax, *r, *z, *p, *q, *_r, *_z, *_q, *_p;
    double b_norm, r_norm;
    double complex rho, rho1 = 1, beta, alpha, omega;
    int iter=0;
    double complex **conj_transpose_left = matrix_transpose_conjugate(array_size, array_size, LeftPart);
    
    p = (double complex *) calloc(array_size, sizeof(double complex));

    double complex * x = (double complex *) calloc(array_size, sizeof(double complex));
    
    Ax = (double complex *) calloc(array_size, sizeof(double complex));

    r = vector_sub_complex(array_size, RightPart, Ax);
    _r = vector_sub_complex(array_size, RightPart, Ax);

    r_norm = vector_complex_norm(array_size,r);
    b_norm = r_norm;

    if(isnan(fabs(b_norm))|| b_norm == 0){
      r_norm = b_norm = 1;  
    }

    while((r_norm/b_norm>cabs(itol)) && (iter<array_size)){
        iter = iter + 1;

        z = solve_diagonal_complex(array_size,LeftPart,r);
        _z = solve_conjugate_diagonal_complex(array_size,LeftPart,_r);

        rho = vector_inner_product_conjugate(array_size, _r, z);

        if(isnan(fabs(rho))){
            break; 
        }
        if(iter == 1){
            p = vector_copy_complex(array_size, z); 
            _p = vector_copy_complex(array_size, _z); 
        }
        else{
            beta = rho/rho1;
            double complex *temp1 = constant_x_vector_complex(array_size, p, beta);
            double complex *temp2 = constant_x_vector_complex(array_size, _p, conj(beta));
            p = vector_add_complex(array_size, z, temp1);
            _p = vector_add_complex(array_size, _z, temp2);
            free(temp1);
            free(temp2);
        }
        rho1 = rho;
        q = matrix_x_vector_complex(array_size, array_size, LeftPart, p);
        _q = matrix_x_vector_complex(array_size, array_size, conj_transpose_left, _p);

        omega = vector_inner_product_conjugate(array_size, _p, q);

        if(isnan(fabs(omega))){
           break; 
        }

        alpha = rho/omega;
        double complex *temp3 = constant_x_vector_complex(array_size, p, alpha);
        double complex *temp4 = x;
        x = vector_add_complex(array_size,x,temp3);
        free(temp3);
        free(temp4);
        temp3 = constant_x_vector_complex(array_size, q, alpha);
        temp4 = r;
        r = vector_sub_complex(array_size,r,temp3);

        r_norm = vector_complex_norm(array_size,r);


        free(temp3);
        free(temp4);
        temp3 = constant_x_vector_complex(array_size, _q, conj(alpha)); 
        temp4 = _r;
        _r = vector_sub_complex(array_size,_r,temp3);
        free(temp3);
        free(temp4);

        free(z);
        free(_z);

    }
    free(p);
    free(_p);
    free(r);
    free(_r);
    free(q);
    free(_q);
    free(Ax);

    return x;
}

double vector_complex_norm(int size,double complex *v){
    double complex res = CMPLX(0,0);
    
    for(int i = 0; i < size; i++){
        res += CMPLX(creal(v[i]),cimag(v[i]))*CMPLX(creal(v[i]),-cimag(v[i]));
    }
    
    return sqrt(creal(res));
}


double complex *solve_diagonal_complex(int m_row,double complex **LeftPart,double complex *RightPart){
    double complex *res = (double complex *) calloc(m_row, sizeof(double complex));

    for(int i = 0; i < m_row; i++){
        if(LeftPart[i][i] != CMPLX(0,0))
            res[i] = RightPart[i]/LeftPart[i][i];
        else
            res[i] = RightPart[i];
    }

    return res;
}

double complex *solve_conjugate_diagonal_complex(int m_row,double complex **LeftPart,double complex *RightPart){
    double complex *res = (double complex *) calloc(m_row, sizeof(double complex));

    for(int i = 0; i < m_row; i++){
        if(LeftPart[i][i] != CMPLX(0,0))
            res[i] = RightPart[i]/conj(LeftPart[i][i]);
        else
            res[i] = RightPart[i];
    }

    return res;
}

double complex *matrix_x_vector_complex(int m_row, int m_col, double complex **m, double complex *v){
  int i, j; 
  double complex *res = (double complex *) calloc(m_row, sizeof(double complex));

  for(i=0; i<m_row; i++){
     for(j=0; j<m_col; j++){
      res[i] += m[i][j] * v[j];
    }
  }
  
  return(res);

}

double complex *constant_x_vector_complex(int m_row, double complex *vector, double complex constant){
  int i; 
  double complex *res = (double complex *) calloc(m_row, sizeof(double complex));

  for(i=0; i<m_row; i++){
      res[i] = constant * vector[i];
  }
  
  return(res);

}

double complex **matrix_add_complex(int m_row,int m_col,double complex **m1,double complex **m2){
    double complex **res = calloc(m_row,sizeof(double complex **));
    for(int i = 0; i < m_row; i++){
        res[i] = calloc(m_col,sizeof(double));
    }

    for(int i = 0; i < m_row; i++){
        for(int j = 0; j < m_col; j++){
            res[i][j] = m1[i][j] + m2[i][j];
        }
    }

    return res;
}

double complex ** matrix_sub_complex(int m_row,int m_col,double complex **m1,double complex **m2){
    double complex**res = calloc(m_row,sizeof(double complex **));
    for(int i = 0; i < m_row; i++){
        res[i] = calloc(m_col,sizeof(double complex));
    }

    for(int i = 0; i < m_row; i++){
        for(int j = 0; j < m_col; j++){
            res[i][j] = m1[i][j] - m2[i][j];
        }
    }

    return res;
}

double complex **matrix_transpose_conjugate(int m_row, int m_col, double complex **m){

    double complex **transpose = (double complex **) calloc(m_col, sizeof(double complex *));
    for(int i = 0; i < m_col; i++)
        transpose[i] = (double complex *) calloc(m_row, sizeof(double complex));

    for(int i = 0; i < m_col; i++)
       for(int j = 0; j < m_row; j++){
          transpose[j][i] = conj(m[i][j]);
    }
    
    return(transpose);
}


double complex vector_inner_product_conjugate(int m_row, double complex *first, double complex *second){

  int i; 
  double complex res = CMPLX(0,0);
  
  for(i=0; i<m_row; i++){
     res += conj(first[i]) * second[i];
  }
  
  return(res);
}

double complex *vector_copy_complex(int m_row, double complex *initial){

  int i; 
  double complex *res = (double complex *) calloc(m_row, sizeof(double complex));
  
  for(i=0; i<m_row; i++){
     res[i] = initial[i];
  }
  
  return(res);
}

double complex *vector_add_complex(int m_row, double complex *first, double complex *second){

  int i; 
  double complex *res = (double complex *) calloc(m_row, sizeof(double complex));
  
  for(i=0; i<m_row; i++){
     res[i] = first[i] + second[i];
  }
  
  return(res);
}

double complex *vector_sub_complex(int m_row, double complex *first, double complex *second){

  int i; 
  double complex *res = (double complex *) calloc(m_row, sizeof (double complex));
  
  for(i=0; i<m_row; i++){
     res[i] = first[i] - second[i];
  }
  
  return(res);
}


double complex *AC_sparse_Bi_CG(cs_ci *LeftPart, double complex* RightPart, int nodes, int col2,double itol){
    int array_size = nodes + col2;
    double complex *Ax, *r, *z, *p, *q, *_r, *_z, *_q, *_p;
    double b_norm, r_norm;
    double complex rho, rho1 = 1, beta, alpha, omega;
    int iter=0;

    p = (double complex *) calloc(array_size, sizeof(double complex));

    double complex * x = (double complex *) calloc(array_size, sizeof(double complex));
    
    Ax = (double complex *) calloc(array_size, sizeof(double complex));

    r = vector_sub_complex(array_size, RightPart, Ax);
    _r = vector_sub_complex(array_size, RightPart, Ax);

    r_norm = vector_complex_norm(array_size,r);
    b_norm = r_norm;

    if(isnan(fabs(b_norm))|| b_norm == 0){
      r_norm = b_norm = 1;  
    }

    while((r_norm/b_norm>cabs(itol)) && (iter<array_size)){
        iter = iter + 1;

        z = (double complex *)calloc(array_size,sizeof(double complex));
        _z = (double complex *)calloc(array_size,sizeof(double complex));
        double complex precond_val;

        for (int j = 0 ; j < nodes+col2 ; j++){
            precond_val = CMPLX(0,0);
            for (int p = LeftPart->p[j] ; p < LeftPart->p[j+1] ; p++){
                if (LeftPart->i[p] == j)
                    precond_val = LeftPart->x[p];
            }        
            if(precond_val != CMPLX(0,0))
                z[j] = r[j]/precond_val;
            else
                z[j] = r[j]; 
        }

        for (int j = 0 ; j < nodes+col2 ; j++){
            precond_val = CMPLX(0,0);
            for (int p = LeftPart->p[j] ; p < LeftPart->p[j+1] ; p++){
                if (LeftPart->i[p] == j)
                    precond_val = conj(LeftPart->x[p]);
            }
            if(precond_val != CMPLX(0,0))
                _z[j] = _r[j]/precond_val;
            else
                _z[j] = _r[j];
        }   


        rho = vector_inner_product_conjugate(array_size, _r, z);

        if(isnan(fabs(rho))){
            break; 
        }
        if(iter == 1){
            p = vector_copy_complex(array_size, z); 
            _p = vector_copy_complex(array_size, _z); 
        }
        else{
            beta = rho/rho1;
            double complex *temp1 = constant_x_vector_complex(array_size, p, beta);
            double complex *temp2 = constant_x_vector_complex(array_size, _p, conj(beta));
            p = vector_add_complex(array_size, z, temp1);
            _p = vector_add_complex(array_size, _z, temp2);
            free(temp1);
            free(temp2);
        }
        rho1 = rho;

        q = (double complex *)calloc(array_size,sizeof(double complex));
        _q = (double complex *)calloc(array_size,sizeof(double complex));

        for (int j = 0 ; j < array_size ; j++){
            for (int idx = LeftPart->p[j] ; idx < LeftPart->p[j+1] ; idx++)
                q[LeftPart->i[idx]] = q[LeftPart->i[idx]] + LeftPart->x[idx] * p[j];
        }

        for (int i = 0 ; i < array_size ; i++){
            for (int idx = LeftPart->p[i] ; idx < LeftPart->p[i+1] ; idx++)
                _q[i] = _q[i] + conj(LeftPart->x[idx]) * _p[LeftPart->i[idx]];
        }

        omega = vector_inner_product_conjugate(array_size, _p, q);

        if(isnan(fabs(omega))){
           break; 
        }

        alpha = rho/omega;
        double complex *temp3 = constant_x_vector_complex(array_size, p, alpha);
        double complex *temp4 = x;
        x = vector_add_complex(array_size,x,temp3);
        free(temp3);
        free(temp4);
        temp3 = constant_x_vector_complex(array_size, q, alpha);
        temp4 = r;
        r = vector_sub_complex(array_size,r,temp3);

        r_norm = vector_complex_norm(array_size,r);


        free(temp3);
        free(temp4);
        temp3 = constant_x_vector_complex(array_size, _q, conj(alpha)); 
        temp4 = _r;
        _r = vector_sub_complex(array_size,_r,temp3);
        free(temp3);
        free(temp4);

        free(z);
        free(_z);

    }
    free(p);
    free(_p);
    free(r);
    free(_r);
    free(q);
    free(_q);
    free(Ax);

    return x;
}
//************************************************************************************************************************************************************
