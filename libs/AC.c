#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<ctype.h>
#include<string.h>
#include<math.h>
#include"hash_table.h"
#include"sol_methods.h"
#include"AC.h"
#include"CXSparse/Sparse/CXSparse/Include/cs.h"


gsl_vector_complex *RightPartAC(struct dataInfo *ListHead,int nodes,int col2,char **nodes_sorted_buffer){
    int k = 0;
    int array_size = nodes+col2;
    struct dataInfo *curr;
    gsl_vector_complex *RightPart = gsl_vector_complex_alloc(array_size);
    for(int i = 0; i < array_size; i++){
        gsl_vector_complex_set(RightPart,i,0.0);
    }
    
    for(curr = ListHead; curr != NULL;curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);
        if(strcmp(curr->type,"Current Source") == 0){
            if(pos != -1){
                gsl_vector_complex_set(RightPart,pos,gsl_vector_complex_get(RightPart,pos)+CMPLX(-creal(curr->ac),-cimag(curr->ac)));
            }
            if(neg != -1){
                gsl_vector_complex_set(RightPart,neg,gsl_vector_complex_get(RightPart,neg)+CMPLX(creal(curr->ac),cimag(curr->ac)));
            }
        }        

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(strcmp(curr->type,"Voltage Source") == 0){
                gsl_vector_complex_set(RightPart,nodes+k,gsl_vector_complex_get(RightPart,nodes+k)+CMPLX(creal(curr->ac),cimag(curr->ac)));
            }
            k++;
        }
    }   
    return RightPart;
}

void SetComplexValue(cs_ci *A, int row, int column, double complex value, int *counter, int nz){
    if(*counter >= nz){
        printf("[ERROR] Sparse matrix overflow.\n");
        exit(0);
    }
    A->i[*counter] = row;
    A->p[*counter] = column;
    A->x[*counter] = value;
    (*counter)++;
}

cs_ci *LeftPartACSparse(struct dataInfo *ListHead,int nodes, int col2,double frequency,char **nodes_sorted_buffer){
    cs_ci *A;
    cs_ci *C;
    int nz = 0;
    int k = 0;
    int counter = 0;

    struct dataInfo *curr;
    
    for(curr = ListHead; curr != NULL;curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);
        if(strcmp(curr->type,"Resistance") == 0 || strcmp(curr->type,"Capacitor") == 0){
            if(pos != -1 && neg != -1){
               nz += 4;
            }
            else
                nz += 1;
        }
        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(pos != -1){
                nz += 2;
            }
            if(neg!= -1){
                nz += 2;
            }
            if(strcmp(curr->type,"Inductor") == 0)
                nz += 1;
            k++;
        }
    }

    A = cs_ci_spalloc(nodes+col2,nodes+col2,nz,1,1);
    A->nz = nz;
    k = 0;

    for(struct dataInfo *curr = ListHead; curr != NULL; curr = curr->next){   
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);
        if(strcmp(curr->type,"Resistance") == 0){
            if(pos != -1 && neg != -1){
                SetComplexValue(A, pos, pos, CMPLX((double)(1.0/curr->value),0), &counter, nz);
                SetComplexValue(A, neg, neg, CMPLX((double)(1.0/curr->value),0), &counter, nz);
                SetComplexValue(A, pos, neg, CMPLX((double)(-1.0/curr->value),0), &counter, nz);
                SetComplexValue(A, neg, pos, CMPLX((double)(-1.0/curr->value),0), &counter, nz);
            }
            else if(pos == -1)
                SetComplexValue(A, neg, neg, CMPLX((double)(1.0/curr->value),0), &counter, nz);
            else
                SetComplexValue(A, pos, pos, CMPLX((double)(1.0/curr->value),0), &counter, nz);
        }
        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(pos != -1){
                SetComplexValue(A, nodes+k, pos, CMPLX(1.0,0), &counter, nz);
                SetComplexValue(A, pos, nodes+k, CMPLX(1.0,0), &counter, nz);
            }
            if(neg!= -1){
                SetComplexValue(A, nodes+k, neg, CMPLX(-1.0,0), &counter, nz);
                SetComplexValue(A, neg, nodes+k, CMPLX(-1.0,0), &counter, nz);
            }
            if(strcmp(curr->type,"Inductor") == 0){
                SetComplexValue(A, nodes+k, nodes+k, CMPLX(0,-(2*M_PI*frequency)*curr->value), &counter, nz);
            }
            k++;
        }
        else if(strcmp(curr->type,"Capacitor") == 0){
            if(pos != -1 && neg != -1){
                SetComplexValue(A, pos, pos, CMPLX(0,(2*M_PI*frequency)*curr->value), &counter, nz);
                SetComplexValue(A, neg, neg, CMPLX(0,(2*M_PI*frequency)*curr->value), &counter, nz);
                SetComplexValue(A, pos, neg, CMPLX(0,-(2*M_PI*frequency)*curr->value), &counter, nz);
                SetComplexValue(A, neg, pos, CMPLX(0,-(2*M_PI*frequency)*curr->value), &counter, nz);
            }
            else if(pos == -1){
                SetComplexValue(A, neg, neg, CMPLX(0,(2*M_PI*frequency)*curr->value), &counter, nz);
            }
            else{
                SetComplexValue(A, pos, pos, CMPLX(0,(2*M_PI*frequency)*curr->value), &counter, nz);
            }
        }
    }
    

    C = cs_ci_compress(A);
    cs_ci_spfree(A);
    cs_ci_dupl(C);
    return C;
}

gsl_matrix_complex *LeftPartDC(struct dataInfo *ListHead,int nodes, int col2,char **nodes_sorted_buffer){
    int array_size = nodes+col2;
	gsl_matrix_complex *LeftPart = gsl_matrix_complex_alloc(array_size,array_size);
    struct dataInfo *curr;
    int k = 0;

    for(int i = 0; i < array_size; i++){
        for(int j = 0; j < array_size; j++){
            gsl_matrix_complex_set(LeftPart,i,j,CMPLX(0,0));
        }
    }

	for(curr = ListHead; curr != NULL;curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);
        if(strcmp(curr->type,"Resistance") == 0){
            if(pos != -1 && neg != -1){
                gsl_matrix_complex_set(LeftPart, pos, pos, gsl_matrix_complex_get(LeftPart,pos,pos)+CMPLX(1.0/curr->value,0));
                gsl_matrix_complex_set(LeftPart, neg, neg, gsl_matrix_complex_get(LeftPart,neg,neg)+CMPLX(1.0/curr->value,0));
                gsl_matrix_complex_set(LeftPart, pos, neg, gsl_matrix_complex_get(LeftPart,pos,neg)-CMPLX(1.0/curr->value,0));
                gsl_matrix_complex_set(LeftPart, neg, pos, gsl_matrix_complex_get(LeftPart,neg,pos)-CMPLX(1.0/curr->value,0));
            }
            else if(pos == -1)
                gsl_matrix_complex_set(LeftPart, neg, neg, gsl_matrix_complex_get(LeftPart,neg,neg)+CMPLX(1.0/curr->value,0));
            else
                gsl_matrix_complex_set(LeftPart, pos, pos, gsl_matrix_complex_get(LeftPart,pos,pos)+CMPLX(1.0/curr->value,0));
        }

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
                if(pos != -1){
                    gsl_matrix_complex_set(LeftPart, nodes+k, pos, gsl_matrix_complex_get(LeftPart,nodes+k,pos)+CMPLX(1.0,0));
                    gsl_matrix_complex_set(LeftPart, pos, nodes+k, gsl_matrix_complex_get(LeftPart,pos,nodes+k)+CMPLX(1.0,0));
                }
                if(neg!= -1){
                    gsl_matrix_complex_set(LeftPart, nodes+k, neg, gsl_matrix_complex_get(LeftPart,nodes+k,neg)-CMPLX(1.0,0));
                    gsl_matrix_complex_set(LeftPart, neg, nodes+k, gsl_matrix_complex_get(LeftPart,neg,nodes+k)-CMPLX(1.0,0));
                }
            k++;
        }
    }

    return LeftPart;
}

void LeftPartAC(struct dataInfo *ListHead,gsl_matrix_complex * LeftPart,int nodes,int col2,double frequency,char **nodes_sorted_buffer){
    struct dataInfo *curr;
    int k = 0;
    
    
    for(curr = ListHead; curr != NULL;curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);
        

        if(strcmp(curr->type,"Capacitor") == 0){
            if(pos != -1 && neg != -1){
                gsl_matrix_complex_set(LeftPart, pos, pos, CMPLX(creal(gsl_matrix_complex_get(LeftPart,pos,pos)),cimag(gsl_matrix_complex_get(LeftPart,pos,pos))+((2*M_PI*frequency)*curr->value)));
                gsl_matrix_complex_set(LeftPart, neg, neg, CMPLX(creal(gsl_matrix_complex_get(LeftPart,neg,neg)),cimag(gsl_matrix_complex_get(LeftPart,neg,neg))+(2*M_PI*frequency)*curr->value));
                gsl_matrix_complex_set(LeftPart, pos, neg, CMPLX(creal(gsl_matrix_complex_get(LeftPart,pos,neg)),cimag(gsl_matrix_complex_get(LeftPart,pos,neg))-(2*M_PI*frequency)*curr->value));
                gsl_matrix_complex_set(LeftPart, neg, pos, CMPLX(creal(gsl_matrix_complex_get(LeftPart,neg,pos)),cimag(gsl_matrix_complex_get(LeftPart,neg,pos))-(2*M_PI*frequency)*curr->value));
                
            }
            else if(pos == -1)
                gsl_matrix_complex_set(LeftPart, neg, neg, CMPLX(creal(gsl_matrix_complex_get(LeftPart,neg,neg)),cimag(gsl_matrix_complex_get(LeftPart,neg,neg))+(2*M_PI*frequency)*curr->value));
            else
                gsl_matrix_complex_set(LeftPart, pos, pos, CMPLX(creal(gsl_matrix_complex_get(LeftPart,pos,pos)),cimag(gsl_matrix_complex_get(LeftPart,pos,pos))+((2*M_PI*frequency)*curr->value)));
        }

        
        else if(strcmp(curr->type,"Voltage Source") == 0)
            k++;

        else if(strcmp(curr->type,"Inductor") == 0){
            gsl_matrix_complex_set(LeftPart, nodes+k, nodes+k, CMPLX(creal(gsl_matrix_complex_get(LeftPart,nodes+k,nodes+k)),cimag(gsl_matrix_complex_get(LeftPart,nodes+k,nodes+k))-(2*M_PI*frequency)*curr->value));
            k++;
        }
    }
}

void ACAnalysis(struct dataInfo *ListHead, int nodes, char **nodes_sorted_buffer,int number_of_actions, char **actions_ptr){
    struct dataInfo *curr;    
    int col2 = 0;
    int points;
    double start_freq,end_freq;
    int sweep = 0;
    bool sparse = false;
    bool iter = false;
    double itol = 1e-3;
    int FileCounter = 0;
    for(int i = 0; i < number_of_actions; i++){
        for(int j = 0; j < strlen(actions_ptr[i]); j++){
            if(actions_ptr[i][j] == 'V'){FileCounter++;}
        }
    }
    
    FILE *out[FileCounter];
    char **NodesToPrint = (char **)malloc(FileCounter*sizeof(char *));
    for(int j = 0; j < FileCounter; j++){
        NodesToPrint[j] = (char *)calloc(50,sizeof(char));
    }
    int k = 0;
    for(int i = 0; i < number_of_actions; i++){
        if(strstr(actions_ptr[i],".OPTIONS")){
            if(strstr(actions_ptr[i],"SPARSE")){sparse=true;}
            if(strstr(actions_ptr[i],"ITER")){iter=true;}
            if(strstr(actions_ptr[i],"ITOL=") != NULL){
                char *start = strstr(actions_ptr[i],"ITOL=");
                char *end = start;
                while(*(end+1)!=' ' && *(end+1)!= '\n' && *(end+1) != '\0'){
                    end++;
                }
                size_t length = end - start + 1;
                char tolerance[length + 1];
                memcpy(tolerance, start + strlen("ITOL=")-1, length);
                tolerance[length] = '\0';
                for(int k = 0; tolerance[k]; k++){
                    tolerance[k] = tolower(tolerance[k]);
                }
                itol = atof(tolerance);
            }

        }
        if(strstr(actions_ptr[i],".AC")){
            strtok(actions_ptr[i]," ");
            sweep = !strcmp(strtok(NULL," "),"LIN")?0:1;
            points = atoi(strtok(NULL," "));
            start_freq = atof(strtok(NULL," "));
            end_freq = atof(strtok(NULL," "));
        }

        if(strstr(actions_ptr[i],".PRINT") || strstr(actions_ptr[i],".PLOT")){
            int number_of_nodes = 0;
            printf("%s\n",actions_ptr[i]);
            for(int j = 0; j < strlen(actions_ptr[i]); j++){
                if(actions_ptr[i][j] == 'V'){
                    number_of_nodes++;
                }
            }   

            
            for(int j = 0; j < strlen(actions_ptr[i]); j++){
                if(actions_ptr[i][j] == '('){
                    int counter = 0;
                    j++;
                    while(actions_ptr[i][j] != ')'){
                        NodesToPrint[k][counter++] = actions_ptr[i][j++];
                    }
                    k++;
                }
            }        
        }         
    }
    
    for(int j = 0; j < FileCounter; j++){
        out[j] = fopen(NodesToPrint[j],"w+");
    }

    for(curr = ListHead; curr != NULL;curr = curr->next){
        if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            col2++;
        }        
    }
          
    double step = (end_freq - start_freq) / (points-1);
    double f;  
    int point_num;

    if(sweep == 0)
        point_num = points;
    else{
        int decades = (int)(log10(end_freq) - log10(start_freq));
        point_num = 10*decades+1;
    }

    double sweep_points[point_num];
    if(sweep == 0){
        for(int i = 0; i < point_num; i++)
            sweep_points[i] = start_freq +i*step;
    }
    else{
        for(int i = 0; i < point_num; i++)
            sweep_points[i] = pow(10,log10(start_freq) + i*(log10(end_freq) - log10(start_freq))/(point_num-1));
    }
    
    
    gsl_matrix_complex *LeftPartInit; 
    gsl_vector_complex *RightPart; 
    gsl_vector_complex *x = NULL;
    double complex *y = NULL;

    if(!sparse){
        LeftPartInit = LeftPartDC(ListHead,nodes, col2,nodes_sorted_buffer);
        RightPart = RightPartAC(ListHead,nodes,col2,nodes_sorted_buffer);
        for(int i = 0; i < point_num; i++){
            f = sweep_points[i];

        
            gsl_matrix_complex *LeftPart = gsl_matrix_complex_alloc(nodes+col2,nodes+col2);
            for(int i = 0; i < nodes+col2; i++){
                for(int j = 0; j < nodes+col2; j++){
                    gsl_matrix_complex_set(LeftPart,i,j,gsl_matrix_complex_get(LeftPartInit,i,j));
                }
            }

            LeftPartAC(ListHead,LeftPart,nodes,col2,f,nodes_sorted_buffer);

            printf("**** Frequency %lf ****\n",f);

            if(!iter){
                x = AC_LU_solution(LeftPart, RightPart, nodes, col2);
                for(int j = 0; j < nodes; j++){
                    complex z = gsl_vector_complex_get(x,j);
                    printf("%s %lf +i%lf\n",nodes_sorted_buffer[j],creal(z),cimag(z));
                }
                printf("\n");
            }
        
            else{
            
                double complex **LeftPartTemp = (double complex **)malloc((nodes+col2)*sizeof(double complex *));
                double complex *RightPartTemp = (double complex *)malloc((nodes+col2)*sizeof(double complex));
                for(int i = 0; i < nodes+col2; i++){
                    LeftPartTemp[i] = (double complex *)malloc((nodes+col2)*sizeof(double complex));
                }
                
                for(int i = 0; i < nodes+col2;i++){
                    for(int j = 0; j < nodes+col2;j++){
                        gsl_complex z = gsl_matrix_complex_get(LeftPart,i,j);
                        LeftPartTemp[i][j] = CMPLX(creal(z),cimag(z));
                    }
                }

                for(int i = 0; i < nodes+col2; i++){    
                    gsl_complex z = gsl_vector_complex_get(RightPart,i);
                    RightPartTemp[i] =  CMPLX(creal(z),cimag(z));
                }
                
                y = AC_Bi_CG(LeftPartTemp,RightPartTemp,nodes,col2,itol);

                for(int j = 0; j < nodes;j++){
                    printf("%s %lf +i%lf\n",nodes_sorted_buffer[j],creal(y[j]),cimag(y[j]));
                }

                for(int j = 0; j < nodes+col2;j++){
                    free(LeftPartTemp[j]);
                }
                free(LeftPartTemp);
                free(RightPartTemp);
            }
        
            for(int i = 0; i < FileCounter; i++){
                fseek(out[i],0,SEEK_END);
                int idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,NodesToPrint[i]);
                if(iter){
                    double mag;
                    if(!sweep)
                        mag = cabs(y[idx]);
                    else
                        mag = 20*log10(cabs(y[idx]));
                    fprintf(out[i],"%lf %lf %lf\n",f,mag,carg(y[idx])*180/M_PI);
                }
                else{
                    gsl_complex z = gsl_vector_complex_get(x,idx);
                    double mag;
                    if(!sweep)
                        mag = cabs(z);
                    else
                        mag = 20*log10(cabs(z));
                    fprintf(out[i],"%lf %lf %lf\n",f,mag,carg(z)*180/M_PI);
                }
            }                
        }
    }

    else{ //sparse
        RightPart = RightPartAC(ListHead,nodes,col2,nodes_sorted_buffer);
        double complex *RightPartTemp = (double complex *)malloc((nodes+col2)*sizeof(double complex));               
        for(int i = 0; i < nodes+col2; i++){    
            gsl_complex z = gsl_vector_complex_get(RightPart,i);
            RightPartTemp[i] =  CMPLX(creal(z),cimag(z));
        }
        for(int i = 0; i < point_num; i++){
            f = sweep_points[i];
            printf("**** Frequency %lf ****\n",f);
            cs_ci *LeftPart = LeftPartACSparse(ListHead,nodes,col2,f,nodes_sorted_buffer);  
            
            if(!iter){
                y = (double complex *) calloc(nodes + col2, sizeof(double complex));
                memcpy(y,RightPartTemp,(nodes + col2)* sizeof(double complex));
                cs_ci_lusol(2,LeftPart,y,1);
            }
            else{
                y = AC_sparse_Bi_CG(LeftPart, RightPartTemp, nodes,col2,itol);
            }
            for(int i = 0; i < FileCounter; i++){
                fseek(out[i],0,SEEK_END);
                int idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,NodesToPrint[i]);
                double mag;
                if(!sweep)
                    mag = cabs(y[idx]);
                else
                    mag = 20*log10(cabs(y[idx]));
                fprintf(out[i],"%lf %lf %lf\n",f,mag,carg(y[idx])*180/M_PI);
            }
            for(int i = 0; i <nodes; i++){
                printf("%s %lf +i%lf\n",nodes_sorted_buffer[i], creal(y[i]),cimag(y[i]));
            }
            free(y);
        }
        free(RightPart);
        free(RightPartTemp);
    }
}