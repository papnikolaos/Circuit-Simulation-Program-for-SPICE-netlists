#include"hash_table.h"
#include"sol_methods.h"
#include"csparse.h"
#include"transient.h"
#include"AC.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <stdbool.h>
#include <gsl/gsl_linalg.h>
#include <ctype.h>
#include <complex.h>

double **create_left_part(struct dataInfo *ListHead,int col2){
    int nodes;
    char **nodes_sorted_buffer = NodesIdentification(ListHead, &nodes);
	double **LeftPart = (double **) calloc(nodes + col2, sizeof(double *));
    struct dataInfo *curr;
    int k = 0;

    //Allocate Memory
    for(int i = 0; i < nodes + col2; i++){
        LeftPart[i] = (double *) calloc(nodes + col2, sizeof(double));
    }

    //Create Left Part
	for(curr = ListHead; curr != NULL;curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);
        if(strcmp(curr->type,"Resistance") == 0){
            if(pos != -1 && neg != -1){
                LeftPart[pos][pos] += 1.0/curr->value;
                LeftPart[neg][neg] += 1.0/curr->value;
                LeftPart[pos][neg] -= 1.0/curr->value;
                LeftPart[neg][pos] -= 1.0/curr->value;
            }
            else if(pos == -1)
                LeftPart[neg][neg] += 1.0/curr->value;
            else
                LeftPart[pos][pos] += 1.0/curr->value;
        }

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(pos != -1){
                LeftPart[nodes + k][pos] += 1;
                LeftPart[pos][nodes + k] += 1;
            }
            if(neg!= -1){
                LeftPart[nodes + k][neg] -= 1;
                LeftPart[neg][nodes + k] -= 1;
            }
            k++;
        }
    }

    return LeftPart;
}

double *create_right_part(struct dataInfo *ListHead,int col2){
    int nodes;
    int k = 0;
    struct dataInfo *curr;
    char **nodes_sorted_buffer = NodesIdentification(ListHead, &nodes);
    double *RightPart = (double *) calloc(nodes + col2, sizeof(double));

    for(curr = ListHead; curr != NULL;curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);
        if(strcmp(curr->type,"Resistance") == 0){
            continue;
        }
        else if(strcmp(curr->type,"Current Source") == 0){
            if(pos != -1)
                RightPart[pos] -=  curr->value;
            if(neg != -1)
                RightPart[neg] += curr->value;
        }

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(strcmp(curr->type,"Voltage Source") == 0)
                RightPart[nodes + k] += curr->value;
            k++;
        }
    }

    return RightPart;
}

void SetValue(cs *A, int row, int column, double value, int *counter, int nz){
    if(*counter >= nz){
        printf("[ERROR] Sparse matrix overflow.\n");
        exit(0);
    }
    A->i[*counter] = row;
    A->p[*counter] = column;
    A->x[*counter] = value;
    (*counter)++;
}

cs *create_left_part_sparse(struct dataInfo *ListHead,int col2){
    int n, nz = 0;
    int counter = 0, k = 0;
    cs *A;
    cs *C;

    char **nodes_sorted_buffer = NodesIdentification(ListHead, &n);

    for(struct dataInfo *curr = ListHead; curr != NULL; curr = curr->next){    
        int pos = NodeIndex(nodes_sorted_buffer,0,n - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,n - 1,curr->neg);
        if(strcmp(curr->type,"Resistance") == 0){
            if(pos != -1 && neg != -1)
                nz += 4;
            else
                nz += 1;
        }
        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(pos != -1)
                nz += 2;
            if(neg != -1)
                nz += 2;
        }
    }

    A = cs_spalloc(n+col2,n+col2,nz,1,1);
    A->nz = nz;

    for(struct dataInfo *curr = ListHead; curr != NULL; curr = curr->next){
        
        int pos = NodeIndex(nodes_sorted_buffer,0,n - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,n - 1,curr->neg);
        if(strcmp(curr->type,"Resistance") == 0){
            if(pos != -1 && neg != -1){
                SetValue(A, pos, pos, (double)(1.0/curr->value), &counter, nz);
                SetValue(A, neg, neg, (double)(1.0/curr->value), &counter, nz);
                SetValue(A, pos, neg, (double)(-1.0/curr->value), &counter, nz);
                SetValue(A, neg, pos,(double)(-1.0/curr->value), &counter, nz);
            }
            else if(pos == -1)
                SetValue(A, neg, neg, (double)(1.0/curr->value), &counter, nz);
            else
                SetValue(A, pos, pos, (double)(1.0/curr->value), &counter, nz);
        }

        else if(strcmp(curr->type,"Current Source") == 0){
            continue;
        }

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(pos != -1){
                SetValue(A, n + k, pos, 1.0, &counter, nz);
                SetValue(A, pos, n + k, 1.0, &counter, nz);
            }
            if(neg!= -1){
                SetValue(A, n + k, neg, -1.0, &counter, nz);
                SetValue(A, neg, n + k, -1.0, &counter, nz);
            }            
            k++;
        }
    }
    C = cs_compress(A);
    cs_spfree(A);
    cs_dupl(C);
    return C;
}

void DCAnalysis(struct dataInfo *ListHead,char **actions_ptr,int number_of_actions){
    struct dataInfo *curr;
    bool iter = false;
    bool cholesky = false;
    bool sparse = false;
    bool transient = false;
    bool euler = false;
    double itol;
    double time_step, fin_time;
    int col2 = 0;
    int nodes;
    char **nodes_sorted_buffer = NodesIdentification(ListHead, &nodes);
    
    for(int i = 0; i < number_of_actions; i++){
        if(strstr(actions_ptr[i],".OPTIONS")){
            if(strstr(actions_ptr[i],"SPD")){cholesky=true;}
            
            if(strstr(actions_ptr[i],"ITER")){
                iter=true;
                if(strstr(actions_ptr[i],"ITOL=") != NULL){
                    char *start = strstr(actions_ptr[i],"ITOL=");
                    char *end = start;
                    while(*(end+1)!=' ' && *(end+1)!= '\n' && *(end+1) != '\0'){
                        end++;
                    }
                    size_t length = end - start + 1;
                    char tolerance[length + 1];
                    memcpy(tolerance, start + strlen("ITOL="), length);
                    tolerance[length] = '\0';
                    for(int k = 0; tolerance[k]; k++){
                        tolerance[k] = tolower(tolerance[k]);
                    }
                    itol = atof(tolerance);
                }
            }

            if(strstr(actions_ptr[i],"SPARSE")){
                sparse = true;
            }

            if(strstr(actions_ptr[i],"METHOD=BE")){
                euler = true;
            }
        }

        if(strstr(actions_ptr[i],".TRAN")){
            transient = true;
            strtok(actions_ptr[i]," ");
            time_step = atof(strtok(NULL," "));
            fin_time = atof(strtok(NULL," "));
        }

        if(strstr(actions_ptr[i],".AC")){
            ACAnalysis(ListHead,nodes, nodes_sorted_buffer,number_of_actions,actions_ptr);
            return;
        }

    }

    for(curr = ListHead; curr != NULL;curr = curr->next){
        if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            col2++;
        }        
    }
    
    double ** LeftPart;
    double *RightPart;
    gsl_vector *x;


    if(transient){   
        RightPart = create_right_part(ListHead,col2);
        double *b;
        if(!sparse){  
            LeftPart = create_left_part(ListHead,col2);  
            x = method(LeftPart,RightPart,nodes,nodes_sorted_buffer,col2,itol,iter,cholesky);
            b = (double *)malloc((nodes+col2)*sizeof(double));
            for(int i=0; i < nodes+col2;i++){
                b[i] = gsl_vector_get(x,i);
            }
        }
        else{
            cs *LeftPart = create_left_part_sparse(ListHead,col2);
            b = method_cs(LeftPart,RightPart,cholesky,iter,nodes,col2,itol);
        }

        TransientAnalysis(ListHead, b,RightPart, actions_ptr,nodes_sorted_buffer,nodes,col2,number_of_actions, euler, time_step, fin_time,iter,itol,cholesky,sparse);
    }

    else{
        if(!sparse){
            printf("Creating Matrices...\n");
            //Create Matrices
            printf("-->Creating Left Part...\n");
            LeftPart = create_left_part(ListHead,col2);
            printf("-->Creating Right Part...\n");
            RightPart = create_right_part(ListHead,col2);
            printf("DC Analysis...\n");
            //Solve For Point Of Operation
            x = method(LeftPart,RightPart,nodes,nodes_sorted_buffer,col2,itol,iter,cholesky);
            FILE *point_of_operation_file = fopen("DCAnalysis.txt","w");
            for(int i = 0; i < nodes; i++){
                fprintf(point_of_operation_file,"%s %.3f\n",nodes_sorted_buffer[i],gsl_vector_get(x,i));
            }
            gsl_vector_free(x);
        }

        else{
            printf("Creating Matrices...\n");
            //Create Matrices
            printf("-->Creating Left Part...\n");
            cs *LeftPart = create_left_part_sparse(ListHead,col2);
            printf("-->Creating Right Part...\n");
            RightPart = create_right_part(ListHead,col2);
            printf("DC Analysis for Sparse Implementation...\n");
            //Solve For Point Of Operation
            double *b = method_cs(LeftPart,RightPart,cholesky,iter,nodes,col2,itol);
            printf("Saving Point Of Operation...\n");
            FILE *fp = fopen("DCAnalysis_Sparse.txt","w");
            for(int i = 0; i < nodes; i++){
                fprintf(fp,"%s %.3f\n",nodes_sorted_buffer[i],b[i]);
            }
        }

        //DC Sweep
        for(int i = 0; i < number_of_actions; i++){
            if(strstr(actions_ptr[i],".DC")){
                strtok(actions_ptr[i], " ");
                char *var = strtok(NULL, " ");
                double start = atof(strtok(NULL, " "));
                double end = atof(strtok(NULL, " "));
                double increment = atof(strtok(NULL, " "));
                char type = var[0];
                memmove(var, var+1, strlen(var));
                        
                int number_of_nodes = 0;
                
                for(int action = 0; action < number_of_actions; action++){
                    if(strstr(actions_ptr[action],".PLOT") || strstr(actions_ptr[action],".PRINT")){
                        for(int j = 0; j < strlen(actions_ptr[action]); j++){
                            if(actions_ptr[action][j] == 'V'){
                                number_of_nodes++;
                            }
                        }
                    }
                }
            
                char **inp_nodes = (char **)malloc(number_of_nodes*sizeof(char *));
                for(int j = 0; j < number_of_nodes; j++){
                    inp_nodes[j] = (char *)calloc(50,sizeof(char));
                }

                int k = 0;
                
                for(int action = 0; action < number_of_actions; action++){
                    if(strstr(actions_ptr[action],".PLOT") || strstr(actions_ptr[action],".PRINT")){
                        for(int j = 0; j < strlen(actions_ptr[action]); j++){
                            if(actions_ptr[action][j] == '('){
                                int counter = 0;
                                j++;
                                while(actions_ptr[action][j] != ')'){
                                    inp_nodes[k][counter++] = actions_ptr[action][j++];
                                }
                                k++;
                            }
                        }
                    }
                }

                FILE *out[number_of_nodes];

                for(int j = 0; j < number_of_nodes; j++){
                    out[j] = fopen(inp_nodes[j],"w+");
                }

                for(curr = ListHead; curr != NULL; curr = curr->next){
                    if(strcmp(curr->name,var) != 0)
                        continue;
                    if((strcmp(curr->type,"Current Source") == 0 && type == 'I') || (strcmp(curr->type,"Voltage Source") == 0 && type == 'V'))
                        break;
                    
                }

                double temp = curr->value;

                for(double val = start; val <= end; val+=increment){
                    if(!sparse){
                        curr->value = val;
                        LeftPart = create_left_part(ListHead,col2);
                        RightPart = create_right_part(ListHead,col2);
                        gsl_vector *x = method(LeftPart,RightPart,nodes,nodes_sorted_buffer,col2,itol,iter,cholesky);
                        for(int j = 0; j < number_of_nodes; j++){
                            fseek(out[j],0,SEEK_END);
                            int idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,inp_nodes[j]);
                            fprintf(out[j],"%lf %lf\n",val,gsl_vector_get(x,idx));
                        }
                        gsl_vector_free(x);
                        for(int i = 0; i < nodes + col2; i++)
                            free(LeftPart[i]);
                        free(LeftPart);
                    }

                    else{
                        curr->value = val;
                        cs *LeftPart = create_left_part_sparse(ListHead,col2);
                        RightPart = create_right_part(ListHead,col2);
                        double *b = method_cs(LeftPart,RightPart,cholesky,iter,nodes,col2,itol);
                        for(int j = 0; j < number_of_nodes; j++){
                            fseek(out[j],0,SEEK_END);
                            int idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,inp_nodes[j]);
                            fprintf(out[j],"%lf %lf\n",val,b[idx]);
                        }
                    }



                }
                curr->value = temp;
                break;
            }
        }
    }
    //Free Memory
    if(!transient)
        free(RightPart);
    for(int i = 0; i < nodes; i++)
        free(nodes_sorted_buffer[i]);
    free(nodes_sorted_buffer);
    
    
}
