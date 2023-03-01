#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <gsl/gsl_linalg.h>
#include"hash_table.h"
#include"sol_methods.h"
#include"transient.h"
#include"csparse.h"
#include"DC.h"
#include <complex.h>

/*FUNCTIONS*/
double EXP(double i1,double i2,double td1, double tc1, double td2, double tc2,double t){
    if(t >= 0 && t <= td1){
        return i1;
    }

    else if(t >= td1 && t <= td2){
        return i1 + (i2-i1)*(1 - exp(-(t-td1)/tc1));
    }

    else{
        return i1 + (i2-i1)*(exp(-(t-td2)/tc2) - exp(-(t-td1)/tc1));
    }
}

double SIN(double i1, double ia, double fr,double td, double df, double ph,double t){

    if(t >= 0 && t <= td){
        return i1 + ia*sin(2*M_PI*ph / 360);
    }

    else if(t >= td){
        return i1 + ia*sin(2*M_PI*fr*(t-td)+2*M_PI*(ph/360))*exp(-(t-td)*df);
    }

    return -1;
}

double PULSE(double i1, double i2, double td, double tr, double tf, double pw, double per, double t){
    int k = (int)(t/per);
    t -= k*per;

    if(t <= td)
        return i1;

    else if(t >= td && t<= td+tr)
        return i1 + (i2-i1)*(t-td)/tr;

    else if(t >= td+tr && t <= td+tr+pw)
        return i2;

    else if(t >= td+tr+pw && t <= td+tr+pw+tf)
        return i2 + (i1-i2)*(t-td-tr-pw)/tf;

    else if(t >= td+tr+pw+tf && t <= td+per)
        return i1;

    return -1;
}

double PWL(struct pwl_data *pwl, double t){
    if(t <= pwl->t)
        return pwl->i;
    
    struct pwl_data *temp = pwl;

    while(temp->next != NULL){
        if(t <= temp->next->t)
            return temp->i + (temp->next->i-temp->i)*(t-temp->t)/(temp->next->t-temp->t);
        temp = temp->next;
    }

    return temp->t;   
}
/*--------------------------*/
/*CREATE MATRICES*/
double **create_c_tilde_div_time_step(struct dataInfo *ListHead,int nodes,int col2, double time_step,char **nodes_sorted_buffer){
    double **c_tilde = (double **)calloc(nodes + col2,sizeof(double *));
    struct dataInfo *curr;
    int k = 0;

    for(int i = 0; i < nodes + col2; i++){
        c_tilde[i] = (double *) calloc(nodes + col2, sizeof(double));
    }

    for(curr = ListHead; curr != NULL; curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);

        if(strcmp(curr->type,"Capacitor") == 0){
            if(pos != -1 && neg != -1){
                c_tilde[pos][pos] += curr->value/time_step;
                c_tilde[neg][neg] += curr->value/time_step;
                c_tilde[pos][neg] -= curr->value/time_step;
                c_tilde[neg][pos] -= curr->value/time_step;
            }
            else if(pos == -1)
                c_tilde[neg][neg] += curr->value/time_step;
            else
                c_tilde[pos][pos] += curr->value/time_step;
        }

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(strcmp(curr->type,"Inductor") == 0)
                c_tilde[nodes + k][nodes + k] -= curr->value/time_step;
            k++;
        }

    }

    return c_tilde;
    
}

cs *create_c_tilde_div_time_step_sparse(struct dataInfo *ListHead,int nodes,int col2, char **nodes_sorted_buffer,double time_step){
    int nz = 0;
    int k = 0;
    int counter = 0;

    for(struct dataInfo *curr = ListHead; curr != NULL; curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);

        if(strcmp(curr->type,"Capacitor") == 0){
            if(pos != -1 && neg != -1){
               nz += 4;
            }
            else if(pos == -1)
                nz += 1;
            else
                nz += 1;
        }

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(strcmp(curr->type,"Inductor") == 0)
                nz += 1;
            k++;
        }

    }
    
    cs *A = cs_spalloc(nodes+col2,nodes+col2,nz,1,1);
    A->nz = nz;
    k = 0;
    for(struct dataInfo *curr = ListHead; curr != NULL; curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);

        if(strcmp(curr->type,"Capacitor") == 0){
            if(pos != -1 && neg != -1){
                SetValue(A, pos, pos, (double)(curr->value/time_step), &counter, nz);
                SetValue(A, neg, neg, (double)(curr->value/time_step), &counter, nz);
                SetValue(A, pos, neg, (double)(-curr->value/time_step), &counter, nz);
                SetValue(A, neg, pos, (double)(-curr->value/time_step), &counter, nz);
            }
            else if(pos == -1)
                SetValue(A, neg, neg, (double)(curr->value/time_step), &counter, nz);
            else
                SetValue(A, pos, pos, (double)(curr->value/time_step), &counter, nz);
        }

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(strcmp(curr->type,"Inductor") == 0)
                SetValue(A, nodes+k, nodes+k, (double)(-curr->value/time_step), &counter, nz);
            k++;
        }

    }

    cs *C = cs_compress(A);
    cs_spfree(A);
    cs_dupl(C);
    return C;
}

double *calc_e(struct dataInfo *ListHead,int nodes,int col2,double t,char **nodes_sorted_buffer){
    int k = 0;
    double *RightPart = (double *) calloc(nodes + col2, sizeof(double));

    for(struct dataInfo *curr = ListHead; curr != NULL;curr = curr->next){
        int pos = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->pos);
        int neg = NodeIndex(nodes_sorted_buffer,0,nodes - 1,curr->neg);
        if(strcmp(curr->type,"Resistance") == 0){
            continue;
        }
        else if(strcmp(curr->type,"Current Source") == 0){
            if(pos != -1){
                if(curr->transient_spec == NULL)
                    RightPart[pos] -=  curr->value;
                else{
                    if(strcmp(curr->transient_spec->type,"EXP") == 0)
                        RightPart[pos] -= EXP(curr->transient_spec->exp->i1,curr->transient_spec->exp->i2,curr->transient_spec->exp->td1,curr->transient_spec->exp->tc1,
                        curr->transient_spec->exp->td2, curr->transient_spec->exp->tc2,t);
                    else if(strcmp(curr->transient_spec->type,"SIN") == 0)
                        RightPart[pos] -= SIN(curr->transient_spec->exp->i1,curr->transient_spec->exp->i2,curr->transient_spec->exp->td1,curr->transient_spec->exp->tc1,
                        curr->transient_spec->exp->td2, curr->transient_spec->exp->tc2,t);
                    else if(strcmp(curr->transient_spec->type,"PULSE") == 0)
                        RightPart[pos] -= PULSE(curr->transient_spec->pulse->i1,curr->transient_spec->pulse->i2,curr->transient_spec->pulse->td,curr->transient_spec->pulse->tr,
                        curr->transient_spec->pulse->tf,curr->transient_spec->pulse->pw,curr->transient_spec->pulse->per,t);
                    else if(strcmp(curr->transient_spec->type,"PWL") == 0)
                        RightPart[pos] -= PWL(curr->transient_spec->pwl,t);
                }
            }
                
            if(neg != -1){ 
                if(curr->transient_spec == NULL)
                    RightPart[neg] += curr->value;
                else{
                    if(strcmp(curr->transient_spec->type,"EXP") == 0)
                        RightPart[neg] += EXP(curr->transient_spec->exp->i1,curr->transient_spec->exp->i2,curr->transient_spec->exp->td1,curr->transient_spec->exp->tc1,
                        curr->transient_spec->exp->td2, curr->transient_spec->exp->tc2,t);
                    else if(strcmp(curr->transient_spec->type,"SIN") == 0)
                        RightPart[neg] += SIN(curr->transient_spec->exp->i1,curr->transient_spec->exp->i2,curr->transient_spec->exp->td1,curr->transient_spec->exp->tc1,
                        curr->transient_spec->exp->td2, curr->transient_spec->exp->tc2,t);
                    else if(strcmp(curr->transient_spec->type,"PULSE") == 0)
                        RightPart[neg] += PULSE(curr->transient_spec->pulse->i1,curr->transient_spec->pulse->i2,curr->transient_spec->pulse->td,curr->transient_spec->pulse->tr,
                        curr->transient_spec->pulse->tf,curr->transient_spec->pulse->pw,curr->transient_spec->pulse->per,t);
                    else if(strcmp(curr->transient_spec->type,"PWL") == 0)
                        RightPart[neg] += PWL(curr->transient_spec->pwl,t);
                }
            }
        }

        else if(strcmp(curr->type,"Voltage Source") == 0 || strcmp(curr->type,"Inductor") == 0){
            if(strcmp(curr->type,"Voltage Source") == 0){
                if(curr->transient_spec == NULL)
                    RightPart[nodes + k] += curr->value;
                else{
                    if(strcmp(curr->transient_spec->type,"EXP") == 0)
                        RightPart[nodes + k] += EXP(curr->transient_spec->exp->i1,curr->transient_spec->exp->i2,curr->transient_spec->exp->td1,curr->transient_spec->exp->tc1,
                        curr->transient_spec->exp->td2, curr->transient_spec->exp->tc2,t);
                    else if(strcmp(curr->transient_spec->type,"SIN") == 0)
                        RightPart[nodes + k] += SIN(curr->transient_spec->exp->i1,curr->transient_spec->exp->i2,curr->transient_spec->exp->td1,curr->transient_spec->exp->tc1,
                        curr->transient_spec->exp->td2, curr->transient_spec->exp->tc2,t);
                    else if(strcmp(curr->transient_spec->type,"PULSE") == 0)
                        RightPart[nodes + k] += PULSE(curr->transient_spec->pulse->i1,curr->transient_spec->pulse->i2,curr->transient_spec->pulse->td,curr->transient_spec->pulse->tr,
                        curr->transient_spec->pulse->tf,curr->transient_spec->pulse->pw,curr->transient_spec->pulse->per,t);
                    else if(strcmp(curr->transient_spec->type,"PWL") == 0)
                        RightPart[nodes + k] += PWL(curr->transient_spec->pwl,t);
                }
            }

                
            k++;
        }
    }

    return RightPart;

}

void TransientAnalysis(struct dataInfo *ListHead, double *dc_solution, double *dc_RightPart, char **actions_ptr,char **nodes_sorted_buffer,
int nodes,int col2,int number_of_actions, bool euler, double time_step, double fin_time,bool iter,double itol,bool cholesky,bool sparse){
    
    char **out_nodes;
    int number_of_nodes = 0;
    int counter = 0;
    int idx;

    for(int i = 0; i < number_of_actions; i++){
        if(strstr(actions_ptr[i],".PLOT") || strstr(actions_ptr[i],".PRINT")){
            for(int j = 0; j < strlen(actions_ptr[i]); j++){
                if(actions_ptr[i][j] == '(')
                    number_of_nodes++;
            }
        }
    }

    out_nodes = (char **)malloc(number_of_nodes*sizeof(char *));
    
    for(int i = 0; i < number_of_actions; i++){
        int m = 0,n = 0,len = 0;
        if(strstr(actions_ptr[i],".PLOT") || strstr(actions_ptr[i],".PRINT")){
            for(int j = 0; j < strlen(actions_ptr[i]); j++){
                if(actions_ptr[i][j] == '(')
                    m = j+1;
                if(actions_ptr[i][j] == ')'){
                    n = j-1;
                    len = n-m+1;
                    out_nodes[counter] = (char*)malloc(sizeof(char) * (len + 1));
                    for(int k = 0; k < len; k++)
                        out_nodes[counter][k] = actions_ptr[i][m+k];
                    out_nodes[counter++][len] = '\0';
                }
            }
        }
    }

    FILE *out[number_of_nodes];

    for(int j = 0; j < number_of_nodes; j++){
        out[j] = fopen(out_nodes[j],"w+");
    }

    for(int j = 0; j < number_of_nodes; j++){        
        idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,out_nodes[j]);
        fprintf(out[j],"%lf %lf\n",0.0,dc_solution[idx]);
    }

    if(!sparse){
        if(euler){//method=be
            double *prev_sol = dc_solution;
            double **g_tilde = create_left_part(ListHead,col2);
            double **c_tilde = create_c_tilde_div_time_step(ListHead,nodes,col2,time_step,nodes_sorted_buffer);
            double **LeftPart = matrix_add(nodes+col2,nodes+col2,c_tilde,g_tilde);
            free(g_tilde);

            for(double t = time_step; t < fin_time + time_step/1000; t += time_step){
                double *temp1 = matrix_x_vector(nodes+col2,nodes+col2,c_tilde,prev_sol);
                free(prev_sol);
                double *curr_e = calc_e(ListHead,nodes,col2,t,nodes_sorted_buffer);
                double *RightPart = vector_add(nodes+col2,curr_e,temp1);
                gsl_vector *x = method(LeftPart,RightPart, nodes,nodes_sorted_buffer,col2,itol,iter,cholesky);
                prev_sol = extract_gsl_vector(x,nodes+col2);
                gsl_vector_free(x);
                for(int j = 0; j < number_of_nodes; j++){        
                    idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,out_nodes[j]);
                    fprintf(out[j],"%lf %lf\n",t,prev_sol[idx]);
                }
                free(temp1);
                free(RightPart);
            }
            free(c_tilde);
            free(prev_sol);
            free(LeftPart);
        }

        else{  //method=tr          
            double *prev_sol = dc_solution;
            double *prev_e = dc_RightPart;
            double **g_tilde = create_left_part(ListHead,col2);
            double **c_tilde = create_c_tilde_div_time_step(ListHead,nodes,col2,(double)((time_step/2)),nodes_sorted_buffer);
            double **LeftPart = matrix_add(nodes+col2,nodes+col2,c_tilde,g_tilde);
            double **g_tilde_minus_c_tilde = matrix_sub(nodes+col2,nodes+col2,g_tilde,c_tilde);
            free(g_tilde);
            free(c_tilde);

            for(double t = time_step; t < fin_time + time_step/1000; t += time_step){
                double *temp1 = matrix_x_vector(nodes+col2,nodes+col2,g_tilde_minus_c_tilde,prev_sol);
                free(prev_sol);
                double *temp2 = vector_sub(nodes+col2,prev_e,temp1);
                double *curr_e = calc_e(ListHead,nodes,col2,t,nodes_sorted_buffer);
                double *RightPart = vector_add(nodes+col2,curr_e,temp2);
                gsl_vector *x = method(LeftPart,RightPart, nodes,nodes_sorted_buffer,col2,itol,iter,cholesky);
                prev_sol = extract_gsl_vector(x,nodes+col2);
                gsl_vector_free(x);
                for(int j = 0; j < number_of_nodes; j++){        
                    idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,out_nodes[j]);
                    fprintf(out[j],"%lf %lf\n",t,prev_sol[idx]);
                }
                free(temp1);
                free(temp2);
                free(prev_e);
                free(RightPart);
                prev_e = curr_e;
            }
            free(prev_sol);
            free(LeftPart);
            free(g_tilde_minus_c_tilde);
        }
    }

    else{ 
        if(euler){ //method = be, sparse
            double *prev_sol = dc_solution;
            cs *g_tilde = create_left_part_sparse(ListHead,col2);
            cs *c_tilde = create_c_tilde_div_time_step_sparse(ListHead,nodes,col2, nodes_sorted_buffer,time_step);
            cs *LeftPart = cs_add(g_tilde,c_tilde,1.0,1.0);
            cs_free(g_tilde);

            for(double t = time_step; t < fin_time + time_step/1000; t += time_step){
                double *curr_e = calc_e(ListHead,nodes,col2,t,nodes_sorted_buffer);
                double *placeholder = calloc(nodes+col2,sizeof(double));
                cs_gaxpy(c_tilde,prev_sol,placeholder);
                free(prev_sol); 
                double *RightPart = vector_add(nodes+col2,curr_e,placeholder);
                free(placeholder);
                free(curr_e);
                prev_sol = method_cs(LeftPart,RightPart,cholesky,iter,nodes,col2,itol);
                for(int j = 0; j < number_of_nodes; j++){        
                    idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,out_nodes[j]);
                    fprintf(out[j],"%lf %lf\n",t,prev_sol[idx]);
                }
            }
            cs_free(c_tilde);
            cs_free(LeftPart);
        }   

        else{ //method = tr,sparse
            cs *g_tilde = create_left_part_sparse(ListHead,col2);
            cs *c_tilde = create_c_tilde_div_time_step_sparse(ListHead,nodes,col2,nodes_sorted_buffer,time_step/2);
            cs *LeftPart = cs_add(g_tilde,c_tilde,1.0,1.0);
            cs *g_tilde_minus_c_tilde = cs_add(g_tilde,c_tilde,1.0,-1.0);
            double *prev_sol = dc_solution;
            double *prev_e = dc_RightPart;
            cs_free(g_tilde);
            cs_free(c_tilde);

            for(double t = time_step; t < fin_time + time_step/1000; t += time_step){
                double *curr_e = calc_e(ListHead,nodes,col2,t,nodes_sorted_buffer);
                double *temp = vector_add(nodes+col2,curr_e,prev_e);
                free(prev_e);
                prev_e = curr_e;
                double *placeholder = calloc(nodes+col2,sizeof(double));
                cs_gaxpy(g_tilde_minus_c_tilde,prev_sol,placeholder);
                free(prev_sol);                
                double *RightPart = vector_sub(nodes+col2,temp,placeholder);
                prev_sol = method_cs(LeftPart,RightPart,cholesky,iter,nodes,col2,itol);
                for(int j = 0; j < number_of_nodes; j++){        
                    idx = NodeIndex(nodes_sorted_buffer,0,nodes - 1,out_nodes[j]);
                    fprintf(out[j],"%lf %lf\n",t,prev_sol[idx]);
                }
                free(temp);
                free(placeholder);
            }
            free(prev_sol);
            cs_free(LeftPart);
            cs_free(g_tilde_minus_c_tilde);
        }
    }
}
