#include<stdio.h>
#include<stdlib.h>
#include <stdbool.h>
#include"hash_table.h"
#include"csparse.h"

#ifndef __transient_H_
#define __transient_H_ 1



void TransientAnalysis(struct dataInfo *ListHead, double *dc_solution, double *dc_RightPart, char **actions_ptr,char **nodes_sorted_buffer,
int nodes,int col2,int number_of_actions, bool euler, double time_step, double fin_time,bool iter,double itol,bool cholesky,bool sparse);
double EXP(double i1,double i2,double td1, double tc1, double td2, double tc2,double t);
double SIN(double i1, double ia, double fr,double td, double df, double ph,double t);
double PULSE(double i1, double i2, double td, double tr, double tf, double pw, double per, double t);
double PWL(struct pwl_data *pwl, double t);
#endif
