#include <stdbool.h>
#include"csparse.h"
#ifndef __DC_H_
#define __DC_H_ 1

double **create_left_part(struct dataInfo *ListHead,int col2);
cs *create_left_part_sparse(struct dataInfo *ListHead,int col2);
void DCAnalysis(struct dataInfo *ListHead,char **actions_ptr,int number_of_actions);
void SetValue(cs *A, int row, int column, double value, int *counter, int nz);
#endif