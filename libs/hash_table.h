#ifndef __HASH_TABLE_H_
#define __HASH_TABLE_H_ 1
#include<math.h>
#include<complex.h>

struct dataInfo{
    char *name;
    char *type;
    char *pos;
    char *neg;
    double value;
    double complex ac;
    struct func *transient_spec;
    struct electronicsInfo *electronic_fields;
    struct dataInfo *next;
};

struct electronicsInfo{
    char *model_name;
    double area;
    char *B;
    char *C;
    char *D;
    char *E;
    char *G;
    double L;
    char *S;
    double W;
};

struct exp_sin_data
{
    double i1;
    double i2;
    double td1; 
    double tc1;
    double td2;
    double tc2;
};

struct pulse_data
{
    double i1;
    double i2;
    double td; 
    double tr;
    double tf;
    double pw;
    double per;
};

struct pwl_data
{
    double i;
    double t;
    struct pwl_data *next;
};


struct func{
    char *type;
    struct exp_sin_data *exp;
    struct pulse_data *pulse;
    struct pwl_data *pwl;
};

struct dataInfo *CreateLinkedList(char *res);
void PrintLinkedList(struct dataInfo *head);
void free_hash_table(struct dataInfo *ListHead);
char **NodesIdentification(struct dataInfo *ListHead,int *nodes);
void Swap(char *array[],int idx1, int idx2);
int NodeIndex(char **nodes_sorted_buffer,int l,int r,char *str);
#endif
