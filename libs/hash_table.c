#include"hash_table.h"
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<complex.h>
#include<math.h>

#define MAX_ARGUMENTS 100 
#define MAX_ARGUMENT_LEN 20
#define MAX_OPT 100

void SortArray(char *array[], int n){
    int i, j, min_index;

    for (i = 0; i < n - 1; i++) {
        min_index = i;
        for (j = i + 1; j < n; j++)
            if (strcmp(array[j],array[min_index]) < 0)
                min_index = j;
        if(min_index != i)
            Swap(array,min_index, i);
    }
}

void Swap(char *array[],int idx1, int idx2){

    char *temp;
    temp = array[idx1];
    array[idx1] = array[idx2];
    array[idx2] = temp;    
}

int NodeIndex(char **nodes_sorted_buffer,int l,int r,char *str){
    if (r >= l) {
        int mid = l + (r - l) / 2;
        if(strcmp(nodes_sorted_buffer[mid],str) == 0)
            return mid;

        if (strcmp(nodes_sorted_buffer[mid],str) > 0)
            return NodeIndex(nodes_sorted_buffer,l,mid - 1,str);
 
        return NodeIndex(nodes_sorted_buffer,mid + 1,r,str);        
    }
    return -1;
}

char **NodesIdentification(struct dataInfo *ListHead,int *nodes){
    int index = 0,i;
    struct dataInfo *curr;
    char **nodes_buffer = (char **)malloc(50*sizeof(char *));
    char **nodes_sorted_buffer;
    *nodes = 0;
    
    for(curr = ListHead; curr != NULL; curr = curr->next){
        nodes_buffer[index] = (char *)malloc(strlen(curr->pos)+1);
        strcpy(nodes_buffer[index++],curr->pos);
        nodes_buffer[index] = (char *)malloc(strlen(curr->neg)+1);
        strcpy(nodes_buffer[index++],curr->neg);

        if(index % 50 == 0)
            nodes_buffer = (char **)realloc(nodes_buffer,(50 + index)*sizeof(char *));

    }


    SortArray(nodes_buffer,index);

    nodes_sorted_buffer = (char **)malloc(sizeof(char *));
    
    for(i = 0; i < index - 1; i++){
        if(strcmp(nodes_buffer[i],nodes_buffer[i+1]) != 0){
            nodes_sorted_buffer = (char **)realloc(nodes_sorted_buffer,(*nodes + 1)*sizeof(char *));
            nodes_sorted_buffer[*nodes] = (char *)malloc(sizeof(nodes_buffer[i+1]));
            strcpy(nodes_sorted_buffer[(*nodes)++],nodes_buffer[i+1]);
        }           
    }

    for(i = 0; i < index; i++)
        free(nodes_buffer[i]);
    free(nodes_buffer);

    return nodes_sorted_buffer;

}

double ValueCalc(char *param){
    int i;

    for(i = 0; i < strlen(param); i++){
        if(param[i] == 'E')
            break;
    }
    
    if(i==strlen(param))
        return atof(param);
    else{
        char prefix[i+1];
        memcpy(prefix, param,i);
        return atof(prefix)*pow((double)10,atof(&param[i+1]));
    }
}

struct dataInfo *CreateNode(char *token){
    char parameters[MAX_ARGUMENTS][MAX_ARGUMENT_LEN];
    int i,j = 0,counter = 0;
    struct dataInfo *Node = (struct dataInfo *)calloc(1,sizeof(struct dataInfo));
    struct electronicsInfo *electronic_fields = NULL;

    for(i = 0;i <= (strlen(token)); i++){       
        if(token[i] == ' ' || token[i] == '\0'){
            parameters[counter][j]='\0';
            counter++;  
            j=0;    
        }
        else{
            parameters[counter][j]=token[i];
            j++;
        }
    }

    Node->name = (char *)malloc(sizeof(&parameters[0][1]));
    strcpy(Node->name,&parameters[0][1]);
    Node->ac = CMPLXF(0,0);

    switch(parameters[0][0]){
        case 'V':
            Node->type = "Voltage Source";
            Node->pos = (char *)malloc(sizeof(parameters[1]));
            strcpy(Node->pos,parameters[1]);
            Node->neg = (char *)malloc(sizeof(parameters[2]));
            strcpy(Node->neg,parameters[2]);
            Node->value = ValueCalc(parameters[3]);
            break;
        case 'I':
            Node->type = "Current Source";
            Node->pos = (char *)malloc(sizeof(parameters[1]));
            strcpy(Node->pos,parameters[1]);
            Node->neg = (char *)malloc(sizeof(parameters[2]));
            strcpy(Node->neg,parameters[2]);
            Node->value = ValueCalc(parameters[3]);
            break;
        case 'R':
            Node->type = "Resistance";
            Node->pos = (char *)malloc(sizeof(parameters[1]));
            strcpy(Node->pos,parameters[1]);
            Node->neg = (char *)malloc(sizeof(parameters[2]));
            strcpy(Node->neg,parameters[2]);
            Node->value = ValueCalc(parameters[3]);
            break;
        case 'C':
            Node->type = "Capacitor";
            Node->pos = (char *)malloc(sizeof(parameters[1]));
            strcpy(Node->pos,parameters[1]);
            Node->neg = (char *)malloc(sizeof(parameters[2]));
            strcpy(Node->neg,parameters[2]);
            Node->value = ValueCalc(parameters[3]);
            break;
        case 'L':
            Node->type = "Inductor";
            Node->pos = (char *)malloc(sizeof(parameters[1]));
            strcpy(Node->pos,parameters[1]);
            Node->neg = (char *)malloc(sizeof(parameters[2]));
            strcpy(Node->neg,parameters[2]);
            Node->value = ValueCalc(parameters[3]);
            break;
        case 'D':
            Node->type = "Diode";
            Node->pos = (char *)malloc(sizeof(parameters[1]));
            strcpy(Node->pos,parameters[1]);
            Node->neg = (char *)malloc(sizeof(parameters[2]));
            strcpy(Node->neg,parameters[2]);
            electronic_fields = (struct electronicsInfo *)malloc(sizeof(struct electronicsInfo));
            electronic_fields->model_name = (char *)malloc(sizeof(parameters[3]));
            strcpy(electronic_fields->model_name,parameters[3]);
            electronic_fields->area = ValueCalc(parameters[4]);
            Node->electronic_fields = electronic_fields;
            break;
        case 'M':
            Node->type = "MOS Transistor";
            electronic_fields = (struct electronicsInfo *)malloc(sizeof(struct electronicsInfo));
            electronic_fields->D = (char *)malloc(sizeof(parameters[1]));
            strcpy(electronic_fields->D,parameters[1]);
            electronic_fields->G = (char *)malloc(sizeof(parameters[2]));
            strcpy(electronic_fields->G,parameters[2]);
            electronic_fields->S = (char *)malloc(sizeof(parameters[3]));
            strcpy(electronic_fields->S,parameters[3]);
            electronic_fields->B = (char *)malloc(sizeof(parameters[4]));
            strcpy(electronic_fields->B,parameters[4]);
            electronic_fields->model_name = (char *)malloc(sizeof(parameters[5]));
            strcpy(electronic_fields->model_name,parameters[5]);
            electronic_fields->L = ValueCalc(&parameters[6][2]);
            electronic_fields->W = ValueCalc(&parameters[7][2]);
            Node->electronic_fields = electronic_fields;
            break;
        case 'Q':
            Node->type = "BJT Transistor";
            electronic_fields = (struct electronicsInfo *)malloc(sizeof(struct electronicsInfo));
            electronic_fields->C = (char *)malloc(sizeof(parameters[1]));
            strcpy(electronic_fields->C,parameters[1]);
            electronic_fields->B = (char *)malloc(sizeof(parameters[2]));
            strcpy(electronic_fields->B,parameters[2]);
            electronic_fields->E = (char *)malloc(sizeof(parameters[3]));
            strcpy(electronic_fields->E,parameters[3]);
            electronic_fields->model_name = (char *)malloc(sizeof(parameters[4]));
            strcpy(electronic_fields->model_name,parameters[4]);
            electronic_fields->area = ValueCalc(parameters[5]);
            Node->electronic_fields = electronic_fields;
            break;
        default:
            printf("[ERROR] There is no circuit element named %c.\n",parameters[i][0]);

    }  
    double mag, phase;
    if(strstr(token,"AC")){
        char *s = strstr(token,"AC");
        s = s+3;
        if(s){
            mag = atof(strtok(s, " "));
            phase = atof(strtok(NULL, " "));
            Node->ac = CMPLX(mag * cos(phase*(M_PI/180)),(mag * sin(phase*(M_PI/180))));
        }
    }
    //Nodes for transient analysis
    if(strstr(token,"EXP") || strstr(token,"SIN") || strstr(token,"PULSE") || strstr(token,"PWL")){
        Node->transient_spec = (struct func *)calloc(1,sizeof(struct func));
        if(strstr(token,"EXP") || strstr(token,"SIN")){
            Node->transient_spec->type = (char *)malloc(4*sizeof(char));
            if(strstr(token,"EXP"))
                strcpy(Node->transient_spec->type,"EXP");
            else
                strcpy(Node->transient_spec->type,"SIN");
                
            Node->transient_spec->exp = (struct exp_sin_data *)calloc(1,sizeof(struct exp_sin_data));
            char *s = strstr(token,"(");
            if(s){
                memmove(s, s+1, strlen(s));
                Node->transient_spec->exp->i1 = atof(strtok(s, " "));
                Node->transient_spec->exp->i2 = atof(strtok(NULL, " "));
                Node->transient_spec->exp->td1 = atof(strtok(NULL, " "));
                Node->transient_spec->exp->tc1 = atof(strtok(NULL, " "));
                Node->transient_spec->exp->td2 = atof(strtok(NULL, " "));
                Node->transient_spec->exp->tc2 = atof(strtok(NULL, " "));               
            }
        }
        else if(strstr(token,"PULSE")){
            Node->transient_spec->type = (char *)malloc(6*sizeof(char));
            strcpy(Node->transient_spec->type,"PULSE");
            Node->transient_spec->pulse = (struct pulse_data *)calloc(1,sizeof(struct pulse_data));
            char *s = strstr(token,"(");
            if(s){
                memmove(s, s+1, strlen(s));
                Node->transient_spec->pulse->i1 = atof(strtok(s, " "));
                Node->transient_spec->pulse->i2 = atof(strtok(NULL, " "));
                Node->transient_spec->pulse->td = atof(strtok(NULL, " "));
                Node->transient_spec->pulse->tr = atof(strtok(NULL, " "));
                Node->transient_spec->pulse->tf = atof(strtok(NULL, " "));
                Node->transient_spec->pulse->pw = atof(strtok(NULL, " "));
                Node->transient_spec->pulse->per = atof(strtok(NULL, " "));               
            }
        }
        else if(strstr(token,"PWL")){
            Node->transient_spec->type = (char *)malloc(4*sizeof(char));
            strcpy(Node->transient_spec->type,"PWL");
            int m = 0,n = 0,len = 0;
            struct pwl_data *temp;
            for(int i = 0; i < strlen(token); i++){
                if(token[i] == '('){
                    m = i+1;
                }
                if(token[i] == ')'){
                    n = i;
                    len = n-m;
                    char *dest = (char*)malloc(sizeof(char) * (len + 1));
                    strncpy(dest, (token + m), len);
                    if(Node->transient_spec->pwl == NULL){
                        Node->transient_spec->pwl = (struct pwl_data *)calloc(1,sizeof(struct pwl_data));
                        Node->transient_spec->pwl->t = atof(strtok(dest, " "));
                        Node->transient_spec->pwl->i = atof(strtok(NULL, " "));
                        temp = Node->transient_spec->pwl;
                    }
                    else{
                        temp->next = (struct pwl_data *)calloc(1,sizeof(struct pwl_data));
                        temp->next->t = atof(strtok(dest, " "));
                        temp->next->i = atof(strtok(NULL, " "));
                        temp = temp->next;
                    }
                    free(dest);
                }
            }
        }

    }
    
    

    return Node;
}

struct dataInfo *CreateLinkedList(char *res){
    char *token;
    struct dataInfo *head, *prev, *temp;
    //int is_head = 1;
    int res_len = strlen(res);
    int start, end;

    for(end = 0; res[end+1] != '\n';end++);
    token = (char*)malloc(sizeof(char) * (end + 3));
    strncpy(token, res, end+1);
    token[end+2] = '\0';
    temp = CreateNode(token);
    free(token);
    head = temp;
    prev = temp;

    for(int i = end; i < res_len; i++){
        if(res[i] == '\n' && res[i+1] != '\0' && res[i+1] != '.'){
            start = ++i;
            while(res[i] != '\n')
                i++;
            end = i--;
            token = (char*)malloc(sizeof(char) * (end - start + 2));
            strncpy(token, (res+start), end-start);
            token[end-start+1] = '\0';
            
            temp = CreateNode(token);
            free(token);
            prev->next = temp;
            prev = temp;
        }
               
    }
    
    return head;
}

void PrintLinkedList(struct dataInfo *head){
    struct dataInfo *curr = head;
    int counter = 1;

    while(curr != NULL){
        printf("***** CIRCUIT ELEMENT %d *****\n",counter);
        printf("Type: %s\n",curr->type);
        printf("Name: %s\n",curr->name);

        if(strcmp(curr->type,"Diode") != 0 && strcmp(curr->type,"MOS Transistor") != 0 && strcmp(curr->type,"BJT Transistor") != 0){
            printf("Pos: %s\n",curr->pos);
            printf("Neg: %s\n",curr->neg);
            printf("Value: %lf\n",curr->value);
        }
        else if(strcmp(curr->type,"Diode") == 0){
            printf("Pos: %s\n",curr->pos);
            printf("Neg: %s\n",curr->neg);
            printf("Model Name: %s\n",curr->electronic_fields->model_name);
            printf("Area: %lf\n",curr->electronic_fields->area);
        }
        else if(strcmp(curr->type,"MOS Transistor") == 0){
            printf("D: %s\n",curr->electronic_fields->D);
            printf("G: %s\n",curr->electronic_fields->G);
            printf("S: %s\n",curr->electronic_fields->S);
            printf("B: %s\n",curr->electronic_fields->B);
            printf("Model Name: %s\n",curr->electronic_fields->model_name);
            printf("L: %lf\n",curr->electronic_fields->L);
            printf("W: %lf\n",curr->electronic_fields->W);
        }
        else if(strcmp(curr->type,"BJT Transistor") == 0){
            printf("C: %s\n",curr->electronic_fields->C);
            printf("B: %s\n",curr->electronic_fields->B);
            printf("E: %s\n",curr->electronic_fields->E);
            printf("Model Name: %s\n",curr->electronic_fields->model_name);
            printf("Area: %lf\n",curr->electronic_fields->area);
        }

        if(curr->transient_spec != NULL){
            printf("Transient spec: %s\n",curr->transient_spec->type);
            if(curr->transient_spec->exp != NULL){
                printf("i1: %lf\n",curr->transient_spec->exp->i1);
                printf("i2: %lf\n",curr->transient_spec->exp->i2);
                printf("td1: %lf\n",curr->transient_spec->exp->td1);
                printf("tc1: %lf\n",curr->transient_spec->exp->tc1);
                printf("td2: %lf\n",curr->transient_spec->exp->td2);
                printf("tc2: %lf\n",curr->transient_spec->exp->tc2);
            }
            if(curr->transient_spec->pulse != NULL){
                printf("i1: %lf\n",curr->transient_spec->pulse->i1);
                printf("i2: %lf\n",curr->transient_spec->pulse->i2);
                printf("td: %lf\n",curr->transient_spec->pulse->td);
                printf("tr: %lf\n",curr->transient_spec->pulse->tr);
                printf("tf: %lf\n",curr->transient_spec->pulse->tf);
                printf("pw: %lf\n",curr->transient_spec->pulse->pw);
                printf("per: %lf\n",curr->transient_spec->pulse->per);
            }
            if(curr->transient_spec->pwl != NULL){
                struct pwl_data *temp = curr->transient_spec->pwl;
                while(temp != NULL){
                    printf("t: %lf\n", temp->t);
                    printf("i: %lf\n", temp->i);
                    temp = temp->next;
                }
            }
        }

        if(curr->ac != CMPLX(0,0))
            printf("AC Value: %lf + j%lf\n", creal(curr->ac), cimag(curr->ac));

        curr = curr->next;
        counter++;
    }
}

void free_hash_table(struct dataInfo *ListHead){
    struct dataInfo *curr;
    for(curr = ListHead; curr->next != NULL; curr = curr->next){
        free(curr->next);
    }
    free(ListHead);
}
