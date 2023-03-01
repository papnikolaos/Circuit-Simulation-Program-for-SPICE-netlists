#include"preprocessor.h"
#include"hash_table.h"
#include"transient.h"
#include"DC.h"
#include "AC.h"
#include<string.h>
#include<stdlib.h>
#include <stdbool.h>
#include <complex.h>


char **ExtractingActions(char *res, char **actions_ptr,int *number_of_actions){
    int start,end;

    for(int i = strlen(res) - 1; i >= 0; i--){
        if(res[i] == '\n' && res[i+1] != '\0'){
            if(res[i+1] != '.')
                break;
            start = i+1;
            for(end=start; res[end] != '\n' && res[end] != '\0' && end < strlen(res); end++);
            end--;
            actions_ptr[*number_of_actions] = (char *)calloc(end - start + 2,sizeof(char));
            for(int j = 0; j < end - start + 1; j++)
                actions_ptr[*number_of_actions][j] = res[start + j];
            (*number_of_actions)++;

            if((*number_of_actions) % 10 == 0){
                char **temp = (char **)malloc((*number_of_actions+10)*sizeof(char *));
                for(int j = 0; j < *number_of_actions; j++){
                    temp[j] = (char *)malloc(sizeof(char)*strlen(actions_ptr[j]));
                    strcpy(temp[j],actions_ptr[j]);
                    free(actions_ptr[j]);
                }
                actions_ptr = temp;
            }
        }        
    }

    printf("\n***** ACTIONS REQUIRED ******\n");
    for(int i = 0; i < *number_of_actions; i++)
        printf("%d) %s\n",i+1,actions_ptr[i]);
    printf("\n");

    return actions_ptr;
}

void control_panel(char *arg){
    int number_of_actions = 0;
    char **actions_ptr = (char **)malloc(10*sizeof(char *));
    //Preprocessor
    char *str = Preprocessing(arg);
    char *temp = malloc(strlen(str));
    memcpy(temp, str, strlen(str));
    printf("***** TEXT FILE AFTER PROCESSING *****\n");
    printf("%s", str);
    printf("\nPress any key to continue...\n");
    getchar();
    system("clear");
    //

    //Extract Actions
    actions_ptr = ExtractingActions(str,actions_ptr,&number_of_actions);
    printf("\nPress any key to continue...\n");
    getchar();
    system("clear");
    //

    //Create Hash Table
    struct dataInfo *ListHead = CreateLinkedList(str);
    PrintLinkedList(ListHead);
    printf("\nPress any key to continue...\n");
    getchar();
    system("clear");
    //
    
    DCAnalysis(ListHead,actions_ptr,number_of_actions);
    free_hash_table(ListHead);
    
}
