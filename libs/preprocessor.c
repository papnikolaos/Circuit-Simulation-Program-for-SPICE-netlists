#include"preprocessor.h"
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include <complex.h>


char *readFile(char *filepath){
    FILE *ptr;
    int length;
    char *string;
    char c;
    int i = 0;

    ptr = fopen(filepath,"r");

    if(ptr == NULL){return "No such file";}
    
    fseek(ptr,0,SEEK_END);
    
    length = ftell(ptr);
    
    fseek(ptr, 0, SEEK_SET);
    
    string = (char *)malloc(sizeof(char) * (length+1));

    while((c = fgetc(ptr)) != EOF){
        if(!(c == '\n' && string[i-1] == '\n' )){
            string[i] = c;
            i++;
        }
    }

    string[i] = '\0';

    fclose(ptr);
    
    return string;
}

void CaseSensitive(char *str){
    int length = strlen(str);

    for(int i = 0; i <= length; i++){
        str[i] = toupper(str[i]);
    }
}

void SpaceRemoval(char *str){
    int length = strlen(str);
    int j = 0;
    int spaceFlag = 0;

    for(int i = 0; i < length; i++){
        if((str[i] == '\t' || str[i] == ' ') && spaceFlag == 0){
            spaceFlag = 1;
            str[j++] = ' ';
        }
        else if(spaceFlag == 1 && (str[i] == '\t' || str[i] == ' ')){
            continue;
        }

        else{
            spaceFlag = 0;
            str[j++] = str[i];
        }
    }
    str[j] = '\0';

}

void RemoveComments(char *str){
    int length = strlen(str);
    int starFlag = 0;
    int j = 0;

    for(int i = 0; i < length; i++){
        if(str[i] == '*'){starFlag = 1;}
        if(!starFlag){
            str[j++] = str[i];
        }
        else{
            if(str[i] == '\n'){starFlag = 0;}
        }
    }
    str[j] = '\0';
}

char *Preprocessing(char *filename){
    char *str = readFile(filename);
    CaseSensitive(str);
    RemoveComments(str);
    SpaceRemoval(str);
    return str;
}


