#include"libs/control_panel.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>


int main(int argc, char **argv){
    
    
    if(argc < 2){
        printf("[ERROR] You have to pass a file as an argument.\n");
        return 1;
    }
    //system("clear");
    control_panel(argv[1]);
    return 0;
}
