#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]){
    int a = 0; 
    int left = ((a-1)%2); 
    int right = (a+1)%2; 
    printf("Me myself is: %d, left is %d, right is %d. \n", a, left, right); 
}
