#include <stdio.h>

int add(int a, int b){
    return a+b; 
}

int main(){
    int a = 2, b = 2; 
    int sum = add(a, b); 
    int array[sum]; 
    for(int i=0; i<sum; i++){
	array[i] = 1; 
    }
    for(int i=0; i<sum; i++){
	printf("%d ", array[i]); 
    }
}
