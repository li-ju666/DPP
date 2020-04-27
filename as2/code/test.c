#include "mat.h"

int main(){
    float A[4] = {1, 2, 3, 4}; 
    float B[4] = {2, 3, 4, 5}; 
    float* C = multiply(A, B, 2); 
    vis(A, 2); 
    vis(B, 2); 
    vis(C, 2); 
    free(C); 
}
