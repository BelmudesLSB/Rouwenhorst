#include <iostream>
#include <vector> 
#include "rouw.hpp" 

//copilot?

int main(){

    const int rows = 3;
    const int columns = 3;

    // Create a 2D array. 
    // LHS: This creates a pointer to a pointer. That means that the variable holds and address, when you go to that adress then you have other addresses!
    // RHS: This is an array of pointers. We have row many POINTERS that belong to elements type double.

    double** myArray = new double*[rows];

    for (int i=0; i<rows; i++)
    {
        // We created row many POINTERS of data type double, now we need to fill them!
        // For each of the [row] many addresses, we create a dynamic array of [column] number double type variables.

        myArray[i] = new double[columns];
    }

    // Initialize the array with values:

    myArray[0][0] = 1;
    myArray[0][1] = 2;
    myArray[0][2] = 3;
    myArray[1][0] = 4;
    myArray[1][1] = 5;
    myArray[1][2] = 6;
    myArray[2][0] = 7;
    myArray[2][1] = 8;
    myArray[2][2] = 9;

    // Display the array:

    displayM(myArray, rows, columns);

    // Test fillMatrix:
    fillMatrix(myArray,3,3,2.0);

    // Display the array:

    displayM(myArray, rows, columns);
    
    // Deallocate the memory:

    for (int i = 0; i < rows; i++)

    {   // For each row address, delete everything inside. 
        // Note that if we would delete myArray, we would never be able to find and delete de columns that were created!

        delete[] myArray[i];
    }

    // We still need to delete the double pointer.
    // Note that the command is delete[] "Memory Address": That would delete everything in that address.

    delete[] myArray;

    // Note that there is still a variable called myArray that points towards a well defined (but empty) memory address.
    // Since that memory address might be used again, we make the pointer go to another place.

    myArray = NULL;

    // We want to call rouwen to see the output array.

    double *y_gird;
    double **Prob;
    const double mu = 0.0;
    const double rho = 0.5;
    const double sig = 0.2;
    const int N = 100;
    
    // Call the function rowenhorst with the parameters:    

    rowenhorst(y_gird,Prob, mu, rho, sig, N);

    // Display the output:

    //displayV(y_gird,N);
    //displayM(Prob,N,N);
    sum_by_rowM(Prob,N,N);

    return 0;


}