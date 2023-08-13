// This function prints the matrix P: 
void displayM(double **P,int M,int N);

// This function sums by row the elements of a matrix:
void sum_by_rowM(double **P, int M, int N);

// This function displays a vector:
void displayV(double *v, int M);

// rouwenhorst rutine:
int rowenhorst(double *&x, double **&P, const double mu, const double rho, const double sig, const int N);

// This function takes a matrix and fill all the elements of the matrix, with the value assign to val:
void fillMatrix(double **P, int M, int N, double val);

// Main function:
int rowenhorst(double *&x, double **&P, const double mu, const double rho, const double sig, const int N);