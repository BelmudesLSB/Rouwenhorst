#include<cmath> 
#include<iostream>
#include<cstdlib>

using namespace std; 

// We need to include the dimensions of the matrix, and then this function will output everything.

void displayM(double **P, int M, int N)
{
	int i, j;
	
	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
		{
			cout << P[i][j] << "  ";
		}
		cout << endl;
		cout << endl;
		cout << endl;
	}
	cout << endl;
}

void sum_by_rowM(double **P, int M, int N)
{
	int i, j;
	float temp;	
	for (i=0;i<M;i++){
		temp = 0.0;
		for (j=0;j<N;j++){
			temp = temp + P[i][j];
		}
		cout << "row " << i << " sums: " << temp << endl;
	}
	cout << endl;
}

void displayV(double *v, int M)
{
    for (int i=0;i<M;i++){
        cout << v[i] << endl;
    }
}

void fillMatrix(double **P, int M, int N, double val){

    for (int i=0;i<M;i++){
        for (int j=0; j<N; j++){
            P[i][j] = val;
        }
    }
}

int rowenhorst(double *&x, double **&P, const double mu, const double rho, const double sig, const int N){
	
	// check parameters
	if (abs(rho)<1 && sig>0 && N>=2){

			x = new double[N];
			P = new double*[N];

			//This makes P to be a matrix:

			for (int i=0; i<N; i++){
				P[i] = new double[N];
			}

			double p, psi;

			p = (1 + rho) / 2;
		    psi = sqrt((N - 1) / (1 - pow(rho,2)))*sig;

			// x.setLinSpaced(N,-psi,psi); (double) changes the data type to ensure that the result is a double.
			double dx = 2*psi/(double)(N-1);


			// This constructs the output grid:
			for (int i=0;i<N;i++){
				x[i] = mu/(1-rho) -psi + i*dx;
			}

            // raise x to exp(x):
            
            for (int i=0;i<N;i++){
                x[i] = exp(x[i]);
            }

			/// This constructs the transition matrix:

			// MatrixXd theta[N], temp1, temp2, temp3, temp4, temp5;
			// Start by creating the double pointers:
			
			double **theta0, **theta, **temp1, **temp2, **temp3, **temp4, **temp5;

			// Create a matrix of pointers, that is, every element of the matrix is an address.
			// Note that we are using dynamic memory.

			theta0 = new double*[2];
			for (int i=0;i<2;i++){
				theta0[i] = new double[2];
			}

			// For each address assing the following values:

			theta0[0][0] = p;
			theta0[0][1] = 1-p;
			theta0[1][0] = 1-p;
			theta0[1][1] = p;

			for (int i=2;i<N;i++){
			
				// temp1.setZero(i+1,i+1);
				// temp1.topLeftCorner(i, i) = theta[i-1];

				// Create a matrix of [3x3], we had a [2x2]:
				temp1 = new double*[i+1];

				for (int j=0;j<i+1;j++){
					temp1[j] = new double[i+1];
				}

				// Fill the matrix with zeros:

				fillMatrix(temp1,i+1,i+1, 0);

				// Copy the values of theta0 to temp1 or the top right corner of temp1.

				for (int i1=0;i1<i;i1++){
					for (int j1=0;j1<i;j1++){
						temp1[i1][j1] = theta0[i1][j1];
					}
				}
				
				/*
				cout << "temp1=\n";
				displayM(temp1,i+1,i+1);
				*/
				
				//temp2.setZero(i+1,i+1);
				//temp2.topRightCorner(i, i) = theta[i-1];

				// Create a matrix of [3x3], we had a [2x2]:

				temp2 = new double*[i+1];
				for (int j=0;j<i+1;j++){
					temp2[j] = new double[i+1];
				}

				// Fill the matrix with zeros:

				fillMatrix(temp2,i+1,i+1,0);

				// Copy the values of theta0 to temp2 or the top left corner of temp2.

				for (int i1=0;i1<i;i1++){
					for (int j1=0;j1<i;j1++){
						temp2[i1][j1+1] = theta0[i1][j1];
					}
				}
				
				// temp3.setZero(i+1,i+1);
				// temp3.bottomLeftCorner(i, i) = theta[i-1];
				
				temp3 = new double*[i+1];
				for (int j=0;j<i+1;j++){
					temp3[j] = new double[i+1];
				}

				fillMatrix(temp3,i+1,i+1,0);
				
				// Copy the values of theta0 to temp3 or the bottom right corner of temp3.

				for (int i1=0;i1<i;i1++){
					for (int j1=0;j1<i;j1++){
						temp3[i1+1][j1] = theta0[i1][j1];
					}
				}
				
				// temp4.setZero(i+1,i+1);
				// temp4.bottomRightCorner(i, i) = theta[i-1];

				temp4 = new double*[i+1];

				for (int j=0;j<i+1;j++){
					temp4[j] = new double[i+1];
				}

				fillMatrix(temp4,i+1,i+1,0);

				// Copy the values of theta0 to temp3 or the bottom left corner of temp4.

				for (int i1=0;i1<i;i1++)
				{
					for (int j1=0;j1<i;j1++)
					{
						temp4[i1+1][j1+1] = theta0[i1][j1];
					}
				}
				
				// Complete the whole iteration for the matrix \theta:

				//theta[i] = p*temp1 + (1-p)*temp2 + (1-p)*temp3 + p*temp4  (*);

				// Construct a [3x3] matrix of pointers: Why dynamic memory?

				theta = new double*[i+1];

				for (int j=0;j<i+1;j++){
					theta[j] = new double[i+1];
				}

				// Fill the matrix with the corresponding values following (*):

				for (int i1=0;i1<i+1;i1++)
				{
					for (int j1=0;j1<i+1;j1++)
					{
						theta[i1][j1] = p*temp1[i1][j1] + (1-p)*temp2[i1][j1] + (1-p)*temp3[i1][j1] + p*temp4[i1][j1];
					}
				}
				
				//temp5.setOnes(i+1,i+1);
				//temp5.middleRows(1,i-1).fill(0.5);

				//create a matrix of [3x3]:

				temp5 = new double*[i+1];

				for (int j=0;j<i+1;j++){
					temp5[j] = new double[i+1];
				}

				// Fill the matrix: Remember that up to this point, i=2.

				//// [	1;   1;  0.5]
				//// [	1;   1;  0.5]
				//// [0.5; 0.5;  0.5]

				for (int i1=0;i1<i+1;i1++){
					for (int j1=0;j1<i+1;j1++){
						if (i1==0 || i1==i){
							// This is either the first row or the i^{th} row:
							temp5[i1][j1] = 1;
						}
						else{
							// Everything else takes the value of 0.5:
							temp5[i1][j1] = 0.5;
						}
					}
				}

				// This step ????
				// theta[i] = theta[i].array() * temp5.array(); // dot product

				for (int i1=0;i1<i+1;i1++){
					for (int j1=0;j1<i+1;j1++){
						theta[i1][j1] = theta[i1][j1]*temp5[i1][j1];
					}
				}
				
				// clean up. This deletes all the elements on each row of temp's. This deletes all the columns data.
				for (int j=0;j<i+1;j++){
					delete[] temp1[j];
					delete[] temp2[j];
					delete[] temp3[j];
					delete[] temp4[j];
					delete[] temp5[j]; 
				}

				// This deletes all the rows data.

				delete[] temp1;
				delete[] temp2;
				delete[] temp3;
				delete[] temp4;
				delete[] temp5;
				
				//transfer to theta0, first, we delete all existing data on theta0:

				for (int j=0;j<i;j++){
					delete[] theta0[j]; 
				}

				delete[] theta0;

				// Create a new thetha0, now [3x3]:

				theta0 = new double*[i+1];

				for (int j=0;j<i+1;j++){
					theta0[j] = new double[i+1]; 
				}

				// Fill theta0[3x3] with the values of theta:

				for (int i1=0;i1<i+1;i1++){
					for (int j1=0;j1<i+1;j1++){
						theta0[i1][j1] = theta[i1][j1];
					}
				}

				// Delete theta (not theta0):
				
				for (int j=0;j<i+1;j++){
					delete[] theta[j]; 
				}

				delete[] theta;		
		}

		// After finisishing with all the iteration, P = theta0;

		for (int i=0;i<N;i++)
		{
			for (int j=0;j<N;j++)
			{
				P[i][j] = theta0[i][j];
			}
		}
			
		// clean up
			
		for (int i=0;i<N;i++)
		{
			delete[] theta0[i];
		}
		delete[] theta0;

		return EXIT_SUCCESS;

	} else {

		cout << "Error in rowenhorst: invalid parameters" << endl;

		return EXIT_FAILURE;
	}
}