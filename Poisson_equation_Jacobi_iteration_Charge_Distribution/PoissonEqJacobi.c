//===================================================//
// Solving the Poisson equation on a 2d grid
//===================================================//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <cmath>

using namespace std;


double Estimate_EX(double **grid, int d, double a, int i, int j); // Function for estimating the x field
double Estimate_Ey(double **grid, int d, double a, int i, int j); // Function for estimating the y field

int main ()
{
	// Defining the variables of the system

	int64_t N, i = 0, j=0, step = 0;
    double  prec, var, app, L=1.,a, gamma, E_x, E_y, Energy;
    double** grid;
    bool done =false;
    char fname[] = "GridStep00.dat"; //modify fname[8], fname[9] to change file name

    // Reading and setting the values of the working variables

    
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    cout<< "This program will solve the Poisson equation using a Jacobi relaxation!"<< endl;
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;


    do 
    {
        cout<< "Please Insert the Number of grid points N on each axis (>1)"<< endl;
        cin >> N;
    }
    while (N <= 1);

    grid = new double*[N];

    for(i = 0; i<N; i++)
    {
    	grid[i] = new double[N];
    }

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do {
    		cout<< "Please Insert the value of the Lattice dimension L (>0)"<< endl;
    		cin >> L;
       }
       while (L <= 0.);

       a = L/N;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do {
    		cout<< "Please Insert the value of the gamma for the jacobi iteration (>0 and <2)"<< endl;
    		cin >> gamma;
       }
       while ((gamma <= 0.) || (gamma > 2));



    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do {
    		cout<< "Please Insert the value of the requested percentual precision (>0)"<< endl;
    		cin >> prec;
       }
       while (prec <= 0.);

       var = prec *10;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    // Initializing the array of spins

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The lattice grid will now be initialized !"<< endl;
    
    ofstream outputfile;

    outputfile.open(fname, ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append

    for(i = 0; i<N; i++)
    {
    	for(j = 0; j<N; j++)
    	{
        	grid[i][j] = sqrt(a*i)*(a*j); // Giving non random values to the generic point of the grid 

        	if(i == (N -1)) // Applying first boundary condition when x = xmax
            {
            	grid[i][j] = sqrt(1 - pow((1 - (j +1)*a/L),2));
            }

            if(j == (N -1)) // Applying first boundary condition when j = jmax
            {
            	grid[i][j] = (i+1)*a/L;
            }

            //Writing the initial configuration to file

            if( outputfile.is_open() ) // writing to file the values each 100 iterations
            {
            	
            	outputfile<< setprecision(20)  // precision of following floating point output
        		<< setfill(' ')      // character used to fill the column
        		<< left              // left alignment  -- remember all comes to left of affected object here.
        		<< scientific
        		<< setw(30) << i*a << setw(30) << j*a << setw(30) << grid[i][j] << endl;
			} 
			else 
			{
				cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
			}
    	}
    }

    outputfile.close();

    // Computing the starting value of the energy of the configuration
    
    Energy = 0.;

    for(i = 1; i<N-1; i++)
    {
    	for(j = 1; j<N-1; j++)
    	{
        	E_x = Estimate_EX(grid, N, a, i, j);
        	E_y = Estimate_Ey(grid, N, a, i, j);
        	Energy += pow(E_x,2) + pow(E_y,2);
    	}
    }
  
    Energy *=  pow(a,2);

    cout<< setprecision(20)  // precision of following floating point output
        << setfill(' ')      // character used to fill the column
        << left              // left alignment  -- remember all comes to left of affected object here.
        << scientific
        << "The starting value of the energy is : "<< Energy <<" !"<<endl;

	// Starting the cycle with the implementation of the metropolis algorithm

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The numerical simulation is started ! "<< endl;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    for(step=1; var > prec; step ++)
    {
    	for(i = 0; i<N; i++)
    	{
    		for(j = 0; j<N; j++)
    		{
        		
        		if((i*j != 0) && (i != N-1) && (j != N-1))
        	    {
        	    	grid[i][j] = (0.25*(grid[(i+1)%N][j] + grid[(i - 1 + N)%N][j] + grid[i][(j+1)%N] + grid[i][(j - 1 +N)%N]) - grid[i][j])*gamma + grid[i][j];
        	    }

                //Opening a new file if the step is a multiple of 10 and it was not already opened

                if(((step % 10) == 0) && !(outputfile.is_open()))
                {
                	if((step%100) == 0) 
                    {
                    	fname[8] += 1;
                    	fname[9] = '0';
                    }
                    else
                    {
                    	fname[9] += 1;
                    }
                	
                	outputfile.open(fname, ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append
    				cout<< setprecision(20)  // precision of following floating point output
            		<< setfill(' ')      // character used to fill the column
            		<< left              // left alignment  -- remember all comes to left of affected object here.
            		<< scientific
            		<< "We are at iteration Number : "<<setw(30) << step 
            		<<" and the energy is : " << Energy<<endl;
                }

                if((step%10) == 0)
                {
                	if( outputfile.is_open() ) // writing to file the values each 100 iterations
            		{
            	
            			outputfile<< setprecision(20)  // precision of following floating point output
        				<< setfill(' ')      // character used to fill the column
        				<< left              // left alignment  -- remember all comes to left of affected object here.
        				<< scientific
        				<< setw(30) << i*a << setw(30) << j*a << setw(30) << grid[i][j] << endl;
					} 
					else 
					{
						//cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
					}
                }            	
    		}
    	}

    	if(outputfile.is_open())
    	{

    		outputfile.close();

    	}

        


    	// Estimating the new energy configuration of the setting

    	app = Energy;
    	Energy = 0;

    	for(i = 1; i<N-1; i++)
    	{
    		for(j = 1; j<N-1; j++)
    		{
        		E_x = Estimate_EX(grid, N, a, i, j);
        		E_y = Estimate_Ey(grid, N, a, i, j);
        		Energy += pow(E_x,2) + pow(E_y,2);
    		}
    	}

   		Energy *=  pow(a,2);
   		var = abs(Energy - app);

        // If step > 90 better to stop the cycle

        if(step > 250)
        {
        	var = -1;
        }
    }

    // Printing the final results 

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    if(var == -1)
    {
        cout<< "The simulation couldn't reach the imposed value of precision !"<< endl;
    }
    else
    {
    	cout<< setprecision(20)  // precision of following floating point output
            << setfill(' ')      // character used to fill the column
            << left              // left alignment  -- remember all comes to left of affected object here.
            << scientific
            << "The precision was reached after : "<<setw(30) << step <<setw(30)<<" steps !"<<endl;

        if((step%10) > 0)
        {
        	fname[8] = 'N';
        	fname[9] = 'D';
    		outputfile.open(fname, ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append
    		for(i = 0; i<N; i++)
    		{
    			for(j = 0; j<N; j++)
    			{
        		    //Writing the initial configuration to file

            		if( outputfile.is_open() ) // writing to file the values each 100 iterations
            		{
            	
            			outputfile<< setprecision(20)  // precision of following floating point output
        				<< setfill(' ')      // character used to fill the column
        				<< left              // left alignment  -- remember all comes to left of affected object here.
        				<< scientific
        				<< setw(30) << i*a << setw(30) << j*a << setw(30) << grid[i][j] << endl;
					} 
					else 
					{
					cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
					}
				}
			}
    	}
    }

    if( outputfile.is_open())
        {
            outputfile.close();
        }

}


double Estimate_EX(double **grid, int d, double a, int i, int j) 
{
   /* local variable declaration */
   double result = 0.;
   
   result = (0.5/a)*(grid[(i+1)%d][j] - grid[(i-1 +d)%d][j]);

   return result; 
}

double Estimate_Ey(double **grid, int d, double a, int i, int j) 
{

   /* local variable declaration */
   double result = 0.;
   
   result = (0.5/a)*(grid[i][(j+1)%d] - grid[i][(j - 1 +d)%d]);

   return result; 
}
