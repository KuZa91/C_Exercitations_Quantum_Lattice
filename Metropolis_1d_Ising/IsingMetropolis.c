//===================================================//
// Solution of the Ising Model with Metropolis Alg.
//===================================================//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <cmath>

using namespace std;

double NRG(int *spin, int d, double Beta, double h); // Function for estimating the total energy of the configuration
double SumArray(int *spin, int d); // Function for estimating the sum of spins

int main ()
{
	// Defining the variables of the system

	int64_t N, i = 0, seed = 9;
    double  Beta_tilde, h_tilde, prec, var, Magn_old, Magn_new, energy, nrg_exp, avg_spin, app ;
    int* spins;

    // Reading and setting the values of the working variables

    
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    cout<< "This program will estimate the magnetization in a 1d Ising model, using the Metropolis Algorithm !"<< endl;
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;


    do {
    		cout<< "Please Insert the Number of 1d spins to simulate (>1)"<< endl;
    		cin >> N;
       }
       while (N <= 1);

    spins = new int[N];

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do {
    		cout<< "Please Insert the value of Beta Tilde (>0)"<< endl;
    		cin >> Beta_tilde;
       }
       while (false);

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do {
    		cout<< "Please Insert the value of h Tilde (>0)"<< endl;
    		cin >> h_tilde;
       }
       while (false);

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do {
    		cout<< "Please Insert the value of the requested precision (>0)"<< endl;
    		cin >> prec;
       }
       while (prec <= 0.);

       var = prec *10;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    // Initializing the randomic functions for the program

    mt19937 genrnd(seed); // Use the Mersenne twister algorithm to generate random number

    uniform_real_distribution<double> Unif(0.0,1.0); // Generate a random number between 0 and 1, will be used to initialize randomically the spins to +-1

    // Initializing the array of spins

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The array of spins will now be initialized randomly"<< endl;

    for(i = 0; i<N; i++)
    {
    	if(Unif(genrnd)<0.5)
        {
            spins[i] = -1;
        }
        else
        {
            spins[i] = 1;
        }
    }

    // Printing the starting values of average spin and energy to screen

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    energy = NRG(spins, N, Beta_tilde, h_tilde);
    nrg_exp = log(energy);

    cout    << setprecision(20)  // precision of following floating point output
            << setfill(' ')      // character used to fill the column
            << left              // left alignment  -- remember all comes to left of affected object here.
            << scientific
            << "The starting value of energy of the configuration is : "<<setw(30) << energy<<endl;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    avg_spin = SumArray(spins, N)/N;
    app = avg_spin;

    cout    << setprecision(20)  // precision of following floating point output
            << setfill(' ')      // character used to fill the column
            << left              // left alignment  -- remember all comes to left of affected object here.
            << scientific
            << "The starting value of the average spin is : "<<setw(30) << avg_spin<<endl;


    // Opening the file to save the data

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    
    cout<< "Opening the file for saving the result ! "<< endl;
    
    ofstream outputfile;

    outputfile.open("State_Variables.dat", ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append
                                                         

    

    // Starting the cycle with the implementation of the metropolis algorithm

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The numerical simulation is started ! "<< endl;

    for(i=0; var > prec; i ++)
    {
    	if((nrg_exp - 2*Beta_tilde*(spins[(i+1)%N]*spins[i%N] + spins[(i+1)%N]*spins[i%N]) - 2*h_tilde*spins[i%N]) < nrg_exp) // Checking the acceptance/rejection condition !
    	{
    		nrg_exp =  nrg_exp - 2*Beta_tilde*(spins[(i+1)%N]*spins[i%N] + spins[(i+1)%N]*spins[i%N]) - 2*h_tilde*spins[i%N]; // new energy with the new setup will be given by
    		spins[i] *= -1; // If accepted, flip spin
    	}

    	if(i%100 == 0)
    	{
            avg_spin = SumArray(spins, N)/N;
            energy = exp(nrg_exp);

            if( outputfile.is_open() ) // writing to file the values each 100 iterations
            {
            	
    		  	outputfile<< setprecision(20)  // precision of following floating point output
              	<< setfill(' ')      // character used to fill the column
              	<< left              // left alignment  -- remember all comes to left of affected object here.
              	<< scientific
              	<< setw(30) << i << setw(30) << avg_spin << setw(30) << energy << endl;

  			} 
  			else 
  			{

    			cerr << "Cannot write to file 'State_Variables.dat' " << endl; //cerr use no buffer, direct message regardless of stack turn

  			}

  			if(i>(N/5)) // adaptive condition for the block of the cycle
  			{
  				var = abs((avg_spin - app)/(0.5*(avg_spin+ app)));
  				app = avg_spin;
  			}

  			if(i >= 5*N)
  			{
  				app = -1;
  				var = -1;

  			}

    	}


    }

    // Printing the final results 

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    if(app == -1)
    {
        cout<< "The simulation couldn't reach the imposed value of precision !"<< endl;
    }
    else
    {
    	cout<< setprecision(20)  // precision of following floating point output
            << setfill(' ')      // character used to fill the column
            << left              // left alignment  -- remember all comes to left of affected object here.
            << scientific
            << "The precision was reached after : "<<setw(30) <<i<<setw(30)<<" steps !"<<endl;
    }
    

}

double NRG(int *spin, int d, double Beta, double h) 
{

   /* local variable declaration */
   double result = 0, kin = 0, field = 0;
   int k;
 
   for (k = 0; k<d; k++)
   {
       kin += spin[(k+1)%d]*spin[(k)%d];
       field += spin[k];
   }
   
   result = exp(Beta*kin + h* field);
   return result; 
}

double SumArray(int *spin, int d) 
{

   /* local variable declaration */
   double result = 0;
   int k;
 
   for(k = 0; k<d; k++)
   {
       result += spin[k];
   }

   return result; 
}