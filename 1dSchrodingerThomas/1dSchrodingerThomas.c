//===================================================//
// Solving the Schrodinger equation using Thomas Alg.
//===================================================//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>

using namespace std;


void solve(complex<double> *a, complex<double> *b, complex<double> *c, complex<double> *d, int64_t n); // Function used to solve a tridiagonal matrix taken from https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm


int main ()
{
	// Defining the variables of the system

	int64_t N,cnt, step = 1;
    double  T, T_end, Delta_t,  omega, H_const, m,L,a, pi = 2 * acos(0.0), Norm = 0.;
    double* x;
    complex<double> app;
    complex<double> *waveform, *a_diag, *b_diag, *c_diag, *app_a, *app_b, *app_c, *d_vec, *app_d;
    char fname[] = "WaveEvol00.dat"; //modify fname[8], fname[9] to change file name

    //Defining the complex unit i

    complex<double> i = complex<double>(0,1.);
    

    // Reading and setting the values of the working variables

    
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    cout<< "This program will solve the Schrodinger eqt. using a Thomas Algorith for the tridiagonal matrix!"<< endl;
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    
    do 
    {
   		cout<< "Please Insert the value of the armonic constant omega squared (>0)"<< endl;
   		cin >> omega;
    }
    while (omega <= 0.);   

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
   		cout<< "Please Insert the value of the normalized mass of the particle (>0 sugg 1.)"<< endl;
   		cin >> m;
    }
    while (m <= 0.);   

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;


    do 
    {
   		cout<< "Please Insert the value of the total normalized lenght of the waveform (>0)"<< endl;
   		cin >> L;
    }
    while (L <= 0.);   

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        cout<< "Please Insert the Number of grid points for the waveform Lenght (>L No scientific notation)"<< endl;
        cin >> N ;
    }
    while (N <= int64_t(L));

    a = L/N;

    // Initializing all the vectors to the input dimension

    waveform = new complex<double>[N];
    x = new double[N];
    a_diag = new complex<double>[N];  
    b_diag = new complex<double>[N];
    c_diag = new complex<double>[N];
    app_a = new complex<double>[N];
    app_b = new complex<double>[N];
    app_c = new complex<double>[N];
    d_vec = new complex<double>[N];
    app_d = new complex<double>[N];

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do {
    		cout<< "Please Insert the ending normalized time T for the evolution (>0)"<< endl;
    		cin >> T_end;
       }
       while (T_end <= 0.);

       T = 0.;


    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do {
    		cout<< "Please Insert the value of the normalized time step interval (>T/100000 and <T)"<< endl;
    		cin >> Delta_t;
       }
       while ((Delta_t <= T_end/100000) || (Delta_t > T_end));

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    // Initializing the starting waveform as a gaussian distribtion and the file where the waveform would be saved

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The Initial waveform will now be initialized !"<< endl;
    
    ofstream outputfile;

    outputfile.open(fname, ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append

    for(cnt = 0; cnt<N; cnt++)
    {
        x[cnt] = cnt*a -(L/2.);
    	waveform[cnt] = /*(1./sqrt(2.*pi))*exp(-0.5*pow(x[cnt],2.))*/ (1./L); // Choose between an uniform starting condition and a gaussian starting condition
    	if( outputfile.is_open() ) // writing to file the values each 100 iterations
    	{     	
       		outputfile<< setprecision(20)  // precision of following floating point output
   			<< setfill(' ')      // character used to fill the column
   			<< left              // left alignment  -- remember all comes to left of affected object here.
   			<< scientific
   			<< setw(30) << x[cnt] << setw(30) << real(waveform[cnt]) << setw(30) << imag(waveform[cnt]) << endl;

   			Norm += a*sqrt((pow(real(waveform[cnt]),2.) + pow(imag(waveform[cnt]),2.)));
		} 
		else 
		{
			cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
		}

	}
	cout<< "The starting value for the norm of the waveform is "<< Norm << " ! "<< endl;



	outputfile.close();

	// Starting the cycle with the implementation of the Thomas algorithm for solving the Tridiagonal with PBC matrix

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The numerical simulation has started ! "<< endl;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    // Defining the constant term of the matrix for the simulation

    H_const = (0.25*Delta_t)/(m*pow(a,2.));

    for(T=0; T < T_end; step ++)
    {
    	for(cnt = 0; cnt<N; cnt++)
    	{  

    		// Initializing the vectors needed to solve the Thomas algorithm

            d_vec[cnt] = waveform[cnt] + i*H_const*(waveform[(cnt-1+N)%N] + waveform[(cnt+1)%N] - 2.*waveform[cnt]) - 0.25*i*Delta_t*m*omega*pow(x[cnt],2.)*waveform[cnt];
            b_diag[cnt] = 1. + 2.*H_const*i + 0.25*i*Delta_t*m*omega*pow(x[cnt],2.);
            app_d[cnt] = 0. ;

            if((cnt != (N-1))&&(cnt != 0)) 
            {
            	a_diag[cnt] = -H_const*i;
            	c_diag[cnt] = -H_const*i;
            }
            else
            {
            	if(cnt == (N-1))
            	{
            		a_diag[cnt] = -H_const*i;
            		c_diag[cnt] = 0.;
            	}
            	else
            	{
            		a_diag[cnt] = 0.;
            		c_diag[cnt] = -H_const*i;
            	}
            }

            // We need to setup a copy of our vectors as our solver will change them

            app_a[cnt] = a_diag[cnt];
            app_b[cnt] = b_diag[cnt];
            app_c[cnt] = c_diag[cnt];

    	}

    	// Adding the PBC to the matrix

    	app_d[0] = H_const*i;
    	app_d[N-2] = H_const*i;
    	app = -app_d[0]; // Last row PBC saved for later uses

    	// Solving the system with Thomas algorithm ,first we solve for the (n-1)*(n-1) matrix and then we solve the system for the last column

    	solve(a_diag, b_diag, c_diag, d_vec, (N-1));
    	solve(app_a, app_b, app_c, app_d, (N-1));

    	app_d[N - 1] = (d_vec[N -1] - (app*d_vec[0]) - (app * d_vec[N-2]) )/(b_diag[N - 1] + (app*app_d[0]) + (app * app_d[N -2]));

    	for(cnt=0; cnt < N; cnt ++)
    	{
    		if(cnt != (N-1))
    		{
    			waveform[cnt] = d_vec[cnt] + (app_d[cnt]*app_d[N-1]);
    		}
    		else
    		{
    			waveform[cnt] = app_d[N-1]; 
    		}
    		
    	}

    	//Opening a new file if the step is a multiple of 10 and it was not already opened


        if((step%100) == 0)
        {
          	if((step%1000) == 0) 
            {
               	fname[8] += 1;
               	fname[9] = '0';
            }
            else
            {
               	fname[9] += 1;
            }
    		outputfile.open(fname, ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append

    		Norm = 0.;
    		
    		for(cnt = 0; cnt<N; cnt++)
    		{
    			if( outputfile.is_open() ) // writing to file the values each 100 iterations
    			{     	
       				outputfile<< setprecision(20)  // precision of following floating point output
   					<< setfill(' ')      // character used to fill the column
   					<< left              // left alignment  -- remember all comes to left of affected object here.
   					<< scientific
   					<< setw(30) << x[cnt] << setw(30) << real(waveform[cnt]) << setw(30) << imag(waveform[cnt]) << endl;

   					Norm += a*sqrt((pow(real(waveform[cnt]),2.) + pow(imag(waveform[cnt]),2.)));
				} 
				else 
				{
					cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
				}
    		}

    		cout<< "We are at time "<< T <<" for the simulation and the value for the norm of the waveform is "<< Norm << " ! "<< endl;

    		outputfile.close();

        }

    	//Increasing the time evolution

    	T += Delta_t;

    }

    // Printing the final results 

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

 	// Writing the final results to file

 	fname[8] = 'N';
    fname[9] = 'D';
    outputfile.open(fname, ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append
    Norm = 0;

    for(cnt = 0; cnt<N; cnt++)
    {
    	//Writing the initial configuration to file

       	if( outputfile.is_open() ) // writing to file the values each 100 iterations
  		{
           	
           	outputfile<< setprecision(20)  // precision of following floating point output
   			<< setfill(' ')      // character used to fill the column
   			<< left              // left alignment  -- remember all comes to left of affected object here.
   			<< scientific
   			<< setw(30) << x[cnt] << setw(30) << real(waveform[cnt]) << setw(30) << imag(waveform[cnt]) << endl;

   			Norm += a*sqrt((pow(real(waveform[cnt]),2.) + pow(imag(waveform[cnt]),2.)));
		} 
		else 
		{
			cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
		}
		    	
    }

    cout<< "The ending value for the norm of the waveform is "<< Norm << " ! "<< endl;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    outputfile.close();

    cout<< "The numerical simulation is over and the results may be readed from the file ! "<< endl;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
}

	
void solve(complex<double> *a, complex<double> *b, complex<double> *c, complex<double> *d, int64_t n) 
{
    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
    Finally we have a triangular matrix:
    |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report results
    Written by Keivan Moradi, 2014
    */
    int64_t i,j;
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}