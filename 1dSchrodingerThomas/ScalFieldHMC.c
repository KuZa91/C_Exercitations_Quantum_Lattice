//===================================================//
// Solving the Schrodinger equation using Thomas Alg.
//===================================================//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <random>

using namespace std;


double NRG(double *impulse, double *field, double a_tau, double m, double lambda, int64_t n); // Function used to estimate the energy of the field configuration given the field and the impulses of the hamiltonian

int main ()
{
	// Defining the variables of the system

	int64_t N, Nt, N_upd, N_burn, N_AutoCorr, cnt, upd, acc, i, n_it, check_stuck, seed = 9;
    double lambda, m, dt, tau, a_tau, phi_avg, phi_var, Energy, app;
    bool updated = false;
    double *oldfield, *newfield, *pi, *phi_corr ;
    char fname[] = "CorrConfg0000.dat"; //modify fname[9], fname[10], fname[11], fname[12] to change file name

    // Initializing the seed and the normal distribution

    mt19937 genrnd(seed);

    normal_distribution<> Gauss{0,1};

    uniform_real_distribution<double> Unif(0.0,1.0); // Generate a random number between 0 and 1, will be used to initialize randomically the spins to +-1
    

    // Reading and setting the values of the working variables

    
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    cout<< "This program will simulate the evolution of a 1d scalar field using HMC !"<< endl;
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    
    do 
    {
   		cout<< "Please Insert the value of the coupling constant lambda (>0 sugg 24)"<< endl;
   		cin >> lambda;
    }
    while (lambda <= 0.);   

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
   		cout<< "Please Insert the value of the normalized mass of the field (>0 sugg 1.)"<< endl;
   		cin >> m;
    }
    while (m <= 0.);   

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        cout<< "Please Insert the Number of steps in euclidean time tau (>0 sugg 12)"<< endl;
        cin >> N ;
    }
    while (N <= 0);

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        cout<< "Please insert the final tau for the evolution of the field in euclidean time (>0 sugg 1.)"<< endl;
        cin >> tau;
    }
    while (tau <= 0.);

    a_tau = tau/N;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
   		cout<< "Please insert the delta t interval in MC time (>0 sugg 0.001)"<< endl;
   		cin >> dt;
    }
    while (dt <= 0.);   

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        cout<< "Please Insert the Number of delta t evolution step to evolve the field in MC time(> 0 sugg 1000)"<< endl;
        cin >> Nt ;
    }
    while (N <= 0);


    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        cout<< "Please Insert the Number of requested accepted updates to the field "<< endl;
        cin >> N_upd ;
    }
    while (N_upd <= 0);

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        cout<< "Please Insert the Number of step used for burning time during the simulation (>0 sugg 10000)"<< endl;
        cin >> N_burn ;
    }
    while (N_burn <= 0);

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        cout<< "Please Insert the Number of accepted samples to drop due to autocorrelation (0 for nor autocorr, sugg 20)"<< endl;
        cin >> N_AutoCorr ;
    }
    while (N_AutoCorr < 0);



    // Initializing all the vectors to the input dimension

    oldfield = new double[N];
    newfield = new double[N];
    pi = new double[N];
    phi_corr = new double[N];

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    // Initializing the starting field and momenta from a normal distribution

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The Initial field and momenta will now be initialized using a gaussian variable !"<< endl;
    
    ofstream outputfile;


    for(cnt = 0; cnt<N; cnt++)
    {
        oldfield[cnt] = Gauss(genrnd);
        newfield[cnt] = oldfield[cnt];
        phi_corr[cnt] = 0.;
    	
	}

	// Starting the cycle with the implementation of the HMC algorythm to draw samples of the field

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The numerical simulation has started ! "<< endl;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    // Initializing the counters for the cycle

    upd = 0;
    acc = 0;
    n_it = 0;
    phi_avg = 0;
    phi_var = 0;
    Energy = -1;
    check_stuck = 0;

    do
    {
        //Evolving pi to the middle step in order to apply leap frog algorithm 

        for(cnt = 0; cnt <N; cnt ++)
        {
        	pi[cnt] = Gauss(genrnd);
            pi[cnt] = pi[cnt] - (dt/4.)*((-2.*oldfield[(cnt - 1 + N)%N] +4.*oldfield[cnt] - 2.*oldfield[(cnt + 1)%N])/(a_tau*a_tau) + 2.*oldfield[cnt]*pow(m,2.) +(lambda*lambda/3.)*pow(oldfield[cnt],3.));

        }

        for(i =0; i <Nt; i ++)
        {
            for(cnt = 0; cnt <N; cnt ++)
            {
                newfield[cnt] = newfield[cnt] + dt*pi[cnt];
                pi[cnt] = pi[cnt] - (dt/2.)*((-2.*newfield[(cnt - 1 + N)%N] +4.*newfield[cnt] - 2.*newfield[(cnt + 1)%N])/(a_tau*a_tau) + 2.*newfield[cnt]*pow(m,2.) +(lambda*lambda/3.)*pow(newfield[cnt],3.));
            }
            
        }

        n_it = n_it + 1;

        // Checking if burning time is over to proceed with the analysis

        if(n_it > N_burn)
        {

        	// If that's the first after-burning evolution of the field, save the first value of energy to compare with others

        	if(Energy == -1)
        	{
        		Energy = NRG(pi, newfield, a_tau, m, lambda, N);
        	}

        	// Count the number of iterations to count if stuck

        	check_stuck += 1;

        	// Accept, reject proposal

   
        	if((NRG(pi, newfield, a_tau, m, lambda, N) <= Energy)||(Unif(genrnd) <= exp(-(NRG(pi, newfield, a_tau, m, lambda, N)- Energy) ))||(check_stuck == N_burn*10))
        	{

        		cout<< "We passed the burning time, the value of the energy for this sample is "<< NRG(pi, newfield, a_tau, m, lambda, N) << " !"<<endl;
        		acc +=1;
        		check_stuck = 0;  
        		
        		Energy = NRG(pi, newfield, a_tau, m, lambda, N);
        		
        		

            	for(cnt = 0; cnt<N; cnt++)
    			{
        			oldfield[cnt] = newfield[cnt];
				}

        		// saving the step each N_autocorr to drop correlated samples
        		if((N_AutoCorr == 0)||(acc%N_AutoCorr == 0))
        		{
        			
        			upd += 1;
        			
        			// Updating the name of the saved waveform file, works up to 10000 files with different names

        			if((upd%1000) == 0)
        			{

        				fname[9] += 1;
               	    	fname[10] = '0';
               	    	fname[11] = '0';
               	    	fname[12] = '0';          			
    				}
    				else
    				{
    					if((upd%100) == 0)
    					{
    						fname[10] += 1;
               	    		fname[11] = '0';
               	    		fname[12] = '0'; 
    					}
    					else
    					{
    						if((upd%10) == 0)
    						{
    							fname[11] += 1;
               	    			fname[12] = '0';
    						}
    						else
    						{
    							fname[12] += 1;
    						}
    					}
    				}

    		    	outputfile.open(fname, ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append    		
    		
    				for(cnt = 0; cnt<N; cnt++)
    				{
    					// Updating the values of the correlators 

    					phi_avg += (1./N)*oldfield[cnt];
    					phi_var += (1./N)*pow(oldfield[cnt],2.);

    					if (cnt == 0)
    					{
    						app = pow(oldfield[cnt],2.);
    					}
    					else
    					{
    						app = oldfield[cnt]*oldfield[0];
    					}

    					phi_corr[cnt] += app;

    					if( outputfile.is_open() ) // writing to file the values each 100 iterations
    					{     	
       						outputfile<< setprecision(20)  // precision of following floating point output
   							<< setfill(' ')      // character used to fill the column
   							<< left              // left alignment  -- remember all comes to left of affected object here.
   							<< scientific
   							<< setw(30) << cnt*a_tau << setw(30) << app << endl;

   						} 
						else 
						{
							cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
						}
    				}    		

    				outputfile.close();
    			}

        	}
        	else
        	{
        		for(cnt = 0; cnt<N; cnt++)
    			{
        			newfield[cnt] = oldfield[cnt];
				}
        	}

        }
        else
        {
        	for(cnt = 0; cnt<N; cnt++)
    		{
        		oldfield[cnt] = newfield[cnt];
			}
        }

        if((((n_it- N_burn)%1000) == 0) && (n_it > N_burn))
        {
        	cout<< "After "<< (n_it- N_burn) << " iterations we have "<< upd << " accepted proposals ! "<<endl;
        }


    }while(upd <= N_upd);

    // Estimating and reporting the final results for the measurements

    phi_avg /= upd;

    phi_var /= upd;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    cout<< "The numerical simulation is over and the sampled configurations, as well as the values of the last correlator, may be readed from the files ! "<< endl;

    cout<< "The final obtained value for the average of the fields over the samples is  "<< phi_avg << " ! "<< endl;

    cout<< "The final obtained value for the average of the squared fields over the samples is  "<< phi_var << " ! "<< endl;

    outputfile.open("avg_correlator.dat", ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append

    for(cnt = 0; cnt<N; cnt++)
    {
    	if( outputfile.is_open() ) // writing to file the values each 100 iterations
    	{     	
  			outputfile<< setprecision(20)  // precision of following floating point output
   			<< setfill(' ')      // character used to fill the column
   			<< left              // left alignment  -- remember all comes to left of affected object here.
  			<< scientific
   			<< setw(30) << cnt*a_tau << setw(30) << (1./upd)*phi_corr[cnt]<< endl;

   		} 
		else 
		{
			cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
		}
    } 
    outputfile.close();



    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
}

double NRG(double *impulse, double *field, double a_tau, double m, double lambda, int64_t n)
{

   /* local variable declaration */
   double result = 0;
   int k;
 
   for(k = 0; k<n; k++)
   {
       result += 0.5*pow(impulse[k],2.) + 0.5*a_tau*(pow((field[(k +1)%n] - field[k])/(a_tau),2.) + pow(m*field[k],2.) + (1./12.)*pow(lambda*pow(field[k],2.),2.));
   }

   return result; 
}