
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>

#define PI 3.141592654

// define the number of branches and nodes
#define N 14


// define phi. 
// in this model, we set that one evolutionary event will result in a trait change of a random value from an exponential distribution whose mean is 0.05, see Kutsukake and Innan (2013, Evolution)
float expdev(long *);
double phi_mean = 0.05;

// a function to choose a random value from a poisson distribtion
int poisson_rnd(double);

// a likelihood function
double get_LH(double, double, double);

long int seed;

// phylogeny & phenotype
double phenotype[N];

// set ancestral and descenent nodes for each branch N
int branch_a[N];
int branch_d[N];

double er[N];	//evolutionary rate
double bl[N];	//blanch length



// data (mean and sd) //
// from Toju & Sota 2006 Mol Ecol //

double camelY = 17.70;		//phenotype[3]
double camel_H = 10.70;		//phenotype[4]
double sasan = 3.15;	//phenotype[5]
double dist = 4.41;		//phenotype[6]
double dent = 6.72;		//phenotype[10]
double robus = 5.85;		//phenotype[11]
double hilgen = 5.19;		//phenotype[12]
double sikk = 6.75;		//phenotype[13]


double sd_camelY = 3.06;		//phenotype[3]
double sd_camel_H = 1.38;		//phenotype[4]
double sd_sasan = 0.12;	//phenotype[5]
double sd_dist = 0.60;		//phenotype[6]
double sd_dent = 0.67;		//phenotype[10]
double sd_robus = 0.26;		//phenotype[11]
double sd_hilgen = 0.10;		//phenotype[12]
double sd_sikk = 0.91;		//phenotype[13]

//end




// parameter set - set each parameter and its prior distribution //
double theta;
double MIN_theta = 3.15;
double MAX_theta = 9.;

double mu;
double MIN_mu = 1.;
double MAX_mu = 50.;


double k;
double MIN_k = 0.00001; // should be larger than zero
double MAX_k = 30.;

// end //


double get_LH(double theta, double mu, double k){
	
	int x;
	int mu_plus, mu_minus;
	int i, j;
	
	double ln_lik;
	double lik = 0;
		
	for (i=0; i<N; i++) phenotype[i] = 0;
		
	for (i=0; i<N; i++) {
			
		mu_plus = poisson_rnd( mu*bl[i] );
		mu_minus = poisson_rnd( mu*bl[i] );
		
		
		if (i==2) {
			for (x=0; x<mu_plus*k; x++) {
				phenotype[i] += -log(drand48())*phi_mean;
			}
			for (x=0; x<mu_minus/k; x++) {
				phenotype[i] -= -log(drand48())*phi_mean;
			}
		}
		
		else {
			for (x=0; x<mu_plus; x++) {
				phenotype[i] += -log(drand48())*phi_mean;
			}
			for (x=0; x<mu_minus; x++) {
				phenotype[i] -= -log(drand48())*phi_mean;
			}
		}
		
		for (j=0; j<N; j++) {
			if (branch_d[i]==branch_a[j]) {
				phenotype[j] += phenotype[i];
			}
		}
		
	}
	
	
	for (i=0; i<N; i++) phenotype[i] += theta;

	// calculate likelihood here, assuming a normal distribution for traits
	ln_lik = 
	log((1/(sqrt(2*PI)*sd_camelY)) * (exp( -(phenotype[3]-camelY)*(phenotype[3]-camelY) / (2*sd_camelY*sd_camelY) )) )
	+ log((1/(sqrt(2*PI)*sd_camel_H)) * (exp( -(phenotype[4]-camel_H)*(phenotype[4]-camel_H) / (2*sd_camel_H*sd_camel_H) )) )
	+ log((1/(sqrt(2*PI)*sd_sasan)) * (exp( -(phenotype[5]-sasan)*(phenotype[5]-sasan) / (2*sd_sasan*sd_sasan) ))) 
	+ log((1/(sqrt(2*PI)*sd_dist)) * (exp( -(phenotype[6]-dist)*(phenotype[6]-dist) / (2*sd_dist*sd_dist) )))
	+ log((1/(sqrt(2*PI)*sd_dent)) * (exp( -(phenotype[10]-dent)*(phenotype[10]-dent) / (2*sd_dent*sd_dent) ))) 
	+ log((1/(sqrt(2*PI)*sd_robus)) * (exp( -(phenotype[11]-robus)*(phenotype[11]-robus) / (2*sd_robus*sd_robus) ))) 
	+ log((1/(sqrt(2*PI)*sd_hilgen)) * (exp( -(phenotype[12]-hilgen)*(phenotype[12]-hilgen) / (2*sd_hilgen*sd_hilgen) ))) 
	+ log((1/(sqrt(2*PI)*sd_sikk)) * (exp( -(phenotype[13]-sikk)*(phenotype[13]-sikk) / (2*sd_sikk*sd_sikk) )));	
	
	
	return ln_lik;

 }


// a function of returning a random value from a Poisson distribution whose mean is lambda
int poisson_rnd(double lambda){
	int m;
	
	lambda = lambda+log(drand48());
	m=0;
	while (lambda>0){
		lambda += log(drand48());
		m++;
	}
	return m;
}



int main(void)
{
	srand48(time(NULL));
	
	struct tm result;
	time_t timep;
	struct timeb timeb;
	FILE *fp;
	char fname[80];
	
	ftime(&timeb);
	timep = timeb.time;
	localtime_r(&timep, &result);
	
	// here an output file will be created whose name can be set here.  
	snprintf(fname, 80, "ABC_weevil_sel_%02d%02d%02d_%02d%02d%02d.txt",
			 result.tm_year-100,result.tm_mon+1, result.tm_mday, result.tm_hour, result.tm_min,result.tm_sec);
	fp = fopen(fname, "w");
	
	freopen (fname, "a", fp);		
	
	// a title of this program
	time_t timer;
	struct tm *t_st;
	time(&timer);
	printf("\n*** %s*** PCM ABC weevil model ***\n\n", ctime(&timer));
	//end
	
	// phylogeny data //
	bl[0]=9.62; // (MRCA, 0)
	bl[1]=7.9; // (0, 1)
	bl[2]=6.53; // (1, 2)
	bl[3]=0.25; // (2, 3)
	bl[4]=0.25; // (2, 4)
	bl[5]=6.78; // (1, 5)
	bl[6]=14.68; // (0, 6)
	bl[7]=5.96; // (MRCA, 7)
	bl[8]=4.74; // (7,8)
	bl[9]=3.27; // (8, 9)
	bl[10]=10.35; // (9, 10)
	bl[11]=10.35; // (9, 11)
	bl[12]=13.62; // (8, 12)
	bl[13]=18.36; // (7, 13)
	
	branch_a[0]=-99; branch_d[0]=0; // -99 is a dummy number
	branch_a[1]=0; branch_d[1]=1; 
	branch_a[2]=1; branch_d[2]=2;
	branch_a[3]=2; branch_d[3]=3; 
	branch_a[4]=2; branch_d[4]=4; 
	branch_a[5]=1; branch_d[5]=5;
	branch_a[6]=0; branch_d[6]=6; 
	branch_a[7]=-99; branch_d[7]=7;
	branch_a[8]=7; branch_d[8]=8; 
	branch_a[9]=8; branch_d[9]=9; 
	branch_a[10]=9; branch_d[10]=10; 
	branch_a[11]=9; branch_d[11]=11; 
	branch_a[12]=8; branch_d[12]=12; 
	branch_a[13]=7; branch_d[13]=13; 
	
	
	
	int i,j;
	double temp_LH;
	double pro = 6.;	// this is an adujusting value for efficient computation because model likelihood is far less than a random number from U[0,1]
	double bar;			// by comparing a (log-)likelihood to this random value from U[0,1], acceptance rate will be proportional to a likelihood.
	int count = 0;		// # accepted data
	
	// using this value, we check an unrealistic case that the phenotypic value becomes negative by simulation. 
	int negative_check;

	
	printf("theta\tk\tmu\tP[0]\tP[1]\tloglik\tbar\n");
	fprintf(fp, "theta\tk\tmu\tP[0]\tP[1]\tloglik\tbar\n");
	

	for (i=0;; i++) {
		
		// propose parameters from a prior distribution for each parameter
		theta = (MAX_theta-MIN_theta)*drand48()+MIN_theta;
		k = (MAX_k - MIN_k)*drand48()+MIN_k;
		mu = (MAX_mu - MIN_mu) * drand48() + MIN_mu;
		
		
		bar = drand48();
		
		// calculate a likelihood under a simualted dataset
		temp_LH = get_LH(theta, mu, k);

		
		// initialization of negative_check
		negative_check = 0;

		// check here whether the simulated data include a negative value
		for (j=0; j<N; j++) {
			if (phenotype[j] < 0 || phenotype[j]==0) negative_check++;
		}

		// accept if all simulated data is positive and likelihood is larger than bar
		if (negative_check==0 & log(bar)-pro<temp_LH) {
			printf("%2.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t%.2lf\t\n", theta, k, mu, phenotype[0], phenotype[1], temp_LH, bar);
			fprintf(fp, "%2.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t%.2lf\t\n", theta, k, mu, phenotype[0], phenotype[1], temp_LH, bar);
			count++;
			
		}
		
		// end simulation when an accepted data reaches to 100. In this example, an acceptance of one parameter takes about one hour. 
		if (count==100) break;
	}
	
	fclose (fp);
	
		
}
