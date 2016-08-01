using namespace std;
// Select the potential you want:

#define hm_fac 20.73553	// MeV
#define     a0 0.67	// 
#define     r0 1.27	// fm
#define    Vls 0.22	//
// #define      q 1.439978	// e^2
#define      q 1.

// Variables
double  wBox;	// Box width  [fm]
double wWell;	// Well width [fm] (don't use yet)
double     h;	// Mesh width [fm]
double    V0;	// Well depth [MeV] (don't use yet)
double  eMin;	// Min limit  [MeV]
double  eMax;	// Max limit  [MeV]
double eStep;	// Step between different energies in brute force approach, needs getting rid of

double nProton;	// no. of protons for calculating Woods-Saxon
double nNeutron;	// no. of neutrons for calculating Woods-Saxon

double wfStep1;
double wfStep2;

double convEng;	// Convergence energy

double normFac;

// Variable used to check if the wavefunction has crossed x axis
// between different energies
double wfPrev, wfLast, wfTmp;

// Variables for reading in from files
vector<double> variable;
double var;

// Function declarations
double numerovAlgorithm(double E, double f_x, double f_x_h, double r, int isoSpin, int L, int spin);
double converge(double eLo, double eHi, double wfLo, double wfHi, int isoSpin, int L, int spin);
double normalise(double eigenEng, int isoSpin, int L, int spin);
double V(double r, int isoSpin, int L, int spin);
double woodsSaxon(double r);
double spinOrbit(double r, int l, double spin);
double centrifugal(double r, double l);
double coulomb(double r);
double totalMatterDensity(vector<double> density);