using namespace std;
// Select the potential you want:
int selectFunc;		// 1. Infinite square well
			// 2. Finite square well
			// 3. Woods-Saxon potential

// Variables
double  wBox;	// Box width  [fm]
double wWell;	// Well width [fm] (don't use yet)
double     h;	// Mesh width [fm]
double    V0;	// Well depth [MeV] (don't use yet)
double  eMin;	// Min limit  [MeV]
double  eMax;	// Max limit  [MeV]
double eStep;	// Step between different energies in brute force approach, needs getting rid of

int nProton;	// no. of protons for calculating Woods-Saxon
int nNeutron;	// no. of neutrons for calculating Woods-Saxon

double wfStep1;
double wfStep2;

double convEng;	// Convergence energy
//double convEng = 0.5;	// Convergence energy

double hm_fac = 20.75;	// units? MeV?

// Variable used to check if the wavefunction has crossed x axis
// between different energies
double wfPrev, wfLast, wfTmp;

// Variables for reading in from files
vector<double> variable;
double var;

// Function declarations
double numerovAlgorithm(double E, double f_x, double f_x_h, double x);
double converge(double eLo, double eHi, double wfLo, double wfHi);
double normalise(double eigenEng);
double V(double x);
double infSW();
double finSW(double x);
double woodsSaxon(double x);
double spinOrbit(double r, int l, double s);
double centrifugal(double x, double l);
double coulomb(double x);
