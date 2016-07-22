#include <iostream>

using namespace std;

int main()
{
  for (double r = -1; r < 1; r += 0.1)
    cout << r << "\t" << radialWaveFunction(1,0,r) << endl; 
  return 0;
}