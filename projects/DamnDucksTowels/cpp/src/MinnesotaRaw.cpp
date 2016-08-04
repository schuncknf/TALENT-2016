#include <iostream>
#include <iomanip>
#include <fstream>

#include "global.h"
#include "quadrature.h"
#include "MinnesotaRaw.h"

MinnesotaRaw::MinnesotaRaw(FullSpBasis &_basis, int _nParticleTypes, std::string _filename) :
  Interaction( _basis, _nParticleTypes ),
  TBME(_basis.size, _basis.size),
  nParticleTypes(_nParticleTypes)
{
///////////////////////////////////////////////
  for (int state0 = 0; state0 < _basis.size;    state0++)
    for (int state1 = 0; state1 < _basis.size;    state1++)
    {
      TBME(state0, state1) = arma::zeros(_basis.size, _basis.size);
    }

/////////////////////////////////////////////////
  int nMax = _basis.nMax;
  int lMax = _basis.lMax;
  int size_nljmt = 0;

  for (int n = 0; n <= nMax; n++)
    for (int l = 0; l <= lMax; l++)
      for (int _2j = 2 * l + 1; (_2j >= 2 * l - 1) && (_2j > 0); _2j -= 2)
        for (int _2t = -1; _2t <= 1; _2t += 2)
          for (int _2m = -_2j; _2m <= _2j; _2m += 2)
            size_nljmt++;

  arma::ivec lst(size_nljmt);
  arma::ivec vec_l(size_nljmt + 1);
  arma::ivec vec_2j(size_nljmt + 1);
  arma::ivec vec_2m(size_nljmt + 1);
  arma::ivec vec_2t(size_nljmt + 1);
  int i_nljmt = 0;
  int i_nlj = 0;

  for (int n = 0; n <= nMax; n++)
    for (int l = 0; l <= lMax; l++)
      for (int _2j = 2 * l + 1; (_2j >= 2 * l - 1) && (_2j > 0); _2j -= 2)
      {
        i_nlj++;
        {
          for (int _2t = -1; _2t <= 1; _2t += 2)
            for (int _2m = -_2j; _2m <= _2j; _2m += 2)
            {
              i_nljmt++;
              lst(i_nljmt - 1) = i_nlj - 1;
              vec_l(i_nljmt) = l;
              vec_2j(i_nljmt) = _2j;
              vec_2m(i_nljmt) = _2m;
              vec_2t(i_nljmt) = _2t;
            }
        }
      }

  //int a, b, c, d;
  std::ifstream input(_filename.c_str());

  if (!input)
  {
    throw std::runtime_error("The file doesn't exist !!!!!!!!!!!!!!!!");
  }

  std::string line;
  std::getline(input, line);
  std::getline(input, line);

  while (std::getline(input, line))
  {
    int a = std::stoi(line.substr(1, 3));
    if (vec_2t(a) != 1) continue;
    int c = std::stoi(line.substr(9, 3));
    if (vec_2t(c) != 1) continue;
    if (vec_l(a) != vec_l(c)) continue;
    if (vec_2j(a) != vec_2j(c)) continue;
    if (vec_2m(a) != vec_2m(c)) continue;
    int b = std::stoi(line.substr(5, 3));
    if (vec_2t(b) != 1) continue;
    int d = std::stoi(line.substr(13, 3));
    if (vec_2t(d) != 1) continue;
    if (vec_l(b) != vec_l(d)) continue;
    if (vec_2j(b) != vec_2j(d)) continue;
    if (vec_2m(b) != vec_2m(d)) continue;
    double V = std::stod(line.substr(21));
    TBME(lst(a - 1), lst(b - 1))(lst(c - 1), lst(d - 1)) += V / (vec_2j(a) + 1) / (vec_2j(b) + 1);
    
  }
}

MinnesotaRaw::MinnesotaRaw(TotalFullSpBasis &_basis, int _nParticleTypes, std::string _filename) :
  Interaction( _basis, _nParticleTypes ),
  TBME(_basis.size, _basis.size),
  nParticleTypes(_nParticleTypes)
{
///////////////////////////////////////////////
  for (int state0 = 0; state0 < _basis.size;    state0++)
    for (int state1 = 0; state1 < _basis.size;    state1++)
    {
      TBME(state0, state1) = arma::zeros(_basis.size, _basis.size);
    }
/////////////////////////////////////////////////
  int nMax = _basis.nMax;
  int lMax = _basis.lMax;
  int size_nljmt = 0;

  for (int n = 0; n <= nMax; n++)
    for (int l = 0; l <= lMax; l++)
      for (int _2j = 2 * l + 1; (_2j >= 2 * l - 1) && (_2j > 0); _2j -= 2)
        for (int _2t = -1; _2t <= 1; _2t += 2)
          for (int _2m = -_2j; _2m <= _2j; _2m += 2)
            size_nljmt++;

  arma::ivec lst(size_nljmt);
  arma::ivec vec_l(size_nljmt + 1);
  arma::ivec vec_2j(size_nljmt + 1);
  arma::ivec vec_2m(size_nljmt + 1);
  arma::ivec vec_2t(size_nljmt + 1);
  int i_nljmt = 0;
  int i_nljm = 0;
  for (int n = 0; n <= nMax; n++)
    for (int l = 0; l <= lMax; l++)
      for (int _2j = 2 * l + 1; (_2j >= 2 * l - 1) && (_2j > 0); _2j -= 2)
        for (int _2t = -1; _2t <= 1; _2t += 2)
	  for (int _2m = -_2j; _2m <= _2j; _2m += 2)
	  {
	    if (_2t == 1)
	      i_nljm++;
	    i_nljmt++;
	    lst(i_nljmt - 1) = i_nljm - 1;
	    vec_l(i_nljmt) = l;
	    vec_2j(i_nljmt) = _2j;
	    vec_2m(i_nljmt) = _2m;
	    vec_2t(i_nljmt) = _2t;
	  }

  int a, b, c, d;
  double V;
  std::ifstream input(_filename.c_str());
  if (!input)
  {
    throw std::runtime_error("The file doesn't exist !!!!!!!!!!!!!!!!");
  }
  std::string line;
  getline(input, line);
  getline(input, line);
  //std::cout << "Reading..." << std::endl;

  while (input >> a >> b >> c >> d >> V)
  {

    if ((vec_2t(a) == 1) && (vec_2t(b) == 1) && (vec_2t(c) == 1) && (vec_2t(d) == 1))
    {
      TBME(lst(a-1), lst(b-1))(lst(c-1), lst(d-1)) = V;
    }
  }
}

MinnesotaRaw::~MinnesotaRaw()
{
}

double MinnesotaRaw::get(arma::field<arma::mat> &R, arma::ivec &bType, arma::ivec &bId, arma::ivec &kType, arma::ivec &kId)
{
  (void) R;

  if (bType.n_elem != 2 || kType.n_elem != 2 || bId.n_elem != 2 || kId.n_elem != 2)
  {
    return 0.;
  }

  return TBME(bId(0), bId(1))(kId(0), kId(1));
}

std::string MinnesotaRaw::info ()
{
  return std::string();
}


std::string MinnesotaRaw::toString ()
{
  return std::string();
}
