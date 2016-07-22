#include <iostream>
#include "Basis.h"

Basis::Basis(std::string _type, std::vector<std::string> _qNames)
{
  type = _type;
  qNames = _qNames;
  qNumSize = qNames.size();
}

Basis::~Basis()
{
}

void Basis::printBasis()
{
  //Printing the names of quantum numbers
  for (int j = 0; j < qNumSize; j++)
      std::cout << "\t" << qNames[j];
  std::cout << "\n";
  //Printing quantum number for each basis state
  for (int i = 0; i < basisSize; i++)
  {
    std::cout << i+1;
    for (int j = 0; j < qNumSize; j++)
    {
      std::cout << "\t" << qNumbers(i,j);
    }
    std::cout << "\n";
  }
  //Printing basis size calculated independently
  std::cout << "Basis size:\t" << basisSize << std::endl;
}