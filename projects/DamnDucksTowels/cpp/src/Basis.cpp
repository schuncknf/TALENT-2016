#include <sstream>
#include <string>
#include <vector>
#include "Basis.h"

Basis::Basis(std::string _type, std::vector<std::string> _qNames) : type(_type), qNames(_qNames), qNumSize(_qNames.size())
{
}

Basis::~Basis()
{
}

std::string Basis::info()
{
  std::string info("");
  return info;
}

std::string Basis::toString()
{
  std::stringstream info;

  //Printing the names of quantum numbers
  for (int j = 0; j < qNumSize; j++)
    info << "\t" << qNames[j];

  info << "\n";

  //Printing quantum number for each basis state
  for (int i = 0; i < size; i++)
  {
    info << i + 1;

    for (int j = 0; j < qNumSize; j++)
    {
      info << "\t" << qNumbers(i, j);
    }

    info << "\n";
  }

  //Printing basis size calculated independently
  info << "Basis size:\t" << size << std::endl;
  return info.str();
}
