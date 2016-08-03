%module shf

%include "typemaps.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_vectora.i"
namespace std
{
  %template(IntVector) vector<int>;
  %template(StrVector) vector<string>;
}

%include "exception.i"
%exception
{
  try
  {
    $action
  }
  catch (const std::runtime_error& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
  catch (const std::logic_error& e) {
    SWIG_exception(SWIG_IndexError, e.what());
  }
}

%{
#define SWIG_FILE_WITH_INIT

#include "utils.h"
#include "Basis.h"
#include "SpBasis.h"
#include "FullSpBasis.h"
#include "ReducedSpBasis.h"
#include "System.h"
#include "NeutronDrop.h"
#include "Nucleus.h"
#include "Solver.h"
#include "HartreeFock.h"
#include "HartreeFockBogo.h"
#include "Interaction.h"
#include "RawInteraction.h"
#include "MinnesotaS0.h"
#include "MinnesotaRaw.h"
%}

%include "armanpy.i"

%include "utils.h"
%include "Basis.h"
%include "SpBasis.h"
%include "FullSpBasis.h"
%include "ReducedSpBasis.h"
%include "System.h"
%include "NeutronDrop.h"
%include "Nucleus.h"
%include "Solver.h"
%include "HartreeFock.h"
%include "HartreeFockBogo.h"
%include "Interaction.h"
%include "RawInteraction.h"
%include "MinnesotaS0.h"
%include "MinnesotaRaw.h"
