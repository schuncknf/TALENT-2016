#include "HartreeFock.h"

HartreeFock::HartreeFock(System &_system) : Solver(_system), D(_system.particleNumbers.n_rows), occ(_system.particleNumbers.n_rows)
{
  int basisSize = _system.basis->size;

  for (unsigned int pType = 0; pType < _system.particleNumbers.n_rows; pType++)
  {
    D(pType).eye(basisSize, basisSize);
    occ(pType) = arma::zeros<arma::vec>(basisSize);

    for (int state = 0; state < _system.particleNumbers(pType); state++)
    {
      occ(pType)(state) = 1;
    }
  }
}

HartreeFock::~HartreeFock()
{
}

void HartreeFock::run()
{
  int nb_state = system->basis->size;
  arma::field<arma::mat> &R = system->R;

  for (unsigned int pType = 0; pType < system->particleNumbers.n_rows; pType++)
  {
    arma::mat H = system->Kinetic(pType) + system->Gamma(pType);
    // Temporary vectors and matrices to store eigenvecs and energies
    arma::vec old_indivE = indivEnergies.row(pType).t();
    arma::vec new_indivE = arma::zeros<arma::vec>(nb_state);

    // Hamiltonian diagonalization to extract D and e
    if ((system->basis->type == "FullSpBasis") && (dynamic_cast<MinnesotaRaw *>(system->inter) != NULL))
    {
      FullSpBasis &fullSpBasis = static_cast<FullSpBasis &>(*system->basis);

      for (int l = 0; l < fullSpBasis.lMax + 1; l++)
      {
        for (int _2j = 2 * l + 1; (_2j >= 2 * l - 1 ) && (_2j > 0); _2j -= 2)
        {
          std::vector<arma::uword> v;

          for (unsigned int i = 0; i < fullSpBasis.qNumbers.n_rows; i++ )
          {
            if (fullSpBasis.qNumbers(i, 1) != l)
              continue;

            if (fullSpBasis.qNumbers(i, 2) != _2j)
              continue;

            v.push_back(i);
          }

          arma::uvec vv(v);
          arma::mat subH = H( vv, vv);
          arma::mat subD ;
          arma::vec subE ;
          arma::eig_sym(subE, subD, subH);
          new_indivE.elem(vv) = subE;
          D(pType)(vv, vv) = subD;
        }
      }
    }
    else
    {
      arma::eig_sym(new_indivE, D(pType), H);
    }

    // Extraction of eigenenergies' sequence
    arma::uvec sorted_ind = arma::sort_index(new_indivE);
    sorted_ind = sorted_ind.head(system->particleNumbers(pType));
    occ(pType).zeros();

    if (system->basis->type == "FullSpBasis")
    {
      FullSpBasis &fullSpBasis = static_cast<FullSpBasis &>(*system->basis);
      arma::vec vecocc = arma::zeros<arma::vec>(sorted_ind.n_elem);
      int accu = 0;

      for (unsigned int i = 0; i < sorted_ind.n_elem; i++)
      {
        int state = sorted_ind(i);
        int _2j = fullSpBasis.qNumbers(state, 2);

        if (accu >= system->particleNumbers(pType))
          break;

        if (accu + _2j + 1 >= system->particleNumbers(pType))
        {
          vecocc(i) = system->particleNumbers(pType) - accu;
          accu += system->particleNumbers(pType) - accu;
          break;
        }

        accu += _2j + 1;
        vecocc(i) = _2j + 1;
      }

      occ(pType)(sorted_ind) = vecocc;
    }
    else
    {
      occ(pType)(sorted_ind) = arma::ones<arma::vec>(sorted_ind.n_rows);
    }

    // New derivation of rho
    R(0, pType) = D(pType) * arma::diagmat(occ(pType)) * D(pType).t();
    indivEnergies.row(pType) = new_indivE.t();
    // Convergence check
    cvg = arma::accu(arma::abs(new_indivE - old_indivE)) / nb_state;
  }
}

void HartreeFock::calcH()
{
}

