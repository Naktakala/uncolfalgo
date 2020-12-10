#include "uncolfalgo.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/** Computes point-wise values.*/
void UncolFAlgo::Solver::ComputePWLDTransformations()
{
  chi_log.Log() << "Computing PWLD Transformations";
  for (auto& cell : grid->local_cells)
  {
    if (cell.local_id==0) continue;

    auto cell_pwl_view = pwl->MapFeViewL(cell.local_id);
    int Nn = cell_pwl_view->dofs;

    MatDbl A = Ak[cell.local_id];

    MatDbl Ainv(Nn,VecDbl(Nn,0.0));
    for (int i=0; i<Nn; ++i) Ainv[i][i] = 1.0;

    if ( std::fabs(chi_math::Determinant(A)) > 1.0e-8 )
      Ainv = chi_math::Inverse(A);

    VecDbl b(cell_pwl_view->dofs,0.0);

    for (int m=0; m<num_moments; ++m)
    {
      for (int g=0; g<num_groups; ++g)
      {
        //Assemble b
        for (int i=0; i<cell_pwl_view->dofs; ++i)
        {
          int ir = pwl->MapDFEMDOFLocal(&cell,i,&uk_man,m,g);
          b[i] = phi_old_local[ir];
        }

        //Compute x
        VecDbl x = chi_math::MatMul(Ainv,b);

        //Unpack x
        for (int i=0; i<cell_pwl_view->dofs; ++i)
        {
          int ir = pwl->MapDFEMDOFLocal(&cell,i,&uk_man,m,g);
          phi_new_local[ir] += x[i];
        }

      }//for g
    }//for m
  }//for cell

  for (auto& cell : grid->local_cells)
  {

    auto cell_pwl_view = pwl->MapFeViewL(cell.local_id);

    MatDbl A = cell_pwl_view->IntV_shapeI_shapeJ;

    MatDbl Ainv = chi_math::Inverse(A);

    VecDbl b(cell_pwl_view->dofs,0.0);

    for (int m=0; m<num_moments; ++m)
    {
      for (int g=0; g<num_groups; ++g)
      {
        //Assemble b
        for (int i=0; i<cell_pwl_view->dofs; ++i)
        {
          int ir = pwl->MapDFEMDOFLocal(&cell,i,&uk_man,m,g);
          b[i] = phi_old_local_src[ir]*cell_pwl_view->IntV_shapeI[i];
        }

        //Compute x
        VecDbl x = chi_math::MatMul(Ainv,b);

        //Unpack x
        for (int i=0; i<cell_pwl_view->dofs; ++i)
        {
          int ir = pwl->MapDFEMDOFLocal(&cell,i,&uk_man,m,g);
          phi_new_local[ir] += x[i];
        }

      }//for g
    }//for m
  }//for cell

}