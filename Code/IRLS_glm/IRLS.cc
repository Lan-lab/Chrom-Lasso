/** \file IRLS.cc
*
* `IRLS' is a C++ implementation of the IRLS algorithm for GLM
* Copyright (C) 2013 Xioaquan Wen, Timothee Flutre
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstring>
#include <iostream>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "IRLS.h"
#include "LogLink.h"

using namespace std;

IRLS::IRLS(const char * link_type)
{
  if(strcmp(link_type,"log-link") == 0)
    link = new LogLink();
  link->quasi = false;
  bv = 0;
  VB = 0;
}

// Xv should not contain the intercept
// if offv is empty, offset 0 will be used
void IRLS::load_data(const vector<double> & yv,
		     const vector<vector<double> > &Xv,
		     const vector<double> & offv)
{
  free_data = true;
  
  n = yv.size();
  p = 1 + Xv.size();
  
  y = gsl_vector_calloc(n);
  X = gsl_matrix_calloc(n, p);
  offset = gsl_vector_calloc(n);
  
  for(size_t i = 0; i < n; ++i){
    gsl_vector_set(y, i, yv[i]);
    gsl_matrix_set(X, i, 0, 1.0); // intercept
    for(size_t j = 1; j < p; ++j)
      gsl_matrix_set(X, i, j, Xv[j-1][i]);
  }
  if(! offv.empty())
    for(size_t i = 0; i < n; ++i)
      gsl_vector_set(offset, i, offv[i]);
}

// Xv should contain the intercept
void IRLS::set_data(gsl_vector * yv,
		    gsl_matrix * Xv,
		    gsl_vector * offv)
{
  free_data = false;
  
  n = yv->size;
  p = Xv->size2;
  
  y = yv;
  X = Xv;
  if(offv != NULL)
    offset = offv;
  else{
    offset = gsl_vector_calloc(n);
  }
}

void IRLS::fit_model()
{
  gsl_vector * mv = gsl_vector_calloc(n);
  ////std::cout<<"1"<<endl;
  link->init_mv(y, mv);
  ////std::cout<<"2"<<endl;
  
  gsl_vector * z = gsl_vector_calloc(n);
  //std::cout<<"3"<<endl;
  
  gsl_vector * w = gsl_vector_calloc(n);
  //std::cout<<"4"<<endl;
  
  bv = gsl_vector_alloc(p);
  //std::cout<<"5"<<endl;
  
  gsl_matrix * cov = gsl_matrix_alloc(p, p);
  //std::cout<<"6"<<endl;
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (n, p);
  //std::cout<<"7"<<endl;
  double old_chisq = -1, chisq;
  
  while(true){
    link->compute_z(y, mv, offset, z);
  //std::cout<<"8"<<endl;
    link->compute_weights(mv, w);
  //std::cout<<"9"<<endl;
    
    // weighted least square fitting
    gsl_multifit_wlinear_svd(X, w, z, GSL_DBL_EPSILON, &rank, bv, cov, &chisq, work);
  //std::cout<<"10"<<endl;
    
    if(fabs(chisq - old_chisq) < 1e-6){ // check convergence
      psi = link->compute_dispersion(y, X, bv, offset, mv, rank, link->quasi);
  //std::cout<<"11"<<endl;
      compute_variance(w);
  //std::cout<<"12"<<endl;
      break;
    }
    
    old_chisq = chisq;
    link->compute_mv(bv, X, offset, mv);
  //std::cout<<"13"<<endl;
  }
  
  gsl_vector_free(mv);
  gsl_vector_free(z);
  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(work);
}

// if quasi-likelihood, link->compute_dispersion() must be called before
void IRLS::compute_variance(gsl_vector * w)
{
  if(VB != 0)
    gsl_matrix_free(VB);
  
  VB = gsl_matrix_calloc(p, p);
  gsl_matrix * t2 = gsl_matrix_calloc(p, p);
  
  for (size_t i=0; i<p; i++) {
    for (size_t j=i; j<p; j++) {
      double val=0.0;
      for (size_t k=0; k<n; k++) {
        val+=gsl_matrix_get(X, k, i)*gsl_matrix_get(X, k, j)*gsl_vector_get(w, k);
      }
      gsl_matrix_set (t2, i, j, val);
      gsl_matrix_set (t2, j, i, val);
    }
  }

  //gsl_matrix * W = gsl_matrix_calloc(n, n);
  //for(size_t i = 0; i < n; ++i)
  //  gsl_matrix_set(W, i, i, gsl_vector_get(w, i));
  //
  //gsl_matrix * t1 = gsl_matrix_calloc(p, n);
  //gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, X, W, 0, t1);
  //
  //gsl_matrix * t2 = gsl_matrix_calloc(p, p);
  //gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, t1, X, 0, t2);
  
  // invert t2
  int ss;
  gsl_permutation * pp = gsl_permutation_alloc(p);
  gsl_linalg_LU_decomp(t2, pp, &ss);
  gsl_linalg_LU_invert(t2, pp, VB);
  
  gsl_matrix_scale(VB, psi); // cf. quasi-likelihood
  
  //gsl_matrix_free(W);
  //gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_permutation_free(pp);
}

vector<double> IRLS::get_coef()
{
  vector<double> coev;
  for(size_t i=0; i < p; ++i)
    coev.push_back(gsl_vector_get(bv, i));
  return coev;
}

vector<double> IRLS::get_stderr()
{
  vector<double> sev;
  for(size_t i = 0; i < p; ++i)
    sev.push_back(sqrt(gsl_matrix_get(VB, i, i)));
  return sev;
}

IRLS::~IRLS()
{
  delete link;
  
  if(free_data){
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(offset);
  }
  
  if(bv !=0)
    gsl_vector_free(bv);
  
  if(VB != 0)
    gsl_matrix_free(VB);
}
