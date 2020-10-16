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
#include <ctime>

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
  p = Xv.size();
  //p = 1 + Xv.size();
  
  y = gsl_vector_calloc(n);
  X = gsl_matrix_calloc(n, p);
  offset = gsl_vector_calloc(n);
  
  for(size_t i = 0; i < n; ++i){
    gsl_vector_set(y, i, yv[i]);
    //gsl_matrix_set(X, i, 0, 1.0); // intercept
    //for(size_t j = 1; j < p; ++j)
      //gsl_matrix_set(X, i, j, Xv[j-1][i]);
    for(size_t j = 0; j < p; ++j)
      gsl_matrix_set(X, i, j, Xv[j][i]);
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
  clock_t time_prev, time_now;
  double executeTime = 0.0;
  gsl_vector * mv = gsl_vector_calloc(n);
  link->init_mv(y, mv);

  gsl_vector * z = gsl_vector_calloc(n);
  //sleep(30);
  
  gsl_vector * w = gsl_vector_calloc(n);
  gsl_vector * w_prev = gsl_vector_calloc(n);
  
  bv = gsl_vector_alloc(p);
  gsl_vector* bv_prev = gsl_vector_calloc(p);
  //gsl_matrix* VB_prev = gsl_matrix_calloc(p, p);
  
  gsl_matrix * cov = gsl_matrix_alloc(p, p);
  //time_prev = clock();
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (n, p);
  //time_now = clock();
  //executeTime = double(time_now-time_prev)/CLOCKS_PER_SEC; 
  //printf("workspace %f\n", executeTime);
  double old_chisq = -1, chisq;
  int itNum = 0;
 
  while(true){
    link->compute_z(y, mv, offset, z);
    link->compute_weights(mv, w);
    
    // weighted least square fitting
    //time_prev = clock();
    gsl_multifit_wlinear_svd(X, w, z, GSL_DBL_EPSILON, &rank, bv, cov, &chisq, work);
    //time_now = clock();
    //executeTime = double(time_now-time_prev)/CLOCKS_PER_SEC; 
    //printf("svd %f\n", executeTime);

    ++itNum;
    //printf("%d\t%f\t%f\t%f\n", itNum, chisq, old_chisq, chisq - old_chisq);
    if(isnan(chisq))
    {
      //printf("chisq is nan\n");
      for (int i=0; i<n; ++i)
      {
        gsl_vector_set(w, i, gsl_vector_get(w_prev, i));
      }
      //for (int i=0; i<p; ++i)
      //{
      //  for (int j=0; j<p; ++j)
      //  {
      //    gsl_matrix_set(VB, i, j, gsl_matrix_get(VB_prev, i, j));
      //  }
      //}
      //printf("set back bv\n");
      for (int i=0; i<p; ++i)
      {
        gsl_vector_set(bv, i, gsl_vector_get(bv_prev, i));
      }
      //printf("set back done\n");

      //psi = link->compute_dispersion(y, X, bv, offset, mv, rank, link->quasi);
      psi = 1.0;
      //time_prev = clock();
      compute_variance(w);
      //time_now = clock();
      //executeTime = double(time_now-time_prev)/CLOCKS_PER_SEC; 
      //printf("compute variance %f\n", executeTime);
      break;
    }
    else if(fabs(chisq - old_chisq) < EPSILON){ // check convergence
      //psi = link->compute_dispersion(y, X, bv, offset, mv, rank, link->quasi);
      psi=1.0;
      //time_prev = clock();
      compute_variance(w);
      //time_now = clock();
      //executeTime = double(time_now-time_prev)/CLOCKS_PER_SEC; 
      //printf("compute variance %f\n", executeTime);
      break;
    }
    
    old_chisq = chisq;
    //printf("set w\n");
    for (int i=0; i<n; ++i)
    {
      gsl_vector_set(w_prev, i, gsl_vector_get(w, i));
    }
    //for (int i=0; i<p; ++i)
    //{
    //  for (int j=0; j<p; ++j)
    //  {
    //    gsl_matrix_set(VB_prev, i, j, gsl_matrix_get(VB, i, j));
    //  }
    //}
    //printf("set bv\n");
    for (int i=0; i<p; ++i)
    {
      gsl_vector_set(bv_prev, i, gsl_vector_get(bv, i));
    }
    //printf("set values done\n");
    //time_prev = clock();
    link->compute_mv(bv, X, offset, mv);
    //time_now = clock();
    //executeTime = double(time_now-time_prev)/CLOCKS_PER_SEC; 
    //printf("compute mv %f\n", executeTime);
  }
  
  gsl_vector_free(bv_prev);
  gsl_vector_free(w_prev);
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
  
  // invert t2
  int ss;
  gsl_permutation * pp = gsl_permutation_alloc(p);
  gsl_linalg_LU_decomp(t2, pp, &ss);
  gsl_linalg_LU_invert(t2, pp, VB);
  
  gsl_matrix_scale(VB, psi); // cf. quasi-likelihood
  
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
