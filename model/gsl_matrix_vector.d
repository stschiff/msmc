/* Copyright (c) 2012,2013 Genome Research Ltd.
 *
 * Author: Stephan Schiffels <stephan.schiffels@sanger.ac.uk>
 *
 * This file is part of msmc.
 * msmc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
module model.gsl_matrix_vector;
import std.stdio;
import std.exception;
import std.conv;

extern (C) {
  
  struct gsl_block {
    size_t size;
    double* data;
  }
       
  struct gsl_matrix {
    size_t size1;
    size_t size2;
    size_t tda;
    double* data;
    gsl_block* block;
    int owner;
  }
  
  struct gsl_vector {
    size_t size;
    size_t stride;
    double * data;
    gsl_block * block;
    int owner;
  }
  
  struct gsl_vector_view {
    gsl_vector vector;
  }
  
  alias const(gsl_vector_view) gsl_vector_const_view;
  
  enum CBLAS_TRANSPOSE_t {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113}
  
  gsl_matrix* gsl_matrix_alloc(size_t n1, size_t n2);
  void gsl_matrix_free (gsl_matrix* m);
  double gsl_matrix_get (const gsl_matrix* m, size_t i, size_t j);
  void gsl_matrix_set (gsl_matrix* m, size_t i, size_t j, double x);
  void gsl_matrix_set_zero (gsl_matrix* m);
  
  gsl_vector* gsl_vector_alloc(size_t n1);
  void gsl_vector_free (gsl_vector * m);
  void gsl_vector_set (gsl_vector* m, size_t i, double x);
  void gsl_vector_set_all (gsl_vector * v, double x);
  void gsl_vector_set_zero (gsl_vector* m);
  double gsl_vector_get (const gsl_vector* m, size_t i);
  int gsl_vector_memcpy (gsl_vector * dest, const gsl_vector * src);
  int gsl_vector_scale (gsl_vector * a, const double x);
  int gsl_vector_mul (gsl_vector * a, const gsl_vector * b);
  int gsl_vector_add (gsl_vector * a, const gsl_vector * b);
  
  gsl_vector_view gsl_matrix_row(gsl_matrix *m, size_t i);
  gsl_vector_const_view gsl_matrix_const_row(const gsl_matrix *m, size_t i);

  int gsl_blas_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha,
      const gsl_matrix* A, const gsl_matrix* B, double beta, gsl_matrix* C);
  int gsl_blas_dgemv(CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A,
      const gsl_vector * x, double beta, gsl_vector * y);
  

}

void gsl_blas_dgemm_checked(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha,
    const gsl_matrix* A, const gsl_matrix* B, double beta, gsl_matrix* C)
{
  assert((*A).size1 == (*B).size2);
  assert((*A).size2 == (*B).size1);
  assert((*A).size1 == (*C).size1);
  assert((*A).size2 == (*C).size2);
  gsl_blas_dgemm(TransA, TransB, alpha, A, B, beta, C);
}

void gsl_blas_dgemv_checked(CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A,
    const gsl_vector * x, double beta, gsl_vector * y)
{
  assert((*A).size1 == (*y).size);
  assert((*A).size2 == (*x).size);
  gsl_blas_dgemv(TransA, alpha, A, x, beta, y);
}

double gsl_vector_norm(const gsl_vector* v) {
  auto sum = 0.0;
  foreach(i; 0 .. (*v).size)
    sum += gsl_vector_get(v, i);
  return sum;
}

unittest {
  writeln("test gsl_matrix");
  auto mat = gsl_matrix_alloc(10, 10);
  gsl_matrix_set_zero(mat);
  gsl_matrix_set(mat, 1, 2, 3.4);
  gsl_matrix_set(mat, 2, 3, 5.4);
  gsl_matrix_set(mat, 2, 6, 8.4);
  
  gsl_vector_set_zero(&(gsl_matrix_row(mat, 2).vector));
  
  assert(gsl_matrix_get(mat, 2, 3) == 0.0);
  assert(gsl_matrix_get(mat, 1, 2) == 3.4);
  
}

unittest {
  writeln("test gsl_matrix");
  auto mat = gsl_matrix_alloc(10, 10);
  gsl_matrix_set_zero(mat);
  gsl_matrix_set(mat, 1, 2, 3.4);
  gsl_matrix_set(mat, 2, 3, 5.4);
  gsl_matrix_set(mat, 2, 6, 8.4);
  
  gsl_vector_set_zero(&(gsl_matrix_row(mat, 2).vector));
  
  assert(gsl_matrix_get(mat, 2, 3) == 0.0);
  assert(gsl_matrix_get(mat, 1, 2) == 3.4);
  
}
