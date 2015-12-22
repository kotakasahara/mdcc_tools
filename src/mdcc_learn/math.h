#ifndef __MATH_H__
#define __MATH_H__

#include <iostream>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
namespace ublas = boost::numeric::ublas;

#include "define.h"

//using namespace std;
template <class M> double determinant(const M& m);
template <class M> double trace(const M& m);
template <class M, class MI> void invert(const M& m, MI& mi);
const ublas::matrix<double> vectorToMatrixTranspose(const ublas::vector<double>& v);
const ublas::matrix<double> vectorToMatrix(const ublas::vector<double>& v);
const ublas::matrix<double> v2m(const ublas::vector<double>& v);
const ublas::matrix<double> v2mT(const ublas::vector<double>& v);
double calGaussian(const ublas::vector<double>& gX,
		   const ublas::vector<double>& gMu,
		   const ublas::matrix<double>& gSigma);
double calGaussianInv(const ublas::vector<double>& gX,
		      const ublas::vector<double>& gMu,
		      const ublas::matrix<double>& gSigmaInv,
		      const double det);
double calMaharanobisDistInv(const ublas::vector<double>& gX,
			    const ublas::vector<double>& gMu,
			    const ublas::matrix<double>& gSigmaInv);
double digamma(double z);
//  double psi(double x);
const ublas::vector<double> zero_vector(int size);
const ublas::matrix<double> zero_matrix(int size1,int size2);
const ublas::matrix<double> id_matrix(int size);

template <class M> double determinant(const M& m){
  BOOST_UBLAS_CHECK(m.size1()==m.size2(), ublas::external_logic());
  ublas::matrix<double> lu(m);
  ublas::permutation_matrix<> pm(m.size1());
  ublas::lu_factorize(lu,pm);
  double det=1;
  typedef ublas::permutation_matrix<>::size_type size_type;
  for(size_type i=0; i<pm.size(); ++i){
    det *= (i==pm(i)) ? +lu(i,i) : -lu(i,i);
  }
  return det;
}


template <class M> double trace(const M& m){
  BOOST_UBLAS_CHECK(m.size1()==m.size2(), ublas::external_logic());
  double tr=0;
  ublas::matrix<double> mtx(m);
  for(int i=0; i<(int)m.size1(); i++){
    tr += mtx(i,i);
  }
  return tr;
}
template <class M, class MI> void invert(const M& m, MI& mi){
  //    std::cout << "invert 1"<<std::endl;
  BOOST_UBLAS_CHECK(m.size1()==m.size2(),ublas::external_logic());
  //    std::cout << "invert 2"<<std::endl;
  ublas::matrix<double> lhs(m);
  ublas::matrix<double> rhs(ublas::identity_matrix<double>(m.size1()));
  ublas::permutation_matrix<> pm(m.size1());
  //    std::cout << "invert 3"<<std::endl;
  ublas::lu_factorize(lhs,pm);
  //    std::cout << "invert 4"<<std::endl;
  ublas::lu_substitute(lhs,pm,rhs);
  //    std::cout << "invert 5"<<std::endl;
  BOOST_UBLAS_CHECK(rhs.size1()==m.size1(), ublas::internal_logic());
  BOOST_UBLAS_CHECK(rhs.size2()==m.size2(), ublas::internal_logic());
  //    std::cout << "invert 6"<<std::endl;
  /*
    #if BOOST_UBLAS_TYPE_CHECK
  BOOST_UBLAS_CHECK
    (ublas::detail::expression_type_check
     (ublas::prod(m,rhs),
      ublas::identity_matrix<typename M::value_type>(m.size1())
      ),
     ublas::internal_logic()
     );
#endif
  */
  //  std::cout << "invert 7"<<std::endl;
  mi.resize(rhs.size1(),rhs.size2(),false);
  mi.assign(rhs);
}


#endif
