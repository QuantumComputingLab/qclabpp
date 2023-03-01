#include "qclab/qgates/Identity.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void Identity< T >::apply( Op op , const int nbQubits ,
                             std::vector< T >& vector ,
                             const int offset ) const {
    // do nothing
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void Identity< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                    const int offset ) const {
    // do nothing
  }
#endif

  // apply
  template <typename T>
  void Identity< T >::apply( Side side , Op op , const int nbQubits ,
                             qclab::dense::SquareMatrix< T >& matrix ,
                             const int offset ) const {
    QGate1< T >::apply( side , op , nbQubits , matrix , offset ) ;
  }

  template class Identity< float > ;
  template class Identity< double > ;
  template class Identity< std::complex< float > > ;
  template class Identity< std::complex< double > > ;

} // namespace qclab::qgates

