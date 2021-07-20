//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_QAdjustable_hpp
#define qclab_QAdjustable_hpp

namespace qclab {

  /**
   * \class QAdjustable
   * \brief Class for representing an adjustable quantum object.
   *
   * This class allows to adjust a quantum object from fixed to variable, and
   * vice versa.
   */
  class QAdjustable
  {

    public:
      /// Constructs a quantum adjustable object with the given flag `fixed`.
      QAdjustable( const bool fixed = false )
      : fixed_( fixed )
      { } // QAdjustable

      /// Checks if this quantum adjustable object is fixed.
      inline bool fixed() const { return fixed_ ; }

      /// Checks if this quantum adjustable object is variable.
      inline bool variable() const { return !fixed_ ; }

      /// Makes this quantum adjustable object fixed.
      inline void makeFixed() { fixed_ = true ; }

      /// Makes this quantum adjustable object variable.
      inline void makeVariable() { fixed_ = false ; }

    private:
      bool fixed_ ;

  } ; // clase QAdjustable

} // namespace qclab

#endif

