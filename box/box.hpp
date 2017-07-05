#ifndef __BOX_HPP__
#define __BOX_HPP__

/**
  * Class for handling boxes
  */

#include <string.h>

namespace snowgoose {
    
template <class FT> struct Box {
  Box(int dim = 0)
  {
    mDim = dim;
    mA = new FT[dim];
    mB = new FT[dim];
  }

  ~Box(){
      delete[] mA;
      delete[] mB;
  }

  Box(const Box& b){
      mDim = b.mDim;
      mA = new FT[mDim];
      mB = new FT[mDim];
      memcpy(mA, b.mA, sizeof(FT) * mDim);
      memcpy(mB, b.mB, sizeof(FT) * mDim);
  }
  
  Box(Box&& b){
      mA = b.mA;
      mB = b.mB;
      mDim = b.mDim;   
      b.mA = nullptr;
      b.mB = nullptr;
      b.mDim = 0;
  }

 Box& operator=(Box&& b){return *this;};

 Box& operator=(const Box& b) {return *this;}

  /**
    * Box number of dimensions
    */
  int mDim;

 /**
   * "Left" box bounds
   */
  FT* mA;

  /**
   * "Right" box bounds
   */
  FT* mB;
};

}

#endif
