#ifndef __BOX_HPP__
#define __BOX_HPP__

#include <cstring>
#include <utility>
/**
 * Class for handling boxes
 */

namespace snowgoose {

    template <class FT> struct Box {

        Box(int dim = 0) {
            mDim = dim;
            mA = new FT[dim];
            mB = new FT[dim];
        }

        ~Box() {
            delete[] mA;
            delete[] mB;
        }

        Box(const Box& b) {
            mDim = b.mDim;
            mA = new FT[mDim];
            mB = new FT[mDim];
            std::memcpy(mA, b.mA, sizeof (FT) * mDim);
            std::memcpy(mB, b.mB, sizeof (FT) * mDim);
        }

        Box(Box&& b) {
            mA = b.mA;
            mB = b.mB;
            mDim = b.mDim;

            b.mA = nullptr;
            b.mB = nullptr;
        }

        /**
         * Move assignment
         * @param b
         * @return 
         */
        Box& operator=(Box&& b) {
            delete [] mA;
            delete [] mB;
            mA = b.mA;
            mB = b.mB;
            mDim = b.mDim;

            b.mA = nullptr;
            b.mB = nullptr;
        }

        /**
         * Move assignment
         * @param b
         * @return 
         *
        Box& operator=(Box&& b) {
            std::swap(mA, b.mA);
            std::swap(mB, b.mB);
            std::swap(mDim, b.mDim);
        }*/

        Box& operator=(const Box& b) = delete;

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
