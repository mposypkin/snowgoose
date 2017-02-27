/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   randpointgen.hpp
 * Author: mposypkin
 *
 * Created on February 22, 2017, 1:32 PM
 */

#ifndef RANDPOINTGEN_HPP
#define RANDPOINTGEN_HPP


#include <random>
#include <ctime>
#include <box/box.hpp>

#include "pointgenerator.hpp"

namespace snowgoose {

    template <class FT> class RandomPointGenerator : public PointGenerator <FT> {
    public:

        /**
         * Constructor
         * @param box input box
         * @param seed seed (0 - not use seed)
         */
        RandomPointGenerator(const snowgoose::Box<FT>& box, int seed = 0) :
        mBox(box),
        mSeed(seed) {
            mGen.seed(time(0));
            reset();
        }

        bool getPoint(FT* x) {
            const int n = mBox.mDim;
            for (int i = 0; i < n; i++) {
                std::uniform_real_distribution<> urd(mBox.mA[i], mBox.mB[i]);
                x[i] = urd(mGen);
            }
            return true;
        }
        bool getPoint(std::vector<FT> &x) {
            const int n = mBox.mDim;
            for (int i = 0; i < n; i++) {
                std::uniform_real_distribution<> urd(mBox.mA[i], mBox.mB[i]);
                x[i] = urd(mGen);
            }
            return true;
        }


        void reset() {
            if (mSeed != 0) {
                mGen.seed(mSeed);
            }
        }


    private:
        const snowgoose::Box<FT>& mBox;
        const int mSeed;
        std::mt19937 mGen;
    };
}

#endif /* RANDPOINTGEN_HPP */

