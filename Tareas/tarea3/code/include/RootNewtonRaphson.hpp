/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 * @author: Andre Herrera
 * @date  : 28.08.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {

    /**
     * Find the roots of the function funct looking by means of the
     * Newton-Raphson method
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xi initial root guess
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootNewtonRaphson(const std::function<T(T)> &funct, T xi, const T eps) {
        const int ITMAX = 20;                     // set iteration max to 20
        T x = xi;                                 // initial x value
        T xp = x - T(1);                          // value near xi
        T h = eps / T(2);                         // tolerance
        for (int i = 0; i < ITMAX; ++i) {
            T fxm = funct(xp - h);                // centered derivative
            T fxp = funct(xp + h);
            T dfxl = fxp - fxm;
            T dev = dfxl / (2 * h);               // derivative
            x = xp;
            T fx = funct(x);                      // evaluated point x
            T dfx = dev;
            xp = x - (fx / dfx);                  // expected root
            if (std::abs(funct(xp)) < eps || funct(xp) == T(0))
                return xp;
        }
        return std::numeric_limits<T>::quiet_NaN();
    }

}

#endif
