/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 * @author: Andre Herrera
 * @date  : 27.08.2018
 * @cite  : Numerical Recipes
 * The code contains changes to pass boost unit test and to improve functionality.
 *
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP

namespace anpi {

    /**
     * Find a root of the function funct looking for it starting at xi
     * by means of the secant method.
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xi initial position
     * @param xii second initial position
     *
     * @return root found, or NaN if no root could be found
     */
    template<typename T>
    T rootSecant(const std::function<T(T)> &funct, T xi, T xii, const T eps) {

        if (xi > xii) throw (anpi::Exception("inverted interval"));
        T fl = funct(xi);                   // set value to point f(xi)
        T fu = funct(xii);                  // set value to point f(xii)
        if (fl * fu > T(0))                 // evaluate if its a unenclosed root
            throw (anpi::Exception("unenclosed root"));
        const int MAXIT = 50;               // set iteration max to 50
        T xl, rts;                          // define variables
        if (std::abs(fl) < std::abs(fu)) {  // evaluate if have to swap values if its necessary this evaluation can be a exception
            rts = xi;                       // swap values
            xl = xii;
            std::swap(fl, fu);
        } else {
            xl = xi;                        // swap values
            rts = xii;
        }
        for (int j = 0; j < MAXIT; j++) {
            T dx = (xl - rts) * fu / (fu - fl); // calculates next value
            xl = rts;                           // change the highest value to the actual root
            fl = fu;                            // change lower with highest
            rts += dx;                          // calculate new approx root
            fu = funct(rts);                    // evaluate the new root in the function
            if (std::abs(funct(rts)) < eps || fu == T(0)) // if the root pass the two evaluations return actual root
                return rts;
        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }

}

#endif

