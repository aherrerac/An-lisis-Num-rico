/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 04.08.2018
 * @author: Andre Herrera
 * @date  :26.08.2018
 * @cite  : Numerical Recipes Third Edition
 * The code contains changes to pass boost unit test and to improve functionality.
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_RIDDER_HPP
#define ANPI_ROOT_RIDDER_HPP

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
    T rootRidder(const std::function<T(T)> &funct, T xi, T xii, const T eps) {
        if (xi > xii)                                // evaluate  exception for inverted interval
            throw (anpi::Exception("inverted interval"));
        T fl = funct(xi);                           // set value to point f(xl) lower
        T fu = funct(xii);                          // set value to point f(fu) upper
        if ((fl*fu < T(0))) {                       // evaluate exception for unenclosed root
            const int MAXIT = 60;                   // set iteration max to 60
            T xl = xi;                              // set lower to xi
            T xh = xii;                             // set highest to xi
            T ans = -9.99e99;                       // value tu simplify logic
            for (int i = 0; i < MAXIT; ++i) {
                T xm = T(0.5) * (xl + xh);          // calculate the next point("root") value
                T fm = funct(xm);                   // calculate the value of the function in xm
                T s = std::sqrt(fm * fm - fl * fu); // calculate sqrt for s
                if (s == T(0)) return ans;          // evaluate if the sqrt its zero
                T xnew = xm + (xm - xl) * ((fl >= fu ? T(1.0) : T(-1.0)) * fm / s); // update the xnew formula
                if (std::abs(funct(ans)) <= eps) return ans;   // return the value if the evaluation point is lower than tolerance
                ans = xnew;                         // change the actual result
                T fnew = funct(ans);                // calculate the evaluation of f(ans)
                if (fnew == T(0))return ans;        // return if its a root
                if (fnew > T(0) ? std::abs(fm) : -std::abs(fm) != fm) { //sign fm and fnew
                    xl = xm;                        // changing values to find ans on next iteration
                    fl = fm;
                    xh = ans;
                    fu = fnew;
                } else if (fnew > T(0) ? std::abs(fl) : -std::abs(fl) != fl) { //sign fl and fnew
                    xh = ans;                       // changing values to find ans on next iteration
                    fu = fnew;
                } else if (fnew > T(0) ? std::abs(fu) : -std::abs(fu) != fu) { //sign fu and fnew
                    xl = ans;                       // changing values to find ans on next iteration
                    fl = fnew;
                } else throw (anpi::Exception("Error")); // return error
            }
            // Return NaN if no root was found
            return std::numeric_limits<T>::quiet_NaN();
        } else {
            if (fl == T(0)) return xi;              //evaluate if xi its a root
            if (fu == T(0)) return xii;             //evaluate if xii its a root
            throw (anpi::Exception("unenclosed root"));
        }
    }
}

#endif

