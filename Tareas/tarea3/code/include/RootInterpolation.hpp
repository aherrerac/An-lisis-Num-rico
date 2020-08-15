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
 * @cite  : Numerical Recipes Third Edition
 * The code contains changes to pass boost unit test and to improve functionality..
 *
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_INTERPOLATION_HPP
#define ANPI_ROOT_INTERPOLATION_HPP

namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], by means of the interpolation method.
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xl lower interval limit
     * @param xu upper interval limit
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootInterpolation(const std::function<T(T)> &funct, T xl, T xu, const T eps) {
      if (xl > xu)                                // evaluate  exception for inverted interval
        throw (anpi::Exception("inverted interval"));
      T fl = funct(xl);
      T fh = funct(xu);
      if (fl * fh > T(0))                         // evaluate exception for unenclosed root
        throw (anpi::Exception("unenclosed root"));
      const int MAXIT = 30;                       // set iteration max to 30
      T xh;                                       // set x higher
      if (fl < T(0)) {                            // evaluate if fl its lower than zero
        xh = xu;                                // change the value if not pass
      } else {
        xh = xl;                                // swap the xu and xl values
        xl = xu;
        std::swap(fl, fh);                      // swap initial evaluation points
      }
      T dx = xh - xl;                             // initialize next point
      for (int i = 0; i < MAXIT; ++i) {
        T rtf = xl + dx * fl / (fl - fh);       // calculate the final root throw iterations
        T f = funct(rtf);                       // calculate the function value at rtf point
        if (f < T(0)) {                         // evaluate if the previous value its lower than zero
          xl = rtf;                           // change the lowest value to previous point
          fl = f;                             // change the lowest evaluation function
        } else {
          xh = rtf;                           // change the highest to new root point
          fh = f;                             // change the highest evaluation function
        }
        dx = xh - xl;
        if (std::abs(funct(rtf)) < eps || f == T(0)) // final evaluation
          return rtf;
      }
      // Return NaN if no root was found
      return std::numeric_limits<T>::quiet_NaN();
    }
}

#endif