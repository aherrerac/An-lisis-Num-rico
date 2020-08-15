/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 * @author: Andre Herrera
 * @update: 27.08.2018
 * @cite  : Numerical Recipes Third Edition
 * The code contains changes to pass boost unit test and to improve functionality.
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], using the bisection method.
     *
     * @param funct a std::function of the form "T funct(T x)"
     * @param xl lower interval limit
     * @param xu upper interval limit
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootBisection(const std::function<T(T)> &funct, T xl, T xu, const T eps) {
      if (xl > xu)
        throw (anpi::Exception("inverted interval"));//evaluate  exception for inverted interval
      T fl = funct(xl);                       // set value to point f(xl)
      T fu = funct(xu);                       // set value to point f(xu)
      if (fl * fu > T(0))                     // evaluate exception for unenclosed root
        throw (anpi::Exception("unenclosed root"));
      const int MAXIT = 50;                   // set iteration max to 50
      T xr;                                   // initialize root final variable
      for (int i = 0; i < MAXIT; ++i) {
        xr = (xl + xu) / T(2);              // centered root
        T fr = funct(xr);                   // shadow for function
        T cond = fl * fr;                   // negative if the points are opposites

        if (cond < T(0)) {                  // evaluate previous condition
          xu = xr;                        // negative -> continue with left side
        } else if (cond > T(0)) {
          xl = xr;                        // positive -> continue with right side
          fl = fr;                        // update fl point
        } else {
          xr = (std::abs(fl) < std::numeric_limits<T>::epsilon())
               ? xl : xr;                 // some border its zero
        }
        if (std::abs(funct(xr)) < eps)      // Evaluate if the current xr its lower than the tolerance
          return xr;
      }
      // Return NaN if no root was found
      return std::numeric_limits<T>::quiet_NaN();
    }

}

#endif

