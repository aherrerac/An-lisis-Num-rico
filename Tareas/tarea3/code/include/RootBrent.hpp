/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 * @author: Andre Herrera
 * @date:   26.08.2018
 * @cite https://en.wikipedia.org/wiki/Brent%27s_method (pseudo code)
 * The code contains changes to pass boost unit test
 */

#include <cmath>
#include <limits>
#include <functional>


#include "Exception.hpp"

#ifndef ANPI_ROOT_BRENT_HPP
#define ANPI_ROOT_BRENT_HPP

namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], using the Brent's method.
     *
     * @param funct a std::function of the form "T funct(T x)"
     * @param xl lower interval limit
     * @param xu upper interval limit
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootBrent(const std::function<T(T)> &funct, T xl, T xu, const T eps) {
      const int MAXIT = 100;
      T fl = funct(xl);                        // set value to point f(xl)
      T fu = funct(xu);                        // set value to point f(xu)
      T fs = 0;                                // initialize
      if (fl * fu > T(0)) {                    // evaluate  exception for unenclosed root
        throw (anpi::Exception("unenclosed root"));
      }
      T c = xl;                                // c equal lower x point
      T fc = fl;                               // set fc equal to the lower f(x) point
      bool flag = true;                        // flag used to evaluate
      T s = T(0);                              // Our Root that will be returned
      T d = T(0);                              // Only used if mflag is unset (mflag == false)

      for (int i = 0; i< MAXIT; ++i) {
        // stop if root less than tolerance eps
        if (std::abs(funct(s)) < eps) {
          return s;                        // return result
        }

        if (fl != fc && fu != fc) {
          // inverse quadratic interopolation method
          s = (xl * fu * fc / ((fl - fu) * (fl - fc)))
              + (xu * fl * fc / ((fu - fl) * (fu - fc)))
              + (c * fl * fu / ((fc - fl) * (fc - fu)));
        } else {
          // secant method
          s = xu - fu * (xu - xl) / (fu - fl);
        }

        // run converging quadratic, secant methods or bisection
        if (((s < (T(3) * xl + xu) * T(0.25)) || (s > xu)) ||
            (flag && (std::abs(s - xu) >= (std::abs(xu - c) * T(0.5)))) ||
            (!flag && (std::abs(s - xu) >= (std::abs(c - d) * T(0.5)))) ||
            (flag && (std::abs(xu - c) < eps)) ||
            (!flag && (std::abs(c - d) < eps))) {
          // bisection method
          s = (xl + xu) * T(0.5);

          flag = true;
        } else {
          flag = false;
        }

        fs = funct(s);                          // calculate fs
        d = c;                              // use d
        c = xu;                             // set c equal to upper point
        fc = fu;                            // set f(c) = f(b)

        if (fl * fs < T(0))                 // evaluate if fa and fs have opposite signs
        {
          xu = s;
          fu = fs;                        // set f(b) = f(s)
        } else {
          xl = s;
          fl = fs;                        // set f(a) = f(s)
        }

        if (std::abs(fl) < std::abs(fu))        // if magnitude of fa is less than magnitude of fb
        {
          std::swap(xl, xu);              // swap a and b
          std::swap(fl, fu);              // swap f(a) and f(b)
        }

      }                                       // end for

      // Return NaN if no root was found
      return std::numeric_limits<T>::quiet_NaN();
    }
}

#endif
