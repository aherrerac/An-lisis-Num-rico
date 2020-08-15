/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   05.08.2018
 * @author Andre Herrera
 * @date   27.08.2018
 */
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <PlotPy.hpp>
#include <Matrix.hpp>

#include "Exception.hpp"

/**
 * Load the root finders themselves
 */
#include "RootBisection.hpp"
#include "RootInterpolation.hpp"
#include "RootSecant.hpp"
#include "RootBrent.hpp"
#include "RootNewtonRaphson.hpp"
#include "RootRidder.hpp"

#include "Allocator.hpp"

namespace anpi {
    // namespace encapsulating benchmarking functionality
    namespace bm {

        /// Square of a number
        template<typename T>
        inline T sqr(const T x) { return x * x; }

        /// Cube of a number
        template<typename T>
        inline T cube(const T x) { return x * x * x; }

        /// First testing function for roots |x|=e^(-x)
        template<typename T>
        T t1(const T x) { return std::abs(x) - std::exp(-x); }

        /// Second testing function for roots e^(-x²) = e^(-(x-3)²/3 )
        template<typename T>
        T t2(const T x) { return std::exp(-x * x) - std::exp(-sqr(x - T(3)) / T(3)); }

        /// Third testing function for roots x² = atan(x)
        template<typename T>
        T t3(const T x) { return x * x - std::atan(x); }

        /// Fourth testing function for roots x² = atan(x)
        template<typename T>
        T t4(const T x) {
            const T x0 = x - T(2);
            return cube(x0) + T(0.01) * x0;
        }

        // Global variable to create new instances of Plot2d class
        int plotNumber = 0;


        /**
         * Wrapper class to count function calls
         *
         * This wrapper fulfills the requirements to act as a
         * std::function<T(T)>, and it simply counts the number
         * of calls made to the operator(), before calling
         * the functor provided at construction time.
         */
        template<typename T>
        class CallCounter {
        protected:
            /// Maximum allowed size for the square matrices
            mutable size_t _counter;

            std::function<T(T)> _f;
        public:
            /// Construct
            CallCounter(std::function<T(T)> f) : _counter(0u), _f(f) {}

            /// Access the counter
            inline size_t counter() const { return _counter; }

            /// Reset the counter
            inline void reset() { _counter = 0u; }

            /// Call the function
            T operator()(const T x) {
                ++_counter;
                return _f(x);
            }

            /// Call the function
            T operator()(const T x) const {
                ++_counter;
                return _f(x);
            }
        };

        /**
         * Initialize Plot2d class,sets the plots for the vectors
         * @param esp Vector with the tolerance values
         * @param f1  Vector with the callCounter values for function 1
         * @param f2  Vector with the callCounter values for function 2
         * @param f3  Vector with the callCounter values for function 3
         * @param f4  Vector with the callCounter values for function 4
         *
         * */
        template<typename T>
        void plot(const std::vector<T> &esp, const std::vector<T> &f1,
                  const std::vector<T> &f2,
                  const std::vector<T> &f3,
                  const std::vector<T> &f4,
                  const std::string &method) {
            static anpi::Plot2d<T> plot2d;    // create a static anpi::Plot2d<T> object to show graphs
            plot2d.initialize(plotNumber);    // initialize plot2d object
            plotNumber++;                     // increment the figure id
            plot2d.setTitle(method + "" + "<" + typeid(T).name() +
                            ">"); // set title to previous graph and set the type of typename
            plot2d.setXLabel("Tolerance");    // set the label to x
            plot2d.setYLabel("Function Calls"); // set the label to y
            plot2d.plot(esp, f1, "Function 1", "r"); // create a plot for function 1
            plot2d.plot(esp, f2, "Function 2", "g"); // create a plot for function 2
            plot2d.plot(esp, f3, "Function 3", "b"); // create a plot for function 3
            plot2d.plot(esp, f4, "Function 4", "m"); // create a plot for function 4

        }

        /**
         *  Show the actual plot, needs a template<typename>T to instantiate
         *  ::anpi::bm::show()<T>;
         * */
        template<typename T>
        void show() {
            static anpi::Plot2d<T> plot2d; // get the static object created in plot
            plot2d.show();                 // show the graph with the given plots
        }

        /**
             * Test the given _closed_ root finder
             *
             * The solver must be itself a std::function expecting another
             * std::function (the one whose roots are being looked for), the
             * two limits of the interval enclosing the root, and the
             * tolerance.
             *
             * The tolerances will start from "start", then progressing with
             *   eps = eps*factor
             * until the end value is reached.
             */
        template<typename T>
        void rootBench(const std::function<T(const std::function<T(T)> &,
                                             T,
                                             T,
                                             const T)> &solver,
                       const T start,
                       const T end,
                       const T factor,
                       const std::string &method) {

            if ((factor >= static_cast<T>(1)) &&
                (factor < static_cast<T>(0))) {
                throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
            }

            // Alias of the function type, for which the roots are being looked for.
            typedef std::function<T(T)> f_type;

            //Initialize Error Vector
            std::vector<T> epsVector;
            //Initialize CallCounter Vector per function
            std::vector<T> f1Vector;
            std::vector<T> f2Vector;
            std::vector<T> f3Vector;
            std::vector<T> f4Vector;
            // Try a series of tolerances
            for (T eps = start; eps > end; eps *= factor) {
                epsVector.push_back(eps);
                std::cout << "eps=" << eps << "; ";

                // Create an std::function instance, which wraps the function
                // t1 with the function counter
                f_type c1(CallCounter<T>(t1<T>));
                solver(c1, T(0), T(2), eps);
                //Intercept and save CallCounter
                int counter = c1.template target<CallCounter<T> >()->counter();
                std::cout << counter << "; ";
                //Fill vector for function t1
                f1Vector.push_back(counter);

                // now the same with function t2
                f_type c2(CallCounter<T>(t2<T>));
                solver(c2, T(0), T(2), eps);
                //Intercept and save CallCounter
                counter = c2.template target<CallCounter<T> >()->counter();
                std::cout << counter << "; ";
                //Fill vector for function t2
                f2Vector.push_back(counter);

                // now the same with function t3
                f_type c3(CallCounter<T>(t3<T>));
                solver(c3, T(0), T(0.5), eps);
                //Intercept and save CallCounter
                counter = c3.template target<CallCounter<T> >()->counter();
                std::cout << counter << "; ";
                //Fill vector for function t3
                f3Vector.push_back(counter);

                // now the same with function t4
                f_type c4(CallCounter<T>(t4<T>));
                solver(c4, T(1), T(3), eps);
                //Intercept and save CallCounter
                counter = c4.template target<CallCounter<T> >()->counter();
                std::cout << counter << "; " << std::endl;
                //Fill vector for function t4
                f4Vector.push_back(counter);
            }
            // Set the values to create a new plot
            plot(epsVector, f1Vector, f2Vector, f3Vector, f4Vector, method);
        }

        /**
         * Test the given _open_ root finder
         *
         * The solver must be itself a std::function expecting another
         * std::function (the one whose roots are being looked for), the
         * starting root guess, and the tolerance.
         */
        template<typename T>
        void rootBench(const std::function<T(const std::function<T(T)> &,
                                             T,
                                             const T)> &solver,
                       const T start,
                       const T end,
                       const T factor,
                       const std::string &method) {

            if ((factor >= static_cast<T>(1)) &&
                (factor < static_cast<T>(0))) {
                throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
            }

            // Alias of the function type, for which the roots are being looked for.
            typedef std::function<T(T)> f_type;

            //Initialize Error Vector
            std::vector<T> epsVector;
            //Initialize CallCounter Vector per function
            std::vector<T> f1Vector;
            std::vector<T> f2Vector;
            std::vector<T> f3Vector;
            std::vector<T> f4Vector;
            // Try a series of tolerances
            for (T eps = start; eps > end; eps *= factor) {
                std::cout << "eps=" << eps << "; ";
                epsVector.push_back(eps);
                // Create an std::function instance, which wraps the function
                // t1 with the function counter
                f_type c1(CallCounter<T>(t1<T>));
                solver(c1, T(0), eps);
                //Intercept and save CallCounter
                int counter = c1.template target<CallCounter<T> >()->counter();
                std::cout << counter << "; ";
                //Fill vector for function t1
                f1Vector.push_back(counter);

                // now the same with function t2
                f_type c2(CallCounter<T>(t2<T>));
                solver(c2, T(2), eps);
                //Intercept and save CallCounter
                counter = c2.template target<CallCounter<T> >()->counter();
                std::cout << counter << "; ";
                //Fill vector for function t2
                f2Vector.push_back(counter);

                // now the same with function t3
                f_type c3(CallCounter<T>(t3<T>));
                solver(c3, T(0), eps);
                //Intercept and save CallCounter
                counter = c3.template target<CallCounter<T> >()->counter();
                std::cout << counter << "; ";
                //Fill vector for function t3
                f3Vector.push_back(counter);

                // now the same with function t4
                f_type c4(CallCounter<T>(t4<T>));
                solver(c4, T(1), eps);
                //Intercept and save CallCounter
                counter = c4.template target<CallCounter<T> >()->counter();
                std::cout << counter << "; " << std::endl;
                //Fill vector for function t4
                f4Vector.push_back(counter);
            }
            // Set the values to create a new plot
            plot(epsVector, f1Vector, f2Vector, f3Vector, f4Vector, method);
        }

        /**
         * Benchmark all solvers using a range of tolerances geometrically changing
         * multiplying from the start point until the end with the given factor
         */
        template<typename T>
        void allSolvers(const T start, const T end, const T factor) {

            std::cout << "Bisection" << std::endl;
            anpi::bm::rootBench<T>(anpi::rootBisection<T>, start, end, factor, "Bisection");

            std::cout << "Interpolation" << std::endl;
            anpi::bm::rootBench<T>(anpi::rootInterpolation<T>, start, end, factor, "Interpolation");

            std::cout << "Secant" << std::endl;
            anpi::bm::rootBench<T>(anpi::rootSecant<T>, start, end, factor, "Secant");

            std::cout << "NewtonRaphson" << std::endl;
            anpi::bm::rootBench<T>(anpi::rootNewtonRaphson<T>, start, end, factor, "NewtonRaphson");

            std::cout << "Brent" << std::endl;
            anpi::bm::rootBench<T>(anpi::rootBrent<T>, start, end, factor, "Brent");

            std::cout << "Ridder" << std::endl;
            anpi::bm::rootBench<T>(anpi::rootRidder<T>, start, end, factor, "Ridder");
        }
    } // bm
}  // anpi

BOOST_AUTO_TEST_SUITE(RootFinders)

/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE(RootFinders) {

        // Benchmark the solvers using float
        std::cout << "<float>" << std::endl;
        anpi::bm::allSolvers<float>(0.1f, 1.e-7f, 0.125f);
        // Show the plot created in allSolvers function with float precision
        ::anpi::bm::show<float>();

        // Wait until the user closes all the previous plots

        // Benchmark the solvers using double
        std::cout << "<double>" << std::endl;
        // Show the plot created in allSolvers function with double precision
        anpi::bm::allSolvers<double>(0.1f, 1.e-15f, 0.125f);
        ::anpi::bm::show<double>();
    }

BOOST_AUTO_TEST_SUITE_END()
