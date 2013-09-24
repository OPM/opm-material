// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Felix Bode                                        *
 *   Copyright (C) 2011-2012 by Philipp Nuske                                *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Opm::N2
 */
#ifndef OPM_N2_HH
#define OPM_N2_HH

#include <opm/material/IdealGas.hpp>

#include "Component.hpp"

#include <cmath>

namespace Opm
{

/*!
 * \ingroup Components
 *
 * \brief Properties of pure molecular nitrogen \f$N_2\f$.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class N2 : public Component<Scalar, N2<Scalar> >
{
    typedef Opm::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for nitrogen.
     */
    static const char *name()
    { return "N2"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of molecular nitrogen.
     */
    static Scalar molarMass()
    { return 28.0134e-3;}

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of molecular nitrogen
     */
    static Scalar criticalTemperature()
    { return 126.192; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of molecular nitrogen.
     */
    static Scalar criticalPressure()
    { return 3.39858e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at molecular nitrogen's triple point.
     */
    static Scalar tripleTemperature()
    { return 63.151; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at molecular nitrogen's triple point.
     */
    static Scalar triplePressure()
    { return 12.523e3; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure molecular nitrogen
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * Taken from:
     *
     * R. Span, E.W. Lemmon, et al.: "A Reference Equation of State
     * for the Thermodynamic Properties of Nitrogen for Temperatures
     * from 63.151 to 1000 K and Pressures to 2200 MPa", Journal of
     * Physical and Chemical Refefence Data, Vol. 29, No. 6,
     * pp. 1361-1433
     */
    static Scalar vaporPressure(Scalar T)
    {
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // N2 is solid: We don't take sublimation into
                      // account

        // note: this is the ancillary equation given on page 1368
        Scalar sigma = Scalar(1.0) - T/criticalTemperature();
        Scalar sqrtSigma = std::sqrt(sigma);
        const Scalar N1 = -6.12445284;
        const Scalar N2 = 1.26327220;
        const Scalar N3 = -0.765910082;
        const Scalar N4 = -1.77570564;
        return
            criticalPressure() *
            std::exp(criticalTemperature()/T*
                     (sigma*(N1 +
                             sqrtSigma*N2 +
                             sigma*(sqrtSigma*N3 +
                                    sigma*sigma*sigma*N4))));
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of \f$N_2\f$ gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

     /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure of gaseous \f$N_2\f$ in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure nitrogen gas.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 154, 657, 665
     */
    static const Scalar gasEnthalpy(Scalar T,
                                    Scalar pressure)
    {
        // method of Joback
        const Scalar cpVapA = 31.15;
        const Scalar cpVapB = -0.01357;
        const Scalar cpVapC = 2.680e-5;
        const Scalar cpVapD = -1.168e-8;

        // calculate: \int_0^T c_p dT
        return
            1/molarMass()* // conversion from [J/(mol K)] to [J/(kg K)]

            T*(cpVapA + T*
               (cpVapB/2 + T*
                (cpVapC/3 + T*
                 (cpVapD/4))));

//#warning NIST DATA STUPID INTERPOLATION
//        Scalar T2 = 300.;
//        Scalar T1 = 285.;
//        Scalar h2 = 311200.;
//        Scalar h1 = 295580.;
//        Scalar h = h1+ (h2-h1) / (T2-T1) * (T-T1);
//        return h ;
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure nitrogen gas.
     *
     *        Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     *
     *        Rearranging for internal energy yields: \f$u = h - pv\f$.
     *
     *        Exploiting the Ideal Gas assumption (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     *        The universal gas constant can only be used in the case of molar formulations.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            1/molarMass()* // conversion from [J/(mol K)] to [J/(kg K)]
            IdealGas::R*temperature; // = pressure * spec. volume for an ideal gas
    }

    /*!
     * \brief Specific isobaric heat capacity \f$[J/(kg K)]\f$ of pure
     *        nitrogen gas.
     *
     * This is equivalent to the partial derivative of the specific
     * enthalpy to the temperature.
     */
    static const Scalar gasHeatCapacity(Scalar T,
                                        Scalar pressure)
    {
        // method of Joback
        const Scalar cpVapA = 31.15;
        const Scalar cpVapB = -0.01357;
        const Scalar cpVapC = 2.680e-5;
        const Scalar cpVapD = -1.168e-8;

        return
            1/molarMass()* // conversion from [J/(mol K)] to [J/(kg K)]

            cpVapA + T*
            (cpVapB + T*
             (cpVapC + T*
              (cpVapD)));
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of \f$N_2\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids,
     * 4th edition, McGraw-Hill, 1987, pp 396-397,
     * 5th edition, McGraw-Hill, 2001  pp 9.7-9.8 (omega and V_c taken from p. A.19)
     *
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 90.1; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.037; // accentric factor
        const Scalar M = molarMass() * 1e3; // molar mas [g/mol]
        const Scalar dipole = 0.0; // dipole moment [debye]

        Scalar mu_r4 = 131.3 * dipole / std::sqrt(Vc * Tc);
        mu_r4 *= mu_r4;
        mu_r4 *= mu_r4;

        Scalar Fc = 1 - 0.2756*omega + 0.059035*mu_r4;
        Scalar Tstar = 1.2593 * temperature/Tc;
        Scalar Omega_v =
            1.16145*std::pow(Tstar, -0.14874) +
            0.52487*std::exp(- 0.77320*Tstar) +
            2.16178*std::exp(- 2.43787*Tstar);
        Scalar mu = 40.785*Fc*std::sqrt(M*temperature)/(std::pow(Vc, 2./3)*Omega_v);

        // convertion from micro poise to Pa s
        return mu/1e6 / 10;
    }

    /*!
     * \brief Specific heat conductivity of steam \f$\mathrm{[W/(m K)]}\f$.
     *
     * Isobaric Properties for Nitrogen in: NIST Standard Reference
     * Database Number 69, Eds. P.J. Linstrom and W.G. Mallard
     * evaluated at p=.1 MPa, T=8°C, does not change dramatically with
     * p,T
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasThermalConductivity(Scalar temperature,
                                                  Scalar pressure)
    { return 0.024572; }
};

} // end namepace

#endif
