// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::BrineCo2Pvt
 */
#ifndef OPM_BRINE_CO2_PVT_HPP
#define OPM_BRINE_CO2_PVT_HPP

#include <vector>

namespace Opm {

template<class Scalar> class SimpleHuDuanH2O;
template<class Scalar, class H2O> class Brine;


class EclipseState;
class Schedule;

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the liquid phase
 * for a CO2-Brine system
 */
template <class Scalar>
class BrineCo2Pvt
{
    static constexpr bool extrapolate = true;

public:
    explicit BrineCo2Pvt() = default;
    BrineCo2Pvt(const std::vector<Scalar>& brineReferenceDensity,
                const std::vector<Scalar>& co2ReferenceDensity,
                const std::vector<Scalar>& salinity);

    BrineCo2Pvt(const std::vector<Scalar>& salinity,
                Scalar T_ref = 288.71, //(273.15 + 15.56)
                Scalar P_ref = 101325);

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for Brine-CO2 system using an ECL deck.
     *
     */
    void initFromState(const EclipseState& eclState, const Schedule&);
#endif

    void setNumRegions(size_t numRegions);

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar rhoRefBrine,
                               Scalar rhoRefCO2,
                               Scalar /*rhoRefWater*/);

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    {
    }

    /*!
     * \brief Specify whether the PVT model should consider that the CO2 component can
     *        dissolve in the brine phase
     *
     * By default, dissolved co2 is considered.
     */
    void setEnableDissolvedGas(bool yesno)
    { enableDissolution_ = yesno; }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return brineReferenceDensity_.size(); }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const Evaluation& Rs) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& /*Rs*/) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas at given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned /*regionIdx*/,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const;

    /*!
     * \brief Returns the formation volume factor [-] of brine saturated with CO2 at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const;

    /*!
     * \brief Returns the saturation pressure of the brine phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param Rs
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rs*/) const;

    /*!
     * \brief Returns the gas dissoluiton factor \f$R_s\f$ [m^3/m^3] of the liquid phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& /*oilSaturation*/,
                                             const Evaluation& /*maxOilSaturation*/) const;

    /*!
     * \brief Returns thegas dissoluiton factor  \f$R_s\f$ [m^3/m^3] of the liquid phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const;

    const Scalar oilReferenceDensity(unsigned regionIdx) const
    { return brineReferenceDensity_[regionIdx]; }

    const Scalar gasReferenceDensity(unsigned regionIdx) const
    { return co2ReferenceDensity_[regionIdx]; }

    const Scalar salinity(unsigned regionIdx) const
    { return salinity_[regionIdx]; }

    bool operator==(const BrineCo2Pvt<Scalar>& data) const;

    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& temperature,
                                    const Evaluation& pressure,
                                    unsigned /*compIdx*/) const;

    static Scalar brineMolarMass();
    static Scalar co2MolarMass();

private:
    std::vector<Scalar> brineReferenceDensity_;
    std::vector<Scalar> co2ReferenceDensity_;
    std::vector<Scalar> salinity_;
    bool enableDissolution_ = true;

    template <class LhsEval>
    LhsEval density_(unsigned regionIdx,
                     const LhsEval& temperature,
                     const LhsEval& pressure,
                     const LhsEval& Rs) const;

    template <class LhsEval>
    LhsEval liquidDensity_(const LhsEval& T,
                           const LhsEval& pl,
                           const LhsEval& xlCO2) const;

    template <class LhsEval>
    LhsEval liquidDensityWaterCO2_(const LhsEval& temperature,
                                          const LhsEval& pl,
                                          const LhsEval& xlCO2) const;

    /*!
     * \brief Convert a gas dissolution factor to the the corresponding mass fraction
     *        of the gas component in the oil phase.
     */
    template <class LhsEval>
    LhsEval convertRsToXoG_(const LhsEval& Rs, unsigned regionIdx) const;

    /*!
     * \brief Convert a gas mass fraction in the oil phase the corresponding mole fraction.
     */
    template <class LhsEval>
    LhsEval convertXoGToxoG_(const LhsEval& XoG) const;

    /*!
     * \brief Convert a gas mole fraction in the oil phase the corresponding mass fraction.
     */
    template <class LhsEval>
    LhsEval convertxoGToXoG(const LhsEval& xoG) const;

    /*!
     * \brief Convert the mass fraction of the gas component in the oil phase to the
     *        corresponding gas dissolution factor.
     */
    template <class LhsEval>
    LhsEval convertXoGToRs(const LhsEval& XoG, unsigned regionIdx) const;

    template <class LhsEval>
    LhsEval rsSat_(unsigned regionIdx,
                   const LhsEval& temperature,
                   const LhsEval& pressure) const;

    template <class LhsEval>
    static LhsEval liquidEnthalpyBrineCO2_(const LhsEval& T,
                                           const LhsEval& p,
                                           Scalar S, // salinity
                                           const LhsEval& X_CO2_w);
};

} // namespace Opm

#ifndef OPM_USE_PRIVATE_TEMPLATES
#include <opm/material/fluidsystems/blackoilpvt/BrineCo2Pvt_impl.hpp>
#endif

#endif
