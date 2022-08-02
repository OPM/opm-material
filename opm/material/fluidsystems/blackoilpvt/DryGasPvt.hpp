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
 * \copydoc Opm::DryGasPvt
 */
#ifndef OPM_DRY_GAS_PVT_HPP
#define OPM_DRY_GAS_PVT_HPP

#include <opm/material/common/Tabulated1DFunction.hpp>

#include <vector>

namespace Opm {

class EclipseState;
class Schedule;

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 *        without vaporized oil.
 */
template <class Scalar>
class DryGasPvt
{
    using SamplingPoints = std::vector<std::pair<Scalar, Scalar>>;

public:
    using TabulatedOneDFunction = Tabulated1DFunction<Scalar>;

    explicit DryGasPvt() = default;
    DryGasPvt(const std::vector<Scalar>& gasReferenceDensity,
              const std::vector<TabulatedOneDFunction>& inverseGasB,
              const std::vector<TabulatedOneDFunction>& gasMu,
              const std::vector<TabulatedOneDFunction>& inverseGasBMu);

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for dry gas using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVDG keywords.
     */
    void initFromState(const EclipseState& eclState, const Schedule&);
#endif

    void setNumRegions(size_t numRegions);


    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar /*rhoRefOil*/,
                               Scalar rhoRefGas,
                               Scalar /*rhoRefWater*/);

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setMolarMasses(unsigned /*regionIdx*/,
                        Scalar /*MOil*/,
                        Scalar /*MGas*/,
                        Scalar /*MWater*/)
    { }

    /*!
     * \brief Initialize the viscosity of the gas phase.
     *
     * This is a function of \f$(p_g)\f$...
     */
    void setGasViscosity(unsigned regionIdx, const TabulatedOneDFunction& mug);

    /*!
     * \brief Initialize the function for the formation volume factor of dry gas
     *
     * \param samplePoints A container of \f$(p_g, B_g)\f$ values
     */
    void setGasFormationVolumeFactor(unsigned regionIdx, const SamplingPoints& samplePoints);

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd();

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return gasReferenceDensity_.size(); }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned,
                        const Evaluation&,
                        const Evaluation&,
                        const Evaluation&) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& /*Rv*/,
                         const Evaluation& /*Rvw*/) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas at given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& pressure) const;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& /*Rv*/,
                                            const Evaluation& /*Rvw*/) const;

    /*!
     * \brief Returns the formation volume factor [-] of oil saturated gas at given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& /*temperature*/,
                                                     const Evaluation& pressure) const;

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param Rv The surface volume of oil component dissolved in what will yield one cubic meter of gas at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rv*/) const;

    /*!
     * \brief Returns the water vaporization factor \f$R_v\f$ [m^3/m^3] of the water phase.
     */
    template <class Evaluation>
    Evaluation saturatedWaterVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/) const;

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/,
                                              const Evaluation& /*oilSaturation*/,
                                              const Evaluation& /*maxOilSaturation*/) const;

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/) const;

    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& /*temperature*/,
                                    const Evaluation& /*pressure*/,
                                    unsigned /*compIdx*/) const;

    const Scalar gasReferenceDensity(unsigned regionIdx) const
    { return gasReferenceDensity_[regionIdx]; }

    const std::vector<TabulatedOneDFunction>& inverseGasB() const
    { return inverseGasB_; }

    const std::vector<TabulatedOneDFunction>& gasMu() const
    { return gasMu_; }

    const std::vector<TabulatedOneDFunction> inverseGasBMu() const
    { return inverseGasBMu_; }

    bool operator==(const DryGasPvt<Scalar>& data) const;

private:
    std::vector<Scalar> gasReferenceDensity_;
    std::vector<TabulatedOneDFunction> inverseGasB_;
    std::vector<TabulatedOneDFunction> gasMu_;
    std::vector<TabulatedOneDFunction> inverseGasBMu_;
};

} // namespace Opm

#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt_impl.hpp>

#endif
