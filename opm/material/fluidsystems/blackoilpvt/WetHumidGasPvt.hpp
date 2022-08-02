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
 * \copydoc Opm::WetHumidGasPvt
 */
#ifndef OPM_WET_HUMID_GAS_PVT_HPP
#define OPM_WET_HUMID_GAS_PVT_HPP

#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

namespace Opm {

class EclipseState;
class Schedule;
class SimpleTable;

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 *        with vaporized oil and vaporized water.
 */
template <class Scalar>
class WetHumidGasPvt
{
    using SamplingPoints = std::vector<std::pair<Scalar, Scalar>>;

public:
    using TabulatedTwoDFunction = UniformXTabulated2DFunction<Scalar>;
    using TabulatedOneDFunction = Tabulated1DFunction<Scalar>;

    WetHumidGasPvt()
    {
        vapPar1_ = 0.0;
    }

    WetHumidGasPvt(const std::vector<Scalar>& gasReferenceDensity,
                   const std::vector<Scalar>& oilReferenceDensity,
                   const std::vector<Scalar>& waterReferenceDensity,
                   const std::vector<TabulatedTwoDFunction>& inverseGasBRvwSat,
                   const std::vector<TabulatedTwoDFunction>& inverseGasBRvSat,
                   const std::vector<TabulatedOneDFunction>& inverseSaturatedGasB,
                   const std::vector<TabulatedTwoDFunction>& gasMuRvwSat,
                   const std::vector<TabulatedTwoDFunction>& gasMuRvSat,
                   const std::vector<TabulatedTwoDFunction>& inverseGasBMuRvwSat,
                   const std::vector<TabulatedTwoDFunction>& inverseGasBMuRvSat,
                   const std::vector<TabulatedOneDFunction>& inverseSaturatedGasBMu,
                   const std::vector<TabulatedOneDFunction>& saturatedWaterVaporizationFactorTable,
                   const std::vector<TabulatedOneDFunction>& saturatedOilVaporizationFactorTable,
                   const std::vector<TabulatedOneDFunction>& saturationPressure,
                   Scalar vapPar1);

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for wet gas using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVTG keywords.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule);

private:
    void extendPvtgwTable_(unsigned regionIdx,
                          unsigned xIdx,
                          const SimpleTable& curTable,
                          const SimpleTable& masterTable);

    void extendPvtgTable_(unsigned regionIdx,
                          unsigned xIdx,
                          const SimpleTable& curTable,
                          const SimpleTable& masterTable);

public:
#endif // HAVE_ECL_INPUT

    void setNumRegions(size_t numRegions);

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar rhoRefOil,
                               Scalar rhoRefGas,
                               Scalar rhoRefWater);

    /*!
     * \brief Initialize the function for the water vaporization factor \f$R_v\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedGasWaterVaporizationFactor(unsigned regionIdx, const SamplingPoints& samplePoints);

    /*!
     * \brief Initialize the function for the oil vaporization factor \f$R_v\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedGasOilVaporizationFactor(unsigned regionIdx, const SamplingPoints& samplePoints);


    /*!
     * \brief Finish initializing the gas phase PVT properties.
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
                         const Evaluation& /*temperature*/,
                         const Evaluation& pressure,
                         const Evaluation& Rv,
                         const Evaluation& Rvw) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& pressure) const;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    // template <class Evaluation>
    // Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
    //                                         const Evaluation& /*temperature*/,
    //                                         const Evaluation& pressure,
    //                                         const Evaluation& Rw) const
    // { return inverseGasB_[regionIdx].eval(pressure, Rw, /*extrapolate=*/true); }

    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& /*temperature*/,
                                            const Evaluation& pressure,
                                            const Evaluation& Rv,
                                            const Evaluation& Rvw) const;

    /*!
     * \brief Returns the formation volume factor [-] of water saturated gas at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& /*temperature*/,
                                                     const Evaluation& pressure) const;

    /*!
     * \brief Returns the water vaporization factor \f$R_vw\f$ [m^3/m^3] of the water phase.
     */
    template <class Evaluation>
    Evaluation saturatedWaterVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure) const;

    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure) const;

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the gas phase.
     *
     * This variant of the method prevents all the oil to be vaporized even if the gas
     * phase is still not saturated. This is physically quite dubious but it corresponds
     * to how the Eclipse 100 simulator handles this. (cf the VAPPARS keyword.)
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure,
                                              const Evaluation& oilSaturation,
                                              Evaluation maxOilSaturation) const;

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the water component
     *
     * \param Rw The surface volume of water component dissolved in what will yield one
     *           cubic meter of gas at the surface [-]
     */
    //PJPE assume dependence on Rv
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation&,
                                  const Evaluation& Rw) const;

    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& /*temperature*/,
                                    const Evaluation& /*pressure*/,
                                    unsigned /*compIdx*/) const;

    const Scalar gasReferenceDensity(unsigned regionIdx) const
    { return gasReferenceDensity_[regionIdx]; }

     const Scalar oilReferenceDensity(unsigned regionIdx) const
    { return oilReferenceDensity_[regionIdx]; }

    const Scalar waterReferenceDensity(unsigned regionIdx) const
    { return waterReferenceDensity_[regionIdx]; }

    const std::vector<TabulatedTwoDFunction>& inverseGasB() const {
        return inverseGasBRvSat_;
    }

    const std::vector<TabulatedOneDFunction>& inverseSaturatedGasB() const {
        return inverseSaturatedGasB_;
    }

    const std::vector<TabulatedTwoDFunction>& gasMu() const {
        return gasMuRvSat_;
    }

    const std::vector<TabulatedTwoDFunction>& inverseGasBMu() const {
        return inverseGasBMuRvSat_;
    }

    const std::vector<TabulatedOneDFunction>& inverseSaturatedGasBMu() const {
        return inverseSaturatedGasBMu_;
    }

    const std::vector<TabulatedOneDFunction>& saturatedWaterVaporizationFactorTable() const {
        return saturatedWaterVaporizationFactorTable_;
    }
     
    const std::vector<TabulatedOneDFunction>& saturatedOilVaporizationFactorTable() const {
        return saturatedOilVaporizationFactorTable_;
    }


    const std::vector<TabulatedOneDFunction>& saturationPressure() const {
        return saturationPressure_;
    }

    Scalar vapPar1() const {
        return vapPar1_;
    }

    bool operator==(const WetHumidGasPvt<Scalar>& data) const;

private:
    void updateSaturationPressure_(unsigned regionIdx);

    std::vector<Scalar> gasReferenceDensity_;
    std::vector<Scalar> oilReferenceDensity_;
    std::vector<Scalar> waterReferenceDensity_;
    std::vector<TabulatedTwoDFunction> inverseGasBRvwSat_;
    std::vector<TabulatedTwoDFunction> inverseGasBRvSat_;
    std::vector<TabulatedOneDFunction> inverseSaturatedGasB_;
    std::vector<TabulatedTwoDFunction> gasMuRvwSat_;
    std::vector<TabulatedTwoDFunction> gasMuRvSat_;
    std::vector<TabulatedTwoDFunction> inverseGasBMuRvwSat_;
    std::vector<TabulatedTwoDFunction> inverseGasBMuRvSat_;
    std::vector<TabulatedOneDFunction> inverseSaturatedGasBMu_;
    std::vector<TabulatedOneDFunction> saturatedWaterVaporizationFactorTable_;
    std::vector<TabulatedOneDFunction> saturatedOilVaporizationFactorTable_;
    std::vector<TabulatedOneDFunction> saturationPressure_;

    Scalar vapPar1_;
};

} // namespace Opm

#include <opm/material/fluidsystems/blackoilpvt/WetHumidGasPvt_impl.hpp>

#endif
