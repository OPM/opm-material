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
 * \copydoc Opm::ConstantCompressibilityWaterPvt
 */
#ifndef OPM_CONSTANT_COMPRESSIBILITY_WATER_PVT_HPP
#define OPM_CONSTANT_COMPRESSIBILITY_WATER_PVT_HPP

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#endif

#include <vector>

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 *        without vaporized oil.
 */
template <class Scalar>
class ConstantCompressibilityWaterPvt
{
public:
    ConstantCompressibilityWaterPvt() = default;
    ConstantCompressibilityWaterPvt(const std::vector<Scalar>& waterReferenceDensity,
                                    const std::vector<Scalar>& waterReferencePressure,
                                    const std::vector<Scalar>& waterReferenceFormationVolumeFactor,
                                    const std::vector<Scalar>& waterCompressibility,
                                    const std::vector<Scalar>& waterViscosity,
                                    const std::vector<Scalar>& waterViscosibility)
        : waterReferenceDensity_(waterReferenceDensity)
        , waterReferencePressure_(waterReferencePressure)
        , waterReferenceFormationVolumeFactor_(waterReferenceFormationVolumeFactor)
        , waterCompressibility_(waterCompressibility)
        , waterViscosity_(waterViscosity)
        , waterViscosibility_(waterViscosibility)
    { }
#if HAVE_ECL_INPUT
    /*!
     * \brief Sets the pressure-dependent water viscosity and density
     *        using a table stemming from the Eclipse PVTW keyword.
     */
    void initFromState(const EclipseState& eclState, const Schedule&)
    {
        const auto& pvtwTable = eclState.getTableManager().getPvtwTable();
        const auto& densityTable = eclState.getTableManager().getDensityTable();

        assert(pvtwTable.size() == densityTable.size());

        size_t numRegions = pvtwTable.size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            waterReferenceDensity_[regionIdx] = densityTable[regionIdx].water;

            waterReferencePressure_[regionIdx] = pvtwTable[regionIdx].reference_pressure;
            waterReferenceFormationVolumeFactor_[regionIdx] = pvtwTable[regionIdx].volume_factor;
            waterCompressibility_[regionIdx] = pvtwTable[regionIdx].compressibility;
            waterViscosity_[regionIdx] = pvtwTable[regionIdx].viscosity;
            waterViscosibility_[regionIdx] = pvtwTable[regionIdx].viscosibility;
        }

        initEnd();
    }
#endif

    void setNumRegions(size_t numRegions)
    {
        waterReferenceDensity_.resize(numRegions);
        waterReferencePressure_.resize(numRegions);
        waterReferenceFormationVolumeFactor_.resize(numRegions);
        waterCompressibility_.resize(numRegions);
        waterViscosity_.resize(numRegions);
        waterViscosibility_.resize(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            setReferenceDensities(regionIdx, 650.0, 1.0, 1000.0);
            setReferenceFormationVolumeFactor(regionIdx, 1.0);
            setReferencePressure(regionIdx, 1e5);
        }
    }

    /*!
     * \brief Set the water reference density [kg / m^3]
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar /*rhoRefOil*/,
                               Scalar /*rhoRefGas*/,
                               Scalar rhoRefWater)
    { waterReferenceDensity_[regionIdx] = rhoRefWater; }

    /*!
     * \brief Set the water reference pressure [Pa]
     */
    void setReferencePressure(unsigned regionIdx, Scalar p)
    { waterReferencePressure_[regionIdx] = p; }

    /*!
     * \brief Set the viscosity and "viscosibility" of the water phase.
     */
    void setViscosity(unsigned regionIdx, Scalar muw, Scalar waterViscosibility = 0.0)
    {
        waterViscosity_[regionIdx] = muw;
        waterViscosibility_[regionIdx] = waterViscosibility;
    }

    /*!
     * \brief Set the compressibility of the water phase.
     */
    void setCompressibility(unsigned regionIdx, Scalar waterCompressibility)
    { waterCompressibility_[regionIdx] = waterCompressibility; }

    /*!
     * \brief Set the water reference formation volume factor [-]
     */
    void setReferenceFormationVolumeFactor(unsigned regionIdx, Scalar BwRef)
    { waterReferenceFormationVolumeFactor_[regionIdx] = BwRef; }

    /*!
     * \brief Set the water "viscosibility" [1/ (Pa s)]
     */
    void setViscosibility(unsigned regionIdx, Scalar muComp)
    { waterViscosibility_[regionIdx] = muComp; }

    /*!
     * \brief Finish initializing the water phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return waterReferenceDensity_.size(); }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of water given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned,
                        const Evaluation&,
                        const Evaluation&,
                        const Evaluation&,
                        const Evaluation&) const
    {
        throw std::runtime_error("Requested the enthalpy of water but the thermal option is not enabled");
    }


    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& saltconcentration) const
    {
        Scalar BwMuwRef = waterViscosity_[regionIdx]*waterReferenceFormationVolumeFactor_[regionIdx];
        const Evaluation& bw = saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure, saltconcentration);

        Scalar pRef = waterReferencePressure_[regionIdx];
        const Evaluation& Y =
            (waterCompressibility_[regionIdx] - waterViscosibility_[regionIdx])
            * (pressure - pRef);
        return BwMuwRef*bw/(1 + Y*(1 + Y/2));
    }
    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rsw,
                         const Evaluation& saltconcentration) const
    {
        Scalar BwMuwRef = waterViscosity_[regionIdx]*waterReferenceFormationVolumeFactor_[regionIdx];
        const Evaluation& bw = inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rsw, saltconcentration);

        Scalar pRef = waterReferencePressure_[regionIdx];
        const Evaluation& Y =
            (waterCompressibility_[regionIdx] - waterViscosibility_[regionIdx])
            * (pressure - pRef);
        return BwMuwRef*bw/(1 + Y*(1 + Y/2));
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                    const Evaluation& temperature,
                                                    const Evaluation& pressure,
                                                    const Evaluation& saltconcentration) const
    {
      Evaluation Rsw = 0.0;
      return inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rsw, saltconcentration);
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& /*temperature*/,
                                            const Evaluation& pressure,
                                            const Evaluation& /*Rsw*/,
                                            const Evaluation& /*saltconcentration*/) const
    {
        // cf. ECLiPSE 2011 technical description, p. 116
        Scalar pRef = waterReferencePressure_[regionIdx];
        const Evaluation& X = waterCompressibility_[regionIdx]*(pressure - pRef);

        Scalar BwRef = waterReferenceFormationVolumeFactor_[regionIdx];

        // TODO (?): consider the salt concentration of the brine
        return (1.0 + X*(1.0 + X/2.0))/BwRef;
    }

    /*!
     * \brief Returns the saturation pressure of the water phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param Rs The surface volume of gas component dissolved in what will yield one cubic meter of oil at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rs*/,
                                  const Evaluation& /*saltconcentration*/) const
    { return 0.0; /* this is dead water, so there isn't any meaningful saturation pressure! */ }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the water phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned /*regionIdx*/,
                                             const Evaluation& /*temperature*/,
                                             const Evaluation& /*pressure*/,
                                             const Evaluation& /*saltconcentration*/) const
    { return 0.0; /* this is dead water! */ }

    const Scalar waterReferenceDensity(unsigned regionIdx) const
    { return waterReferenceDensity_[regionIdx]; }

    const std::vector<Scalar>& waterReferencePressure() const
    { return waterReferencePressure_; }

    const std::vector<Scalar>& waterReferenceFormationVolumeFactor() const
    { return waterReferenceFormationVolumeFactor_; }

    const std::vector<Scalar>& waterCompressibility() const
    { return waterCompressibility_; }

    const std::vector<Scalar>& waterViscosity() const
    { return waterViscosity_; }

    const std::vector<Scalar>& waterViscosibility() const
    { return waterViscosibility_; }

    bool operator==(const ConstantCompressibilityWaterPvt<Scalar>& data) const
    {
        return this->waterReferenceDensity_ == data.waterReferenceDensity_ &&
               this->waterReferencePressure() == data.waterReferencePressure() &&
               this->waterReferenceFormationVolumeFactor() == data.waterReferenceFormationVolumeFactor() &&
               this->waterCompressibility() == data.waterCompressibility() &&
               this->waterViscosity() == data.waterViscosity() &&
               this->waterViscosibility() == data.waterViscosibility();
    }

private:
    std::vector<Scalar> waterReferenceDensity_;
    std::vector<Scalar> waterReferencePressure_;
    std::vector<Scalar> waterReferenceFormationVolumeFactor_;
    std::vector<Scalar> waterCompressibility_;
    std::vector<Scalar> waterViscosity_;
    std::vector<Scalar> waterViscosibility_;
};

} // namespace Opm

#endif
