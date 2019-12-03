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
 * \copydoc Opm::ConstantCompressibilitySaltWaterPvt
 */
#ifndef OPM_CONSTANT_COMPRESSIBILITY_SALTWATER_PVT_HPP
#define OPM_CONSTANT_COMPRESSIBILITY_SALTWATER_PVT_HPP

#include <opm/material/common/Tabulated1DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/Deck/DeckItem.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvtwsaltTable.hpp>
#endif

#include <vector>

namespace Opm {
template <class Scalar, bool enableThermal, bool enableSaltWater>
class WaterPvtMultiplexer;
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 *        without vaporized oil.
 */
template <class Scalar>
class ConstantCompressibilitySaltWaterPvt
{
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef typename Opm::Tabulated1DFunction<Scalar> TabulatedFunction;
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
#if HAVE_ECL_INPUT
    /*!
     * \brief Sets the pressure-dependent water viscosity and density
     *        using a table stemming from the Eclipse PVTWSALT keyword.
     */
    void initFromDeck(const Deck& deck, const EclipseState& eclState)
    {
        const auto& tableManager = eclState.getTableManager();
        size_t numRegions = tableManager.getTabdims().getNumPVTTables();
        const auto& densityKeyword = deck.getKeyword("DENSITY");

        formationVolumeTables_.resize(numRegions);
        compressibilityTables_.resize(numRegions);
        viscosityTables_.resize(numRegions);
        viscosibilityTables_.resize(numRegions);
        referencePressure_.resize(numRegions);

        const auto& pvtwsaltTables = tableManager.getPvtwSaltTables();
        if(!pvtwsaltTables.empty()){
            assert(numRegions == pvtwsaltTables.size());
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                const auto& pvtwsaltTable = pvtwsaltTables[regionIdx];
                const auto& c = pvtwsaltTable.getSaltConcentrationColumn();

                const auto& B = pvtwsaltTable.getFormationVolumeFactorColumn();
                formationVolumeTables_[regionIdx].setXYContainers(c, B);

                const auto& compressibility = pvtwsaltTable.getCompressibilityColumn();
                compressibilityTables_[regionIdx].setXYContainers(c, compressibility);

                const auto& viscositytable = pvtwsaltTable.getViscosityColumn();
                viscosityTables_[regionIdx].setXYContainers(c, viscositytable);

                const auto& viscosibility = pvtwsaltTable.getViscosibilityColumn();
                viscosibilityTables_[regionIdx].setXYContainers(c, viscosibility);
                referencePressure_[regionIdx] = pvtwsaltTable.getReferencePressureValue();
            }
        }
        else {
            throw std::runtime_error("PVTWSALT must be specified in SALTWATER runs\n");
        }


        size_t numPvtwRegions = numRegions;
        setNumRegions(numPvtwRegions);

        for (unsigned regionIdx = 0; regionIdx < numPvtwRegions; ++ regionIdx) {
            auto densityRecord = densityKeyword.getRecord(regionIdx);

            waterReferenceDensity_[regionIdx] =
                densityRecord.getItem("WATER").getSIDouble(0);
        }

        initEnd();
    }
#endif

    void setNumRegions(size_t numRegions)
    {
        waterReferenceDensity_.resize(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            setReferenceDensities(regionIdx, 650.0, 1.0, 1000.0);
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
    Evaluation internalEnergy(unsigned regionIdx OPM_UNUSED,
                        const Evaluation& temperature OPM_UNUSED,
                        const Evaluation& pressure OPM_UNUSED) const
    {
        throw std::runtime_error("Requested the enthalpy of water but the thermal option is not enabled");
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& saltwaterconcentration) const
    {
        // cf. ECLiPSE 2013.2 technical description, p. 114
        Scalar pRef = referencePressure_[regionIdx];
        const Evaluation C = compressibilityTables_[regionIdx].eval(saltwaterconcentration, /*extrapolate=*/true);
        const Evaluation Cv = viscosibilityTables_[regionIdx].eval(saltwaterconcentration, /*extrapolate=*/true);
        const Evaluation BwRef = formationVolumeTables_[regionIdx].eval(saltwaterconcentration, /*extrapolate=*/true);
        const Evaluation Y = (C-Cv)* (pressure - pRef);
        Evaluation MuwRef = viscosityTables_[regionIdx].eval(saltwaterconcentration, /*extrapolate=*/true);

        const Evaluation& bw = inverseFormationVolumeFactor(regionIdx, temperature, pressure, saltwaterconcentration);

        return MuwRef*BwRef*bw/(1 + Y*(1 + Y/2));
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& /*temperature*/,
                                            const Evaluation& pressure,
                                            const Evaluation& saltwaterconcentration) const
    {
        Scalar pRef = referencePressure_[regionIdx];

        const Evaluation BwRef = formationVolumeTables_[regionIdx].eval(saltwaterconcentration, /*extrapolate=*/true);
        const Evaluation C = compressibilityTables_[regionIdx].eval(saltwaterconcentration, /*extrapolate=*/true);
        const Evaluation X = C * (pressure - pRef);

        return (1.0 + X*(1.0 + X/2.0))/BwRef;

    }

private:
    std::vector<TabulatedFunction> formationVolumeTables_;
    std::vector<TabulatedFunction> compressibilityTables_;
    std::vector<TabulatedFunction> viscosityTables_;
    std::vector<TabulatedFunction> viscosibilityTables_;
    std::vector<Scalar> referencePressure_;
    std::vector<Scalar> waterReferenceDensity_;

};

} // namespace Opm

#endif
