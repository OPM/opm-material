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
 * \copydoc Opm::ConstantCompressibilityBrinePvt
 */
#ifndef OPM_CONSTANT_COMPRESSIBILITY_BRINE_PVT_HPP
#define OPM_CONSTANT_COMPRESSIBILITY_BRINE_PVT_HPP

#include <opm/material/common/Tabulated1DFunction.hpp>

#include <vector>

namespace Opm {
template <class Scalar, bool enableThermal, bool enableBrine>
class WaterPvtMultiplexer;

class EclipseState;
class Schedule;

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 *        without vaporized oil.
 */
template <class Scalar>
class ConstantCompressibilityBrinePvt
{
public:
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    ConstantCompressibilityBrinePvt() = default;
    ConstantCompressibilityBrinePvt(const std::vector<Scalar>& waterReferenceDensity,
                                    const std::vector<Scalar>& referencePressure,
                                    const std::vector<TabulatedFunction> formationVolumeTables,
                                    const std::vector<TabulatedFunction> compressibilityTables,
                                    const std::vector<TabulatedFunction> viscosityTables,
                                    const std::vector<TabulatedFunction> viscosibilityTables);

#if HAVE_ECL_INPUT
    /*!
     * \brief Sets the pressure-dependent water viscosity and density
     *        using a table stemming from the Eclipse PVTWSALT keyword.
     */
    void initFromState(const EclipseState& eclState, const Schedule&);
#endif

    void setNumRegions(size_t numRegions);

    /*!
     * \brief Set the water reference density [kg / m^3]
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar /*rhoRefOil*/,
                               Scalar /*rhoRefGas*/,
                               Scalar rhoRefWater);

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
                        const Evaluation&) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& saltconcentration) const;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& /*temperature*/,
                                            const Evaluation& pressure,
                                            const Evaluation& saltconcentration) const;

    const Scalar waterReferenceDensity(unsigned regionIdx) const
    { return waterReferenceDensity_[regionIdx]; }

    const std::vector<Scalar>& referencePressure() const
    { return referencePressure_; }

    const std::vector<TabulatedFunction>& formationVolumeTables() const
    { return formationVolumeTables_; }

    const std::vector<TabulatedFunction>& compressibilityTables() const
    { return compressibilityTables_; }

    const std::vector<TabulatedFunction>& viscosityTables() const
    { return viscosityTables_; }

    const std::vector<TabulatedFunction>& viscosibilityTables() const
    { return viscosibilityTables_; }

    bool operator==(const ConstantCompressibilityBrinePvt<Scalar>& data) const;

private:
    std::vector<Scalar> waterReferenceDensity_;
    std::vector<Scalar> referencePressure_;
    std::vector<TabulatedFunction> formationVolumeTables_;
    std::vector<TabulatedFunction> compressibilityTables_;
    std::vector<TabulatedFunction> viscosityTables_;
    std::vector<TabulatedFunction> viscosibilityTables_;

};

} // namespace Opm

#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityBrinePvt_impl.hpp>

#endif
