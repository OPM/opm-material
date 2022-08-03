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

#include <vector>

namespace Opm {

class EclipseState;
class Schedule;

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
                                    const std::vector<Scalar>& waterViscosibility);

#if HAVE_ECL_INPUT
    /*!
     * \brief Sets the pressure-dependent water viscosity and density
     *        using a table stemming from the Eclipse PVTW keyword.
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
     * \brief Set the water reference pressure [Pa]
     */
    void setReferencePressure(unsigned regionIdx, Scalar p);

    /*!
     * \brief Set the viscosity and "viscosibility" of the water phase.
     */
    void setViscosity(unsigned regionIdx, Scalar muw, Scalar waterViscosibility = 0.0);

    /*!
     * \brief Set the compressibility of the water phase.
     */
    void setCompressibility(unsigned regionIdx, Scalar waterCompressibility);

    /*!
     * \brief Set the water reference formation volume factor [-]
     */
    void setReferenceFormationVolumeFactor(unsigned regionIdx, Scalar BwRef);

    /*!
     * \brief Set the water "viscosibility" [1/ (Pa s)]
     */
    void setViscosibility(unsigned regionIdx, Scalar muComp);

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
                                            const Evaluation& /*saltconcentration*/) const;

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

    bool operator==(const ConstantCompressibilityWaterPvt<Scalar>& data) const;

private:
    std::vector<Scalar> waterReferenceDensity_;
    std::vector<Scalar> waterReferencePressure_;
    std::vector<Scalar> waterReferenceFormationVolumeFactor_;
    std::vector<Scalar> waterCompressibility_;
    std::vector<Scalar> waterViscosity_;
    std::vector<Scalar> waterViscosibility_;
};

} // namespace Opm

#ifndef OPM_USE_PRIVATE_TEMPLATES
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt_impl.hpp>
#endif

#endif
