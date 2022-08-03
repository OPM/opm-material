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
 * \copydoc Opm::WaterPvtMultiplexer
 */
#ifndef OPM_WATER_PVT_MULTIPLEXER_HPP
#define OPM_WATER_PVT_MULTIPLEXER_HPP

#include "ConstantCompressibilityWaterPvt.hpp"
#include "ConstantCompressibilityBrinePvt.hpp"
#include "WaterPvtThermal.hpp"

namespace Opm {

class EclipseState;
class Schedule;

enum class WaterPvtApproach {
    NoWaterPvt,
    ConstantCompressibilityBrinePvt,
    ConstantCompressibilityWaterPvt,
    ThermalWaterPvt
};

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the water
 *        phase in the black-oil model.
 */
template <class Scalar, bool enableThermal = true, bool enableBrine = true>
class WaterPvtMultiplexer
{
public:
    WaterPvtMultiplexer();

    WaterPvtMultiplexer(WaterPvtApproach approach, void* realWaterPvt);

    WaterPvtMultiplexer(const WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>& data);

    ~WaterPvtMultiplexer();

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for water using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVDG keywords.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule);
#endif // HAVE_ECL_INPUT

    void initEnd();

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const;

    /*!
     * \brief Return the reference density which are considered by this PVT-object.
     */
    const Scalar waterReferenceDensity(unsigned regionIdx);

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const Evaluation& saltconcentration) const;

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
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& saltconcentration) const;

    void setApproach(WaterPvtApproach appr);

    /*!
     * \brief Returns the concrete approach for calculating the PVT relations.
     *
     * (This is only determined at runtime.)
     */
    WaterPvtApproach approach() const
    { return approach_; }

    // get the concrete parameter object for the water phase
    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == WaterPvtApproach::ConstantCompressibilityWaterPvt, ConstantCompressibilityWaterPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<ConstantCompressibilityWaterPvt<Scalar>* >(realWaterPvt_);
    }

    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == WaterPvtApproach::ConstantCompressibilityWaterPvt, const ConstantCompressibilityWaterPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<ConstantCompressibilityWaterPvt<Scalar>* >(realWaterPvt_);
    }

    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == WaterPvtApproach::ConstantCompressibilityBrinePvt, ConstantCompressibilityBrinePvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<ConstantCompressibilityBrinePvt<Scalar>* >(realWaterPvt_);
    }

    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == WaterPvtApproach::ConstantCompressibilityBrinePvt, const ConstantCompressibilityBrinePvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<ConstantCompressibilityBrinePvt<Scalar>* >(realWaterPvt_);
    }

    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == WaterPvtApproach::ThermalWaterPvt, WaterPvtThermal<Scalar, enableBrine> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<WaterPvtThermal<Scalar, enableBrine>* >(realWaterPvt_);
    }

    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == WaterPvtApproach::ThermalWaterPvt, const WaterPvtThermal<Scalar, enableBrine> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<WaterPvtThermal<Scalar, enableBrine>* >(realWaterPvt_);
    }

    const void* realWaterPvt() const { return realWaterPvt_; }

    bool operator==(const WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>& data) const;

    WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>& operator=(const WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>& data);

private:
    WaterPvtApproach approach_;
    void* realWaterPvt_;
};

} // namespace Opm

#ifndef OPM_USE_PRIVATE_TEMPLATES
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer_impl.hpp>
#endif

#endif
