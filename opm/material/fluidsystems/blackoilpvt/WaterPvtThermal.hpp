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
 * \copydoc Opm::WaterPvtThermal
 */
#ifndef OPM_WATER_PVT_THERMAL_HPP
#define OPM_WATER_PVT_THERMAL_HPP

#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

namespace Opm {
template <class Scalar, bool enableThermal, bool enableBrine>
class WaterPvtMultiplexer;

class EclipseState;
class Schedule;

/*!
 * \brief This class implements temperature dependence of the PVT properties of water
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar, bool enableBrine>
class WaterPvtThermal
{
public:
    using TabulatedOneDFunction = Tabulated1DFunction<Scalar>;
    using IsothermalPvt = WaterPvtMultiplexer<Scalar, /*enableThermal=*/false, enableBrine>;

    WaterPvtThermal()
    {
        enableThermalDensity_ = false;
        enableThermalViscosity_ = false;
        enableInternalEnergy_ = false;
        isothermalPvt_ = nullptr;
    }

    WaterPvtThermal(IsothermalPvt* isothermalPvt,
                    const std::vector<Scalar>& viscrefPress,
                    const std::vector<Scalar>& watdentRefTemp,
                    const std::vector<Scalar>& watdentCT1,
                    const std::vector<Scalar>& watdentCT2,
                    const std::vector<Scalar>& watJTRefPres,
                    const std::vector<Scalar>& watJTC,
                    const std::vector<Scalar>& pvtwRefPress,
                    const std::vector<Scalar>& pvtwRefB,
                    const std::vector<Scalar>& pvtwCompressibility,
                    const std::vector<Scalar>& pvtwViscosity,
                    const std::vector<Scalar>& pvtwViscosibility,
                    const std::vector<TabulatedOneDFunction>& watvisctCurves,
                    const std::vector<TabulatedOneDFunction>& internalEnergyCurves,
                    bool enableThermalDensity,
                    bool enableJouleThomson,
                    bool enableThermalViscosity,
                    bool enableInternalEnergy);

    WaterPvtThermal(const WaterPvtThermal& data);

    ~WaterPvtThermal();

#if HAVE_ECL_INPUT
    /*!
     * \brief Implement the temperature part of the water PVT properties.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule);
#endif // HAVE_ECL_INPUT

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions);

    /*!
     * \brief Finish initializing the thermal part of the water phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Returns true iff the density of the water phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

     /*!
     * \brief Returns true iff Joule-Thomson effect for the water phase is active.
     */
    bool enableJouleThomson() const
    { return enableJouleThomson_; }

    /*!
     * \brief Returns true iff the viscosity of the water phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    size_t numRegions() const
    { return pvtwRefPress_.size(); }

    /*!
     * \brief Returns the specific internal energy [J/kg] of water given a set of parameters.
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

    const IsothermalPvt* isoThermalPvt() const
    { return isothermalPvt_; }

    const Scalar waterReferenceDensity(unsigned regionIdx) const
    { return isothermalPvt_->waterReferenceDensity(regionIdx); }

    const std::vector<Scalar>& viscrefPress() const
    { return viscrefPress_; }

    const std::vector<Scalar>& watdentRefTemp() const
    { return watdentRefTemp_; }

    const std::vector<Scalar>& watdentCT1() const
    { return watdentCT1_; }

    const std::vector<Scalar>& watdentCT2() const
    { return watdentCT2_; }

    const std::vector<Scalar>& pvtwRefPress() const
    { return pvtwRefPress_; }

    const std::vector<Scalar>& pvtwRefB() const
    { return pvtwRefB_; }

    const std::vector<Scalar>& pvtwCompressibility() const
    { return pvtwCompressibility_; }

    const std::vector<Scalar>& pvtwViscosity() const
    { return pvtwViscosity_; }

    const std::vector<Scalar>& pvtwViscosibility() const
    { return pvtwViscosibility_; }

    const std::vector<TabulatedOneDFunction>& watvisctCurves() const
    { return watvisctCurves_; }

    const std::vector<TabulatedOneDFunction> internalEnergyCurves() const
    { return internalEnergyCurves_; }

    bool enableInternalEnergy() const
    { return enableInternalEnergy_; }

    const std::vector<Scalar>& watJTRefPres() const
    { return  watJTRefPres_; }

     const std::vector<Scalar>&  watJTC() const
    { return watJTC_; }


    bool operator==(const WaterPvtThermal<Scalar,enableBrine>& data) const;

    WaterPvtThermal<Scalar, enableBrine>& operator=(const WaterPvtThermal<Scalar,enableBrine>& data);

private:
    IsothermalPvt* isothermalPvt_;

    // The PVT properties needed for temperature dependence. We need to store one
    // value per PVT region.
    std::vector<Scalar> viscrefPress_;

    std::vector<Scalar> watdentRefTemp_;
    std::vector<Scalar> watdentCT1_;
    std::vector<Scalar> watdentCT2_;

    std::vector<Scalar> watJTRefPres_;
    std::vector<Scalar> watJTC_;

    std::vector<Scalar> pvtwRefPress_;
    std::vector<Scalar> pvtwRefB_;
    std::vector<Scalar> pvtwCompressibility_;
    std::vector<Scalar> pvtwViscosity_;
    std::vector<Scalar> pvtwViscosibility_;

    std::vector<TabulatedOneDFunction> watvisctCurves_;

    // piecewise linear curve representing the internal energy of water
    std::vector<TabulatedOneDFunction> internalEnergyCurves_;

    bool enableThermalDensity_;
    bool enableJouleThomson_;
    bool enableThermalViscosity_;
    bool enableInternalEnergy_;
};

} // namespace Opm

#include <opm/material/fluidsystems/blackoilpvt/WaterPvtThermal_impl.hpp>

#endif
