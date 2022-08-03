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
 * \copydoc Opm::GasPvtThermal
 */
#ifndef OPM_GAS_PVT_THERMAL_HPP
#define OPM_GAS_PVT_THERMAL_HPP

#include <opm/material/common/Tabulated1DFunction.hpp>

namespace Opm {
template <class Scalar, bool enableThermal>
class GasPvtMultiplexer;

class EclipseState;
class Schedule;

/*!
 * \brief This class implements temperature dependence of the PVT properties of gas
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class GasPvtThermal
{
public:
    using IsothermalPvt = GasPvtMultiplexer<Scalar, /*enableThermal=*/false>;
    using TabulatedOneDFunction = Tabulated1DFunction<Scalar>;

    GasPvtThermal()
    {
        enableThermalDensity_ = false;
        enableJouleThomson_ = false;
        enableThermalViscosity_ = false;
        isothermalPvt_ = nullptr;
    }

    GasPvtThermal(IsothermalPvt* isothermalPvt,
                  const std::vector<TabulatedOneDFunction>& gasvisctCurves,
                  const std::vector<Scalar>& gasdentRefTemp,
                  const std::vector<Scalar>& gasdentCT1,
                  const std::vector<Scalar>& gasdentCT2,
                  const std::vector<Scalar>& gasJTRefPres,
                  const std::vector<Scalar>& gasJTC,
                  const std::vector<TabulatedOneDFunction>& internalEnergyCurves,
                  bool enableThermalDensity,
                  bool enableJouleThomson,
                  bool enableThermalViscosity,
                  bool enableInternalEnergy);

    GasPvtThermal(const GasPvtThermal& data);

    ~GasPvtThermal();

#if HAVE_ECL_INPUT
    /*!
     * \brief Implement the temperature part of the gas PVT properties.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule);
#endif // HAVE_ECL_INPUT

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions);

    /*!
     * \brief Finish initializing the thermal part of the gas phase PVT properties.
     */
    void initEnd()
    { }

    size_t numRegions() const
    { return gasvisctCurves_.size(); }

    /*!
     * \brief Returns true iff the density of the gas phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

     /*!
     * \brief Returns true iff Joule-Thomson effect for the gas phase is active.
     */
    bool enableJouleThomson() const
    { return enableJouleThomson_; }

    /*!
     * \brief Returns true iff the viscosity of the gas phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    /*!
     * \brief Returns the specific internal energy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure,
                              const Evaluation& Rv) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rv,
                         const Evaluation& Rvw) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the oil-saturated gas phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const;

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rv,
                                            const Evaluation& /*Rvw*/) const;

    /*!
     * \brief Returns the formation volume factor [-] of oil-saturated gas.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const;

    /*!
     * \brief Returns the water vaporization factor \f$R_v\f$ [m^3/m^3] of the water phase.
     */
    template <class Evaluation>
    Evaluation saturatedWaterVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/) const;
    

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the gas phase.
     *
     * This method implements temperature dependence and requires the gas pressure,
     * temperature and the oil saturation as inputs. Currently it is just a dummy method
     * which passes through the isothermal oil vaporization factor.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure) const;

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the gas phase.
     *
     * This method implements temperature dependence and requires the gas pressure,
     * temperature and the oil saturation as inputs. Currently it is just a dummy method
     * which passes through the isothermal oil vaporization factor.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure,
                                              const Evaluation& oilSaturation,
                                              const Evaluation& maxOilSaturation) const;

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *
     * This method implements temperature dependence and requires isothermal satuation
     * pressure and temperature as inputs. Currently it is just a dummy method which
     * passes through the isothermal saturation pressure.
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const;

    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& temperature,
                                    const Evaluation& pressure,
                                    unsigned compIdx) const;

    const IsothermalPvt* isoThermalPvt() const
    { return isothermalPvt_; }

    const Scalar gasReferenceDensity(unsigned regionIdx) const
    { return isothermalPvt_->gasReferenceDensity(regionIdx); }

    const std::vector<TabulatedOneDFunction>& gasvisctCurves() const
    { return gasvisctCurves_; }

    const std::vector<Scalar>& gasdentRefTemp() const
    { return gasdentRefTemp_; }

    const std::vector<Scalar>& gasdentCT1() const
    { return gasdentCT1_; }

    const std::vector<Scalar>& gasdentCT2() const
    { return gasdentCT2_; }

    const std::vector<TabulatedOneDFunction>& internalEnergyCurves() const
    { return internalEnergyCurves_; }

    bool enableInternalEnergy() const
    { return enableInternalEnergy_; }

    const std::vector<Scalar>& gasJTRefPres() const
    { return  gasJTRefPres_; }

     const std::vector<Scalar>&  gasJTC() const
    { return gasJTC_; }


    bool operator==(const GasPvtThermal<Scalar>& data) const;

    GasPvtThermal<Scalar>& operator=(const GasPvtThermal<Scalar>& data);

private:
    IsothermalPvt* isothermalPvt_;

    // The PVT properties needed for temperature dependence of the viscosity. We need
    // to store one value per PVT region.
    std::vector<TabulatedOneDFunction> gasvisctCurves_;

    std::vector<Scalar> gasdentRefTemp_;
    std::vector<Scalar> gasdentCT1_;
    std::vector<Scalar> gasdentCT2_;

    std::vector<Scalar> gasJTRefPres_;
    std::vector<Scalar> gasJTC_;

    std::vector<Scalar> rhoRefO_;

    // piecewise linear curve representing the internal energy of gas
    std::vector<TabulatedOneDFunction> internalEnergyCurves_;

    bool enableThermalDensity_;
    bool enableJouleThomson_;
    bool enableThermalViscosity_;
    bool enableInternalEnergy_;
};

} // namespace Opm

#ifndef OPM_USE_PRIVATE_TEMPLATES
#include <opm/material/fluidsystems/blackoilpvt/GasPvtThermal_impl.hpp>
#endif

#endif
