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
 * \copydoc Opm::OilPvtThermal
 */
#ifndef OPM_OIL_PVT_THERMAL_HPP
#define OPM_OIL_PVT_THERMAL_HPP

#include <opm/material/common/Tabulated1DFunction.hpp>


namespace Opm {
template <class Scalar, bool enableThermal>
class OilPvtMultiplexer;

class EclipseState;
class Schedule;

/*!
 * \brief This class implements temperature dependence of the PVT properties of oil
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class OilPvtThermal
{
public:
    using TabulatedOneDFunction = Tabulated1DFunction<Scalar>;
    using IsothermalPvt = OilPvtMultiplexer<Scalar, /*enableThermal=*/false>;

    OilPvtThermal()
    {
        enableThermalDensity_ = false;
        enableJouleThomson_ = false;
        enableThermalViscosity_ = false;
        enableInternalEnergy_ = false;
        isothermalPvt_ = nullptr;
    }

    OilPvtThermal(IsothermalPvt* isothermalPvt,
                  const std::vector<TabulatedOneDFunction>& oilvisctCurves,
                  const std::vector<Scalar>& viscrefPress,
                  const std::vector<Scalar>& viscrefRs,
                  const std::vector<Scalar>& viscRef,
                  const std::vector<Scalar>& oildentRefTemp,
                  const std::vector<Scalar>& oildentCT1,
                  const std::vector<Scalar>& oildentCT2,
                  const std::vector<Scalar>& oilJTRefPres,
                  const std::vector<Scalar>& oilJTC,
                  const std::vector<TabulatedOneDFunction>& internalEnergyCurves,
                  bool enableThermalDensity,
                  bool enableJouleThomson,
                  bool enableThermalViscosity,
                  bool enableInternalEnergy);

    OilPvtThermal(const OilPvtThermal& data);

    ~OilPvtThermal();

#if HAVE_ECL_INPUT
    /*!
     * \brief Implement the temperature part of the oil PVT properties.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule);
#endif // HAVE_ECL_INPUT

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions);

    /*!
     * \brief Finish initializing the thermal part of the oil phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Returns true iff the density of the oil phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

    /*!
     * \brief Returns true iff Joule-Thomson effect for the oil phase is active.
     */
    bool enableJouleThomson() const
    { return enableJouleThomson_; }

    /*!
     * \brief Returns true iff the viscosity of the oil phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    size_t numRegions() const
    { return viscrefRs_.size(); }

    /*!
     * \brief Returns the specific internal energy [J/kg] of oil given a set of parameters.
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
                         const Evaluation& Rs) const;

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
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
                                            const Evaluation& Rs) const;

    /*!
     * \brief Returns the formation volume factor [-] of gas-saturated oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const;

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     *
     * This method implements temperature dependence and requires the isothermal gas
     * dissolution factor for gas saturated oil and temperature as inputs. Currently it
     * is just a dummy method which passes through the isothermal gas dissolution factor.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const;

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     *
     * This method implements temperature dependence and requires the isothermal gas
     * dissolution factor for gas saturated oil and temperature as inputs. Currently it
     * is just a dummy method which passes through the isothermal gas dissolution factor.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& oilSaturation,
                                             const Evaluation& maxOilSaturation) const;

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
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

    const Scalar oilReferenceDensity(unsigned regionIdx) const
    { return isothermalPvt_->oilReferenceDensity(regionIdx); }

    const std::vector<TabulatedOneDFunction>& oilvisctCurves() const
    { return oilvisctCurves_; }

    const std::vector<Scalar>& viscrefPress() const
    { return viscrefPress_; }

    const std::vector<Scalar>& viscrefRs() const
    { return viscrefRs_; }

    const std::vector<Scalar>& viscRef() const
    { return viscRef_; }

    const std::vector<Scalar>& oildentRefTemp() const
    { return oildentRefTemp_; }

    const std::vector<Scalar>& oildentCT1() const
    { return oildentCT1_; }

    const std::vector<Scalar>& oildentCT2() const
    { return oildentCT2_; }

    const std::vector<TabulatedOneDFunction> internalEnergyCurves() const
    { return internalEnergyCurves_; }

    bool enableInternalEnergy() const
    { return enableInternalEnergy_; }

    const std::vector<Scalar>& oilJTRefPres() const
    { return  oilJTRefPres_; }

     const std::vector<Scalar>&  oilJTC() const
    { return oilJTC_; }

    bool operator==(const OilPvtThermal<Scalar>& data) const;

    OilPvtThermal<Scalar>& operator=(const OilPvtThermal<Scalar>& data);

private:
    IsothermalPvt* isothermalPvt_;

    // The PVT properties needed for temperature dependence of the viscosity. We need
    // to store one value per PVT region.
    std::vector<TabulatedOneDFunction> oilvisctCurves_;
    std::vector<Scalar> viscrefPress_;
    std::vector<Scalar> viscrefRs_;
    std::vector<Scalar> viscRef_;

    // The PVT properties needed for temperature dependence of the density.
    std::vector<Scalar> oildentRefTemp_;
    std::vector<Scalar> oildentCT1_;
    std::vector<Scalar> oildentCT2_;

    std::vector<Scalar> oilJTRefPres_;
    std::vector<Scalar> oilJTC_;

    std::vector<Scalar> rhoRefG_;

    // piecewise linear curve representing the internal energy of oil
    std::vector<TabulatedOneDFunction> internalEnergyCurves_;

    bool enableThermalDensity_;
    bool enableJouleThomson_;
    bool enableThermalViscosity_;
    bool enableInternalEnergy_;
};

} // namespace Opm

#ifndef OPM_USE_PRIVATE_TEMPLATES
#include <opm/material/fluidsystems/blackoilpvt/OilPvtThermal_impl.hpp>
#endif

#endif
