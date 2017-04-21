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

#include <opm/material/Constants.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {
template <class Scalar, bool enableThermal>
class OilPvtMultiplexer;

/*!
 * \brief This class implements temperature dependence of the PVT properties of oil
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class OilPvtThermal
{
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef OilPvtMultiplexer<Scalar, /*enableThermal=*/false> IsothermalPvt;

public:
    ~OilPvtThermal()
    { delete isothermalPvt_; }

#if HAVE_OPM_PARSER
    /*!
     * \brief Implement the temperature part of the oil PVT properties.
     */
    void initFromDeck(const EclipseState& eclState)
    {
        //////
        // initialize the isothermal part
        //////
        isothermalPvt_ = new IsothermalPvt;
        isothermalPvt_->initFromDeck(eclState);

        //////
        // initialize the thermal part
        //////
        const auto& tables = eclState.getTableManager();

        enableThermalDensity_ = !tables.get_refs_material().THERMEX1.empty();
        enableThermalViscosity_ = !tables.getViscrefTable().empty();

        unsigned numRegions = isothermalPvt_->numRegions();
        setNumRegions(numRegions);

        // viscosity
        if (!tables.getViscrefTable().empty()) {
            const auto& oilvisctTables = tables.getOilvisctTables();
            const auto& viscrefTable = tables.getViscrefTable();

            assert(oilvisctTables.size() == numRegions);
            assert(viscrefTable.size() == numRegions);

            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& TCol = oilvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
                const auto& muCol = oilvisctTables[regionIdx].getColumn("Viscosity").vectorCopy();
                oilvisctCurves_[regionIdx].setXYContainers(TCol, muCol);

                viscrefPress_[regionIdx] = viscrefTable[regionIdx].reference_pressure;
                viscrefRs_[regionIdx] = viscrefTable[regionIdx].reference_rs;

                // temperature used to calculate the reference viscosity [K]. the
                // value does not really matter if the underlying PVT object really
                // is isothermal...
                Scalar Tref = 273.15 + 20;

                // compute the reference viscosity using the isothermal PVT object.
                viscRef_[regionIdx] =
                    isothermalPvt_->viscosity(regionIdx,
                                              Tref,
                                              viscrefPress_[regionIdx],
                                              viscrefRs_[regionIdx]);
            }
        }

        // quantities required for density. note that we just always use the values
        // for the first EOS. (since EOS != PVT region.)
        refTemp_ = 0.0;
        if (enableThermalDensity_) {
            const auto& refs = tables.get_refs_material();
            int oilCompIdx = refs.oil_component_index;

            // always use the values of the first EOS
            refTemp_ = refs.TREF.at(oilCompIdx);
            refPress_ = refs.PREF.at(oilCompIdx);
            refC_ = refs.CREF.at(oilCompIdx);
            thermex1_ = refs.THERMEX1.at(oilCompIdx);
        }
    }
#endif // HAVE_OPM_PARSER

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions)
    {
        oilvisctCurves_.resize(numRegions);
        viscrefPress_.resize(numRegions);
        viscrefRs_.resize(numRegions);
        viscRef_.resize(numRegions);
    }

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
     * \brief Returns true iff the viscosity of the oil phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    size_t numRegions() const
    { return viscrefRs_.size(); }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rs) const
    {
        const auto& isothermalMu = isothermalPvt_->viscosity(regionIdx, temperature, pressure, Rs);
        if (!enableThermalViscosity())
            return isothermalMu;

        // compute the viscosity deviation due to temperature
        const auto& muOilvisct = oilvisctCurves_[regionIdx].eval(temperature);
        return muOilvisct/viscRef_[regionIdx]*isothermalMu;
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        const auto& isothermalMu = isothermalPvt_->saturatedViscosity(regionIdx, temperature, pressure);
        if (!enableThermalViscosity())
            return isothermalMu;

        // compute the viscosity deviation due to temperature
        const auto& muOilvisct = oilvisctCurves_[regionIdx].eval(temperature);
        return muOilvisct/viscRef_[regionIdx]*isothermalMu;
    }


    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const
    {
        const auto& b =
            isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs);
        if (!enableThermalDensity())
            return b;

        // we use equation (3.208) from the Eclipse 2011.1 Reference Manual, but we
        // calculate rho_ref using the isothermal keyword instead of using the value for
        // the components, so the oil compressibility is already dealt with there. Note
        // that we only do the part for the oil component here, the part for dissolved
        // gas is ignored so far.
        const auto& alpha = 1.0/(1 + thermex1_*(temperature - refTemp_));
        return alpha*b;
    }

    /*!
     * \brief Returns the formation volume factor [-] of gas-saturated oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    {
        const auto& b =
            isothermalPvt_->saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure);
        if (!enableThermalDensity())
            return b;

        // we use equation (3.208) from the Eclipse 2011.1 Reference Manual, but we
        // calculate rho_ref using the isothermal keyword instead of using the value for
        // the components, so the oil compressibility is already dealt with there. Note
        // that we only do the part for the oil component here, the part for dissolved
        // gas is ignored so far.
        const auto& alpha = 1.0/(1 + thermex1_*(temperature - refTemp_));
        return alpha*b;
    }

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
                                             const Evaluation& pressure) const
    { return isothermalPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure); }

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
                                             Scalar maxOilSaturation) const
    { return isothermalPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation); }

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
                                  const Evaluation& pressure) const
    { return isothermalPvt_->saturationPressure(regionIdx, temperature, pressure); }

private:
    IsothermalPvt* isothermalPvt_;

    // The PVT properties needed for temperature dependence of the viscosity. We need
    // to store one value per PVT region.
    std::vector<TabulatedOneDFunction> oilvisctCurves_;
    std::vector<Scalar> viscrefPress_;
    std::vector<Scalar> viscrefRs_;
    std::vector<Scalar> viscRef_;

    // The PVT properties needed for temperature dependence of the density. This is
    // specified as one value per EOS in the manual, but we unconditionally use the
    // expansion coefficient of the first EOS...
    Scalar refTemp_;
    Scalar refPress_;
    Scalar refC_;
    Scalar thermex1_;

    bool enableThermalDensity_;
    bool enableThermalViscosity_;
};

} // namespace Opm

#endif
