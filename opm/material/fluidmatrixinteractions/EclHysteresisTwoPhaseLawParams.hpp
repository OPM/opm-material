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
 * \copydoc Opm::EclHysteresisTwoPhaseLawParams
 */
#ifndef OPM_ECL_HYSTERESIS_TWO_PHASE_LAW_PARAMS_HPP
#define OPM_ECL_HYSTERESIS_TWO_PHASE_LAW_PARAMS_HPP

#include "EclHysteresisConfig.hpp"
#include "EclEpsScalingPoints.hpp"

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

#include <string>
#include <memory>
#include <cassert>
#include <algorithm>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief A default implementation of the parameters for the material law which
 *        implements the ECL relative permeability and capillary pressure hysteresis
 */
template <class EffLawT>
class EclHysteresisTwoPhaseLawParams
{
    typedef typename EffLawT::Params EffLawParams;
    typedef typename EffLawParams::Traits::Scalar Scalar;

public:
    typedef typename EffLawParams::Traits Traits;

    EclHysteresisTwoPhaseLawParams()
    {
        pcSwMdc_ = 2.0;
        krnSwMdc_ = 2.0;
        // krwSwMdc_ = 2.0;

        deltaSwImbKrn_ = 0.0;
        // deltaSwImbKrw_ = 0.0;
        deltaSwImbPc_ = 0.0;

#ifndef NDEBUG
        finalized_ = false;
#endif
    }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        if (config().enableHysteresis()) {
            //C_ = 1.0/(Sncri_ - Sncrd_) + 1.0/(Snmaxd_ - Sncrd_);

            updateDynamicParams_();
        }

#ifndef NDEBUG
        finalized_ = true;
#endif
    }

    /*!
     * \brief Set the endpoint scaling configuration object.
     */
    void setConfig(std::shared_ptr<EclHysteresisConfig> value)
    { config_ = value; }

    /*!
     * \brief Returns the endpoint scaling configuration object.
     */
    const EclHysteresisConfig& config() const
    { return *config_; }

    /*!
     * \brief Sets the parameters used for the drainage curve
     */
    void setDrainageParams(std::shared_ptr<EffLawParams> value,
                           const EclEpsScalingPointsInfo<Scalar>& /* info */,
                           EclTwoPhaseSystemType /* twoPhaseSystem */)

    {
        drainageParams_ = *value;
    }

    /*!
     * \brief Returns the parameters used for the drainage curve
     */
    const EffLawParams& drainageParams() const
    { return drainageParams_; }

    EffLawParams& drainageParams()
    { return drainageParams_; }

    /*!
     * \brief Sets the parameters used for the imbibition curve
     */
    void setImbibitionParams(std::shared_ptr<EffLawParams> value,
                             const EclEpsScalingPointsInfo<Scalar>& /* info */,
                             EclTwoPhaseSystemType /* twoPhaseSystem */)
    {
        imbibitionParams_ = *value;

/*
        if (twoPhaseSystem == EclGasOilSystem) {
            Sncri_ = info.Sgcr;
        }
        else {
            assert(twoPhaseSystem == EclOilWaterSystem);
            Sncri_ = info.Sowcr;
        }
*/
    }

    /*!
     * \brief Returns the parameters used for the imbibition curve
     */
    const EffLawParams& imbibitionParams() const
    { return imbibitionParams_; }

    EffLawParams& imbibitionParams()
    { return imbibitionParams_; }

    /*!
     * \brief Set the saturation of the wetting phase where the last switch from the main
     *        drainage curve (MDC) to imbibition happend on the capillary pressure curve.
     */
    void setPcSwMdc(Scalar value)
    { pcSwMdc_ = value; }

    /*!
     * \brief Set the saturation of the wetting phase where the last switch from the main
     *        drainage curve to imbibition happend on the capillary pressure curve.
     */
    Scalar pcSwMdc() const
    { return pcSwMdc_; }

    /*!
     * \brief Set the saturation of the wetting phase where the last switch from the main
     *        drainage curve (MDC) to imbibition happend on the relperm curve for the
     *        wetting phase.
     */
    void setKrwSwMdc(Scalar /* value */)
    {}
    //    { krwSwMdc_ = value; };

    /*!
     * \brief Set the saturation of the wetting phase where the last switch from the main
     *        drainage curve to imbibition happend on the relperm curve for the
     *        wetting phase.
     */
    Scalar krwSwMdc() const
    { return 2.0; }
    //    { return krwSwMdc_; };

    /*!
     * \brief Set the saturation of the wetting phase where the last switch from the main
     *        drainage curve (MDC) to imbibition happend on the relperm curve for the
     *        non-wetting phase.
     */
    void setKrnSwMdc(Scalar value)
    { krnSwMdc_ = value; }

    /*!
     * \brief Set the saturation of the wetting phase where the last switch from the main
     *        drainage curve to imbibition happend on the relperm curve for the
     *        non-wetting phase.
     */
    Scalar krnSwMdc() const
    { return krnSwMdc_; }

    /*!
     * \brief Sets the saturation value which must be added if krw is calculated using
     *        the imbibition curve.
     *
     * This means that krw(Sw) = krw_drainage(Sw) if Sw < SwMdc and
     * krw(Sw) = krw_imbibition(Sw + Sw_shift,krw) else
     */
    void setDeltaSwImbKrw(Scalar value)
    { deltaSwImbKrw_ = value; }

    /*!
     * \brief Returns the saturation value which must be added if krw is calculated using
     *        the imbibition curve.
     *
     * This means that krw(Sw) = krw_drainage(Sw) if Sw < SwMdc and
     * krw(Sw) = krw_imbibition(Sw + Sw_shift,krw) else
     */
    Scalar deltaSwImbKrw() const
    { return deltaSwImbKrw_; }

    /*!
     * \brief Sets the saturation value which must be added if krn is calculated using
     *        the imbibition curve.
     *
     * This means that krn(Sw) = krn_drainage(Sw) if Sw < SwMdc and
     * krn(Sw) = krn_imbibition(Sw + Sw_shift,krn) else
     */
    void setDeltaSwImbKrn(Scalar value)
    { deltaSwImbKrn_ = value; }

    /*!
     * \brief Returns the saturation value which must be added if krn is calculated using
     *        the imbibition curve.
     *
     * This means that krn(Sw) = krn_drainage(Sw) if Sw < SwMdc and
     * krn(Sw) = krn_imbibition(Sw + Sw_shift,krn) else
     */
    Scalar deltaSwImbKrn() const
    { return deltaSwImbKrn_; }

    /*!
     * \brief Sets the saturation value which must be added if the capillary pressure is
     *        calculated using the imbibition curve.
     *
     * This means that p_cnw(Sw) = pc_cnw,drainage(Sw) if Sw < SwMdc and pcnw(Sw) =
     * pcnw_imbibition(Sw + Sw_shift,pc) else
     */
    void setDeltaSwImbPc(Scalar value)
    { deltaSwImbPc_ = value; }

    /*!
     * \brief Returns the saturation value which must be added if the capillary pressure is
     *        calculated using the imbibition curve.
     *
     * This means that p_cnw(Sw) = pc_cnw,drainage(Sw) if Sw < SwMdc and pcnw(Sw) =
     * pcnw_imbibition(Sw + Sw_shift,pc) else
     */
    Scalar deltaSwImbPc() const
    { return deltaSwImbPc_; }

    /*!
     * \brief Notify the hysteresis law that a given wetting-phase saturation has been seen
     *
     * This updates the scanning curves and the imbibition<->drainage reversal points as
     * appropriate.
     */
    void update(Scalar pcSw, Scalar krwSw, Scalar krnSw)
    {
        bool updateParams = false;
        if (pcSw < pcSwMdc_) {
            pcSwMdc_ = pcSw;
            updateParams = true;
        }

        if (krwSw < krwSwMdc_) {
            krwSwMdc_ = krwSw;
            updateParams = true;
        }

        if (krnSw < krnSwMdc_) {
            krnSwMdc_ = krnSw;
            updateParams = true;
        }

        if (updateParams)
            updateDynamicParams_();
    }

private:

#ifndef NDEBUG
    void assertFinalized_() const
    { assert(finalized_); }

    bool finalized_;
#else
    void assertFinalized_() const
    { }
#endif

    void updateDynamicParams_()
    {
        // calculate the saturation deltas for the relative permeabilities
        Scalar krwMdcDrainage = EffLawT::twoPhaseSatKrw(drainageParams(), krwSwMdc_);
        Scalar SwKrwMdcImbibition = EffLawT::twoPhaseSatKrwInv(imbibitionParams(), krwMdcDrainage);
        deltaSwImbKrw_ = SwKrwMdcImbibition - krwSwMdc_;

        Scalar krnMdcDrainage = EffLawT::twoPhaseSatKrn(drainageParams(), krnSwMdc_);
        Scalar SwKrnMdcImbibition = EffLawT::twoPhaseSatKrnInv(imbibitionParams(), krnMdcDrainage);
        deltaSwImbKrn_ = SwKrnMdcImbibition - krnSwMdc_;

        Scalar pcMdcDrainage = EffLawT::twoPhaseSatPcnw(drainageParams(), pcSwMdc_);
        Scalar SwPcMdcImbibition = EffLawT::twoPhaseSatPcnwInv(imbibitionParams(), pcMdcDrainage);
        deltaSwImbPc_ = SwPcMdcImbibition - pcSwMdc_;

        assert(std::abs(EffLawT::twoPhaseSatPcnw(imbibitionParams(), pcSwMdc_ + deltaSwImbPc_)
                        - EffLawT::twoPhaseSatPcnw(drainageParams(), pcSwMdc_)) < 1e-8);
        assert(std::abs(EffLawT::twoPhaseSatKrn(imbibitionParams(), krnSwMdc_ + deltaSwImbKrn_)
                        - EffLawT::twoPhaseSatKrn(drainageParams(), krnSwMdc_)) < 1e-8);
        assert(std::abs(EffLawT::twoPhaseSatKrw(imbibitionParams(), krwSwMdc_ + deltaSwImbKrw_)
                        - EffLawT::twoPhaseSatKrw(drainageParams(), krwSwMdc_)) < 1e-8);
    }

    std::shared_ptr<EclHysteresisConfig> config_;
    EffLawParams imbibitionParams_;
    EffLawParams drainageParams_;

    // largest wettinging phase saturation which is on the main-drainage curve. These are
    // three different values because the sourounding code can choose to use different
    // definitions for the saturations for different quantities
    Scalar krwSwMdc_;
    Scalar krnSwMdc_;
    Scalar pcSwMdc_;

    // offsets added to wetting phase saturation uf using the imbibition curves need to
    // be used to calculate the wetting phase relperm, the non-wetting phase relperm and
    // the capillary pressure
    Scalar deltaSwImbKrw_;
    Scalar deltaSwImbKrn_;
    Scalar deltaSwImbPc_;
};

} // namespace Opm

#endif
