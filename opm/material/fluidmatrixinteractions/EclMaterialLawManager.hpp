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
 * \copydoc Opm::EclMaterialLawManager
 */
#if ! HAVE_OPM_PARSER
#error "The opm-parser module is required to use the ECL material manager!"
#endif

#ifndef OPM_ECL_MATERIAL_LAW_MANAGER_HPP
#define OPM_ECL_MATERIAL_LAW_MANAGER_HPP

#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <algorithm>


namespace Opm {

/*!
 * \ingroup fluidmatrixinteractions
 *
 * \brief Provides an simple way to create and manage the material law objects
 *        for a complete ECL deck.
 */
template <class TraitsT>
class EclMaterialLawManager
{
private:
    typedef TraitsT Traits;
    typedef typename Traits::Scalar Scalar;
    enum { waterPhaseIdx = Traits::wettingPhaseIdx };
    enum { oilPhaseIdx = Traits::nonWettingPhaseIdx };
    enum { gasPhaseIdx = Traits::gasPhaseIdx };
    enum { numPhases = Traits::numPhases };

    typedef TwoPhaseMaterialTraits<Scalar, oilPhaseIdx, gasPhaseIdx> GasOilTraits;
    typedef TwoPhaseMaterialTraits<Scalar, waterPhaseIdx, oilPhaseIdx> OilWaterTraits;

    // the two-phase material law which is defined on effective (unscaled) saturations
    typedef PiecewiseLinearTwoPhaseMaterial<GasOilTraits> GasOilUnscaledTwoPhaseLaw;
    typedef PiecewiseLinearTwoPhaseMaterial<OilWaterTraits> OilWaterUnscaledTwoPhaseLaw;
    typedef typename GasOilUnscaledTwoPhaseLaw::Params GasOilEffectiveTwoPhaseParams;
    typedef typename OilWaterUnscaledTwoPhaseLaw::Params OilWaterEffectiveTwoPhaseParams;

    // the two-phase material law which is defined on absolute (scaled) saturations
    typedef EclEpsTwoPhaseLaw<GasOilUnscaledTwoPhaseLaw> GasOilScaledTwoPhaseLaw;
    typedef EclEpsTwoPhaseLaw<OilWaterUnscaledTwoPhaseLaw> OilWaterScaledTwoPhaseLaw;
    typedef typename GasOilScaledTwoPhaseLaw::Params GasOilEpsTwoPhaseParams;
    typedef typename OilWaterScaledTwoPhaseLaw::Params OilWaterEpsTwoPhaseParams;

    // the scaled two-phase material laws with hystersis
    typedef EclHysteresisTwoPhaseLaw<GasOilScaledTwoPhaseLaw> GasOilTwoPhaseLaw;
    typedef EclHysteresisTwoPhaseLaw<OilWaterScaledTwoPhaseLaw> OilWaterTwoPhaseLaw;
    typedef typename GasOilTwoPhaseLaw::Params GasOilTwoPhaseParams;
    typedef typename OilWaterTwoPhaseLaw::Params OilWaterTwoPhaseParams;

public:
    // the three-phase material law used by the simulation
    typedef EclMultiplexerMaterial<Traits, GasOilTwoPhaseLaw, OilWaterTwoPhaseLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

private:
    // internal typedefs
    typedef std::vector<std::shared_ptr<GasOilEffectiveTwoPhaseParams> > GasOilEffectiveParamVector;
    typedef std::vector<std::shared_ptr<OilWaterEffectiveTwoPhaseParams> > OilWaterEffectiveParamVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar> > > GasOilScalingPointsVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar> > > OilWaterScalingPointsVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPointsInfo<Scalar> > > GasOilScalingInfoVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPointsInfo<Scalar> > > OilWaterScalingInfoVector;
    typedef std::vector<std::shared_ptr<GasOilTwoPhaseParams> > GasOilParamVector;
    typedef std::vector<std::shared_ptr<OilWaterTwoPhaseParams> > OilWaterParamVector;
    typedef std::vector<std::shared_ptr<MaterialLawParams> > MaterialLawParamsVector;

public:
    EclMaterialLawManager()
    {}

    void initFromDeck(const Opm::Deck& deck,
                      const Opm::EclipseState& eclState,
                      const std::vector<int>& compressedToCartesianElemIdx)
    {
        // get the number of saturation regions and the number of cells in the deck
        const size_t numSatRegions = eclState.runspec().tabdims().getNumSatTables();
        size_t numCompressedElems = compressedToCartesianElemIdx.size();

        // copy the SATNUM grid property. in some cases this is not necessary, but it
        // should not require much memory anyway...
        std::vector<int> satnumRegionArray(numCompressedElems);
        if (eclState.get3DProperties().hasDeckIntGridProperty("SATNUM")) {
            const auto& satnumRawData = eclState.get3DProperties().getIntGridProperty("SATNUM").getData();
            for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                unsigned cartesianElemIdx = static_cast<unsigned>(compressedToCartesianElemIdx[elemIdx]);
                satnumRegionArray[elemIdx] = satnumRawData[cartesianElemIdx] - 1;
            }
        }
        else
            std::fill(satnumRegionArray.begin(), satnumRegionArray.end(), 0);

        readGlobalEpsOptions_(deck, eclState);
        readGlobalHysteresisOptions_(deck);
        readGlobalThreePhaseOptions_(deck);

        unscaledEpsInfo_.resize(numSatRegions);
        for (unsigned satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx)
            unscaledEpsInfo_[satnumIdx].extractUnscaled(deck, eclState, satnumIdx);

        initParamsForElements_(deck, eclState, compressedToCartesianElemIdx, satnumRegionArray);
    }

    /*!
     * \brief Modify the initial condition according to the SWATINIT keyword.
     *
     * The method returns the water saturation which yields a givenn capillary
     * pressure. The reason this method is not folded directly into initFromDeck() is
     * that the capillary pressure given depends on the particuars of how the simulator
     * calculates its initial condition.
     */
    void applySwatinit(unsigned elemIdx,
                       Scalar pcow,
                       Scalar Sw)
    {
        // apply SWATINIT to the drainage parameters
        auto& elemScaledEpsInfoDrain = *oilWaterScaledEpsInfoDrainage_[elemIdx];
        auto& elemEclEpsScalingPointsDrain = oilWaterScaledEpsPointsDrainage(elemIdx);
        auto& matParamsDrain = oilWaterParams_[elemIdx]->drainageParams();

        applySwatinitToBranch_<OilWaterScaledTwoPhaseLaw>(pcow,
                                                          Sw,
                                                          matParamsDrain,
                                                          elemScaledEpsInfoDrain,
                                                          elemEclEpsScalingPointsDrain);

        if (enableHysteresis()) {
            // apply SWATINIT to the imbibition parameters
            auto& elemScaledEpsInfoImb = *oilWaterScaledEpsInfoImbibition_[elemIdx];
            auto& elemEclEpsScalingPointsImb = oilWaterScaledEpsPointsImbibition(elemIdx);
            auto& matParamsImb = oilWaterParams_[elemIdx]->imbibitionParams();

            applySwatinitToBranch_<OilWaterScaledTwoPhaseLaw>(pcow,
                                                              Sw,
                                                              matParamsImb,
                                                              elemScaledEpsInfoImb,
                                                              elemEclEpsScalingPointsImb);
        }
    }

protected:
    template <class MatLaw, class MatLawParams>
    void applySwatinitToBranch_(Scalar pcow,
                                Scalar Sw,
                                const MatLawParams& matLawParams,
                                Opm::EclEpsScalingPointsInfo<Scalar>& scaledEpsInfo,
                                EclEpsScalingPoints<Scalar>& scaledEpsScalingPoints)
    {
        // compute the oil-water capillary pressure at the water-oil contact
        Scalar pcowWoc = MatLaw::twoPhaseSatPcnw(matLawParams, 1.0);
        if (pcow <= pcowWoc) {
            // the cell is below the water-oil contact depth.
            Sw = 1.0;
            pcow = pcowWoc;
        }

        Scalar pcowAtSw = MatLaw::twoPhaseSatPcnw(matLawParams, Sw);
        if (std::abs(pcowAtSw) > 1e-5) {
            // TODO (?): maybe it would be better to add an offset to the capillary
            // pressure instead of scaling it.
            scaledEpsInfo.maxPcow *= pcow/pcowAtSw;
            scaledEpsScalingPoints.init(scaledEpsInfo, *oilWaterEclEpsConfig_, Opm::EclOilWaterSystem);
        }
    }

public:
    bool enableEndPointScaling() const
    { return enableEndPointScaling_; }

    bool enableHysteresis() const
    { return hysteresisConfig_->enableHysteresis(); }

    MaterialLawParams& materialLawParams(unsigned elemIdx)
    {
        assert(0 <= elemIdx && elemIdx <  materialLawParams_.size());
        return *materialLawParams_[elemIdx];
    }

    const MaterialLawParams& materialLawParams(unsigned elemIdx) const
    {
        assert(0 <= elemIdx && elemIdx <  materialLawParams_.size());
        return *materialLawParams_[elemIdx];
    }

    std::shared_ptr<MaterialLawParams>& materialLawParamsPointerReferenceHack(unsigned elemIdx)
    {
        assert(0 <= elemIdx && elemIdx <  materialLawParams_.size());
        return materialLawParams_[elemIdx];
    }

    template <class FluidState>
    void updateHysteresis(const FluidState& fluidState, unsigned elemIdx)
    {
        if (!enableHysteresis())
            return;

        auto threePhaseParams = materialLawParams_[elemIdx];
        MaterialLaw::updateHysteresis(*threePhaseParams, fluidState);
    }

    EclEpsScalingPoints<Scalar>& oilWaterScaledEpsPointsDrainage(unsigned elemIdx)
    {
        auto& materialParams = *materialLawParams_[elemIdx];
        switch (materialParams.approach()) {
        case EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone1Approach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone2Approach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclDefaultApproach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclTwoPhaseApproach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        default:
            OPM_THROW(std::logic_error, "Enum value for material approach unknown!");
        }
    }

    EclEpsScalingPoints<Scalar>& oilWaterScaledEpsPointsImbibition(unsigned elemIdx)
    {
        auto& materialParams = *materialLawParams_[elemIdx];
        switch (materialParams.approach()) {
        case EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone1Approach>();
            return realParams.oilWaterParams().imbibitionParams().scaledPoints();
        }

        case EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone2Approach>();
            return realParams.oilWaterParams().imbibitionParams().scaledPoints();
        }

        case EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclDefaultApproach>();
            return realParams.oilWaterParams().imbibitionParams().scaledPoints();
        }

        case EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclTwoPhaseApproach>();
            return realParams.oilWaterParams().imbibitionParams().scaledPoints();
        }

        default:
            OPM_THROW(std::logic_error, "Enum value for material approach unknown!");
        }
    }

    const Opm::EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoDrainage(size_t elemIdx) const
    {
        return *oilWaterScaledEpsInfoDrainage_[elemIdx];
    }

    std::shared_ptr<EclEpsScalingPointsInfo<Scalar> >& oilWaterScaledEpsInfoDrainagePointerReferenceHack(unsigned elemIdx)
    {
        return oilWaterScaledEpsInfoDrainage_[elemIdx];
    }

    const Opm::EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoImbibition(size_t elemIdx) const
    {
        return *oilWaterScaledEpsInfoImbibition_[elemIdx];
    }

private:
    void readGlobalEpsOptions_(const Opm::Deck& deck, const Opm::EclipseState& eclState)
    {
        oilWaterEclEpsConfig_ = std::make_shared<Opm::EclEpsConfig>();
        oilWaterEclEpsConfig_-> initFromDeck(deck, eclState, Opm::EclOilWaterSystem);

        enableEndPointScaling_ = deck.hasKeyword("ENDSCALE");
    }

    void readGlobalHysteresisOptions_(const Opm::Deck& deck)
    {
        hysteresisConfig_ = std::make_shared<Opm::EclHysteresisConfig>();
        hysteresisConfig_->initFromDeck(deck);
    }

    void readGlobalThreePhaseOptions_(const Opm::Deck& deck)
    {
        bool gasEnabled = deck.hasKeyword("GAS");
        bool oilEnabled = deck.hasKeyword("OIL");
        bool waterEnabled = deck.hasKeyword("WATER");

        int numEnabled =
            (gasEnabled?1:0)
            + (oilEnabled?1:0)
            + (waterEnabled?1:0);

        if (numEnabled < 2)
            OPM_THROW(std::runtime_error,
                      "At least two fluid phases must be enabled. (Is: " << numEnabled << ")");

        if (numEnabled == 2) {
            threePhaseApproach_ = Opm::EclTwoPhaseApproach;
            if (!gasEnabled)
                twoPhaseApproach_ = Opm::EclTwoPhaseOilWater;
            else if (!oilEnabled)
                twoPhaseApproach_ = Opm::EclTwoPhaseGasWater;
            else if (!waterEnabled)
                twoPhaseApproach_ = Opm::EclTwoPhaseGasOil;
        }
        else {
            assert(numEnabled == 3);

            threePhaseApproach_ = Opm::EclDefaultApproach;
            if (deck.hasKeyword("STONE") || deck.hasKeyword("STONE2"))
                threePhaseApproach_ = Opm::EclStone2Approach;
            else if (deck.hasKeyword("STONE1"))
                threePhaseApproach_ = Opm::EclStone1Approach;
        }
    }

    void initParamsForElements_(const Deck& deck, const EclipseState& eclState,
                                const std::vector<int>& compressedToCartesianElemIdx,
                                const std::vector<int>& satnumRegionArray)
    {
        const size_t numSatRegions = eclState.runspec().tabdims().getNumSatTables();
        unsigned numCompressedElems = static_cast<unsigned>(compressedToCartesianElemIdx.size());

        // read the end point scaling configuration. this needs to be done only once per
        // deck.
        auto gasOilConfig = std::make_shared<Opm::EclEpsConfig>();
        auto oilWaterConfig = std::make_shared<Opm::EclEpsConfig>();
        gasOilConfig->initFromDeck(deck, eclState, Opm::EclGasOilSystem);
        oilWaterConfig->initFromDeck(deck, eclState, Opm::EclOilWaterSystem);

        // read the saturation region specific parameters from the deck
        GasOilScalingPointsVector gasOilUnscaledPointsVector(numSatRegions);
        OilWaterScalingPointsVector oilWaterUnscaledPointsVector(numSatRegions);
        GasOilEffectiveParamVector gasOilEffectiveParamVector(numSatRegions);
        OilWaterEffectiveParamVector oilWaterEffectiveParamVector(numSatRegions);
        for (unsigned satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
            // unscaled points for end-point scaling
            readGasOilUnscaledPoints_(gasOilUnscaledPointsVector, gasOilConfig, deck, eclState, satnumIdx);
            readOilWaterUnscaledPoints_(oilWaterUnscaledPointsVector, oilWaterConfig, deck, eclState, satnumIdx);

            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector, deck, eclState, satnumIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector, deck, eclState, satnumIdx);

            // read the end point scaling info for the saturation region
            unscaledEpsInfo_[satnumIdx].extractUnscaled(deck, eclState, satnumIdx);

        }

        // read the scaled end point scaling parameters which are specific for each
        // element
        GasOilScalingInfoVector gasOilScaledInfoVector(numCompressedElems);
        oilWaterScaledEpsInfoDrainage_.resize(numCompressedElems);
        GasOilScalingInfoVector gasOilScaledImbInfoVector;

        GasOilScalingPointsVector gasOilScaledPointsVector(numCompressedElems);
        GasOilScalingPointsVector oilWaterScaledEpsPointsDrainage(numCompressedElems);
        GasOilScalingPointsVector gasOilScaledImbPointsVector;
        OilWaterScalingPointsVector oilWaterScaledImbPointsVector;

        if (enableHysteresis()) {
            gasOilScaledImbInfoVector.resize(numCompressedElems);
            gasOilScaledImbPointsVector.resize(numCompressedElems);
            oilWaterScaledEpsInfoImbibition_.resize(numCompressedElems);
            oilWaterScaledImbPointsVector.resize(numCompressedElems);
        }

        EclEpsGridProperties epsGridProperties, epsImbGridProperties;
        epsGridProperties.initFromDeck(deck, eclState, /*imbibition=*/false);
        if (enableHysteresis())
            epsImbGridProperties.initFromDeck(deck, eclState, /*imbibition=*/true);
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            unsigned cartElemIdx = static_cast<unsigned>(compressedToCartesianElemIdx[elemIdx]);
            readGasOilScaledPoints_(gasOilScaledInfoVector,
                                    gasOilScaledPointsVector,
                                    gasOilConfig,
                                    deck,
                                    eclState,
                                    epsGridProperties,
                                    elemIdx,
                                    cartElemIdx);
            readOilWaterScaledPoints_(oilWaterScaledEpsInfoDrainage_,
                                      oilWaterScaledEpsPointsDrainage,
                                      oilWaterConfig,
                                      deck,
                                      eclState,
                                      epsGridProperties,
                                      elemIdx,
                                      cartElemIdx);

            if (enableHysteresis()) {
                readGasOilScaledPoints_(gasOilScaledImbInfoVector,
                                        gasOilScaledImbPointsVector,
                                        gasOilConfig,
                                        deck,
                                        eclState,
                                        epsImbGridProperties,
                                        elemIdx,
                                        cartElemIdx);
                readOilWaterScaledPoints_(oilWaterScaledEpsInfoImbibition_,
                                          oilWaterScaledImbPointsVector,
                                          oilWaterConfig,
                                          deck,
                                          eclState,
                                          epsImbGridProperties,
                                          elemIdx,
                                          cartElemIdx);
            }
        }

        // create the parameter objects for the two-phase laws
        gasOilParams_.resize(numCompressedElems);
        oilWaterParams_.resize(numCompressedElems);

        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");
        bool hasWater = deck.hasKeyword("WATER");

        const auto& imbnumData = eclState.get3DProperties().getIntGridProperty("IMBNUM").getData();
        assert(numCompressedElems == satnumRegionArray.size());
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            unsigned satnumIdx = static_cast<unsigned>(satnumRegionArray[elemIdx]);

            gasOilParams_[elemIdx] = std::make_shared<GasOilTwoPhaseParams>();
            oilWaterParams_[elemIdx] = std::make_shared<OilWaterTwoPhaseParams>();

            gasOilParams_[elemIdx]->setConfig(hysteresisConfig_);
            oilWaterParams_[elemIdx]->setConfig(hysteresisConfig_);

            if (hasGas && hasOil) {
                auto gasOilDrainParams = std::make_shared<GasOilEpsTwoPhaseParams>();
                gasOilDrainParams->setConfig(gasOilConfig);
                gasOilDrainParams->setUnscaledPoints(gasOilUnscaledPointsVector[satnumIdx]);
                gasOilDrainParams->setScaledPoints(gasOilScaledPointsVector[elemIdx]);
                gasOilDrainParams->setEffectiveLawParams(gasOilEffectiveParamVector[satnumIdx]);
                gasOilDrainParams->finalize();

                gasOilParams_[elemIdx]->setDrainageParams(gasOilDrainParams,
                                                          *gasOilScaledInfoVector[elemIdx],
                                                          EclGasOilSystem);
            }

            if (hasOil && hasWater) {
                auto oilWaterDrainParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
                oilWaterDrainParams->setConfig(oilWaterConfig);
                oilWaterDrainParams->setUnscaledPoints(oilWaterUnscaledPointsVector[satnumIdx]);
                oilWaterDrainParams->setScaledPoints(oilWaterScaledEpsPointsDrainage[elemIdx]);
                oilWaterDrainParams->setEffectiveLawParams(oilWaterEffectiveParamVector[satnumIdx]);
                oilWaterDrainParams->finalize();

                oilWaterParams_[elemIdx]->setDrainageParams(oilWaterDrainParams,
                                                            *oilWaterScaledEpsInfoDrainage_[elemIdx],
                                                            EclOilWaterSystem);
            }

            if (enableHysteresis()) {
                unsigned imbRegionIdx = static_cast<unsigned>(imbnumData[elemIdx]) - 1;

                if (hasGas && hasOil) {
                    auto gasOilImbParams = std::make_shared<GasOilEpsTwoPhaseParams>();
                    gasOilImbParams->setConfig(gasOilConfig);
                    gasOilImbParams->setUnscaledPoints(gasOilUnscaledPointsVector[imbRegionIdx]);
                    gasOilImbParams->setScaledPoints(gasOilScaledImbPointsVector[elemIdx]);
                    gasOilImbParams->setEffectiveLawParams(gasOilEffectiveParamVector[imbRegionIdx]);
                    gasOilImbParams->finalize();

                    gasOilParams_[elemIdx]->setImbibitionParams(gasOilImbParams,
                                                                *gasOilScaledImbInfoVector[elemIdx],
                                                                EclGasOilSystem);
                }

                if (hasOil && hasWater) {
                    auto oilWaterImbParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
                    oilWaterImbParams->setConfig(oilWaterConfig);
                    oilWaterImbParams->setUnscaledPoints(oilWaterUnscaledPointsVector[imbRegionIdx]);
                    oilWaterImbParams->setScaledPoints(oilWaterScaledImbPointsVector[elemIdx]);
                    oilWaterImbParams->setEffectiveLawParams(oilWaterEffectiveParamVector[imbRegionIdx]);
                    oilWaterImbParams->finalize();

                    oilWaterParams_[elemIdx]->setImbibitionParams(oilWaterImbParams,
                                                                  *gasOilScaledImbInfoVector[elemIdx],
                                                                  EclGasOilSystem);
                }
            }

            if (hasGas && hasOil)
                gasOilParams_[elemIdx]->finalize();

            if (hasOil && hasWater)
                oilWaterParams_[elemIdx]->finalize();
        }

        // create the parameter objects for the three-phase law
        materialLawParams_.resize(numCompressedElems);
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            materialLawParams_[elemIdx] = std::make_shared<MaterialLawParams>();
            unsigned satnumIdx = static_cast<unsigned>(satnumRegionArray[elemIdx]);

            initThreePhaseParams_(deck,
                                  eclState,
                                  *materialLawParams_[elemIdx],
                                  satnumIdx,
                                  *oilWaterScaledEpsInfoDrainage_[elemIdx],
                                  oilWaterParams_[elemIdx],
                                  gasOilParams_[elemIdx]);

            materialLawParams_[elemIdx]->finalize();
        }
    }

    // The saturation function family.
    // If SWOF and SGOF are specified in the deck it return FamilyI
    // If SWFN, SGFN and SOF3 are specified in the deck it return FamilyII
    // If keywords are missing or mixed, an error is given.
    enum SaturationFunctionFamily {
        noFamily,
        FamilyI,
        FamilyII
    };

    SaturationFunctionFamily getSaturationFunctionFamily(const Opm::Deck& deck, const Opm::EclipseState& eclState) const
    {
        const auto& tableManager = eclState.getTableManager();
        const TableContainer& swofTables = tableManager.getSwofTables();
        const TableContainer& slgofTables= tableManager.getSlgofTables();
        const TableContainer& sgofTables = tableManager.getSgofTables();
        const TableContainer& swfnTables = tableManager.getSwfnTables();
        const TableContainer& sgfnTables = tableManager.getSgfnTables();
        const TableContainer& sof3Tables = tableManager.getSof3Tables();

        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");
        bool hasWater = deck.hasKeyword("WATER");

        bool family1 = false;
        bool family2 = false;
        if (!hasGas) {
            // oil-water case
            family1 = !swofTables.empty();
            family2 = !swfnTables.empty() && !sof3Tables.empty();
        }
        else if (!hasWater) {
            // oil-gas case
            throw std::runtime_error("oil-gas two-phase simulations are currently not supported");
        }
        else if (!hasOil) {
            // water-gas case
            throw std::runtime_error("water-gas two-phase simulations are currently not supported");
        }
        else {
            // three-phase case
            family1 = (!sgofTables.empty() || !slgofTables.empty()) && !swofTables.empty();
            family2 = !swfnTables.empty() && !sgfnTables.empty() && !sof3Tables.empty();
        }

        if (family1 && family2) {
            throw std::invalid_argument("Saturation families should not be mixed \n"
                                        "Use either SGOF and SWOF or SGFN, SWFN and SOF3");
        }

        if (!family1 && !family2) {
            throw std::invalid_argument("Saturations function must be specified using either "
                                        "family 1 or family 2 keywords \n"
                                        "Use either SGOF and SWOF or SGFN, SWFN and SOF3" );
        }

        if (family1 && !family2)
            return SaturationFunctionFamily::FamilyI;
        else if (family2 && !family1)
            return SaturationFunctionFamily::FamilyII;
        return SaturationFunctionFamily::noFamily; // no family or two families
    }

    template <class Container>
    void readGasOilEffectiveParameters_(Container& dest,
                                        const Opm::Deck& deck,
                                        const Opm::EclipseState& eclState,
                                        unsigned satnumIdx)
    {
        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasGas || !hasOil)
            // we don't read anything if either the gas or the oil phase is not active
            return;

        dest[satnumIdx] = std::make_shared<GasOilEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satnumIdx];

        // the situation for the gas phase is complicated that all saturations are
        // shifted by the connate water saturation.
        Scalar Swco = unscaledEpsInfo_[satnumIdx].Swl;
        const auto& tableManager = eclState.getTableManager();

        switch (getSaturationFunctionFamily(deck, eclState)) {
        case FamilyI:
        {
            const TableContainer& sgofTables = tableManager.getSgofTables();
            const TableContainer& slgofTables = tableManager.getSlgofTables();
            if (!sgofTables.empty())
                readGasOilEffectiveParametersSgof_(effParams,
                                                   Swco,
                                                   sgofTables.getTable<SgofTable>(satnumIdx));
            else if (!slgofTables.empty())
                readGasOilEffectiveParametersSlgof_(effParams,
                                                    Swco,
                                                    slgofTables.getTable<SlgofTable>(satnumIdx));
            break;
        }

        case FamilyII:
        {
            const Sof3Table& sof3Table = tableManager.getSof3Tables().getTable<Sof3Table>( satnumIdx );
            const SgfnTable& sgfnTable = tableManager.getSgfnTables().getTable<SgfnTable>( satnumIdx );
            readGasOilEffectiveParametersFamily2_(effParams,
                                                  Swco,
                                                  sof3Table,
                                                  sgfnTable);
            break;
        }

        //default:
        case noFamily:
            throw std::domain_error("No valid saturation keyword family specified");

        }
    }

    void readGasOilEffectiveParametersSgof_(GasOilEffectiveTwoPhaseParams& effParams,
                                            Scalar Swco,
                                            const Opm::SgofTable& sgofTable)
    {
        // convert the saturations of the SGOF keyword from gas to oil saturations
        std::vector<double> SoSamples(sgofTable.numRows());
        std::vector<double> SoKroSamples(sgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < sgofTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = 1 - sgofTable.get("SG", sampleIdx);
            SoKroSamples[sampleIdx] = SoSamples[sampleIdx] - Swco;
        }

        effParams.setKrwSamples(SoKroSamples, sgofTable.getColumn("KROG").vectorCopy());
        effParams.setKrnSamples(SoSamples, sgofTable.getColumn("KRG").vectorCopy());
        effParams.setPcnwSamples(SoSamples, sgofTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersSlgof_(GasOilEffectiveTwoPhaseParams& effParams,
                                             Scalar Swco,
                                             const Opm::SlgofTable& slgofTable)
    {
        // convert the saturations of the SLGOF keyword from "liquid" to oil saturations
        std::vector<double> SoSamples(slgofTable.numRows());
        std::vector<double> SoKroSamples(slgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < slgofTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = slgofTable.get("SL", sampleIdx);
            SoKroSamples[sampleIdx] = slgofTable.get("SL", sampleIdx) - Swco;
        }

        effParams.setKrwSamples(SoKroSamples, slgofTable.getColumn("KROG").vectorCopy());
        effParams.setKrnSamples(SoSamples, slgofTable.getColumn("KRG").vectorCopy());
        effParams.setPcnwSamples(SoSamples, slgofTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               Scalar /* Swco */,
                                               const Opm::Sof3Table& sof3Table,
                                               const Opm::SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        std::vector<double> SoColumn = sof3Table.getColumn("SO").vectorCopy();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = 1 - sgfnTable.get("SG", sampleIdx);
        }

        effParams.setKrwSamples(SoColumn, sof3Table.getColumn("KROG").vectorCopy());
        effParams.setKrnSamples(SoSamples, sgfnTable.getColumn("KRG").vectorCopy());
        effParams.setPcnwSamples(SoSamples, sgfnTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    template <class Container>
    void readOilWaterEffectiveParameters_(Container& dest,
                                          const Opm::Deck& deck,
                                          const Opm::EclipseState& eclState,
                                          unsigned satnumIdx)
    {
        bool hasWater = deck.hasKeyword("WATER");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasOil || !hasWater)
            // we don't read anything if either the water or the oil phase is not active
            return;

        dest[satnumIdx] = std::make_shared<OilWaterEffectiveTwoPhaseParams>();

        const auto& tableManager = eclState.getTableManager();
        auto& effParams = *dest[satnumIdx];

        switch (getSaturationFunctionFamily(deck, eclState)) {
        case FamilyI: {
            const auto& swofTable = tableManager.getSwofTables().getTable<SwofTable>(satnumIdx);
            std::vector<double> SwColumn = swofTable.getColumn("SW").vectorCopy();

            effParams.setKrwSamples(SwColumn, swofTable.getColumn("KRW").vectorCopy());
            effParams.setKrnSamples(SwColumn, swofTable.getColumn("KROW").vectorCopy());
            effParams.setPcnwSamples(SwColumn, swofTable.getColumn("PCOW").vectorCopy());
            effParams.finalize();
            break;
        }
        case FamilyII:
        {
            const auto& swfnTable = tableManager.getSwfnTables().getTable<SwfnTable>(satnumIdx);
            const auto& sof3Table = tableManager.getSof3Tables().getTable<Sof3Table>(satnumIdx);
            std::vector<double> SwColumn = swfnTable.getColumn("SW").vectorCopy();

            // convert the saturations of the SOF3 keyword from oil to water saturations
            std::vector<double> SwSamples(sof3Table.numRows());
            for (size_t sampleIdx = 0; sampleIdx < sof3Table.numRows(); ++ sampleIdx)
                SwSamples[sampleIdx] = 1 - sof3Table.get("SO", sampleIdx);

            effParams.setKrwSamples(SwColumn, swfnTable.getColumn("KRW").vectorCopy());
            effParams.setKrnSamples(SwSamples, sof3Table.getColumn("KROW").vectorCopy());
            effParams.setPcnwSamples(SwColumn, swfnTable.getColumn("PCOW").vectorCopy());
            effParams.finalize();
            break;
        }

        case noFamily:
        //default:
            throw std::domain_error("No valid saturation keyword family specified");

        }
    }


    template <class Container>
    void readGasOilUnscaledPoints_(Container &dest,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const Opm::Deck& deck,
                                   const Opm::EclipseState& /* eclState */,
                                   unsigned satnumIdx)
    {
        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasGas || !hasOil)
            // we don't read anything if either the gas or the oil phase is not active
            return;

        dest[satnumIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satnumIdx]->init(unscaledEpsInfo_[satnumIdx], *config, EclGasOilSystem);
    }

    template <class Container>
    void readOilWaterUnscaledPoints_(Container &dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     const Opm::Deck& deck,
                                     const Opm::EclipseState& /* eclState */,
                                     unsigned satnumIdx)
    {
        bool hasWater = deck.hasKeyword("WATER");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasOil || !hasWater)
            // we don't read anything if either the water or the oil phase is not active
            return;

        dest[satnumIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satnumIdx]->init(unscaledEpsInfo_[satnumIdx], *config, EclOilWaterSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readGasOilScaledPoints_(InfoContainer& destInfo,
                                 PointsContainer& destPoints,
                                 std::shared_ptr<EclEpsConfig> config,
                                 const Opm::Deck& deck,
                                 const Opm::EclipseState& eclState,
                                 const EclEpsGridProperties& epsGridProperties,
                                 unsigned elemIdx,
                                 unsigned cartElemIdx)
    {
        unsigned satnumIdx = static_cast<unsigned>((*epsGridProperties.satnum)[cartElemIdx]) - 1; // ECL uses Fortran indices!

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satnumIdx]);
        destInfo[elemIdx]->extractScaled(deck, eclState, epsGridProperties, cartElemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclGasOilSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readOilWaterScaledPoints_(InfoContainer& destInfo,
                                   PointsContainer& destPoints,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const Opm::Deck& deck,
                                   const Opm::EclipseState& eclState,
                                   const EclEpsGridProperties& epsGridProperties,
                                   unsigned elemIdx,
                                   unsigned cartElemIdx)
    {
        unsigned satnumIdx = static_cast<unsigned>((*epsGridProperties.satnum)[cartElemIdx]) - 1; // ECL uses Fortran indices!

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satnumIdx]);
        destInfo[elemIdx]->extractScaled(deck, eclState, epsGridProperties, cartElemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclOilWaterSystem);
    }

    void initThreePhaseParams_(const Opm::Deck& deck,
                               const Opm::EclipseState& /* eclState */,
                               MaterialLawParams& materialParams,
                               unsigned satnumIdx,
                               const EclEpsScalingPointsInfo<Scalar>& epsInfo,
                               std::shared_ptr<OilWaterTwoPhaseParams> oilWaterParams,
                               std::shared_ptr<GasOilTwoPhaseParams> gasOilParams)
    {
        materialParams.setApproach(threePhaseApproach_);

        switch (materialParams.approach()) {
        case EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone1Approach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);

            if (deck.hasKeyword("STONE1EX")) {
                Scalar eta =
                    deck.getKeyword("STONE1EX").getRecord(satnumIdx).getItem(0).getSIDouble(0);
                realParams.setEta(eta);
            }
            else
                realParams.setEta(1.0);
            realParams.finalize();
            break;
        }

        case EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone2Approach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);
            realParams.finalize();
            break;
        }

        case EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclDefaultApproach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);
            realParams.finalize();
            break;
        }

        case EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclTwoPhaseApproach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setApproach(twoPhaseApproach_);
            realParams.finalize();
            break;
        }
        }
    }

    bool enableEndPointScaling_;
    std::shared_ptr<EclHysteresisConfig> hysteresisConfig_;

    std::shared_ptr<EclEpsConfig> oilWaterEclEpsConfig_;
    std::vector<Opm::EclEpsScalingPointsInfo<Scalar>> unscaledEpsInfo_;
    OilWaterScalingInfoVector oilWaterScaledEpsInfoDrainage_;
    OilWaterScalingInfoVector oilWaterScaledEpsInfoImbibition_;

    OilWaterParamVector oilWaterParams_;
    GasOilParamVector gasOilParams_;

    Opm::EclMultiplexerApproach threePhaseApproach_;

    // this attribute only makes sense for twophase simulations!
    enum EclTwoPhaseApproach twoPhaseApproach_;

    std::vector<std::shared_ptr<MaterialLawParams> > materialLawParams_;
};
} // namespace Opm

#endif
