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
#if ! HAVE_ECL_INPUT
#error "Eclipse input support in opm-common is required to use the ECL material manager!"
#endif

#ifndef OPM_ECL_MATERIAL_LAW_MANAGER_HPP
#define OPM_ECL_MATERIAL_LAW_MANAGER_HPP

#include <opm/material/fluidmatrixinteractions/SatCurveMultiplexer.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

#if HAVE_OPM_COMMON
#include <opm/common/OpmLog/OpmLog.hpp>
#endif

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableColumn.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>

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
    using Traits = TraitsT;
    using Scalar = typename Traits::Scalar;
    enum { waterPhaseIdx = Traits::wettingPhaseIdx };
    enum { oilPhaseIdx = Traits::nonWettingPhaseIdx };
    enum { gasPhaseIdx = Traits::gasPhaseIdx };
    enum { numPhases = Traits::numPhases };

    using GasOilTraits = TwoPhaseMaterialTraits<Scalar, oilPhaseIdx, gasPhaseIdx>;
    using OilWaterTraits = TwoPhaseMaterialTraits<Scalar, waterPhaseIdx, oilPhaseIdx>;
    using GasWaterTraits = TwoPhaseMaterialTraits<Scalar, waterPhaseIdx, gasPhaseIdx>;

    // the two-phase material law which is defined on effective (unscaled) saturations
    using GasOilEffectiveTwoPhaseLaw = SatCurveMultiplexer<GasOilTraits>;
    using OilWaterEffectiveTwoPhaseLaw = SatCurveMultiplexer<OilWaterTraits>;
    using GasWaterEffectiveTwoPhaseLaw = SatCurveMultiplexer<GasWaterTraits>;

    using GasOilEffectiveTwoPhaseParams = typename GasOilEffectiveTwoPhaseLaw::Params;
    using OilWaterEffectiveTwoPhaseParams = typename OilWaterEffectiveTwoPhaseLaw::Params;
    using GasWaterEffectiveTwoPhaseParams = typename GasWaterEffectiveTwoPhaseLaw::Params;

    // the two-phase material law which is defined on absolute (scaled) saturations
    using GasOilEpsTwoPhaseLaw = EclEpsTwoPhaseLaw<GasOilEffectiveTwoPhaseLaw>;
    using OilWaterEpsTwoPhaseLaw = EclEpsTwoPhaseLaw<OilWaterEffectiveTwoPhaseLaw>;
    using GasWaterEpsTwoPhaseLaw = EclEpsTwoPhaseLaw<GasWaterEffectiveTwoPhaseLaw>;
    using GasOilEpsTwoPhaseParams = typename GasOilEpsTwoPhaseLaw::Params;
    using OilWaterEpsTwoPhaseParams = typename OilWaterEpsTwoPhaseLaw::Params;
    using GasWaterEpsTwoPhaseParams = typename GasWaterEpsTwoPhaseLaw::Params;

    // the scaled two-phase material laws with hystersis
    using GasOilTwoPhaseLaw = EclHysteresisTwoPhaseLaw<GasOilEpsTwoPhaseLaw>;
    using OilWaterTwoPhaseLaw = EclHysteresisTwoPhaseLaw<OilWaterEpsTwoPhaseLaw>;
    using GasWaterTwoPhaseLaw = EclHysteresisTwoPhaseLaw<GasWaterEpsTwoPhaseLaw>;
    using GasOilTwoPhaseHystParams = typename GasOilTwoPhaseLaw::Params;
    using OilWaterTwoPhaseHystParams = typename OilWaterTwoPhaseLaw::Params;
    using GasWaterTwoPhaseHystParams = typename GasWaterTwoPhaseLaw::Params;

public:
    // the three-phase material law used by the simulation
    using MaterialLaw = EclMultiplexerMaterial<Traits, GasOilTwoPhaseLaw, OilWaterTwoPhaseLaw, GasWaterTwoPhaseLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;

private:
    // internal typedefs
    using GasOilEffectiveParamVector = std::vector<std::shared_ptr<GasOilEffectiveTwoPhaseParams>>;
    using OilWaterEffectiveParamVector = std::vector<std::shared_ptr<OilWaterEffectiveTwoPhaseParams>>;
    using GasWaterEffectiveParamVector = std::vector<std::shared_ptr<GasWaterEffectiveTwoPhaseParams>>;

    using GasOilScalingPointsVector = std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar>>>;
    using OilWaterScalingPointsVector = std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar>>>;
    using GasWaterScalingPointsVector = std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar>>>;
    using OilWaterScalingInfoVector = std::vector<EclEpsScalingPointsInfo<Scalar>>;
    using GasOilParamVector = std::vector<std::shared_ptr<GasOilTwoPhaseHystParams>>;
    using OilWaterParamVector = std::vector<std::shared_ptr<OilWaterTwoPhaseHystParams>>;
    using GasWaterParamVector = std::vector<std::shared_ptr<GasWaterTwoPhaseHystParams>>;
    using MaterialLawParamsVector = std::vector<std::shared_ptr<MaterialLawParams>>;

public:
    EclMaterialLawManager()
    {}

    void initFromState(const EclipseState& eclState)
    {
        // get the number of saturation regions and the number of cells in the deck
        const auto&  runspec       = eclState.runspec();
        const size_t numSatRegions = runspec.tabdims().getNumSatTables();

        const auto& ph = runspec.phases();
        this->hasGas = ph.active(Phase::GAS);
        this->hasOil = ph.active(Phase::OIL);
        this->hasWater = ph.active(Phase::WATER);

        readGlobalEpsOptions_(eclState);
        readGlobalHysteresisOptions_(eclState);
        readGlobalThreePhaseOptions_(runspec);

        // Read the end point scaling configuration (once per run).
        gasOilConfig = std::make_shared<EclEpsConfig>();
        oilWaterConfig = std::make_shared<EclEpsConfig>();
        gasWaterConfig = std::make_shared<EclEpsConfig>();
        gasOilConfig->initFromState(eclState, EclGasOilSystem);
        oilWaterConfig->initFromState(eclState, EclOilWaterSystem);
        gasWaterConfig->initFromState(eclState, EclGasWaterSystem);


        const auto& tables = eclState.getTableManager();

        {
            const auto& stone1exTables = tables.getStone1exTable();

            if (! stone1exTables.empty()) {
                stoneEtas.clear();
                stoneEtas.reserve(numSatRegions);

                for (const auto& table : stone1exTables) {
                    stoneEtas.push_back(table.eta);
                }
            }
        }

        this->unscaledEpsInfo_.resize(numSatRegions);

        if (this->hasGas + this->hasOil + this->hasWater == 1) {
            // Single-phase simulation.  Special case.  Nothing to do here.
            return;
        }

        // Multiphase simulation.  Common case.
        const auto tolcrit = runspec.saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto rtep  = satfunc::getRawTableEndpoints(tables, ph, tolcrit);
        const auto rfunc = satfunc::getRawFunctionValues(tables, ph, rtep);

        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            this->unscaledEpsInfo_[satRegionIdx]
                .extractUnscaled(rtep, rfunc, satRegionIdx);
        }
    }

    void initParamsForElements(const EclipseState& eclState, size_t numCompressedElems)
    {
        // get the number of saturation regions
        const size_t numSatRegions = eclState.runspec().tabdims().getNumSatTables();

        // setup the saturation region specific parameters
        gasOilUnscaledPointsVector_.resize(numSatRegions);
        oilWaterUnscaledPointsVector_.resize(numSatRegions);
        gasWaterUnscaledPointsVector_.resize(numSatRegions);

        gasOilEffectiveParamVector_.resize(numSatRegions);
        oilWaterEffectiveParamVector_.resize(numSatRegions);
        gasWaterEffectiveParamVector_.resize(numSatRegions);
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            // unscaled points for end-point scaling
            readGasOilUnscaledPoints_(gasOilUnscaledPointsVector_, gasOilConfig, eclState, satRegionIdx);
            readOilWaterUnscaledPoints_(oilWaterUnscaledPointsVector_, oilWaterConfig, eclState, satRegionIdx);
            readGasWaterUnscaledPoints_(gasWaterUnscaledPointsVector_, gasWaterConfig, eclState, satRegionIdx);

            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector_, eclState, satRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector_, eclState, satRegionIdx);
            readGasWaterEffectiveParameters_(gasWaterEffectiveParamVector_, eclState, satRegionIdx);
        }

        // copy the SATNUM grid property. in some cases this is not necessary, but it
        // should not require much memory anyway...
        satnumRegionArray_.resize(numCompressedElems);
        if (eclState.fieldProps().has_int("SATNUM")) {
            const auto& satnumRawData = eclState.fieldProps().get_int("SATNUM");
            for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                satnumRegionArray_[elemIdx] = satnumRawData[elemIdx] - 1;
            }
        }
        else {
            std::fill(satnumRegionArray_.begin(), satnumRegionArray_.end(), 0);
        }
        auto copy_krnum = [&eclState, numCompressedElems](std::vector<int>& dest, const std::string keyword) {
            if (eclState.fieldProps().has_int(keyword)) {
                dest.resize(numCompressedElems);
                const auto& satnumRawData = eclState.fieldProps().get_int(keyword);
                for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                    dest[elemIdx] = satnumRawData[elemIdx] - 1;
                }
            }
        };
        copy_krnum(krnumXArray_, "KRNUMX");
        copy_krnum(krnumYArray_, "KRNUMY");
        copy_krnum(krnumZArray_, "KRNUMZ");

        // create the information for the imbibition region (IMBNUM). By default this is
        // the same as the saturation region (SATNUM)
        imbnumRegionArray_ = satnumRegionArray_;
        if (eclState.fieldProps().has_int("IMBNUM")) {
            const auto& imbnumRawData = eclState.fieldProps().get_int("IMBNUM");
            for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                imbnumRegionArray_[elemIdx] = imbnumRawData[elemIdx] - 1;
            }
        }

        assert(numCompressedElems == satnumRegionArray_.size());
        assert(!enableHysteresis() || numCompressedElems == imbnumRegionArray_.size());

        // read the scaled end point scaling parameters which are specific for each
        // element
        oilWaterScaledEpsInfoDrainage_.resize(numCompressedElems);

        std::unique_ptr<EclEpsGridProperties> epsImbGridProperties;

        if (enableHysteresis()) {
            epsImbGridProperties = std::make_unique<EclEpsGridProperties>(eclState, true);
        }

        EclEpsGridProperties epsGridProperties(eclState, false);
        materialLawParams_.resize(numCompressedElems);

        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            unsigned satRegionIdx = static_cast<unsigned>(satnumRegionArray_[elemIdx]);
            auto gasOilParams = std::make_shared<GasOilTwoPhaseHystParams>();
            auto oilWaterParams = std::make_shared<OilWaterTwoPhaseHystParams>();
            auto gasWaterParams = std::make_shared<GasWaterTwoPhaseHystParams>();
            gasOilParams->setConfig(hysteresisConfig_);
            oilWaterParams->setConfig(hysteresisConfig_);
            gasWaterParams->setConfig(hysteresisConfig_);

            auto [gasOilScaledInfo, gasOilScaledPoint] =
                readScaledPoints_(*gasOilConfig,
                                  eclState,
                                  epsGridProperties,
                                  elemIdx,
                                  EclGasOilSystem);

            auto [owinfo, oilWaterScaledEpsPointDrainage] =
                readScaledPoints_(*oilWaterConfig,
                                  eclState,
                                  epsGridProperties,
                                  elemIdx,
                                  EclOilWaterSystem);
            oilWaterScaledEpsInfoDrainage_[elemIdx] = owinfo;

            auto [gasWaterScaledInfo, gasWaterScaledPoint] =
                readScaledPoints_(*gasWaterConfig,
                                  eclState,
                                  epsGridProperties,
                                  elemIdx,
                                  EclGasWaterSystem);

            if (hasGas && hasOil) {
                GasOilEpsTwoPhaseParams gasOilDrainParams;
                gasOilDrainParams.setConfig(gasOilConfig);
                gasOilDrainParams.setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
                gasOilDrainParams.setScaledPoints(gasOilScaledPoint);
                gasOilDrainParams.setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
                gasOilDrainParams.finalize();

                gasOilParams->setDrainageParams(gasOilDrainParams,
                                                gasOilScaledInfo,
                                                EclGasOilSystem);
            }

            if (hasOil && hasWater) {
                OilWaterEpsTwoPhaseParams oilWaterDrainParams;
                oilWaterDrainParams.setConfig(oilWaterConfig);
                oilWaterDrainParams.setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
                oilWaterDrainParams.setScaledPoints(oilWaterScaledEpsPointDrainage);
                oilWaterDrainParams.setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
                oilWaterDrainParams.finalize();

                oilWaterParams->setDrainageParams(oilWaterDrainParams,
                                                  owinfo,
                                                  EclOilWaterSystem);
            }

            if (hasGas && hasWater && !hasOil) {
                GasWaterEpsTwoPhaseParams gasWaterDrainParams;
                gasWaterDrainParams.setConfig(gasWaterConfig);
                gasWaterDrainParams.setUnscaledPoints(gasWaterUnscaledPointsVector_[satRegionIdx]);
                gasWaterDrainParams.setScaledPoints(gasWaterScaledPoint);
                gasWaterDrainParams.setEffectiveLawParams(gasWaterEffectiveParamVector_[satRegionIdx]);
                gasWaterDrainParams.finalize();

                gasWaterParams->setDrainageParams(gasWaterDrainParams,
                                                  gasWaterScaledInfo,
                                                  EclGasWaterSystem);
            }

            if (enableHysteresis()) {
                auto [gasOilScaledImbInfo, gasOilScaledImbPoint] =
                    readScaledPoints_(*gasOilConfig,
                                      eclState,
                                      *epsImbGridProperties,
                                      elemIdx,
                                      EclGasOilSystem);

                auto [oilWaterScaledImbInfo, oilWaterScaledImbPoint] =
                    readScaledPoints_(*oilWaterConfig,
                                      eclState,
                                      *epsImbGridProperties,
                                      elemIdx,
                                      EclOilWaterSystem);

                auto [gasWaterScaledImbInfo, gasWaterScaledImbPoint] =
                    readScaledPoints_(*gasWaterConfig,
                                      eclState,
                                      *epsImbGridProperties,
                                      elemIdx,
                                      EclGasWaterSystem);

                unsigned imbRegionIdx = imbnumRegionArray_[elemIdx];
                if (hasGas && hasOil) {
                    GasOilEpsTwoPhaseParams gasOilImbParamsHyst;
                    gasOilImbParamsHyst.setConfig(gasOilConfig);
                    gasOilImbParamsHyst.setUnscaledPoints(gasOilUnscaledPointsVector_[imbRegionIdx]);
                    gasOilImbParamsHyst.setScaledPoints(gasOilScaledImbPoint);
                    gasOilImbParamsHyst.setEffectiveLawParams(gasOilEffectiveParamVector_[imbRegionIdx]);
                    gasOilImbParamsHyst.finalize();

                    gasOilParams->setImbibitionParams(gasOilImbParamsHyst,
                                                      gasOilScaledImbInfo,
                                                      EclGasOilSystem);
                }

                if (hasOil && hasWater) {
                    OilWaterEpsTwoPhaseParams oilWaterImbParamsHyst;
                    oilWaterImbParamsHyst.setConfig(oilWaterConfig);
                    oilWaterImbParamsHyst.setUnscaledPoints(oilWaterUnscaledPointsVector_[imbRegionIdx]);
                    oilWaterImbParamsHyst.setScaledPoints(oilWaterScaledImbPoint);
                    oilWaterImbParamsHyst.setEffectiveLawParams(oilWaterEffectiveParamVector_[imbRegionIdx]);
                    oilWaterImbParamsHyst.finalize();

                    oilWaterParams->setImbibitionParams(oilWaterImbParamsHyst,
                                                        oilWaterScaledImbInfo,
                                                        EclOilWaterSystem);
                }

                if (hasGas && hasWater && !hasOil) {
                    GasWaterEpsTwoPhaseParams gasWaterImbParamsHyst;
                    gasWaterImbParamsHyst.setConfig(gasWaterConfig);
                    gasWaterImbParamsHyst.setUnscaledPoints(gasWaterUnscaledPointsVector_[imbRegionIdx]);
                    gasWaterImbParamsHyst.setScaledPoints(gasWaterScaledImbPoint);
                    gasWaterImbParamsHyst.setEffectiveLawParams(gasWaterEffectiveParamVector_[imbRegionIdx]);
                    gasWaterImbParamsHyst.finalize();

                    gasWaterParams->setImbibitionParams(gasWaterImbParamsHyst,
                                                        gasWaterScaledImbInfo,
                                                        EclGasWaterSystem);
                }
            }

            if (hasGas && hasOil)
                gasOilParams->finalize();

            if (hasOil && hasWater)
                oilWaterParams->finalize();

            if (hasGas && hasWater && !hasOil)
                gasWaterParams->finalize();

            initThreePhaseParams_(eclState,
                                  materialLawParams_[elemIdx],
                                  satRegionIdx,
                                  oilWaterScaledEpsInfoDrainage_[elemIdx],
                                  oilWaterParams,
                                  gasOilParams,
                                  gasWaterParams);

            materialLawParams_[elemIdx].finalize();
        }
    }


    /*!
     * \brief Modify the initial condition according to the SWATINIT keyword.
     *
     * The method returns the water saturation which yields a givenn capillary
     * pressure. The reason this method is not folded directly into initFromState() is
     * that the capillary pressure given depends on the particuars of how the simulator
     * calculates its initial condition.
     */
    Scalar applySwatinit(unsigned elemIdx,
                         Scalar pcow,
                         Scalar Sw)
    {
        auto& elemScaledEpsInfo = oilWaterScaledEpsInfoDrainage_[elemIdx];

        // TODO: Mixed wettability systems - see ecl kw OPTIONS switch 74

        if (pcow < 0.0)
            Sw = elemScaledEpsInfo.Swu;
        else {

            if (Sw <= elemScaledEpsInfo.Swl)
                Sw = elemScaledEpsInfo.Swl;

            // specify a fluid state which only stores the saturations
            using FluidState = SimpleModularFluidState<Scalar,
                                                       numPhases,
                                                       /*numComponents=*/0,
                                                       /*FluidSystem=*/void, /* -> don't care */
                                                       /*storePressure=*/false,
                                                       /*storeTemperature=*/false,
                                                       /*storeComposition=*/false,
                                                       /*storeFugacity=*/false,
                                                       /*storeSaturation=*/true,
                                                       /*storeDensity=*/false,
                                                       /*storeViscosity=*/false,
                                                       /*storeEnthalpy=*/false>;
            FluidState fs;
            fs.setSaturation(waterPhaseIdx, Sw);
            fs.setSaturation(gasPhaseIdx, 0);
            fs.setSaturation(oilPhaseIdx, 0);
            std::array<Scalar, numPhases> pc = { 0 };
            MaterialLaw::capillaryPressures(pc, materialLawParams(elemIdx), fs);

            Scalar pcowAtSw = pc[oilPhaseIdx] - pc[waterPhaseIdx];
            constexpr const Scalar pcowAtSwThreshold = 1.0; //Pascal
            // avoid divison by very small number
            if (std::abs(pcowAtSw) > pcowAtSwThreshold) {
                elemScaledEpsInfo.maxPcow *= pcow/pcowAtSw;
                auto& elemEclEpsScalingPoints = oilWaterScaledEpsPointsDrainage(elemIdx);
                elemEclEpsScalingPoints.init(elemScaledEpsInfo, *oilWaterEclEpsConfig_, EclOilWaterSystem);
            }
        }

        return Sw;
    }

    bool enableEndPointScaling() const
    { return enableEndPointScaling_; }

    bool enableHysteresis() const
    { return hysteresisConfig_->enableHysteresis(); }

    MaterialLawParams& materialLawParams(unsigned elemIdx)
    {
        assert(elemIdx <  materialLawParams_.size());
        return materialLawParams_[elemIdx];
    }

    const MaterialLawParams& materialLawParams(unsigned elemIdx) const
    {
        assert(elemIdx <  materialLawParams_.size());
        return materialLawParams_[elemIdx];
    }

    /*!
     * \brief Returns a material parameter object for a given element and saturation region.
     *
     * This method changes the saturation table idx in the original material law parameter object.
     * In the context of ECL reservoir simulators, this is required to properly handle
     * wells with its own saturation table idx. In order to reset the saturation table idx
     * in the materialLawparams_ call the method with the cells satRegionIdx
     */
    const MaterialLawParams& connectionMaterialLawParams(unsigned satRegionIdx, unsigned elemIdx) const
    {
        MaterialLawParams& mlp = const_cast<MaterialLawParams&>(materialLawParams_[elemIdx]);

#if HAVE_OPM_COMMON
        if (enableHysteresis())
            OpmLog::warning("Warning: Using non-default satnum regions for connection is not tested in combination with hysteresis");
#endif
        // Currently we don't support COMPIMP. I.e. use the same table lookup for the hysteresis curves.
        // unsigned impRegionIdx = satRegionIdx;

        // change the sat table it points to.
        switch (mlp.approach()) {
        case EclMultiplexerApproach::EclStone1Approach: {
            auto& realParams = mlp.template getRealParams<EclMultiplexerApproach::EclStone1Approach>();

            realParams.oilWaterParams().drainageParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
            realParams.oilWaterParams().drainageParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
//            if (enableHysteresis()) {
//                realParams.oilWaterParams().imbibitionParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[impRegionIdx]);
//                realParams.oilWaterParams().imbibitionParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setUnscaledPoints(gasOilUnscaledPointsVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setEffectiveLawParams(gasOilEffectiveParamVector_[impRegionIdx]);
//            }
        }
            break;

        case EclMultiplexerApproach::EclStone2Approach: {
            auto& realParams = mlp.template getRealParams<EclMultiplexerApproach::EclStone2Approach>();
            realParams.oilWaterParams().drainageParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
            realParams.oilWaterParams().drainageParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
//            if (enableHysteresis()) {
//                realParams.oilWaterParams().imbibitionParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[impRegionIdx]);
//                realParams.oilWaterParams().imbibitionParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setUnscaledPoints(gasOilUnscaledPointsVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setEffectiveLawParams(gasOilEffectiveParamVector_[impRegionIdx]);
//            }
        }
            break;

        case EclMultiplexerApproach::EclDefaultApproach: {
            auto& realParams = mlp.template getRealParams<EclMultiplexerApproach::EclDefaultApproach>();
            realParams.oilWaterParams().drainageParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
            realParams.oilWaterParams().drainageParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
//            if (enableHysteresis()) {
//                realParams.oilWaterParams().imbibitionParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[impRegionIdx]);
//                realParams.oilWaterParams().imbibitionParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setUnscaledPoints(gasOilUnscaledPointsVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setEffectiveLawParams(gasOilEffectiveParamVector_[impRegionIdx]);
//            }
        }
            break;

        case EclMultiplexerApproach::EclTwoPhaseApproach: {
            auto& realParams = mlp.template getRealParams<EclMultiplexerApproach::EclTwoPhaseApproach>();
            realParams.oilWaterParams().drainageParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
            realParams.oilWaterParams().drainageParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
//            if (enableHysteresis()) {
//                realParams.oilWaterParams().imbibitionParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[impRegionIdx]);
//                realParams.oilWaterParams().imbibitionParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setUnscaledPoints(gasOilUnscaledPointsVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setEffectiveLawParams(gasOilEffectiveParamVector_[impRegionIdx]);
//            }
        }
            break;

        default:
            throw std::logic_error("Enum value for material approach unknown!");
        }

        return mlp;
    }

    int satnumRegionIdx(unsigned elemIdx) const
    { return satnumRegionArray_[elemIdx]; }

    int getKrnumSatIdx(unsigned elemIdx, FaceDir::DirEnum facedir) const {
        using Dir = FaceDir::DirEnum;
        const std::vector<int>* array = nullptr;
        switch(facedir) {
            case Dir::XPlus:
                array = &krnumXArray_;
                break;
            case Dir::YPlus:
                array = &krnumYArray_;
                break;
            case Dir::ZPlus:
                array = &krnumZArray_;
                break;
            default:
                throw std::runtime_error("Unknown face direction");
        }
        if (array->size() > 0) {
            return (*array)[elemIdx];
        }
        else {
            return satnumRegionArray_[elemIdx];
        }
    }
    bool hasDirectionalRelperms() const {
        if (krnumXArray_.size() > 0 || krnumYArray_.size() > 0 || krnumZArray_.size() > 0) {
            return true;
        }
        return false;
    }
    int imbnumRegionIdx(unsigned elemIdx) const
    { return imbnumRegionArray_[elemIdx]; }

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

        MaterialLaw::updateHysteresis(materialLawParams_[elemIdx], fluidState);
    }

    void oilWaterHysteresisParams(Scalar& pcSwMdc,
                                  Scalar& krnSwMdc,
                                  unsigned elemIdx) const
    {
        if (!enableHysteresis())
            throw std::runtime_error("Cannot get hysteresis parameters if hysteresis not enabled.");

        const auto& params = materialLawParams(elemIdx);
        MaterialLaw::oilWaterHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void setOilWaterHysteresisParams(const Scalar& pcSwMdc,
                                     const Scalar& krnSwMdc,
                                     unsigned elemIdx)
    {
        if (!enableHysteresis())
            throw std::runtime_error("Cannot set hysteresis parameters if hysteresis not enabled.");

        auto& params = materialLawParams(elemIdx);
        MaterialLaw::setOilWaterHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void gasOilHysteresisParams(Scalar& pcSwMdc,
                                Scalar& krnSwMdc,
                                unsigned elemIdx) const
    {
        if (!enableHysteresis())
            throw std::runtime_error("Cannot get hysteresis parameters if hysteresis not enabled.");

        const auto& params = materialLawParams(elemIdx);
        MaterialLaw::gasOilHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void setGasOilHysteresisParams(const Scalar& pcSwMdc,
                                   const Scalar& krnSwMdc,
                                   unsigned elemIdx)
    {
        if (!enableHysteresis())
            throw std::runtime_error("Cannot set hysteresis parameters if hysteresis not enabled.");

        auto& params = materialLawParams(elemIdx);
        MaterialLaw::setGasOilHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    EclEpsScalingPoints<Scalar>& oilWaterScaledEpsPointsDrainage(unsigned elemIdx)
    {
        auto& materialParams = materialLawParams_[elemIdx];
        switch (materialParams.approach()) {
        case EclMultiplexerApproach::EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclStone1Approach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclMultiplexerApproach::EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclStone2Approach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclMultiplexerApproach::EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclDefaultApproach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclMultiplexerApproach::EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclTwoPhaseApproach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }
        default:
            throw std::logic_error("Enum value for material approach unknown!");
        }
    }

    const EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoDrainage(size_t elemIdx) const
    { return oilWaterScaledEpsInfoDrainage_[elemIdx]; }

private:
    void readGlobalEpsOptions_(const EclipseState& eclState)
    {
        oilWaterEclEpsConfig_ = std::make_shared<EclEpsConfig>();
        oilWaterEclEpsConfig_->initFromState(eclState, EclOilWaterSystem);

        enableEndPointScaling_ = eclState.getTableManager().hasTables("ENKRVD");
    }

    void readGlobalHysteresisOptions_(const EclipseState& state)
    {
        hysteresisConfig_ = std::make_shared<EclHysteresisConfig>();
        hysteresisConfig_->initFromState(state.runspec());
    }

    void readGlobalThreePhaseOptions_(const Runspec& runspec)
    {
        bool gasEnabled = runspec.phases().active(Phase::GAS);
        bool oilEnabled = runspec.phases().active(Phase::OIL);
        bool waterEnabled = runspec.phases().active(Phase::WATER);

        int numEnabled =
            (gasEnabled?1:0)
            + (oilEnabled?1:0)
            + (waterEnabled?1:0);

        if (numEnabled == 0) {
            throw std::runtime_error("At least one fluid phase must be enabled. (Is: "+std::to_string(numEnabled)+")");
        } else if (numEnabled == 1) {
            threePhaseApproach_ = EclMultiplexerApproach::EclOnePhaseApproach;
        } else if ( numEnabled == 2) {
            threePhaseApproach_ = EclMultiplexerApproach::EclTwoPhaseApproach;
            if (!gasEnabled)
                twoPhaseApproach_ = EclTwoPhaseApproach::EclTwoPhaseOilWater;
            else if (!oilEnabled)
                twoPhaseApproach_ = EclTwoPhaseApproach::EclTwoPhaseGasWater;
            else if (!waterEnabled)
                twoPhaseApproach_ = EclTwoPhaseApproach::EclTwoPhaseGasOil;
        }
        else {
            assert(numEnabled == 3);

            threePhaseApproach_ = EclMultiplexerApproach::EclDefaultApproach;
            const auto& satctrls = runspec.saturationFunctionControls();
            if (satctrls.krModel() == SatFuncControls::ThreePhaseOilKrModel::Stone2)
                threePhaseApproach_ = EclMultiplexerApproach::EclStone2Approach;
            else if (satctrls.krModel() == SatFuncControls::ThreePhaseOilKrModel::Stone1)
                threePhaseApproach_ = EclMultiplexerApproach::EclStone1Approach;
        }
    }

    template <class Container>
    void readGasOilEffectiveParameters_(Container& dest,
                                        const EclipseState& eclState,
                                        unsigned satRegionIdx)
    {
        if (!hasGas || !hasOil)
            // we don't read anything if either the gas or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<GasOilEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satRegionIdx];

        // the situation for the gas phase is complicated that all saturations are
        // shifted by the connate water saturation.
        const Scalar Swco = unscaledEpsInfo_[satRegionIdx].Swl;
        const auto tolcrit = eclState.runspec().saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto& tableManager = eclState.getTableManager();

        switch (eclState.runspec().saturationFunctionControls().family()) {
        case SatFuncControls::KeywordFamily::Family_I:
        {
            const TableContainer& sgofTables = tableManager.getSgofTables();
            const TableContainer& slgofTables = tableManager.getSlgofTables();
            if (!sgofTables.empty())
                readGasOilEffectiveParametersSgof_(effParams, Swco, tolcrit,
                                                   sgofTables.getTable<SgofTable>(satRegionIdx));
            else if (!slgofTables.empty())
                readGasOilEffectiveParametersSlgof_(effParams, Swco, tolcrit,
                                                    slgofTables.getTable<SlgofTable>(satRegionIdx));
            else if ( !tableManager.getSgofletTable().empty() ) {
                const auto& letSgofTab = tableManager.getSgofletTable()[satRegionIdx];
                const std::vector<Scalar> dum; // dummy arg to comform with existing interface

                effParams.setApproach(SatCurveMultiplexerApproach::LETApproach);
                auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::LETApproach>();

                // S=(So-Sogcr)/(1-Sogcr-Sgcr-Swco),  krog = Krt*S^L/[S^L+E*(1.0-S)^T]
                const Scalar s_min_w = letSgofTab.s2_critical;
                const Scalar s_max_w = 1.0-letSgofTab.s1_critical-Swco;
                const std::vector<Scalar>& letCoeffsOil = {s_min_w, s_max_w,
                                                           static_cast<Scalar>(letSgofTab.l2_relperm),
                                                           static_cast<Scalar>(letSgofTab.e2_relperm),
                                                           static_cast<Scalar>(letSgofTab.t2_relperm),
                                                           static_cast<Scalar>(letSgofTab.krt2_relperm)};
                realParams.setKrwSamples(letCoeffsOil, dum);

                // S=(1-So-Sgcr-Swco)/(1-Sogcr-Sgcr-Swco), krg = Krt*S^L/[S^L+E*(1.0-S)^T]
                const Scalar s_min_nw = letSgofTab.s1_critical+Swco;
                const Scalar s_max_nw = 1.0-letSgofTab.s2_critical;
                const std::vector<Scalar>& letCoeffsGas = {s_min_nw, s_max_nw,
                                                           static_cast<Scalar>(letSgofTab.l1_relperm),
                                                           static_cast<Scalar>(letSgofTab.e1_relperm),
                                                           static_cast<Scalar>(letSgofTab.t1_relperm),
                                                           static_cast<Scalar>(letSgofTab.krt1_relperm)};
                realParams.setKrnSamples(letCoeffsGas, dum);

                // S=(So-Sorg)/(1-Sorg-Sgl-Swco), Pc = Pct + (pcir_pc-Pct)*(1-S)^L/[(1-S)^L+E*S^T]
                const std::vector<Scalar>& letCoeffsPc = {static_cast<Scalar>(letSgofTab.s2_residual),
                                                          static_cast<Scalar>(letSgofTab.s1_residual+Swco),
                                                          static_cast<Scalar>(letSgofTab.l_pc),
                                                          static_cast<Scalar>(letSgofTab.e_pc),
                                                          static_cast<Scalar>(letSgofTab.t_pc),
                                                          static_cast<Scalar>(letSgofTab.pcir_pc),
                                                          static_cast<Scalar>(letSgofTab.pct_pc)};
                realParams.setPcnwSamples(letCoeffsPc, dum);

                realParams.finalize();
            }
            break;
        }

        case SatFuncControls::KeywordFamily::Family_II:
        {
            const SgfnTable& sgfnTable = tableManager.getSgfnTables().getTable<SgfnTable>( satRegionIdx );
            bool co2store = eclState.runspec().co2Storage();
            if (co2store) {
                const SwfnTable& swfnTable = tableManager.getSwfnTables().getTable<SwfnTable>( satRegionIdx );
                readGasOilEffectiveParametersFamily2_(effParams, Swco, tolcrit, swfnTable, sgfnTable);
            }
            else if (!hasWater) {
                // oil and gas case
                const Sof2Table& sof2Table = tableManager.getSof2Tables().getTable<Sof2Table>( satRegionIdx );
                readGasOilEffectiveParametersFamily2_(effParams, Swco, tolcrit, sof2Table, sgfnTable);
            }
            else {
                const Sof3Table& sof3Table = tableManager.getSof3Tables().getTable<Sof3Table>( satRegionIdx );
                readGasOilEffectiveParametersFamily2_(effParams, Swco, tolcrit, sof3Table, sgfnTable);
            }
            break;
        }

        case SatFuncControls::KeywordFamily::Undefined:
            throw std::domain_error("No valid saturation keyword family specified");
        }
    }

    void readGasOilEffectiveParametersSgof_(GasOilEffectiveTwoPhaseParams& effParams,
                                            const Scalar Swco,
                                            const double tolcrit,
                                            const SgofTable& sgofTable)
    {
        // convert the saturations of the SGOF keyword from gas to oil saturations
        std::vector<double> SoSamples(sgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < sgofTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = (1.0 - Swco) - sgofTable.get("SG", sampleIdx);
        }

        effParams.setApproach(SatCurveMultiplexerApproach::PiecewiseLinearApproach);
        auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::PiecewiseLinearApproach>();

        realParams.setKrwSamples(SoSamples, normalizeKrValues_(tolcrit, sgofTable.getColumn("KROG")));
        realParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, sgofTable.getColumn("KRG")));
        realParams.setPcnwSamples(SoSamples, sgofTable.getColumn("PCOG").vectorCopy());
        realParams.finalize();
    }

    void readGasOilEffectiveParametersSlgof_(GasOilEffectiveTwoPhaseParams& effParams,
                                             const Scalar Swco,
                                             const double tolcrit,
                                             const SlgofTable& slgofTable)
    {
        // convert the saturations of the SLGOF keyword from "liquid" to oil saturations
        std::vector<double> SoSamples(slgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < slgofTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = slgofTable.get("SL", sampleIdx) - Swco;
        }

        effParams.setApproach(SatCurveMultiplexerApproach::PiecewiseLinearApproach);
        auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::PiecewiseLinearApproach>();

        realParams.setKrwSamples(SoSamples, normalizeKrValues_(tolcrit, slgofTable.getColumn("KROG")));
        realParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, slgofTable.getColumn("KRG")));
        realParams.setPcnwSamples(SoSamples, slgofTable.getColumn("PCOG").vectorCopy());
        realParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               const Scalar Swco,
                                               const double tolcrit,
                                               const Sof3Table& sof3Table,
                                               const SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        std::vector<double> SoColumn = sof3Table.getColumn("SO").vectorCopy();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = (1.0 - Swco) - sgfnTable.get("SG", sampleIdx);
        }

        effParams.setApproach(SatCurveMultiplexerApproach::PiecewiseLinearApproach);
        auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::PiecewiseLinearApproach>();

        realParams.setKrwSamples(SoColumn, normalizeKrValues_(tolcrit, sof3Table.getColumn("KROG")));
        realParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, sgfnTable.getColumn("KRG")));
        realParams.setPcnwSamples(SoSamples, sgfnTable.getColumn("PCOG").vectorCopy());
        realParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               const Scalar Swco,
                                               const double tolcrit,
                                               const Sof2Table& sof2Table,
                                               const SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        std::vector<double> SoColumn = sof2Table.getColumn("SO").vectorCopy();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = (1.0 - Swco) - sgfnTable.get("SG", sampleIdx);
        }

        effParams.setApproach(SatCurveMultiplexerApproach::PiecewiseLinearApproach);
        auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::PiecewiseLinearApproach>();

        realParams.setKrwSamples(SoColumn, normalizeKrValues_(tolcrit, sof2Table.getColumn("KRO")));
        realParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, sgfnTable.getColumn("KRG")));
        realParams.setPcnwSamples(SoSamples, sgfnTable.getColumn("PCOG").vectorCopy());
        realParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               const Scalar Swco,
                                               const double tolcrit,
                                               const SwfnTable& swfnTable,
                                               const SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        std::vector<double> SoColumn = swfnTable.getColumn("SW").vectorCopy();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = (1.0 - Swco) - sgfnTable.get("SG", sampleIdx);
        }

        effParams.setApproach(SatCurveMultiplexerApproach::PiecewiseLinearApproach);
        auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::PiecewiseLinearApproach>();

        realParams.setKrwSamples(SoColumn, normalizeKrValues_(tolcrit, swfnTable.getColumn("KRW")));
        realParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, sgfnTable.getColumn("KRG")));
        realParams.setPcnwSamples(SoSamples, sgfnTable.getColumn("PCOG").vectorCopy());
        realParams.finalize();
    }

    template <class Container>
    void readOilWaterEffectiveParameters_(Container& dest,
                                          const EclipseState& eclState,
                                          unsigned satRegionIdx)
    {
        if (!hasOil || !hasWater)
            // we don't read anything if either the water or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<OilWaterEffectiveTwoPhaseParams>();

        const auto tolcrit = eclState.runspec().saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto& tableManager = eclState.getTableManager();
        auto& effParams = *dest[satRegionIdx];

        switch (eclState.runspec().saturationFunctionControls().family()) {
        case SatFuncControls::KeywordFamily::Family_I:
        {
            if (tableManager.hasTables("SWOF")) {
                const auto& swofTable = tableManager.getSwofTables().getTable<SwofTable>(satRegionIdx);
                const std::vector<double> SwColumn = swofTable.getColumn("SW").vectorCopy();

                effParams.setApproach(SatCurveMultiplexerApproach::PiecewiseLinearApproach);
                auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::PiecewiseLinearApproach>();

                realParams.setKrwSamples(SwColumn, normalizeKrValues_(tolcrit, swofTable.getColumn("KRW")));
                realParams.setKrnSamples(SwColumn, normalizeKrValues_(tolcrit, swofTable.getColumn("KROW")));
                realParams.setPcnwSamples(SwColumn, swofTable.getColumn("PCOW").vectorCopy());
                realParams.finalize();
            }
            else if ( !tableManager.getSwofletTable().empty() ) {
                const auto& letTab = tableManager.getSwofletTable()[satRegionIdx];
                const std::vector<Scalar> dum; // dummy arg to conform with existing interface

                effParams.setApproach(SatCurveMultiplexerApproach::LETApproach);
                auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::LETApproach>();

                // S=(Sw-Swcr)/(1-Sowcr-Swcr),  krw = Krt*S^L/[S^L+E*(1.0-S)^T]
                const Scalar s_min_w = letTab.s1_critical;
                const Scalar s_max_w = 1.0-letTab.s2_critical;
                const std::vector<Scalar>& letCoeffsWat = {s_min_w, s_max_w,
                                                           static_cast<Scalar>(letTab.l1_relperm),
                                                           static_cast<Scalar>(letTab.e1_relperm),
                                                           static_cast<Scalar>(letTab.t1_relperm),
                                                           static_cast<Scalar>(letTab.krt1_relperm)};
                realParams.setKrwSamples(letCoeffsWat, dum);

                // S=(So-Sowcr)/(1-Sowcr-Swcr), krow = Krt*S^L/[S^L+E*(1.0-S)^T]
                const Scalar s_min_nw = letTab.s2_critical;
                const Scalar s_max_nw = 1.0-letTab.s1_critical;
                const std::vector<Scalar>& letCoeffsOil = {s_min_nw, s_max_nw,
                                                           static_cast<Scalar>(letTab.l2_relperm),
                                                           static_cast<Scalar>(letTab.e2_relperm),
                                                           static_cast<Scalar>(letTab.t2_relperm),
                                                           static_cast<Scalar>(letTab.krt2_relperm)};
                realParams.setKrnSamples(letCoeffsOil, dum);

                // S=(Sw-Swco)/(1-Swco-Sorw), Pc = Pct + (Pcir-Pct)*(1-S)^L/[(1-S)^L+E*S^T]
                const std::vector<Scalar>& letCoeffsPc = {static_cast<Scalar>(letTab.s1_residual),
                                                          static_cast<Scalar>(letTab.s2_residual),
                                                          static_cast<Scalar>(letTab.l_pc),
                                                          static_cast<Scalar>(letTab.e_pc),
                                                          static_cast<Scalar>(letTab.t_pc),
                                                          static_cast<Scalar>(letTab.pcir_pc),
                                                          static_cast<Scalar>(letTab.pct_pc)};
                realParams.setPcnwSamples(letCoeffsPc, dum);

                realParams.finalize();
            }
            break;
        }

        case SatFuncControls::KeywordFamily::Family_II:
        {
            const auto& swfnTable = tableManager.getSwfnTables().getTable<SwfnTable>(satRegionIdx);
            const std::vector<double> SwColumn = swfnTable.getColumn("SW").vectorCopy();

            effParams.setApproach(SatCurveMultiplexerApproach::PiecewiseLinearApproach);
            auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::PiecewiseLinearApproach>();

            realParams.setKrwSamples(SwColumn, normalizeKrValues_(tolcrit, swfnTable.getColumn("KRW")));
            realParams.setPcnwSamples(SwColumn, swfnTable.getColumn("PCOW").vectorCopy());

            if (!hasGas) {
                const auto& sof2Table = tableManager.getSof2Tables().getTable<Sof2Table>(satRegionIdx);
                // convert the saturations of the SOF2 keyword from oil to water saturations
                std::vector<double> SwSamples(sof2Table.numRows());
                for (size_t sampleIdx = 0; sampleIdx < sof2Table.numRows(); ++ sampleIdx)
                    SwSamples[sampleIdx] = 1 - sof2Table.get("SO", sampleIdx);

                realParams.setKrnSamples(SwSamples, normalizeKrValues_(tolcrit, sof2Table.getColumn("KRO")));
            } else {
                const auto& sof3Table = tableManager.getSof3Tables().getTable<Sof3Table>(satRegionIdx);
                // convert the saturations of the SOF3 keyword from oil to water saturations
                std::vector<double> SwSamples(sof3Table.numRows());
                for (size_t sampleIdx = 0; sampleIdx < sof3Table.numRows(); ++ sampleIdx)
                    SwSamples[sampleIdx] = 1 - sof3Table.get("SO", sampleIdx);

                realParams.setKrnSamples(SwSamples, normalizeKrValues_(tolcrit, sof3Table.getColumn("KROW")));
            }
            realParams.finalize();
            break;
        }

        case SatFuncControls::KeywordFamily::Undefined:
            throw std::domain_error("No valid saturation keyword family specified");
        }
    }

    template <class Container>
    void readGasWaterEffectiveParameters_(Container& dest,
                                        const EclipseState& eclState,
                                        unsigned satRegionIdx)
    {
        if (!hasGas || !hasWater || hasOil)
            // we don't read anything if either the gas or the water phase is not active or if oil is present
            return;

        dest[satRegionIdx] = std::make_shared<GasWaterEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satRegionIdx];

        const auto tolcrit = eclState.runspec().saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto& tableManager = eclState.getTableManager();

        switch (eclState.runspec().saturationFunctionControls().family()) {
        case SatFuncControls::KeywordFamily::Family_I:
        {
            throw std::domain_error("Saturation keyword family I is not applicable for a gas-water system");
        }

        case SatFuncControls::KeywordFamily::Family_II:
        {
            //Todo: allow also for Sgwfn table input as alternative to Sgfn and Swfn table input
            const SgfnTable& sgfnTable = tableManager.getSgfnTables().getTable<SgfnTable>( satRegionIdx );
            const SwfnTable& swfnTable = tableManager.getSwfnTables().getTable<SwfnTable>( satRegionIdx );

            effParams.setApproach(SatCurveMultiplexerApproach::PiecewiseLinearApproach);
            auto& realParams = effParams.template getRealParams<SatCurveMultiplexerApproach::PiecewiseLinearApproach>();

            std::vector<double> SwColumn = swfnTable.getColumn("SW").vectorCopy();

            realParams.setKrwSamples(SwColumn, normalizeKrValues_(tolcrit, swfnTable.getColumn("KRW")));
            std::vector<double> SwSamples(sgfnTable.numRows());
            for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx)
                SwSamples[sampleIdx] = 1 - sgfnTable.get("SG", sampleIdx);
            realParams.setKrnSamples(SwSamples, normalizeKrValues_(tolcrit, sgfnTable.getColumn("KRG")));
            //Capillary pressure is read from SWFN. 
            //For gas-water system the capillary pressure column values are set to 0 in SGFN
            realParams.setPcnwSamples(SwColumn, swfnTable.getColumn("PCOW").vectorCopy());
            realParams.finalize();
                       
            break;
        }

        case SatFuncControls::KeywordFamily::Undefined:
            throw std::domain_error("No valid saturation keyword family specified");
        }
    }

    template <class Container>
    void readGasOilUnscaledPoints_(Container& dest,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const EclipseState& /* eclState */,
                                   unsigned satRegionIdx)
    {
        if (!hasGas || !hasOil)
            // we don't read anything if either the gas or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclGasOilSystem);
    }

    template <class Container>
    void readOilWaterUnscaledPoints_(Container& dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     const EclipseState& /* eclState */,
                                     unsigned satRegionIdx)
    {
        if (!hasOil || !hasWater)
            // we don't read anything if either the water or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclOilWaterSystem);
    }

    template <class Container>
    void readGasWaterUnscaledPoints_(Container& dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     const EclipseState& /* eclState */,
                                     unsigned satRegionIdx)
    {
        if (hasOil)
            // we don't read anything if oil phase is active
            return;

        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclGasWaterSystem);
    }

    std::tuple<EclEpsScalingPointsInfo<Scalar>,
               EclEpsScalingPoints<Scalar>>
    readScaledPoints_(const EclEpsConfig& config,
                      const EclipseState& eclState,
                      const EclEpsGridProperties& epsGridProperties,
                      unsigned elemIdx,
                      EclTwoPhaseSystemType type)
    {
        unsigned satRegionIdx = epsGridProperties.satRegion( elemIdx );

        EclEpsScalingPointsInfo<Scalar> destInfo(unscaledEpsInfo_[satRegionIdx]);
        destInfo.extractScaled(eclState, epsGridProperties, elemIdx);

        EclEpsScalingPoints<Scalar> destPoint;
        destPoint.init(destInfo, config, type);

        return {destInfo, destPoint};
    }

    void initThreePhaseParams_(const EclipseState& /* eclState */,
                               MaterialLawParams& materialParams,
                               unsigned satRegionIdx,
                               const EclEpsScalingPointsInfo<Scalar>& epsInfo,
                               std::shared_ptr<OilWaterTwoPhaseHystParams> oilWaterParams,
                               std::shared_ptr<GasOilTwoPhaseHystParams> gasOilParams,
                               std::shared_ptr<GasWaterTwoPhaseHystParams> gasWaterParams)
    {
        materialParams.setApproach(threePhaseApproach_);

        switch (materialParams.approach()) {
        case EclMultiplexerApproach::EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclStone1Approach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);

            if (!stoneEtas.empty()) {
                realParams.setEta(stoneEtas[satRegionIdx]);
            }
            else
                realParams.setEta(1.0);
            realParams.finalize();
            break;
        }

        case EclMultiplexerApproach::EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclStone2Approach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);
            realParams.finalize();
            break;
        }

        case EclMultiplexerApproach::EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclDefaultApproach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);
            realParams.finalize();
            break;
        }

        case EclMultiplexerApproach::EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclTwoPhaseApproach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setGasWaterParams(gasWaterParams);
            realParams.setApproach(twoPhaseApproach_);
            realParams.finalize();
            break;
        }

        case EclMultiplexerApproach::EclOnePhaseApproach: {
            // Nothing to do, no parameters.
            break;
        }
        }
    }

    // Relative permeability values not strictly greater than 'tolcrit' treated as zero.
    std::vector<double> normalizeKrValues_(const double tolcrit,
                                           const TableColumn& krValues) const
    {
        auto kr = krValues.vectorCopy();
        std::transform(kr.begin(), kr.end(), kr.begin(),
            [tolcrit](const double kri)
        {
            return (kri > tolcrit) ? kri : 0.0;
        });

        return kr;
    }

    bool enableEndPointScaling_;
    std::shared_ptr<EclHysteresisConfig> hysteresisConfig_;

    std::shared_ptr<EclEpsConfig> oilWaterEclEpsConfig_;
    std::vector<EclEpsScalingPointsInfo<Scalar>> unscaledEpsInfo_;
    OilWaterScalingInfoVector oilWaterScaledEpsInfoDrainage_;

    std::shared_ptr<EclEpsConfig> gasWaterEclEpsConfig_;

    GasOilScalingPointsVector gasOilUnscaledPointsVector_;
    OilWaterScalingPointsVector oilWaterUnscaledPointsVector_;
    GasWaterScalingPointsVector gasWaterUnscaledPointsVector_;

    GasOilEffectiveParamVector gasOilEffectiveParamVector_;
    OilWaterEffectiveParamVector oilWaterEffectiveParamVector_;
    GasWaterEffectiveParamVector gasWaterEffectiveParamVector_;

    EclMultiplexerApproach threePhaseApproach_ = EclMultiplexerApproach::EclDefaultApproach;
    // this attribute only makes sense for twophase simulations!
    enum EclTwoPhaseApproach twoPhaseApproach_ = EclTwoPhaseApproach::EclTwoPhaseGasOil;

    std::vector<MaterialLawParams> materialLawParams_;

    std::vector<int> satnumRegionArray_;
    std::vector<int> krnumXArray_;
    std::vector<int> krnumYArray_;
    std::vector<int> krnumZArray_;
    std::vector<int> imbnumRegionArray_;
    std::vector<Scalar> stoneEtas;

    bool hasGas;
    bool hasOil;
    bool hasWater;

    std::shared_ptr<EclEpsConfig> gasOilConfig;
    std::shared_ptr<EclEpsConfig> oilWaterConfig;
    std::shared_ptr<EclEpsConfig> gasWaterConfig;
};
} // namespace Opm

#endif
