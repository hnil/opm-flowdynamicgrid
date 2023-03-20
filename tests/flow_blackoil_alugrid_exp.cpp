/*
  Copyright 2020, NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/common/transfluxmodule.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <dune/alugrid/grid.hh>
#include <ebos/eclproblem.hh>
#include <opm/flowdynamicgrid/eclproblemdynamic.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/fvbaselinearizer.hh>
#include <opm/models/discretization/common/fvbaseintensivequantities.hh>
//#include <opm/material/fluidmatrixinteractions/EclMaterialLawManagerTable.hpp>
namespace Opm{
    template<typename TypeTag>
    class EclProblemNew: public EclProblem<TypeTag>{
        using Parent = EclProblem<TypeTag>;
        using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        static constexpr bool waterEnabled = Indices::waterEnabled;
        static constexpr bool gasEnabled = Indices::gasEnabled;
        static constexpr bool oilEnabled = Indices::oilEnabled;
        using DirectionalMobilityPtr = ::Opm::Utility::CopyablePtr<DirectionalMobility<TypeTag, Evaluation>>;
        using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
        using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
        enum { numPhases = FluidSystem::numPhases };
public:
        EclProblemNew(Simulator& simulator): EclProblem<TypeTag>(simulator){
        }
        template <class FluidState>
        void updateRelperms(
            std::array<Evaluation,numPhases> &mobility,
            DirectionalMobilityPtr &/*dirMob*/,
            FluidState &fluidState,
            unsigned globalSpaceIdx) const
        {
            OPM_TIMEBLOCK_LOCAL(updateRelperms);
            {
                // calculate relative permeabilities. note that we store the result into the
                // mobility_ class attribute. the division by the phase viscosity happens later.
                const auto& materialParams = this->materialLawParams(globalSpaceIdx);
                MaterialLaw::relativePermeabilities(mobility, materialParams, fluidState);
                Valgrind::CheckDefined(mobility);
            }
        };
};
       template<typename TypeTag>
    class BlackOilModelDynamic: public BlackOilModel<TypeTag>{
        using Parent = BlackOilModel<TypeTag>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Element = typename GridView::template Codim<0>::Entity;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
        using ElementIterator = typename GridView::template Codim<0>::Iterator;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
        enum {
        historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>(),
        };
    public:
        BlackOilModelDynamic(Simulator& simulator): BlackOilModel<TypeTag>(simulator){
        }
    void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx) const
    {
       OPM_TIMEBLOCK_LOCAL(invalidateAndUpdateIntensiveQuantities);
       Parent::invalidateIntensiveQuantitiesCache(timeIdx);
    }
// Overwriting that function to avoid throwing error when having dune-fem
      void addAuxiliaryModule(BaseAuxiliaryModule<TypeTag>* auxMod)
    {
        OPM_TIMEBLOCK_LOCAL(addAuxiliaryModule);
        auxMod->setDofOffset(this->numTotalDof());
        this->auxEqModules_.push_back(auxMod);


        size_t numDof = this->numTotalDof();
        for (unsigned timeIdx = 0; timeIdx < this->historySize; ++timeIdx)
            this->solution(timeIdx).resize(numDof);

        auxMod->applyInitial();
    }
    };
}
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemAlugrid {
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}
 template<class TypeTag>
 struct Stencil<TypeTag, TTag::EclFlowProblemAlugrid>
 {
 private:
     using Scalar = GetPropType<TypeTag, Properties::Scalar>;
     using GridView = GetPropType<TypeTag, Properties::GridView>;

 public:
     using type = EcfvStencil<Scalar,
                              GridView,
                              /*needIntegrationPos=*/true,
                              /*needIntegrationPos=*/true>;
 };
    template<class TypeTag>
    struct Model<TypeTag, TTag::EclFlowProblemAlugrid> {
        using type = BlackOilModelDynamic<TypeTag>;
    };
    template<class TypeTag>
    struct SpatialDiscretizationSplice<TypeTag, TTag::EclFlowProblemAlugrid> {
        using type = TTag::EcfvDiscretization;
    };
   template<class TypeTag>
    struct Grid<TypeTag, TTag::EclFlowProblemAlugrid> {
        static const int dim = 3;
        using type = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming >;
    };
    template<class TypeTag>
    struct Problem<TypeTag, TTag::EclFlowProblemAlugrid> {
        using type = EclProblemDynamic<TypeTag>;
    };
    
template<class TypeTag>
struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemAlugrid> {
    static constexpr bool value = false;
};
// use automatic differentiation for this simulator
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::EclFlowProblemAlugrid> { using type = TTag::AutoDiffLocalLinearizer; };
// By default, include the intrinsic permeability tensor to the VTK output files
template<class TypeTag>
struct VtkWriteIntrinsicPermeabilities<TypeTag, TTag::EclFlowProblemAlugrid> { static constexpr bool value = true; };

// enable the storage cache by default for this problem
template<class TypeTag>
struct EnableStorageCache<TypeTag, TTag::EclFlowProblemAlugrid> { static constexpr bool value = true; };

// enable the cache for intensive quantities by default for this problem
template<class TypeTag>
struct EnableIntensiveQuantityCache<TypeTag, TTag::EclFlowProblemAlugrid> { static constexpr bool value = true; };
// this problem works fine if the linear solver uses single precision scalars
template<class TypeTag>
struct LinearSolverScalar<TypeTag, TTag::EclFlowProblemAlugrid> { using type = float; };
template <class TypeTag>
struct FluxModule<TypeTag, TTag::EclFlowProblemAlugrid> {
    using type = TransFluxModule<TypeTag>;
};
}
}
int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblemAlugrid;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
