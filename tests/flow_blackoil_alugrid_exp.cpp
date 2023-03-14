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
namespace Opm{
    template<typename TypeTag>
    class EclProblemNew: public EclProblem<TypeTag>{
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
        using WaterMeaning = typename PrimaryVariables::WaterMeaning;
        using PressureMeaning = typename PrimaryVariables::PressureMeaning;
        using GasMeaning = typename PrimaryVariables::GasMeaning;
        enum class PrimaryVarsMeaning {
        WaterMeaning,  //Sw, Rvw, Rsw, disabled; (Water Meaning)
        PressureMeaning, // Po, Pg, Pw, disable; (Pressure Meaning)
        GasMeaning, // Rg, Rs, Rvm disabled; (Gas Meaning)
        Undef, // The primary variable is not used
        };

        struct ProblemContainer {
               PrimaryVarsMeaning pm;
               int origGridIndex;
              };
    using GlobalContainer = Dune::PersistentContainer< Grid, ProblemContainer>;
    using RestrictProlongOperator = CopyRestrictProlong< Grid, GlobalContainer >;
//private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
        enum { numPhases = FluidSystem::numPhases };
    GlobalContainer container_;
public:
        EclProblemNew(Simulator& simulator): EclProblem<TypeTag>(simulator), container_(simulator.vanguard().grid(), 0 ){
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
}
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemAlugrid {
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}
   template<class TypeTag>
    struct Grid<TypeTag, TTag::EclFlowProblemAlugrid> {
        static const int dim = 3;
        using type = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming >;
    };
    template<class TypeTag>
    struct Problem<TypeTag, TTag::EclFlowProblemAlugrid> {
        using type = EclProblemNew<TypeTag>;
    };
    
template<class TypeTag>
struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemAlugrid> {
    static constexpr bool value = false;
};
//template <class TypeTag>
//struct FluxModule<TypeTag, TTag::EclFlowProblemAlugrid> {
//    using type = Opm::TransFluxModule<TypeTag>;
//};
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::EclFlowProblemAlugrid> { using type = TTag::AutoDiffLocalLinearizer; };
// use the element centered finite volume spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::EclFlowProblemAlugrid> { using type = TTag::EcfvDiscretization; };
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
}
}
int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblemAlugrid;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
