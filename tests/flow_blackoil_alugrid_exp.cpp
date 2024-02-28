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
#define DISABLE_ALUGRID_SFC_ORDERING 1
#include "config.h"
#include <opm/grid/cpgrid/GridHelpers.hpp>
#include <dune/alugrid/grid.hh>
#include <ebos/eclalugridvanguard.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/common/transfluxmodule.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <dune/alugrid/grid.hh>
#include <opm/grid/CpGrid.hpp>
#include <ebos/eclproblem.hh>
#include <opm/flowdynamicgrid/eclproblemdynamic.hh>
#include <opm/models/discretization/common/fvbasegradientcalculator.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/fvbaselinearizer.hh>
#include <opm/models/discretization/common/fvbaseintensivequantities.hh>
#include <opm/flowdynamicgrid/blackoilintensivequantitiesdynamic.hh>
// these are not explicitly instanced in library
#include <ebos/collecttoiorank_impl.hh>
#include <ebos/eclgenericproblem_impl.hh>
#include <ebos/eclgenericthresholdpressure_impl.hh>
#include <ebos/eclgenerictracermodel_impl.hh>
#include <ebos/ecltransmissibility_impl.hh>
#include <ebos/eclgenericwriter_impl.hh>
#include <ebos/equil/initstateequil_impl.hh>
#include <opm/models/discretization/common/fvbasediscretizationfemadapt.hh>
//#include <opm/material/fluidmatrixinteractions/EclMaterialLawManagerTable.hpp>
namespace Opm{
    template<typename TypeTag>
    class EclProblemNew: public EclProblem<TypeTag>{
        using Parent = EclProblem<TypeTag>;
        using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
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
        //template <class Context, class FluidState>
       //void updateRelperms( const Context& context,
      //      std::array<Evaluation,numPhases> &mobility,
      //      DirectionalMobilityPtr &dirMob,
     //       FluidState &fluidState,
     //       unsigned dofIdx, unsigned timeIdx) const
     //   {
     //       OPM_TIMEBLOCK_LOCAL(updateRelperms);
     //       Parent::updateRelperms(context, mobility, dirMob, fluidState, dofIdx, timeIdx);
     //   };
};
       template<typename TypeTag>
    class BlackOilModelDynamic: public FIBlackOilModel<TypeTag>{
        using Parent = BlackOilModel<TypeTag>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Element = typename GridView::template Codim<0>::Entity;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
        using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;
        using ElementIterator = typename GridView::template Codim<0>::Iterator;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        public:
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;

        using Indices = GetPropType<TypeTag, Properties::Indices>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        static constexpr int numEq = Indices::numEq;
        using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
        using BVector = Dune::BlockVector<VectorBlockType>;
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        enum {
        historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>(),
        };
    public:
        BlackOilModelDynamic(Simulator& simulator): FIBlackOilModel<TypeTag>(simulator){
        }
 //   void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx) const
//    {
//       OPM_TIMEBLOCK_LOCAL(invalidateAndUpdateIntensiveQuantities);
//       Parent::invalidateIntensiveQuantitiesCache(timeIdx);
//    }
// Overwriting that function to avoid throwing error when having dune-fem

    void updateCartesianToCompressedMapping_()
    {
        OPM_TIMEBLOCK_LOCAL(updateCartesianToCompressedMapping_);
        std::size_t num_cells = this->asImp_().grid().leafGridView().size(0);
        this->is_interior_.resize(num_cells);

        Opm::Properties::ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());
        for (const auto& element : elements(this->gridView()))
        {
            const auto elemIdx = elemMapper.index(element);
            unsigned cartesianCellIdx = this->cartesianIndex(elemIdx);
            this->cartesianToCompressed_[cartesianCellIdx] = elemIdx;
            if (element.partitionType() == Dune::InteriorEntity)
            {
                this->is_interior_[elemIdx] = 1;
            }
            else
            {
                this->is_interior_[elemIdx] = 0;
            }
        }
    }
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

        /*!
     * \brief Called by the update() method when the grid should be refined.
     */
  //  void adaptGrid()
  //  {

//    }

    void updateSolution(const BVector& dx)
    {
         OPM_TIMEBLOCK(updateSolution);
         auto& ebosNewtonMethod = this->simulator_.model().newtonMethod();
         SolutionVector& solution = this->simulator_.model().solution(/*timeIdx=*/0);

         ebosNewtonMethod.update_(/*nextSolution=*/solution,
                                  /*curSolution=*/solution,
                                  /*update=*/dx,
                                  /*resid=*/dx); // the update routines of the blac
                                                 // oil model do not care about the
                                                 // residual

         // if the solution is updated, the intensive quantities need to be recalculated
         {
            OPM_TIMEBLOCK_LOCAL(invalidateAndUpdateIntensiveQuantities);
            this->simulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
            //ebosSimulator_.problem().eclWriter()->mutableEclOutputModule().invalidateLocalData();
         }
     }
   //  void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
   // {
   //     OPM_TIMEBLOCK_LOCAL(update);
   //     ParentType::update(elemCtx, dofIdx, timeIdx);
   //     OPM_TIMEBLOCK_LOCAL(blackoilIntensiveQuanititiesUpdate);
  //  }

 //          void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx) const
 //          {
 //      OPM_TIMEBLOCK_LOCAL(invalidateAndUpdateIntensiveQuantities);
 //      Parent::invalidateAndUpdateIntensiveQuantities(timeIdx);
 //   }
        //     IntensiveQuantities intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const{
        //     OPM_TIMEBLOCK_LOCAL(intensiveQuantitiesNoCache);
        //     const auto& primaryVar = this->solution(timeIdx)[globalIdx];
        //     const auto& problem = this->simulator_.problem();
        //     //IntensiveQuantities* intQuants = &(this->intensiveQuantityCache_[timeIdx][globalIdx]);
        //     if (!(this->enableIntensiveQuantityCache_) ||
        //         !(this->intensiveQuantityCacheUpToDate_[timeIdx][globalIdx])){
        //         IntensiveQuantities intQuants;
        //         intQuants.update(problem,primaryVar, globalIdx, timeIdx);
        //         return intQuants;// reqiored for updating extrution factor
        //     }else{
        //         IntensiveQuantities intQuants = (this->intensiveQuantityCache_[timeIdx][globalIdx]);
        //         return intQuants;
        //     }

        // }

       };
}
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemAlugrid {
    using InheritsFrom = std::tuple<FlowProblem>;
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
    // alugrid need cp grid as equilgrid
    template<class TypeTag>
    struct EquilGrid<TypeTag, TTag::EclFlowProblemAlugrid> {
    using type = Dune::CpGrid;
    };
    template<class TypeTag>
    struct Problem<TypeTag, TTag::EclFlowProblemAlugrid> {
        using type = EclProblemDynamic<TypeTag>;
    };
    template<class TypeTag>
    struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemAlugrid> {
    using type = BlackOilIntensiveQuantitiesDynamic<TypeTag>;
    };
    template<class TypeTag>
    struct Vanguard<TypeTag, TTag::EclFlowProblemAlugrid> {
        using type = Opm::EclAluGridVanguard<TypeTag>;
    };
// template<class TypeTag>
// struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemAlugrid> {
//     static constexpr bool value = false;
// };
//template<class TypeTag>
 //   struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemAlugrid> {
//    using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
//    };


  //  template<class TypeTag>
 //   struct LocalResidual<TypeTag, TTag::EclFlowProblemTest> { using type = BlackOilLocalResidualTPFA<TypeTag>; };
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

template<class TypeTag>
struct GradientCalculator<TypeTag, TTag::EclFlowProblemAlugrid> {
    using type = FvBaseGradientCalculator<TypeTag>;
};

template<class TypeTag>
struct BaseDiscretizationType<TypeTag,TTag::EclFlowProblemAlugrid> {
    using type = FvBaseDiscretizationFemAdapt<TypeTag>;
};
template<class TypeTag>
struct GridPart<TypeTag, TTag::EclFlowProblemAlugrid>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using type = Dune::Fem::AdaptiveLeafGridPart<Grid>;
};

template<class TypeTag>
struct GridView<TypeTag, TTag::EclFlowProblemAlugrid> { using type = typename GetPropType<TypeTag, Properties::GridPart>::GridViewType; };

template<class TypeTag>
struct DiscreteFunctionSpace<TypeTag, TTag::EclFlowProblemAlugrid>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>  ;
    using GridPart = GetPropType<TypeTag, Properties::GridPart>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using FunctionSpace = Dune::Fem::FunctionSpace<typename GridPart::GridType::ctype,
                                                   Scalar,
                                                   GridPart::GridType::dimensionworld,
                                                   numEq>;
    using type = Dune::Fem::FiniteVolumeSpace< FunctionSpace, GridPart, 0 >;

};

template<class TypeTag>
struct DiscreteFunction<TypeTag, TTag::EclFlowProblemAlugrid> {
    using DiscreteFunctionSpace  = GetPropType<TypeTag, Properties::DiscreteFunctionSpace>;
    using PrimaryVariables  = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using type = Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace, PrimaryVariables>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::EclFlowProblemAlugrid>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::BlackOilFluidSystem<Scalar>;
};


}
}
int main(int argc, char** argv)
{
  //  OPM_TIMEBLOCK(fullSimulation);
  //  using TypeTag = Opm::Properties::TTag::EclFlowProblemAlugrid;

  //  auto mainObject = std::make_unique<Opm::Main>(argc, argv);
  //  auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
  //  mainObject.reset();

    //Opm::Main mainObject(argc, argv);
    //auto ret = mainObject.runStatic<TypeTag>();
   // return ret;

    using TypeTag = Opm::Properties::TTag::EclFlowProblemAlugrid;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
