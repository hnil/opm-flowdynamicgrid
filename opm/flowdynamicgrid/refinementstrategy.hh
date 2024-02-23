#ifndef OPM_REFINEMENTSTRATEGY_HEADER_INCLUDED
#define OPM_REFINEMENTSTRATEGY_HEADER_INCLUDED
#include <array>
#include <cstddef>
#include <vector>
#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm{
    class RefinementStrategy
    {
    public:
        RefinementStrategy(){};
        RefinementStrategy(const Opm::PropertyTree& ptree);
        bool isInitialRefined(std::size_t initialIdx);
        bool isInitialCoarsened(std::size_t initialIdx);
        bool shouldBeRefined(double indicator,int level);
        bool shouldBeCoarsened(bool hasSamePrimaryVarsMeaning, double indicator,int level);
        const PropertyTree& ptree(){
            return ptree_;
        }
    private:
        Opm::PropertyTree ptree_;
        std::vector<int> intialrefined_;
    };
}
#endif
