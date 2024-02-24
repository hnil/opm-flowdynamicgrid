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
        bool isInitialRefined(std::size_t initialIdx, int level);
        bool isInitialCoarsened(std::size_t initialIdx, int level);
        bool shouldBeRefined(double indicator,int level);
        bool shouldBeCoarsened(bool hasSamePrimaryVarsMeaning, double indicator,int level);
        const PropertyTree& ptree(){
            return ptree_;
        }
        int minNumMarked();
    private:
        Opm::PropertyTree ptree_;
        std::vector<double> initind_;
    };
}
#endif
