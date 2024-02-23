#include "refinementstrategy.hh"


namespace Opm{
    RefinementStrategy::RefinementStrategy(const PropertyTree& ptree):
        ptree_(ptree)
    {
    }

    bool RefinementStrategy::shouldBeRefined(double indicator,int level){
        return indicator > 0.3 && level < 2;
    }
    bool RefinementStrategy::shouldBeCoarsened(bool hasSamePrimaryVarsMeaning, double indicator,int level){
        return hasSamePrimaryVarsMeaning && indicator < 0.025 && level > 0 ;
    }
    bool RefinementStrategy::isInitialRefined(size_t initialIdx){
        return true;
    }
    bool RefinementStrategy::isInitialCoarsened(size_t initialIdx){
        return false;
    }
}
