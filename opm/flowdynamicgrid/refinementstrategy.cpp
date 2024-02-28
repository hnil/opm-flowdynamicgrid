#include "refinementstrategy.hh"
#include <string>
#include <iostream>
#include <fstream>
namespace Opm{
    RefinementStrategy::RefinementStrategy(const PropertyTree& ptree):
        ptree_(ptree)
    {
        //read initial indicator file
        bool initialIndicator = ptree.get<bool>("UseInitialFile");
        if(initialIndicator){
            std::cout << "Using indicatorfile" << std::endl;
            std::string filename(ptree.get<std::string>("InitIndicatorFile"));
            std::ifstream myfile(filename);
            if ( myfile.is_open() ){
                while(myfile){
                    double indicator;
                    myfile >> indicator;
                    initind_.push_back(indicator);
                }
            }else{
                std::cout << "No indicator file" << std::endl;
            }
        }
    }

    bool RefinementStrategy::shouldBeRefined(double indicator,int level){
        return (indicator > ptree_.get<double>("RefIndicator")) && (level < std::min(ptree_.get<int>("MaxLevel"),5));
    }
    bool RefinementStrategy::shouldBeCoarsened(bool hasSamePrimaryVarsMeaning, double indicator,int level){
        if(ptree_.get<bool>("ShouldCoarsen")){
            return hasSamePrimaryVarsMeaning && (indicator < ptree_.get<double>("CoarsenIndicator")) && (level > ptree_.get<int>("MinLevel")) ;
        }else{
            return false;
        }
    }
    bool RefinementStrategy::isInitialRefined(size_t initialIdx,int level){
        if(initind_.size()>0){
            return (initind_[initialIdx] > ptree_.get<double>("InitRefIndicator")) && (level < 3);
        }else{
            return true && (level < 4);
        }
    }
    bool RefinementStrategy::isInitialCoarsened(size_t initialIdx, int level){
        return false;
    }
    int RefinementStrategy::minNumMarked(){
        return ptree_.get<int>("MinNumMarked");
    }
}
