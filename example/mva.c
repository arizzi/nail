#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class MVAWrapper {
 public: 
    MVAWrapper() {
        reader = new TMVA::Reader("Silent");
        reader->AddVariable("ll_mass",(float *)0);   
        reader->AddVariable("MqqLog",(float *)0);
        reader->AddVariable("mumujj_pt",(float *)0);
        reader->AddVariable("DeltaEtaQQ",(float *)0);
        reader->AddVariable("softActivityEWK_njets5",(float *)0);

        reader->AddVariable("ll_zstar",(float *)0);
        reader->AddVariable("ll_pt",(float *)0);
        reader->AddVariable("theta2",(float *)0); 
        reader->AddVariable("impulsoZ",(float *)0);
        reader->AddVariable("maxAbsEta",(float *)0);
	
        reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/octoberTraining/Classification_BDTG_1704SeedGen.weights_30000Event_100Tree_2Deep_mll_MqqLog_mumujjPt_DEtajj_softN5_llZstar_llPt_cosMuMuJ2_mumujjPz_maxAbsEta.xml");
}

float eval(std::vector<float> f){
        return  reader->EvaluateMVA(f,"BDTG");
}

TMVA::Reader * reader;
};

MVAWrapper w;
int main(int argv,const char ** argc){
        std::cout <<  w.eval({0,0,0,0,0,0,0,0,0,0}) << std::endl;
}
