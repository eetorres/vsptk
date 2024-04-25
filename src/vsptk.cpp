#include <cmodel.h>
#include <algorithm>
#include <string>
#include <iostream>

class InputParser{
    public:
        //
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        // Get the argument of the option
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }
        // Chekc if the option exist
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

int main(int argc, char **argv){
  // Initialize the model class
  CModel *md = new CModel();
  md->ReadStructureFile("./","./POSCAR",INPUT_FILE_TYPE_VSP);
  
  InputParser input(argc, argv);
  if(input.cmdOptionExists("-h")){
    std::cout<<"help"<<std::endl;
  }
  // Apply strain in X
  const std::string &strainx = input.getCmdOption("-sx");
  if (!strainx.empty()){
    std::cout<<"strain in X = "<<strainx<<std::endl;
    md->vasp.is_dynamics=false;
    md->vasp.StrainCell(std::stod(strainx),0);
    //md->SaveStructureFile("./","./POSCAR-RES",INPUT_FILE_TYPE_VSP);
  }
  // Apply strain in Y
  const std::string &strainy = input.getCmdOption("-sy");
  if (!strainy.empty()){
    std::cout<<"strain in Y = "<<strainy<<std::endl;
    md->vasp.is_dynamics=false;
    md->vasp.StrainCell(std::stod(strainy),1);
    //md->SaveStructureFile("./","./POSCAR-RES",INPUT_FILE_TYPE_VSP);
  }
  // Apply strain in Z
  const std::string &strainz = input.getCmdOption("-sz");
  if (!strainz.empty()){
    std::cout<<"strain in Z = "<<strainx<<std::endl;
    md->vasp.is_dynamics=false;
    md->vasp.StrainCell(std::stod(strainz),2);
    //md->SaveStructureFile("./","./POSCAR-RES",INPUT_FILE_TYPE_VSP);
  }
  
  md->SaveStructureFile("./","./POSCAR-VTK",INPUT_FILE_TYPE_VSP);

  return 0;
}

