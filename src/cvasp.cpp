/*
 * cvasp.h
 *
 * Copyright 2018 Edmanuel <eetorres@gmail.com>
 *
 * MIT License
 *
 * Copyright (c) 2018 Edmanuel Torres
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 *
 */

#include <cvasp.h>
#include <ctype.h>
#include <algorithm>

#define _VASP_DEBUG_MESSAGES_
#define _VASP_DEBUG_DATA_
//#define _POSCAR_SHOW_DATA_

CVasp::CVasp()
{
    is_periodic = false;
}

bool CVasp::GetFormat(std::string d,std::string f)
{
    int res1;
    double c[9];
    std::string s;
    is_periodic = false;
    is_extended = false;
    static int i1 = 0, i2 = 0, i3 = 0;
    double a1 = 0, a2 = 0, a3 = 0;
    char s1, s2, s3;
    char * pch;
    char * pchar;
    //char s1[2],s2[2],s3[2];
    char buff[256];
    std::string text_line;
    total_atoms = 0;
    header_lines=0;
    std::ifstream poscar;
    v_table_index.resize(0);
    v_table_atoms.resize(0);
    try
    {
        poscar.open(f.c_str());
        if(!poscar.is_open())
        {
#ifdef _POSCAR_DEBUG_MESSAGES_
            setvbuf(stdout, NULL, _IONBF, 0);
            DEBUG_VALUE("CPOSCAR: The POSCAR file: ",f)
            DEBUG_MESSAGE("CPOSCAR: was not found")
            std::cout<< std::flush;
#endif
            return false;
        }
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_MESSAGE("<-- The POSCAR was successfuly opened -->")
        DEBUG_MESSAGE("<--       Analysing the POSCAR        -->")
#endif
        while((res1<3 || header_lines < 8) && !poscar.eof())
        {
            //poscar.ignore(256,'\n');
            std::getline(poscar,text_line,'\n');
            strcpy (header_info[header_lines],text_line.c_str());
#ifdef _POSCAR_SHOW_DATA_
            std::cout<<"POSCAR: [header "<<header_lines<<"] "<<text_line<<std::endl;
#endif
            res1 = std::sscanf((const char*)text_line.c_str(),"%lf %lf %lf %c %c %c",&a1,&a2,&a3,&s1,&s2,&s3);
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_VALUE(" res1 =",res1)
#endif
            header_lines++;
        }
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_VALUE("POSCAR: total header lines ",header_lines)
#endif
        ////////////////////////////
        // Check if extended POSCAR
        ////////////////////////////
        strcpy(buff,header_info[5]);
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_VALUE(" Checking Extended = ",buff)
#endif
        res1 = std::sscanf(buff,"%d %d %d",&i1,&i2,&i3);
        //std::cout <<" res1 = "<<res1<<std::endl;
        if(res1 == 0)
            is_extended = true;
        pch = strtok(buff," ");
        while (pch != NULL)
        {
            if(sizeof(pch) != 0)
            {
                if(is_extended)
                {
#ifdef _VASP_DEBUG_MESSAGES_
                    DEBUG_VALUE(" Index = ",GetAtomIndex(pch))
#endif
                    v_table_index.push_back(GetAtomIndex(pch));
                }
                else
                {
#ifdef _VASP_DEBUG_MESSAGES_
                    DEBUG_VALUE(" Atoms = ",pch)
#endif
                    v_table_atoms.push_back(atoi(pch));
                }
            }
            pch = strtok (NULL," ");
        }
        if(is_extended)
        {
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_MESSAGE(" Checking Composition")
#endif
            strcpy(buff,header_info[6]);
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_VALUE(" Composition = ",buff)
#endif
            pch = strtok(buff," ");
            while (pch != NULL)
            {
                v_table_atoms.push_back(atoi(pch));
                pch = strtok (NULL," ");
            }
        }
        else
        {
            ReadPotcar(d);
        }
#ifdef _VASP_DEBUG_MESSAGES_
        if(is_extended)
            DEBUG_MESSAGE(" Extended POSCAR")
            else
                DEBUG_MESSAGE(" No Extended POSCAR")
#endif
                //////////////////////////
                // check for dynamic type
                //////////////////////////
                if(is_extended)
                    strcpy(buff,header_info[7]);
                else
                    strcpy(buff,header_info[6]);
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_VALUE(" Dynamic = ",buff)
#endif
        if((buff[0]=='S') || (buff[0]=='s'))
        {
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_MESSAGE(" Selective Dynamics")
#endif
            is_dynamics = true;
        }
        else
        {
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_MESSAGE(" No Selective Dynamics")
#endif
            is_dynamics = false;
        }
        /////////////////////////////
        // Check for coordinate type
        /////////////////////////////
        if(is_extended && is_dynamics)
            strcpy(buff,header_info[8]);
        else if(is_extended || is_dynamics)
            strcpy(buff,header_info[7]);
        else
            strcpy(buff,header_info[6]);
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_VALUE(" Coordinates = ",buff)
#endif
        if((buff[0]=='D') || (buff[0]=='d'))
        {
            is_direct = true;
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_MESSAGE(" Direct coordinates")
#endif
        }
        else
        {
            is_direct = false;
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_MESSAGE(" Cartesian coordinates")
#endif
        }
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_MESSAGE(" Read Cell")
#endif
        strcpy(buff,header_info[1]);
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_VALUE(" Cell parameter = ",buff)
#endif
        parameter = atof(buff);
        for(int i=2; i<5; i++)
        {
            strcpy(buff,header_info[i]);
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_VALUE(" Lattice vector = ",buff)
            // a1 = 0, a2 = 0, a3 = 0;
#endif
            std::sscanf((const char*)buff,"%lf %lf %lf",&a1,&a2,&a3);
            cell_xyz[i-2] = parameter*glm::dvec3(a1,a2,a3);
        }
        is_periodic = true;
        for (std::vector<uint>::iterator it = v_table_atoms.begin(); it != v_table_atoms.end(); ++it)
        {
            total_atoms += *it;
        }
        //////////////////////////////////////////////////////////////////////
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_MESSAGE(" =======Summary======")
        print_vector("v_table_index",v_table_index);
        print_vector("v_table_atoms",v_table_atoms);
        DEBUG_VALUE(" cell parameter = ",parameter)
        DEBUG_VEC3("cell=",glm::dvec3(cell_xyz[0]));
        DEBUG_VEC3("cell=",glm::dvec3(cell_xyz[1]));
        DEBUG_VEC3("cell=",glm::dvec3(cell_xyz[2]));
        DEBUG_MESSAGE(" =======Summary======")
#endif
        //////////////////////////////////////////////////////////////////////
        if(poscar.is_open())
        {
            poscar.close();
        }
        if((v_table_index.size() == v_table_atoms.size()) && parameter > 0)
            return true;
        else
            return false;
    }
    catch(...)
    {
        return false;
    }
}

// load the POTCAR file
bool CVasp::ReadPotcar(std::string d)
{
    char tmp_buffer[256], *pch;
    std::string atom_symbol;
    std::ifstream potcar;
    // Reading atomic species
    std::string potcar_file = d+"/POTCAR";
    potcar.open(potcar_file.c_str());
    if(!potcar.is_open())
    {
        return false;
    }
    while(!potcar.eof())
    {
        potcar.getline(tmp_buffer,256);
        pch=strstr(tmp_buffer,"VRHFIN");
        if(pch!=NULL)
        {
            pch=strrchr(tmp_buffer,'=');
            pch = strtok (pch," =:");
            atom_symbol = pch;
#ifdef _VASP_DEBUG_MESSAGES_
            DEBUG_VALUE(" VRHFIN = ",atom_symbol)
#endif
            v_table_index.push_back(GetAtomIndex(pch));
            //__potcar_species++;
            //v_atomic_symbols.push_back(atom_symbol);
        }
    }
    if(!potcar.is_open())
    {
        return false;
    }
    potcar.close();
    return false;
}

// load the VASP file
bool CVasp::ReadFile(std::string d,std::string f,std::vector<glm::dvec3>& xyz, std::vector<uint>& idx)
{
    std::ifstream poscar;
    double a1, a2, a3;
    char s1, s2, s3;
    std::string s;
    char vchr[3];
    bool compare = false;
    bool res = true;
    std::string str;
    std::string text_line;
    glm::dvec3 fxyz;
    //total_atoms = 0;
    res = GetFormat(d,f);
#ifdef _VASP_DEBUG_MESSAGES_
    DEBUG_VALUE("CVASP: RES: ",res)
#endif
    if(res)
    {
        poscar.open(f.c_str());
        std::getline(poscar,str,'\n');
        std::sscanf((const char*)str.c_str(),"%i",&total_atoms);
        //header_lines=1;
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_VALUE("CVASP: Total atoms: ",total_atoms)
#endif
        for(uint i=0; i<header_lines-2; i++)
        {
            //poscar.ignore(1024,'\n');
            std::getline(poscar,str,'\n');
        }
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_MESSAGE("CVASP: loading atomic coordinates")
#endif
        xyz.resize(total_atoms);
        idx.resize(total_atoms);
        v_atom_type.resize(total_atoms);
        v_fract_xyz.resize(total_atoms);
        xyz.resize(total_atoms);
        if(is_dynamics)
            v_dynamics.resize(total_atoms);
        uint cont=0;
        if(v_table_index.size() > 0)
        {
            for (uint j=0; j<v_table_atoms.size(); j++)
            {
                for (uint i=0; i<v_table_atoms[j]; i++)
                {
                    idx[cont]=v_table_index[j];
                    v_atom_type[cont]=j;
                    cont++;
                }
            }
        }
        for(uint f=0; f<total_atoms; f++)
        {
            std::getline(poscar,str,'\n');
            if(is_dynamics)
                std::sscanf((const char*)str.c_str(),"%lf %lf %lf %c %c %c",&a1,&a2,&a3,&s1,&s2,&s3);
            else
                std::sscanf((const char*)str.c_str(),"%lf %lf %lf",&a1,&a2,&a3);
            fxyz.x=a1;
            fxyz.y=a2;
            fxyz.z=a3;
            if(is_direct)
            {
                v_fract_xyz[f] = fxyz;
                xyz[f] = cell_xyz*fxyz;
            }
            else
            {
                xyz[f] = fxyz;
            }
            if(is_dynamics)
            {
                v_dynamics[f].d[0] = s1;
                v_dynamics[f].d[1] = s2;
                v_dynamics[f].d[2] = s3;
            }
            //idx[f] = 0;
#ifdef _VASP_DEBUG_DATA_
            //DEBUG_VALUE(" XYZ = ",str)
            std::cout<<idx[f]<<" "<<xyz[f].x<<" "<<xyz[f].y<<" "<<xyz[f].z<<" "<<v_atom_type[f]<<std::endl;
            if(is_dynamics) std::cout<<" "<<s1<<" "<<s2<<" "<<s3<<std::endl;
#endif
            if(compare)
            {
                SetMinMax(xyz[f]);
            }
            else
            {
                min_coord = xyz[f];
                max_coord = xyz[f];
                compare = true;
            }
        }
        if(poscar.is_open())
        {
            poscar.close();
        }
        //std::cout<<"min = "<<glm::to_string(min_coord)<<std::endl;
        //std::cout<<"max = "<<glm::to_string(max_coord)<<std::endl;
#ifdef _VASP_DEBUG_MESSAGES_
        DEBUG_MESSAGE("CVASP: The VASP was successful closed!")
#endif
        return true;
    }
    else
    {
        return false;
    }
}

// load the VASP file
bool CVasp::SaveFile(std::string d, std::string f, std::vector<glm::dvec3>& xyz, std::vector<uint>& idx)
{
    std::ofstream oxyz;
    uint tidx;
    bool is_new;
    glm::mat3 inv_cell = glm::inverse(cell_xyz);
    std::cout<<glm::to_string(cell_xyz)<<std::endl;
    std::cout<<glm::to_string(inv_cell)<<std::endl;
    long double ldtmp;
#ifdef _VASP_DEBUG_DATA_
    print_vector("idx",idx);
#endif
    v_table_index.clear();
    v_table_index.push_back(idx[0]);
    for( uint i=0; i<idx.size(); i++)
    {
        is_new = true;
        for(uint j=0; j<v_table_index.size(); j++)
        {
            tidx = v_table_index[j];
            if( tidx == idx[i] )
            {
                v_atom_type[i]=j;
                is_new = false;
#ifdef _VASP_DEBUG_MESSAGES_
                DEBUG_VALUE("idx found: ",tidx)
#endif
            }
        }
        if(is_new)
        {
            v_table_index.push_back(idx[i]);
            v_atom_type[i]=v_table_index.size()-1;
        }
    }
    v_table_atoms.resize(v_table_index.size());
    for(uint j=0; j<v_table_index.size(); j++)
    {
        tidx = v_table_index[j];
        v_table_atoms[j]=0;
        for( uint i=0; i<idx.size(); i++)
        {
            if( tidx == idx[i] )
            {
                v_table_atoms[j]++;
            }
        }
    }
#ifdef _VASP_DEBUG_DATA_
    print_vector("v_table_index",v_table_index);
    print_vector("v_table_atoms",v_table_atoms);
    print_vector("v_atom_type",v_atom_type);
#endif
    try
    {
        oxyz.open(f.c_str());
        //oxyz<<xyz.size()<<std::endl<<std::endl;
        oxyz.precision(15);
        //for(uint i=0; i<5; i++)
        //{
            //oxyz<<header_info[i]<<std::endl;
        //}
        oxyz<<header_info[0]<<std::endl;
        oxyz<<parameter<<std::endl;
        for(int i=0; i<3; i++)
        {
                glm::dvec3 cxyz = cell_xyz[i];
                oxyz<<"    ";
                //oxyz<<oxyz.width(4);
                oxyz<<std::fixed<<std::right<<cxyz.x;
                oxyz<<"   ";
                //oxyz<<oxyz.width(20);
                oxyz<<std::fixed<<std::right<<cxyz.y;
                oxyz<<"   ";
                //oxyz<<oxyz.width(4);
                oxyz<<std::fixed<<std::right<<cxyz.z<<std::endl;
        }
        for (uint i=0; i<v_table_index.size(); i++)
        {
            oxyz<<symbol[v_table_index[i]];
            if(i<v_table_index.size()-1) oxyz<<" ";
            else oxyz<<std::endl;
        }
        for (uint i = 0; i<v_table_atoms.size(); i++)
        {
            oxyz<<v_table_atoms[i];
            if(i<v_table_atoms.size()-1) oxyz<<" ";
            else oxyz<<std::endl;
        }
        if(is_dynamics)
        {
            oxyz<<"Selective Dynamics"<<std::endl;
        }
        if(is_direct)
        {
            oxyz<<"Direct coordinates"<<std::endl;
        }
        else
        {
            oxyz<<"Cartesian coordinates"<<std::endl;
        }
        for (uint i=0; i<v_table_index.size(); i++)
        {
            for(uint u=0; u<xyz.size(); u++)
            {
                if(v_atom_type[u] == i)
                {
                    if(is_direct)
                    {
                        glm::dvec3 vxyz = v_fract_xyz[u]; //inv_cell*xyz[u];
                        oxyz.width(20);
                        //oxyz<<std::fixed<<std::right<<v_fract_xyz[u].x;
                        oxyz<<std::fixed<<std::right<<vxyz.x;
                        oxyz.width(20);
                        //oxyz<<std::fixed<<std::right<<v_fract_xyz[u].y;
                        oxyz<<std::fixed<<std::right<<vxyz.y;
                        oxyz.width(20);
                        //oxyz<<std::fixed<<std::right<<v_fract_xyz[u].z;
                        oxyz<<std::fixed<<std::right<<vxyz.z;
                        if(is_dynamics)
                        {
                            oxyz<<"   "<<v_dynamics[u].d[0]<<"   "<<v_dynamics[u].d[1]<<"   "<<v_dynamics[u].d[2];
                        }
                        oxyz<<std::endl;
                    }
                    else
                    {
                        oxyz.width(17);
                        oxyz<<std::fixed<<std::right<<xyz[u].x;
                        oxyz.width(17);
                        oxyz<<std::fixed<<std::right<<xyz[u].y;
                        oxyz.width(17);
                        oxyz<<std::fixed<<std::right<<xyz[u].z;
                        if(is_dynamics)
                        {
                            oxyz<<" "<<v_dynamics[u].d[0]<<" "<<v_dynamics[u].d[1]<<" "<<v_dynamics[u].d[2];
                        }
                        oxyz<<std::endl;
                    }
                }
            }
        }
        if(oxyz.is_open())
        {
            oxyz.close();
        }
#ifdef _VASP_DEBUG_MESSAGES_
        std::cout<<"CVASP: The VASP was writen!"<<std::endl;
#endif
        return true;
    }
    catch(...)
    {
        return false;
    }
}

void CVasp::StrainCell(double s, int i) // Working here
{
    cell_xyz[i]*=s;
};

uint CVasp::GetAtomIndex(std::string s)
{
    uint atnum;
    size_t _sw;
    transform(s.begin(), s.begin()+1, s.begin(), toupper);
    transform(s.begin()+1, s.end(), s.begin()+1, tolower);
    for(uint k=0; k<periodic_table_atoms; k++)
    {
        if(k < 95)
        {
            _sw = 2;
        }
        else
        {
            _sw = 1;
        }
        if(strncmp(s.c_str(),symbol_sorted[k].c_str(),_sw)==0)
        {
            atnum = atom_index[k];
            break;
        }
    }
    return atnum;
}

void CVasp::DeleteAtom(uint u)
{
    //xyz.erase(xyz.begin()+u);
    v_atom_type.erase(v_atom_type.begin()+u);
    if(is_dynamics)
    {
        v_dynamics.erase(v_dynamics.begin()+u);
    }
    if(is_direct)
    {
        v_fract_xyz.erase(v_fract_xyz.begin()+u);
    }
}

inline void CVasp::SetMinMax(glm::dvec3 &v)
{
    min_coord.x = std::min(min_coord.x,v.x);
    min_coord.y = std::min(min_coord.y,v.y);
    min_coord.z = std::min(min_coord.z,v.z);
    max_coord.x = std::max(max_coord.x,v.x);
    max_coord.y = std::max(max_coord.y,v.y);
    max_coord.z = std::max(max_coord.z,v.z);
}

glm::dvec3 CVasp::GetMinCoordinate(void)
{
    return min_coord;
}

glm::dvec3 CVasp::GetMaxCoordinate(void)
{
    return max_coord;
}

// END

