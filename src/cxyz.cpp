/*
 * cxyz.h
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

#include <cxyz.h>
#include <ctype.h>
#include <algorithm>

//#define GLM_ENABLE_EXPERIMENTAL
//#include <glm/gtx/string_cast.hpp>
//#include <glm/gtc/type_ptr.hpp>

//#define _XYZ_DEBUG_MESSAGES_
//#define _XYZ_DEBUG_DATA_

// formats
// 1 standard xyz (i.e., avogadro)
// 2 standard including fragments
// 3 standard without number of atoms

// 4 with fragments without number of atoms

CXyz::CXyz()
{
    is_periodic = false;
}

bool CXyz::GetFormat(std::string d,std::string f)
{
    std::ifstream ixyz;
    double c[9];
    std::string s;
    is_periodic = false;
    char buff[200];
    //header_lines=2;
    try
    {
        ixyz.open(f.c_str());
        std::string str;
        //std::getline(ixyz, str); // the first line is the number of atoms
        std::getline(ixyz,str,'\n');
        std::sscanf((const char*)str.c_str(),"%i",&total_atoms);
#ifdef _XYZ_DEBUG_MESSAGES_
        std::cout<<"CXYZ: Total atoms: "<<total_atoms<<std::endl;
#endif
        std::getline(ixyz, str); // Parameters of the extended format
        std::size_t start_found = str.find("Cell");
        if (start_found!=std::string::npos)
        {
            //std::cout << "found at: " << start_found << '\n';
            start_found = str.find("\"");
            //std::cout << "start: " << ++start_found << '\n';
            std::size_t final_found = str.find_last_of("\"");
            final_found = final_found-start_found;
            //std::cout << "final: " << final_found << '\n';
            std::string str2 = str.substr (++start_found,--final_found);
#ifdef _XYZ_DEBUG_MESSAGES_
            std::cout << "Cell = " << str2 << '\n';
#endif
            strcpy(buff,str2.c_str());
            char * pch;
            pch = strtok (buff," ");
            int cont=0;
            //while (pch != NULL)
            for(int i=0; i<9; i++)
            {
                //printf ("%i %s\n",cont,pch);
                c[i] = atof(pch);
                pch = strtok (NULL, " ");
            }
            cell_xyz = glm::make_mat3(c);
            is_periodic = true;
            //DEBUG_VEC3("cell=",glm::dvec3(cell_xyz[0]));
        }
        if(ixyz.is_open())
        {
            ixyz.close();
        }
        if(total_atoms > 0)
            return true;
        else
            return false;
    }
    catch(...)
    {
        return false;
    }
}

// load the XYZ file
bool CXyz::ReadFile(std::string d,std::string f,std::vector<glm::dvec3>& xyz, std::vector<uint>& idx)
{
    std::ifstream ixyz;
    double a, b, c;
    std::string s;
    char vchr[3];
    bool compare = false;
    bool res = true;
    std::string str;
    std::string text_line;
    total_atoms = 0;
    glm::dvec3 fxyz;
    res = GetFormat(d,f);
#ifdef _XYZ_DEBUG_MESSAGES_
    DEBUG_VALUE("XYZ: RES: ",res)
#endif
    if(!res)
        return false;
    try
    {
        ixyz.open(f.c_str());
        std::getline(ixyz,str,'\n');
        //std::sscanf((const char*)str.c_str(),"%i",&total_atoms);
        header_lines=1;
#ifdef _XYZ_DEBUG_MESSAGES_
        std::cout<<"CXYZ: Total atoms: "<<total_atoms<<std::endl;
#endif
        for(uint i=0; i<header_lines; i++)
        {
            //ixyz.ignore(1024,'\n');
            std::getline(ixyz,str,'\n');
        }
#ifdef _XYZ_DEBUG_MESSAGES_
        std::cout<<"CXYZ: loading atomic coordinates"<<std::endl;
#endif
        xyz.resize(total_atoms);
        idx.resize(total_atoms);
        for(uint f=0; f<total_atoms; f++)
        {
            std::getline(ixyz,str,'\n');
            std::sscanf((const char*)str.c_str(),"%s %lf %lf %lf",vchr,&a,&b,&c);
            fxyz.x=a;
            fxyz.y=b;
            fxyz.z=c;
            xyz[f] = fxyz;
            s = vchr;
#ifdef _XYZ_DEBUG_DATA_
            std::cout<<s<<" "<<xyz[f].x<<" "<<xyz[f].y<<"  "<<xyz[f].z<<std::endl;
#endif
            if(sizeof(s) != 0)
            {
                idx[f] = GetAtomIndex(s);
            }
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
        if(ixyz.is_open())
        {
            ixyz.close();
        }
        //std::cout<<"min = "<<glm::to_string(min_coord)<<std::endl;
        //std::cout<<"max = "<<glm::to_string(max_coord)<<std::endl;
#ifdef _XYZ_DEBUG_MESSAGES_
        std::cout<<"CXYZ: The XYZ was successful closed!"<<std::endl;
#endif
        return true;
    }
    catch(...)
    {
        return false;
    }
}

// load the XYZ file
bool CXyz::SaveFile(std::string d, std::string f, std::vector<glm::dvec3>& xyz, std::vector<uint>& idx)
{
    std::ofstream oxyz;
    try
    {
        oxyz.open(f.c_str());
        oxyz<<xyz.size()<<std::endl<<std::endl;
        oxyz.precision(12);
        for(uint u=0; u<xyz.size(); u++)
        {
            oxyz.width(3);
            oxyz<<std::fixed<<std::left<<symbol[idx[u]];
            oxyz.width(20);
            oxyz<<std::fixed<<std::right<<xyz[u].x;
            oxyz.width(20);
            oxyz<<std::fixed<<std::right<<xyz[u].y;
            oxyz.width(20);
            oxyz<<std::fixed<<std::right<<xyz[u].z<<std::endl;
        }
        if(oxyz.is_open())
        {
            oxyz.close();
        }
#ifdef _XYZ_DEBUG_MESSAGES_
        std::cout<<"CXYZ: The XYZ was writen!"<<std::endl;
#endif
        return true;
    }
    catch(...)
    {
        return false;
    }
}

uint CXyz::GetAtomIndex(std::string s)
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

inline void CXyz::SetMinMax(glm::dvec3 &v)
{
    min_coord.x = std::min(min_coord.x,v.x);
    min_coord.y = std::min(min_coord.y,v.y);
    min_coord.z = std::min(min_coord.z,v.z);
    max_coord.x = std::max(max_coord.x,v.x);
    max_coord.y = std::max(max_coord.y,v.y);
    max_coord.z = std::max(max_coord.z,v.z);
}

uint CXyz::GetTotalAtoms(void)
{
    return total_atoms;
}

glm::dvec3 CXyz::GetMinCoordinate(void)
{
    return min_coord;
}

glm::dvec3 CXyz::GetMaxCoordinate(void)
{
    return max_coord;
};


// END

