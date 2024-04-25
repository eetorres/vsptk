/*
 * cpdb.h
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

#include <cpdb.h>
#include <ctype.h>
#include <algorithm>

//#define GLM_ENABLE_EXPERIMENTAL
//#include <glm/gtx/string_cast.hpp>
//#include <glm/gtc/type_ptr.hpp>

//#define _PDB_DEBUG_MESSAGES_
//#define _PDB_DEBUG_DATA_

CPdb::CPdb()
{
    is_periodic = false;
}

bool CPdb::GetFormat(std::string d,std::string f)
{
    int res1, res2, res3, res4, res5, res6, res7;
    int k, l;
    float e, g, h;
    //float a, b, c, d;
    //float a, b, c, e;
    //float alpha, beta, gamma;
    header_info.clear();
    char str[80], stra[80], strb[80], str2[80];
    std::string text_line;
    std::ifstream pdb;
    str2[0]='X';
    is_periodic = false;
    try
    {
        pdb.open(f);
        if(!pdb.is_open())
        {
#ifdef _PDB_DEBUG_MESSAGES_
            std::cout<<"CPDB: The PDB file: "<<f<<std::endl;
            std::cout<<"CPDB: was not found\n"<<std::endl;
#endif
            return false;
        }
        res1=0;
        res2=0;
        header_lines=0;
        while(strcmp("ATOM",str2) && strcmp("HETATM",str2) && !pdb.eof())
        {
            //pdb.ignore(256,'\n');
            std::getline(pdb,text_line);
            header_info.push_back(text_line);
            res1 = std::sscanf((const char*)text_line.c_str(),"%s",str2);
            if(!strcmp("CRYST1",str2))
            {
#ifdef _PDB_DEBUG_MESSAGES_
                std::cout<<" CRYST1 found"<<std::endl;
                std::cout<<" res1: "<<res1<<std::endl;
#endif
                res1 = std::sscanf((const char*)text_line.c_str(),"%*s %f %f %f %f %f %f",&a,&b,&c,&alpha,&beta,&gamma);
#ifdef _PDB_DEBUG_MESSAGES_
                std::cout<<" (a,b,c) =("<<a<<","<<b<<","<<c<<")"<<std::endl;
                std::cout<<" (alpha,beta,gamma) =("<<alpha<<","<<beta<<","<<gamma<<")"<<std::endl;
#endif
                if ( (a > 0) && (b > 0) && (c > 0) )
                {
                    // Compute lattice vectors coordinates in a cartesian frame.
                    // Vector a is along the x axes, vector b is in the (x, y) plane.
                    // conversion en radian
                    is_periodic = true;
                    float alphar = glm::radians(alpha);
                    //float betar  = beta  * DEG_RAD;
                    float gammar = glm::radians(gamma);
                    // calcul des vecteurs
                    // lattice vector 1 //
                    cell_xyz[0][0] = a;
                    // lattice vector 2 //
                    if(gammar != 90)
                    {
                        cell_xyz[1][0] = b*cos(gammar);
                        cell_xyz[1][1] = b*sin(gammar);
                    }
                    else
                    {
                        cell_xyz[1][1] = b;
                    }
                    // lattice vector 3
                    if(beta != 90 || gamma != 90)
                    {
                        cell_xyz[2][0] = c*cos(alphar)*cos(gammar);
                        cell_xyz[2][1] = c*cos(alphar)*sin(gammar);
                        cell_xyz[2][2] = c*sin(alphar);
                    }
                    else
                    {
                        cell_xyz[2][2] = c;
                    }
                    /*
                    self._vecc = [0.,0.,0.]
                    self._vecc[0] = self._c * cos( betar )
                    cy = ( cos(alphar) - cos(gammar) * cos(betar) ) / sin(gammar)
                    self._vecc[1] = self._c * cy
                    cz = sqrt( ( sin( betar ) )**2 - cy**2 )
                    self._vecc[2] = self._c * cz
                    */
                }
            }
            header_lines++;
        }
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<"PDB: header lines: ["<<header_lines<<"]"<<std::endl;
        DEBUG_VALUE(" Text line =",text_line)
#endif
        res1 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %f %f %f",&k,str,&a,&b,&c);
        res2 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %i %f %f %f",&k,str,&l,&a,&b,&c);
        res3 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %i %f %f %f %f",&k,str,&l,&a,&b,&c,&e);
        res4 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %i %f %f %f %f %f",&k,str,strb,&l,&a,&b,&c,&e,&g);
        res5 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %i %f %f %f %f %f %f",&k,str,stra,strb,&l,&a,&b,&c,&e,&g,&h);
        res6 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %i %f %f %f %f %f %f",&k,str,stra,strb,&l,&a,&b,&c,&e,&g,&h);
        res7 = std::sscanf((const char*)text_line.c_str(),"%*s %i %s %*s %*s %i %f %f %f %f %f %f",&k,str,stra,strb,&l,&a,&b,&c,&e,&g,&h);
        //__header_lines++;
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<" res1: "<<res1<<std::endl;
        std::cout<<" res2: "<<res2<<std::endl;
        std::cout<<" res3: "<<res3<<std::endl;
        std::cout<<" res4: "<<res4<<std::endl;
        std::cout<<" res5: "<<res5<<std::endl;
        std::cout<<" res6: "<<res6<<std::endl;
        std::cout<<" res7: "<<res7<<std::endl;
#endif
        /*    if(res3==6){                 // standard format
              if(res2 == 2 && res3 == 2) pdb_format = 1;
              if(res2 == 6 && res3 == 5) pdb_format = 2;
            }else if(res3==7){
              if(res2 == 6){
                pdb_format = 2;
              }else if(res4==8){
                pdb_format = 4;
                is_charges=true;
              }
            }else if(res3==2 && res2 == 6 && res1 == 5){           // standard format
              pdb_format = 3;
            }else if(res5==9){           // full format
              pdb_format = 5;
        */

        if(res5==9)            // full format
        {
            pdb_format = 5;
        }
        else if(res4==8)
        {
            pdb_format = 4;
            is_charges=true;
        }
        else if(res3==7)
        {
            pdb_format = 3;
        }
        else if(res2 == 6)    // standard format
        {
            pdb_format = 2;
        }
        else if(res1 == 5)
        {
            pdb_format = 1;
        }
        else if(res6 >= 6)
        {
            pdb_format = 6;
        }
        else if(res7 >= 8)
        {
            pdb_format = 7;
        }
        else
        {
            return false;
        }
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<"CPDB: file format: "<<pdb_format<<std::endl;
        std::cout<<"CPDB: counting atoms"<<std::endl;
#endif
        total_atoms=0;
        while((!strcmp("ATOM",str2) || !strcmp("HETATM",str2)) && (res1>0) && !pdb.eof())
        {
            total_atoms++;
            std::getline(pdb,text_line);
            str2[0]='X';
            res1 = std::sscanf((const char*)text_line.c_str(),"%s",str2);
        }
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<" total atoms: "<<total_atoms<<std::endl;
        std::cout<<" res1: "<<res1<<std::endl;
#endif

        if(pdb.is_open())
            pdb.close();
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<"CPDB: The PDB file was successful closed"<<std::endl;
#endif
        if(total_atoms>0)
            return true;
        else
            return false;
    }
    catch(...)
    {
        return false;
    }
}

// load the PDB file
bool CPdb::ReadFile(std::string d,std::string f,std::vector<glm::dvec3>& xyz, std::vector<uint>& idx)
{
    std::ifstream ipdb;
    float x, y, z, e;
    std::string s;
    int k, l;
    int res1;
    char vchr[3];
    char str2[80];
    bool compare = false;
    bool res = true;
    std::string str;
    std::string text_line;
    total_atoms = 0;
    glm::vec3 fxyz;
    res = GetFormat(d,f);
#ifdef _PDB_DEBUG_MESSAGES_
    DEBUG_VALUE("CPDB: RES: ",res)
#endif
    if(!res)
        return false;
    try
    {
        ipdb.open(f.c_str());
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<"CPDB: Total Atoms: "<<total_atoms<<std::endl;
#endif
        if(header_info.size()>1)
        {
            for(uint i=0; i<header_info.size()-1; i++)
            {
                //ipdb.ignore(1024,'\n');
                std::getline(ipdb,str,'\n');
                //std::getline(ipdb,str,'\n');
            }
        }
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<"CPDB: loading atomic coordinates"<<std::endl;
#endif
        idx.resize(total_atoms);
        xyz.resize(total_atoms);
        for(uint f=0; f<total_atoms; f++)
        {
            std::getline(ipdb,str,'\n');
            res1 = std::sscanf((const char*)str.c_str(),"%s",str2);
#ifdef _PDB_DEBUG_MESSAGES_
            DEBUG_VALUE("CPDB: RES1: ",res1)
#endif
            if( (!strcmp("ATOM",str2) || !strcmp("HETATM",str2)) && !ipdb.eof() )
            {
#ifdef _PDB_DEBUG_MESSAGES_
                DEBUG_VALUE(" Text line =",str)
#endif
                //std::sscanf((const char*)str.c_str(),"%s %lf %lf %lf",vchr,&a,&b,&c);
                switch(pdb_format)
                {
                case 1:
                    std::sscanf((const char*)str.c_str(),"%*s %i %s %f %f %f",&k,vchr,&x,&y,&z);
                    break;
                case 2:
                case 4:
                    std::sscanf((const char*)str.c_str(),"%*s %i %s %i %f %f %f %*s",&k,vchr,&l,&x,&y,&z);
                    break;
                //std::sscanf((const char*)str.c_str(),"%*s %i %s %*s %*s %f %f %f %f",&k,vchr,&x,&y,&z,&e);
                //break;
                case 5:
                    std::sscanf((const char*)str.c_str(),"%*s %i %s %*s %*s %*s %f %f %f %f",&k,vchr,&a,&b,&c,&e);
                    break;
                case 6:
                    std::sscanf((const char*)str.c_str(),"%*s %i %s %*s %i %f %f %f %*s",&k,vchr,&l,&x,&y,&z);
                    break;
                case 7:
                    std::sscanf((const char*)str.c_str(),"%*s %i %s %*s %*s %i %f %f %f %*s",&k,vchr,&l,&x,&y,&z);
                    break;
                }
#ifdef _PDB_SHOW_DATA_
                std::cout<<"CPDB: "<<k<<" str="<<str<<" a="<<a<<" b="<<b<<" c="<<c<<" e="<<e<<std::endl;
#endif
                fxyz.x=x;
                fxyz.y=y;
                fxyz.z=z;
                xyz[f] = fxyz;
                s = vchr;
#ifdef _PDB_DEBUG_DATA_
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
        }
        if(ipdb.is_open())
        {
            ipdb.close();
        }
#ifdef _PDB_DEBUG_DATA_
        DEBUG_MESSAGE(" =======PDB Summary======")
        print_vector("index",idx);
        //DEBUG_VEC3("atoms",xyz);
        //DEBUG_VEC3("xyz",xyz);
        //DEBUG_VALUE(" cell parameter = ",parameter)
        DEBUG_VEC3("cell=",cell_xyz[0]);
        DEBUG_VEC3("cell=",cell_xyz[1]);
        DEBUG_VEC3("cell=",cell_xyz[2]);
        std::cout<<"min = "<<glm::to_string(min_coord)<<std::endl;
        std::cout<<"max = "<<glm::to_string(max_coord)<<std::endl;
        DEBUG_MESSAGE(" =======Summary======")
        std::cout<<"CPDB: The PDB was successful closed!"<<std::endl;
#endif
        return true;
    }
    catch(...)
    {
        return false;
    }
}

// load the PDB file
bool CPdb::SaveFile(std::string d,std::string f, std::vector<glm::dvec3>& xyz, std::vector<uint>& idx)
{
    std::ofstream opdb;
    try
    {
        opdb.open(f.c_str());
        if(header_info.size()>1)
        {
            for(uint i=0; i<header_info.size()-1; i++)
            {
                opdb<<header_info[i]<<std::endl;
            }
        }
        else
        {
            std::string stime = Timer.GetDate();
            opdb<<"HEADER PDB structure file generated by XMolView (www.xmol.org) ["<<stime<<"]"<<std::endl;
        }
        for(uint u=0; u<xyz.size(); u++)
        {
            opdb<<"ATOM";
            opdb.width(7);
            opdb<<std::fixed<<std::right<<u+1<<"  ";
            opdb.width(4);
            opdb<<std::fixed<<std::left<<symbol[idx[u]]<<"LIG     1   ";
            opdb.precision(3);
            for(uint c=0; c<3; c++)
            {
                opdb.width(8);
                opdb<<std::fixed<<std::right<<xyz[u][c];;
            }
            //opdb<<"  "; //1.00";
            //if(__is_charges){
            //opdb.width(9);
            //opdb.precision(4);
            //opdb<<std::fixed<<std::right<<v_atomic_charges[f];
            //}
            //else
            opdb<<"  1.00  0.00";
            //opdb<<"  0.00           "<<v_atomic_symbols[f];
            opdb<<std::endl;
        }
        //}
        //}
        if(opdb.is_open())
        {
            opdb.close();
        }
#ifdef _PDB_DEBUG_MESSAGES_
        std::cout<<"CPDB: The PDB was writen!"<<std::endl;
#endif
        return true;
    }
    catch(...)
    {
        return false;
    }
}

uint CPdb::GetAtomIndex(std::string s)
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

inline void CPdb::SetMinMax(glm::dvec3 &v)
{
    min_coord.x = std::min(min_coord.x,v.x);
    min_coord.y = std::min(min_coord.y,v.y);
    min_coord.z = std::min(min_coord.z,v.z);
    max_coord.x = std::max(max_coord.x,v.x);
    max_coord.y = std::max(max_coord.y,v.y);
    max_coord.z = std::max(max_coord.z,v.z);
}

glm::dvec3 CPdb::GetMinCoordinate(void)
{
    return min_coord;
}

glm::dvec3 CPdb::GetMaxCoordinate(void)
{
    return max_coord;
};


// END

