/*
 * mouseh
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

#ifndef _CMODEL_H_
#define _CMODEL_H_

#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <vector>
#include <algorithm>    // std::max
#include <glm/glm.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include <glm/gtx/norm.hpp>
#include "elements.h"
#include "ctimer.h"

#include <cxyz.h>
#include <cpdb.h>
#include <cvasp.h>

// Statements like:
// #pragma message(Reminder "Fix this problem!")
// Which will cause messages like:
// C:\Source\Project\main.cpp(47): Reminder: Fix this problem!
// to show up during compiles. Note that you can NOT use the
// words "error" or "warning" in your reminders, since it will
// make the IDE think it should abort execution. You can double
// click on these messages and jump to the line in question.
//#define Stringize( L )     #L
//#define MakeString( M, L ) M(L)
//#define $Line MakeString( Stringize, __LINE__ )
//#define Reminder __FILE__ "(" $Line ") : Reminder: "

// Neighbouring cells use GLM
const glm::tvec3<int> neighbor_cells[27] =
{
    glm::vec3(  0, 0, 0), // 0
    glm::vec3(  1, 0, 0), // 1
    glm::vec3(  1, 1, 0), // 2
    glm::vec3(  0, 1, 0), // 3
    glm::vec3( -1, 1, 0), // 4
    glm::vec3(  0, 0, 1), // 5
    glm::vec3(  1, 0, 1), // 6
    glm::vec3( -1, 0, 1), // 7
    glm::vec3(  1, 1, 1), // 8
    glm::vec3(  0, 1, 1), // 9
    glm::vec3( -1, 1, 1), // 10
    glm::vec3( -1,-1, 1), // 11
    glm::vec3(  0,-1, 1), // 12
    glm::vec3(  1,-1, 1), // 13
    glm::vec3(  0, 0,-1), // 14
    glm::vec3(  1, 0,-1), // 15
    glm::vec3( -1, 0,-1), // 16
    glm::vec3(  1, 1,-1), // 17
    glm::vec3(  0, 1,-1), // 18
    glm::vec3( -1, 1,-1), // 19
    glm::vec3( -1,-1,-1), // 20
    glm::vec3(  0,-1,-1), // 21
    glm::vec3(  1,-1,-1), // 22
    glm::vec3( -1, 0, 0), // 23
    glm::vec3( -1,-1, 0), // 24
    glm::vec3(  0,-1, 0), // 25
    glm::vec3(  1,-1, 0)  // 26
};

const uint INPUT_FILE_TYPE_UNKNOWN = 100;
const uint INPUT_FILE_TYPE_XYZ = 0;
const uint INPUT_FILE_TYPE_VSP = 1;
const uint INPUT_FILE_TYPE_PDB = 2;
const uint INPUT_FILE_TYPE_GAU = 3;
const uint INPUT_FILE_TYPE_DLP = 4;
const uint INPUT_FILE_TYPE_ZMT = 5;

class CModel
{

public:

    CModel();
    ~CModel() {};
    // Get
    uint GetAtomIndex(std::string s);
    uint GetTotalAtoms(void)
    {
        return v_real_xyz.size();
    };

// protected:
    //
    std::vector<uint> v_table_index;
    std::vector<uint> v_table_atoms;
    uint GetAtomTotal(uint u)
    {
        return v_table_atoms[u];
    };
    uint GetAtomIndex(uint u)
    {
        return v_table_index[u];
    };


// private:
    //
    CXyz  xyz;
    CVasp vasp;
    CPdb  pdb;
    // Variables
    int  file_format;
    int  i_neighbor_cells;
    //
    uint u_cell_number;
    uint u_total_cells;
    uint i_number_of_bonds;
    uint i_number_of_bonds_pbc;
    //
    bool is_linked_cell;
    bool is_bond_eval;
    bool is_updated;
    bool is_periodic;
    bool is_linked_list;
    bool is_centered;
    // bond variables
    float max_radius;
    float r_cut_radius;
    float r_cut_radius_2;
    //float f_atom_bond_delta;
    //
    glm::vec3 v_cell_frac;
    glm::vec3 v_box_size;
    glm::vec3 v_box_middle;
    //
    glm::vec3 min_xyz;
    glm::vec3 max_xyz;
    glm::vec3 center_xyz;
    glm::vec3 pcb_images;
    //
    glm::dvec3 dmin_xyz;
    glm::dvec3 dmax_xyz;
    //
    glm::dmat3 m_cell_box;
    glm::ivec2 v_pbc_x;
    glm::ivec2 v_pbc_y;
    glm::ivec2 v_pbc_z;
    glm::ivec3 v_cell_side;
    //glm::tvec3<int> v_cell_side;
    // Temporal variables
    std::vector<int> v_cell_head;
    std::vector<int> v_cell_list;
    std::vector<glm::tvec3<int>> neighbor_cells_xyz;
    //
//protected:
    //
    std::vector<glm::vec3> v_float_raw_xyz;
    std::vector<glm::vec3> v_bond_position_xyz;
    //
    std::vector<glm::dvec3> v_real_xyz;
    std::vector<uint> v_atom_index;
    //
    std::vector<uint> v_bond_table;
    std::vector<uint> v_bond_number;
    std::vector<uint> v_bond_number_pbc;
    //
    std::vector<float> v_bond_length;
    //
    std::vector<glm::vec3> m_bond_xyz;
    std::vector<glm::vec2> m_bond_ang;
    //
    std::vector<glm::tvec2<uint>> m_bond_indices;
    std::vector<glm::tvec2<uint>> m_bond_boundary_pbc;
    std::vector<glm::tvec2<uint>> m_bond_indices_pbc;

    // Timer class for optimization
    CTimer timer;
    //

    //
// public: methods
    bool ReadStructureFile(std::string d, std::string f, int fmt);
    //bool ReadFileXyz(std::string,std::string);
    //bool ReadFileVasp(std::string,std::string);
    //
    // Set
    //void SetNewCoordinates(void);
    void SetModel(void);
    void SetCells(void);
    void SetInverseCell(void);
    void SetNeighbourList(void);
    void SetCellList(void);
    inline void SetMinMax(glm::vec3 &v);
    void SetElement(uint num,uint idx);
    void SetNewCoordinates(uint u,glm::vec3& v3);
    void SetPbcImages(int x, int y, int z);
    //
    bool SaveStructureFile(std::string,std::string, int fmt);
    //
    void CenterView(void);
    void CenterCell(void);
    void AnalyseCoordinates(void);
    //void EvalCenterOfMass(void);
    //void ShowBonds(bool b);
    //
    // Eval
    void EvalCutRadius(void);
    void EvalModel(void);
    void EvalAtomicBonds(void);
    void EvalLinkedList(void);
    void UpdateView(void);
    void UpdateAtomPositions(void);
    void UpdateBondPositions(void);
    //
    void DeleteAtomIndex(uint u);
    //
    //
    void IsPeriodic(bool b)
    {
        is_periodic = b;
    };
    bool IsPeriodic(void)
    {
        return is_periodic;
    };
    void IsSaved(bool b)
    {
        is_updated = b;
    };
    bool IsSaved(void)
    {
        return is_updated;
    };
    void IsCentered(bool b)
    {
        is_centered = b;
    };
    bool IsCentered(void)
    {
        return is_centered;
    };
    //
    inline int GetLinkedCell(glm::tvec3<int>& p,glm::tvec3<int>& s)
    {
        return int(((p.z*s.y+p.y)*s.x)+p.x);
    };
    // Auxilar functions
    void print_vector(std::string s,std::vector<int>& v)
    {
        std::cout<<s<<" ="<<std::endl;
        for (uint i = 0; i < v.size(); ++i)
        {
            std::cout << v[i] <<" ";
        }
        std::cout<<std::endl;
    };
    void print_glm_vec(std::string s,std::vector<glm::tvec2<uint>>& v)
    {
        std::cout<<s<<" ="<<std::endl;
        for (uint i = 0; i < v.size(); ++i)
        {
            std::cout<<glm::to_string(v[i])<<std::endl;
        }
        std::cout<<std::endl;
    };
    void print_glm_vec3(std::string s,std::vector<glm::vec3>& v)
    {
        std::cout<<s<<" ="<<std::endl;
        for (uint i = 0; i < v.size(); ++i)
        {
            std::cout<<glm::to_string(v[i])<<std::endl;
        }
        std::cout<<std::endl;
    };
};

#endif
