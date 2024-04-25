/*
 * cmodel.cpp
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

#include <cmodel.h>
#include <ctype.h>
#include <algorithm>

//#define _DEBUG_MODEL_DATA_
//#define _DEBUG_LINKED_LIST_
//#define _DEBUG_BONDS_DATA_

CModel::CModel()
{
    is_periodic = false;
    is_linked_list = false;
    is_linked_cell = false;
    is_bond_eval = false;
    is_centered = false;
    //is_unchanged = true;
    IsSaved(true);
    r_cut_radius = 0;
    u_total_cells = 1;
    v_pbc_x = glm::ivec2(0,1);
    v_pbc_y = glm::ivec2(0,1);
    v_pbc_z = glm::ivec2(0,1);
    center_xyz = glm::vec3(0.0f);
    v_table_index.clear();
    v_table_atoms.clear();
}

bool CModel::ReadStructureFile(std::string d, std::string f, int fmt)
{
    file_format = fmt;
    bool res = true;
    IsPeriodic(false);
#ifdef _DEBUG_MODEL_DATA_
    DEBUG_MESSAGE("Read file")
    DEBUG_VALUE("File = ",res)
#endif
    switch(file_format)
    {
    case INPUT_FILE_TYPE_XYZ:
        res = xyz.ReadFile(d,f,v_real_xyz,v_atom_index);
        if(res)
        {
            dmin_xyz = xyz.GetMinCoordinate();
            dmax_xyz = xyz.GetMaxCoordinate();
            if(xyz.is_periodic)
            {
                IsPeriodic(xyz.is_periodic);
                m_cell_box = xyz.cell_xyz;
            }
        }
        break;
    case INPUT_FILE_TYPE_VSP:
        res = vasp.ReadFile(d,f,v_real_xyz,v_atom_index);
        if(res)
        {
            dmin_xyz = vasp.GetMinCoordinate();
            dmax_xyz = vasp.GetMaxCoordinate();
            if(vasp.is_periodic)
            {
                IsPeriodic(vasp.is_periodic);
                m_cell_box = vasp.cell_xyz;
            }
            //std::cout<<vasp.v_table_index.size()<<std::endl;
            //std::cout<<vasp.v_table_atoms.size()<<std::endl;
            v_table_index=vasp.v_table_index;
            v_table_atoms=vasp.v_table_atoms;
        }
        break;
    case INPUT_FILE_TYPE_PDB:
        res = pdb.ReadFile(d,f,v_real_xyz,v_atom_index);
        if(res)
        {
            dmin_xyz = pdb.GetMinCoordinate();
            dmax_xyz = pdb.GetMaxCoordinate();
            if(pdb.is_periodic)
            {
                IsPeriodic(pdb.is_periodic);
                m_cell_box = pdb.cell_xyz;
            }
        }
        break;
    }
    if(res)
    {
        min_xyz=glm::vec3(dmin_xyz);
        max_xyz=glm::vec3(dmax_xyz);
        UpdateAtomPositions();
        SetModel();
        // Working: alpha
        EvalModel();
    }
    return res;
}


void CModel::SetModel(void)
{
    is_linked_cell = false;
    is_linked_list = false;
    IsCentered(false);
    // Center of the geometry
    center_xyz = 0.5f*(max_xyz+min_xyz);
#ifdef _DEBUG_MODEL_DATA_
    DEBUG_VEC3("min_xyz = ",min_xyz)
    DEBUG_VEC3("max_xyz = ",max_xyz)
    DEBUG_VEC3("center_xyz = ",center_xyz)
#endif
    EvalCutRadius();
    //UpdateAtomPositions();
    // To avoid atoms at the boundary in non-periodic systems
    //if(!IsPeriodic()){
      min_xyz -= (1.1*max_radius);
      max_xyz += (1.1*max_radius);
    //}
    // Set a box to enclose the atomic structure
    v_box_size = max_xyz-min_xyz;
    // Define a supercell for non-periodic systems
    if(!IsPeriodic())
    {
        m_cell_box = diagonal3x3(v_box_size);
    }
    EvalAtomicBonds();
    UpdateBondPositions();
    // Do not center before the bonds evaluation
    if(!IsCentered()){
      CenterView();
      IsCentered(true);
    }
#ifdef _DEBUG_MODEL_DATA_
    DEBUG_VEC3("min_xyz2 = ",min_xyz)
    DEBUG_VEC3("max_xyz2 = ",max_xyz)
#endif
}

void CModel::EvalModel(void)
{
    bool is_new = false;
    uint tindex;
    v_table_index.clear();
    v_table_atoms.clear();
    if(v_table_index.size()==0)
    {
        v_table_index.push_back(v_atom_index[0]);
        for (uint i=0; i<v_atom_index.size(); i++)
        {
#ifdef _DEBUG_MODEL_DATA_
            std::cout<<"v_atom_index["<<i<<"] = "<<v_atom_index[i]<<std::endl;
#endif
            is_new = true;
            for (uint j=0; j<v_table_index.size(); j++)
            {
                if(v_atom_index[i]==v_table_index[j])
                {
                    is_new = false;
                }
            }
            if(is_new)
                v_table_index.push_back(v_atom_index[i]);
        }
    }
    v_table_atoms.resize(v_table_index.size());
    for(uint j=0; j<v_table_index.size(); j++)
    {
        tindex = v_table_index[j];
        v_table_atoms[j]=0;
        for( uint i=0; i<v_atom_index.size(); i++)
        {
            if( tindex == v_atom_index[i] )
            {
                v_table_atoms[j]++;
            }
        }
    }
#ifdef _DEBUG_MODEL_DATA_
    for (uint i=0; i<v_table_index.size(); i++)
    {
        std::cout<<"v_table_index = "<<v_table_index[i]<<std::endl;
        std::cout<<"v_table_atoms = "<<v_table_atoms[i]<<std::endl;
    }
#endif
}

bool CModel::SaveStructureFile(std::string d, std::string f, int fmt)
{
    bool res = false;
    switch(file_format)
    {
    case INPUT_FILE_TYPE_XYZ:
        res = xyz.SaveFile(d,f,v_real_xyz,v_atom_index);
        break;
    case INPUT_FILE_TYPE_VSP:
        res = vasp.SaveFile(d,f,v_real_xyz,v_atom_index);
        break;
    case INPUT_FILE_TYPE_PDB:
        res = pdb.SaveFile(d,f,v_real_xyz,v_atom_index);
        break;
    }
    IsSaved(true);
    return res;
}

/* deprecated Sat Feb  1 09:37:14 EST 2020
void CModel::SetNewCoordinates(void)
{
    v_real_xyz.resize(v_float_raw_xyz.size());
    for (int i=0; i<v_real_xyz.size(); i++)
    {
        v_real_xyz[i]=glm::dvec3(v_float_raw_xyz[i]);
    }
}*/

void CModel::SetNewCoordinates(uint u,glm::vec3& v3)
{
    v_float_raw_xyz[u]=v3;
    IsSaved(false);
}

void CModel::SetElement(uint num,uint idx)
{
    v_atom_index[num]=idx;
    IsSaved(false);
    EvalModel();
}

void CModel::DeleteAtomIndex(uint u)
{
#ifdef _DEBUG_MODEL_DATA_
    std::cout<<"delete = "<<u<<std::endl;
#endif
    if(file_format == 1)
        vasp.DeleteAtom(u);
    v_real_xyz.erase(v_real_xyz.begin()+u);
    //v_float_raw_xyz.erase(v_float_raw_xyz.begin()+u);
    v_atom_index.erase(v_atom_index.begin()+u);
    UpdateAtomPositions();
    //
    IsSaved(false);
    EvalModel();
}

uint CModel::GetAtomIndex(std::string s)
{
    return 0;
}

void CModel::AnalyseCoordinates(void)
{
    if(is_periodic)
    {
#ifdef _DEBUG_MODEL_DATA_
        std::cout<<" CModel: is_periodic"<<std::endl;
#endif
    }//else{
    //glm::vec3 min_xyz;
    //glm::vec3 max_xyz;
    //float vacuum = 2.2;
    //center_xyz.resize(total_atoms);
    // set the coordinate reference to (0,0,0)
    //for(uint i=0;i<total_atoms;i++){
    // center_xyz[i]= v_float_raw_xyz[i] - min_xyz;
    //}
    //cell_size=max_xyz + vacuum;
    //}
}

void CModel::SetCellList(void)
{
    v_cell_head.clear();
    v_cell_list.clear();
    v_cell_head.resize(u_cell_number);
    v_cell_list.resize(GetTotalAtoms());
    for(uint _n=0; _n<u_cell_number; _n++)
        v_cell_head[_n] = -1;
    for(uint _n=0; _n<GetTotalAtoms(); _n++)
        v_cell_list[_n] = 0;
}



// Fri Jan 13 16:55:51 MST 2012
// beta version
// find bonds between atoms closer than the sum of their van der Waals radius
void CModel::EvalAtomicBonds(void)
{
    int j;
    int u_icell;
    bool use_pbc;
    float rl, r, r2, rr, ri, rj, rlz, rlxy;
    glm::tvec2<uint> vidx;
    glm::tvec2<float> vang;
    glm::vec3 vi;
    glm::vec3 vj;
    glm::vec3 vij;
    //
    uint k;
    glm::vec3       v_positive_r;
    glm::tvec3<int> v_integer_r;
    glm::tvec3<int> v_neighbor_cell;
    glm::tvec3<int> v_pbc;
    //
    uint max_bonds = 15*GetTotalAtoms();
    v_bond_table.resize(0);
    i_number_of_bonds = 0;
    i_number_of_bonds_pbc = 0;
    // clean
    //v_bond_number.clear();
    m_bond_indices.clear();
    m_bond_xyz.clear();
    m_bond_ang.clear();
    ////////////////////////////////////////
    //v_bond_number.resize(max_bonds);
    v_bond_length.resize(max_bonds);
    m_bond_indices.resize(max_bonds);
    m_bond_xyz.resize(max_bonds);
    m_bond_ang.resize(max_bonds);
    //v_bond_number_pbc.resize(uint(max_bonds/2));
    //m_bond_indices_pbc.resize(uint(max_bonds/2));
    //m_bond_boundary_pbc.resize(uint(max_bonds/2));
    ////////////////////////////////////////
    if(!is_linked_list)
    {
        EvalLinkedList();
    }
    //
    for(int i=0; i<GetTotalAtoms()-1; i++)
    {
        if(strcmp(symbol[v_atom_index[i]].c_str(),"X"))
        {
#ifdef _DEBUG_LINKED_LIST_
            std::cout<<"symbol = "<<symbol[v_atom_index[i]]<<std::endl;
#endif
            vi = v_float_raw_xyz[i];
            ri = atom_rgbs[v_atom_index[i]][3];
            if(is_linked_list && (u_cell_number > 8))
            {
                v_positive_r = vi - min_xyz;
                //v_positive_r = vi + 0.5f*v_box_size;
                v_integer_r  = glm::tvec3<int>(v_positive_r*v_cell_frac);
                for (int _m=0; _m<i_neighbor_cells; _m++)
                {
                    v_neighbor_cell = v_integer_r+neighbor_cells[_m];
                    // Used to apply PCB to cells in each dimension
                    for (uint coord=0; coord<3; coord++)
                    {
                        if(v_neighbor_cell[coord] >= v_cell_side[coord])   // check if  PBC is necessary
                        {
                            v_neighbor_cell[coord] = 0;                      // apply PBC to each  cell
                            //if(is_pbc) use_pbc = false;
                        }
                        else if(v_neighbor_cell[coord] < 0)                // check if  PBC is necessary
                        {
                            v_neighbor_cell[coord] = v_cell_side[coord]-1;   // apply PBC to each cell
                            //if(is_pbc) use_pbc = false;
                        }
                    }
                    u_icell = GetLinkedCell(v_neighbor_cell,v_cell_side);    // head atom index in the the cell
                    if(u_icell>=0 && u_icell < (int)u_cell_number)           // inside of a cells
                    {
                        j = v_cell_head[u_icell];                         // head atom in the actual cell
                    }
                    else                                                // out of the box
                    {
                        j = -1;                                           // outside of a cells
                    }
                    while(1)                                            // over all the particles in the cell
                    {
                        if(j<0) break;                                    // stop searching in the cell
                        //if(_m!=0 || j>i)                                // avoid self-interaction
                        //if((j>i) && (v_ft[j] == v_ft[i]))               // avoid self-interaction and double bond
                        if(j>i)                                           // avoid self-interaction and double bond
                        {
                            r2 = 0;                                         // set distance to cero
                            vj = v_float_raw_xyz[j];
                            rj = atom_rgbs[v_atom_index[j]][3];
                            r = (ri+rj);
                            rr = 1.2*(r*r);
                            use_pbc = false;
                            vij = (vj-vi);
                            /*
                            for(uint coord=0; coord<3; coord++){
                                v_pbc[coord] = 0;
                                if(vij[coord] <= -max_xyz[coord]){
                                    vj += (2.0f*max_xyz[coord]); //2.0*m_bbox[coord];       // PBC
                                    use_pbc = true;
                                    v_pbc[coord] = 1;
                                }else if(vij[coord] > max_xyz[coord]){
                                    vj -= (2.0f*max_xyz[coord]); //2.0*m_bbox[coord];       // PBC
                                    use_pbc = true;
                                    v_pbc[coord] = -1;
                                }
                            }*/
                            //if(use_pbc) vij = (vj-vi);
                            r2 = glm::length2(vij);
                            if( (r2 < rr) && (r2 > 0.1) )                                     // atoms inside de cut radius
                            {
                                vidx = glm::tvec2<uint>(i,j);
                                r = sqrt(r2);
                                vij = (vj-vi);
                                rlz = vij.z;
                                vij.z = 0;
                                rlxy = glm::l2Norm(vij);
                                vang.x = atan2(vij.y,vij.x);       // bond precession
                                vang.y = fabs(atan2(rlxy,rlz));      // bond tilt
                                m_bond_indices[i_number_of_bonds]=vidx;
                                m_bond_xyz[i_number_of_bonds]=(vj+vi)/2.0f;
                                m_bond_ang[i_number_of_bonds]=vang;
                                v_bond_length[i_number_of_bonds]=r;
#ifdef _DEBUG_LINKED_LIST_
                                std::cout<<"r = "<<r<<" - ";
#endif
                                i_number_of_bonds++;
                            }
                        }
                        j = v_cell_list[j]; // cell_list loop!
                    }
#ifdef _DEBUG_LINKED_LIST_
                    std::cout<<std::endl;
#endif
                }
            }
            else
            {
                // the code below can be used for small number of atoms.
                // instead searching with a cell_list
                // the code below can be used for small number of atoms.
                // searching inside the neighbour cells
                //for (int _m=0; _m<27; _m++){
                for(int j=i+1; j<GetTotalAtoms(); j++)
                {
                    if(strcmp(symbol[v_atom_index[j]].c_str(),"X"))
                    {
                        vj = v_float_raw_xyz[j];
                        /*for(uint coord=0; coord<3; coord++){
                          vj += neighbor_cells[_m][coord]*get_uvw_to_xyz(coord); //  2.0*m_bbox[coord];       // PBC
                          v_pbc[coord] = neighbor_cells[_m][coord];
                          }*/
                        rj = atom_rgbs[v_atom_index[j]][3];
                        r = (ri+rj);
                        rr = 1.2*(r*r);
                        use_pbc = false;
                        vij = (vj-vi);
                        r2 = glm::length2(vij);
                        if( (r2 <= rr) &&  (r2 > 0.1) )
                        {
                            vidx = glm::tvec2<uint>(i,j);
                            r = sqrt(r2);
                            rlz = vij.z;
                            vij.z = 0;
                            rlxy = glm::l2Norm(vij);
                            vang.x = atan2(vij.y,vij.x);       // bond precession
                            vang.y = fabs(atan2(rlxy,rlz));      // bond tilt
                            m_bond_indices[i_number_of_bonds]=vidx;
                            m_bond_xyz[i_number_of_bonds]=(vj+vi)/2.0f;
                            m_bond_ang[i_number_of_bonds]=vang;
                            v_bond_length[i_number_of_bonds]=r;
                            //std::cout<<"r = "<<r<<" - ";
                            i_number_of_bonds++;
                        }
                    }
                }

                //}
            }
        }
    }
#ifdef _DEBUG_BONDS_DATA_
    //std::cout<<"m_bond_indices = "<<glm::to_string(m_bond_indices)<<std::endl;
    print_glm_vec("m_bond_indices",m_bond_indices);
#endif
    is_bond_eval = true;
}

void CModel::EvalCutRadius(void)
{
    max_radius = 0;
    r_cut_radius = 0;
    for(uint i=0; i<GetTotalAtoms(); i++)
    {
        //r_cut_radius=std::max(r_cut_radius,atom_rgbs[v_atom_index[i]][3]);
        max_radius=std::max(max_radius,atom_rgbs[v_atom_index[i]][3]);
    }
    r_cut_radius =2.0*max_radius;
    //r_cut_radius *=2.0;
    r_cut_radius_2 = (r_cut_radius*r_cut_radius);
#ifdef _DEBUG_MODEL_DATA_
    std::cout<<" Max Radius = "<<max_radius<<std::endl;
    std::cout<<" Cut Radius = "<<r_cut_radius<<std::endl;
    std::cout<<" Cut Radius2 = "<<r_cut_radius_2<<std::endl;
#endif
}

// Linked and shell cell configuration functions
void CModel::SetCells(void)
{
    float inv_cut_rad = float(1.0/r_cut_radius);
    //v_box_size = max_xyz-min_xyz;
    v_cell_side = glm::tvec3<int>(inv_cut_rad*v_box_size);
    u_cell_number = glm::compMul(v_cell_side);
#ifdef _DEBUG_LINKED_LIST_
    std::cout<<"inv_cut_rad = "<<inv_cut_rad<<std::endl;
    std::cout<<"v_box_size = "<<glm::to_string(v_box_size)<<std::endl;
    std::cout<<"v_cell_side = "<<glm::to_string(v_cell_side)<<std::endl;
    std::cout<<"u_cell_number = "<<u_cell_number<<std::endl;
#endif
    if(u_cell_number > 1)
        SetNeighbourList();
}

void CModel::SetInverseCell(void)
{
    v_cell_frac = glm::vec3(v_cell_side)/v_box_size;
#ifdef _DEBUG_LINKED_LIST_
    std::cout<<"v_cell_frac = "<<glm::to_string(v_cell_frac)<<std::endl;
#endif
}

void CModel::EvalLinkedList(void)
{
    SetCells();
    if(u_cell_number > 1)
    {
        SetInverseCell();
        SetCellList();
        uint _n, u_icell;
        //
        glm::vec3 v_positive_r;
        glm::tvec3<int> v_integer_r;
        //
        for(_n=0; _n< GetTotalAtoms(); _n++)
        {
            v_positive_r=v_float_raw_xyz[_n];
            v_positive_r -= min_xyz;
            v_integer_r = glm::tvec3<int>(v_positive_r*v_cell_frac);
            u_icell = ((v_integer_r.z*v_cell_side.y+v_integer_r.y)*v_cell_side.x)+v_integer_r.x;
#ifdef _DEBUG_LINKED_LIST_
            std::cout<<"v_positive_r1 = "<<glm::to_string(v_positive_r)<<std::endl;
            std::cout<<"v_integer_r = "<<glm::to_string(v_integer_r)<<std::endl;
            std::cout<<" u_icell = "<<u_icell<<" n = "<<_n<<std::endl;
#endif
            v_cell_list[_n] = v_cell_head[u_icell];
            v_cell_head[u_icell] = _n;
        }
#ifdef _DEBUG_LINKED_LIST_
        print_vector("v_cell_list",v_cell_list);
        print_vector("v_cell_head",v_cell_head);
#endif
        is_linked_list = true;
    }
}

void CModel::SetNeighbourList(void)
{
    glm::tvec3<int> v1;
    int xpcb=-2, ypcb=-2, zpcb=-2;
    // set the necesary neighbor cells
    //if(v_cell_side[0]==2) xpcb=-1;
    //if(v_cell_side[1]==2) ypcb=-1;
    //if(v_cell_side[2]==2) zpcb=-1;
    neighbor_cells_xyz.clear();
    neighbor_cells_xyz.resize(0);
    i_neighbor_cells=0;
    for(int k=0; k<27; k++)
    {
        v1.x = neighbor_cells[k][0];
        v1.y = neighbor_cells[k][1];
        v1.z = neighbor_cells[k][2];
        if((v1.x > xpcb) && (v1.y > ypcb) && (v1.z > zpcb))
        {
            neighbor_cells_xyz.push_back(v1);
            i_neighbor_cells++;
        }
    }
    is_linked_cell=true;
    for(int i=0; i<3; i++)
    {
        if(v_cell_side[1]<=2)
        {
            is_linked_cell=false;
        }
    }
}

// Depredated Sun Feb  2 08:25:45 EST 2020
void CModel::CenterCell(void)
{
    if(!IsCentered())
    {
        //xyz.CenterCoordinates();
        glm::dvec3 vcenter = 0.5*(dmax_xyz+dmin_xyz);
        for (std::vector<glm::dvec3>::iterator it = v_real_xyz.begin(); it != v_real_xyz.end(); ++it)
        {
#ifdef _DEBUG_MODEL_DATA_
            std::cout<<"max = "<<glm::to_string(*it)<<std::endl;
#endif
            *it -= vcenter;
        }
        center_xyz = glm::vec3(0,0,0);
        IsCentered(true);
    }
    UpdateAtomPositions();
    IsSaved(false);
}

void CModel::UpdateView(void){
  UpdateAtomPositions();
  UpdateBondPositions();
  if(IsCentered())
    CenterView();
}

void CModel::SetPbcImages(int x, int y, int z){
    v_pbc_x.x = -int(ceil(float(x)/2.0)-1.0);
    v_pbc_x.y =  int(floor(float(x)/2.0)+1);
    v_pbc_y.x = -int(ceil(float(y)/2.0)-1.0);
    v_pbc_y.y =  int(floor(float(y)/2.0)+1);
    v_pbc_z.x = -int(ceil(float(z)/2.0)-1.0);
    v_pbc_z.y =  int(floor(float(z)/2.0)+1);
    u_total_cells = uint(x*y*z);
    //DEBUG_VEC3("v_pbc_x = ",v_pbc_x)
    //DEBUG_VEC3("v_pbc_y = ",v_pbc_y)
    //DEBUG_VEC3("v_pbc_z = ",v_pbc_z)

}

void CModel::UpdateAtomPositions(void){
  //uint total_pcb_atoms = (u_total_cells*v_real_xyz.size());
  //DEBUG_VALUE(" total_pcb_atoms = ", total_pcb_atoms)
  if(v_float_raw_xyz.size()!=v_real_xyz.size())
    v_float_raw_xyz.resize(v_real_xyz.size());
  // Add periodic images
  for (int i=0; i<v_real_xyz.size(); i++)
  {
#ifdef _DEBUG_MODEL_DATA_
    std::cout<<"max = "<<glm::to_string(v_real_xyz[i])<<std::endl;
    std::cout<<"v_atom_index = "<<v_atom_index[i]<<std::endl;
#endif
    v_float_raw_xyz[i]=glm::vec3(v_real_xyz[i]);
  }
}

void CModel::UpdateBondPositions(void){
  if(is_bond_eval){
    if(v_bond_position_xyz.size()!=m_bond_xyz.size())
      v_bond_position_xyz.resize(m_bond_xyz.size());
      for(uint u=0; u<i_number_of_bonds; u++)
      {
        v_bond_position_xyz[u]=m_bond_xyz[u];
        //if(IsCentered()){
          //v_bond_position_xyz[u]-= center_xyz;
        //}
      }
  }
}

// deprecated Sat Feb  1 12:01:07 EST 2020
void CModel::CenterView(void)
{
#ifdef _DEBUG_MODEL_DATA_
    std::cout<<"center_xyz = "<<glm::to_string(center_xyz)<<std::endl;
#endif
    //center_xyz = 0.5f*(max_xyz+min_xyz);
    for (std::vector<glm::vec3>::iterator it = v_float_raw_xyz.begin(); it != v_float_raw_xyz.end(); ++it)
    {
#ifdef _DEBUG_MODEL_DATA_
        std::cout<<"xyz = "<<glm::to_string(*it)<<std::endl;
#endif
        *it -= center_xyz;
    }
    if(is_bond_eval){
      for(uint u=0; u<i_number_of_bonds; u++)
      {
        v_bond_position_xyz[u]-= center_xyz;
      }
    }
}

inline void CModel::SetMinMax(glm::vec3 &v)
{
    min_xyz.x = std::min(min_xyz.x,v.x);
    min_xyz.y = std::min(min_xyz.y,v.y);
    min_xyz.z = std::min(min_xyz.z,v.z);
    max_xyz.x = std::max(max_xyz.x,v.x);
    max_xyz.y = std::max(max_xyz.y,v.y);
    max_xyz.z = std::max(max_xyz.z,v.z);
}


// END

