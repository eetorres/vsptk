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

#ifndef _CXYZ_H_
#define _CXYZ_H_

#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

#include "elements.h"

class CXyz
{

public:

    CXyz();
    ~CXyz() {};
    //
    uint header_lines;
    uint total_atoms;
    //
    glm::dvec3 min_coord;
    glm::dvec3 max_coord;
    glm::dmat3 cell_xyz;
    //
    bool is_periodic;
    bool GetFormat(std::string d,std::string f);
    bool ReadFile(std::string d,std::string f, std::vector<glm::dvec3>& xyz, std::vector<uint>& idx);
    bool SaveFile(std::string d,std::string f, std::vector<glm::dvec3>& xyz, std::vector<uint>& idx);
    uint GetAtomIndex(std::string s);
    //

//protected:
    //
    uint GetTotalAtoms(void);
    inline void SetMinMax(glm::dvec3 &v);
    glm::dvec3 GetMinCoordinate(void);
    glm::dvec3 GetMaxCoordinate(void);
};

#endif
