/*
 * Atoms.h
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

#ifndef _ATOM_SYMBOL_H_
#define _ATOM_SYMBOL_H_

#include<string>
#include<glm/glm.hpp>

#if defined(_WIN32)
typedef unsigned int uint;
#endif

const unsigned int periodic_table_atoms=110;
#define PTABLE_ATOMS 110

// Simple debug messages
#define DEBUG_MESSAGE(X) std::cout<<X<<std::endl;
#define DEBUG_VALUE(X,N) std::cout<<X<<N<<std::endl;
#define DEBUG_VEC3(X,V) std::cout<<X<<glm::to_string(V)<<std::endl;

// Periodic table
const std::string symbol[periodic_table_atoms] =
{
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "X"
}; // last symbol is for a dummy atom

// Sorted periodic table
const std::string symbol_sorted[periodic_table_atoms] =
{
    "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Cl", "Ar", //   1 -  10
    "Ca", "Sc", "Ti", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", //  11 -  20
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Zr", "Nb", //  21 -  30
    "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", //  31 -  40
    "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", //  41 -  50
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", //  51 -  60
    "Ta", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", //  61 -  70
    "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "Np", "Pu", //  71 -  80
    "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", //  81 -  90
    "Db", "Sg", "Bh", "Hs", "Mt",                               //  91 -  95
    "H",  "B",  "C",  "N",  "O",  "F",  "P",  "S",  "K",  "V",  //  96 - 105
    "Y",  "I",  "W",  "U",  "X"                                 // 106 - 110
}; // last symbol is for a dummy atom

// Atomic symbol characters
const size_t symbol_t[periodic_table_atoms] =
{
    1, 2, 2, 2, 1, 1, 1, 1, 1, 2,
    2, 2, 2, 2, 1, 1, 2, 2, 1, 2,
    2, 2, 1, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 1, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 1, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 1, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 1
}; // last symbol is for a dummy atom

// Atom index
enum atomic_number
{
    H=0,He, Li, Be, B,  C,  N,  O,  F,  Ne,
    Na, Mg, Al, Si, P,  S,  Cl, Ar, K,  Ca,
    Sc, Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn,
    Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y,  Zr,
    Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn,
    Sb, Te, I,  Xe, Cs, Ba, La, Ce, Pr, Nd,
    Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb,
    Lu, Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg,
    Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th,
    Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm,
    Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, X
};

// Atom index sorted
const uint atom_index[periodic_table_atoms] =
{
    He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, //   1 -  10
    Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, //  11 -  20
    Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, //  21 -  30
    Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, //  31 -  40
    Te, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, //  41 -  50
    Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, //  51 -  60
    Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, //  61 -  70
    Po, At, Rn, Fr, Ra, Ac, Th, Pa, Np, Pu, //  71 -  80
    Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, //  81 -  90
    Db, Sg, Bh, Hs, Mt,                     //  91 -  95
    H,  B,  C,  N,  O,  F,  P,  S,  K,  V,  //  96 - 105
    Y,  I,  W,  U,  X                       // 106 - 110
}; // last symbol is for a dummy atom

// Atomic colors
// special colors
//copper    = {184, 115, 51 }/184;
const float copper[3]    = { 1.00,   0.63,   0.28 };
//gold      = {255, 215, 0  }/ 255;
const float gold[3]      = { 1.00,   0.84,   0.00 };
//steelblue = {176, 196, 222}/ 255;
const float steelblue[3] = { 0.69,   0.77,   0.87 };
//silver    = {164, 170, 173}/ 255;
const float silver[3]    = { 0.64,   0.67,   0.68 };

// Atomic colors
// radius[Angstrom]  colour rgb
// radius, red, green, blue
const glm::vec4 atom_rgbs[110] =
{
    glm::vec4( 0.80,   0.80,   0.80,   0.46 ), // H:
    glm::vec4( 0.64,   0.67,   0.68,   1.60 ), // He:
    glm::vec4( 0.70,   0.70,   0.70,   0.68 ), // Li:
    glm::vec4( 0.64,   0.67,   0.68,   0.35 ), // Be:
    glm::vec4( 0.90,   0.40,   0.00,   0.82 ), // B:
    glm::vec4( 0.35,   0.35,   0.35,   0.85 ), // C:
    glm::vec4( 0.20,   0.20,   0.80,   0.75 ), // N:
    glm::vec4( 0.80,   0.20,   0.20,   0.74 ), // O:
    glm::vec4( 0.70,   0.85,   0.45,   0.72 ), // F:
    glm::vec4( 0.64,   0.67,   0.68,   1.12 ), // Ne:
    glm::vec4( 0.60,   0.60,   0.60,   1.54 ), // Na:
    glm::vec4( 0.60,   0.60,   0.70,   1.37 ), // Mg:
    glm::vec4( 0.64,   0.67,   0.68,   1.33 ), // Al:
    glm::vec4( 0.69,   0.77,   0.87,   1.17 ), // Si:
    glm::vec4( 0.10,   0.70,   0.30,   1.06 ), // P:
    glm::vec4( 0.95,   0.90,   0.20,   1.02 ), // S:
    glm::vec4( 0.15,   0.50,   0.10,   0.99 ), // Cl:
    glm::vec4( 0.64,   0.67,   0.68,   1.54 ), // Ar:
    glm::vec4( 0.50,   0.50,   0.50,   1.33 ), // K:
    glm::vec4( 0.80,   0.80,   0.70,   1.74 ), // Ca:
    glm::vec4( 0.64,   0.67,   0.68,   1.91 ), // Sc:
    glm::vec4( 0.64,   0.67,   0.68,   1.63 ), // Ti:
    glm::vec4( 0.64,   0.67,   0.68,   1.61 ), // V:
    glm::vec4( 0.00,   0.80,   0.00,   1.55 ), // Cr:
    glm::vec4( 0.64,   0.67,   0.68,   1.42 ), // Mn:
    glm::vec4( 0.52,   0.58,   0.65,   1.37 ), // Fe:
    glm::vec4( 0.64,   0.67,   0.68,   1.55 ), // Co:
    glm::vec4( 0.26,   0.27,   0.27,   1.55 ), // Ni:
    glm::vec4( 0.95,   0.79,   0.02,   1.17 ), // Cu:
    glm::vec4( 0.64,   0.67,   0.68,   1.47 ), // Zn:
    glm::vec4( 0.90,   0.00,   1.00,   1.52 ), // Ga:
    glm::vec4( 0.64,   0.67,   0.68,   1.53 ), // Ge:
    glm::vec4( 1.00,   1.00,   0.30,   1.55 ), // As:
    glm::vec4( 0.64,   0.67,   0.68,   1.46 ), // Se:
    glm::vec4( 0.50,   0.08,   0.12,   1.44 ), // Br:
    glm::vec4( 0.64,   0.67,   0.68,   1.60 ), // Kr:
    glm::vec4( 0.64,   0.67,   0.68,   2.78 ), // Rb:
    glm::vec4( 0.64,   0.67,   0.68,   2.45 ), // Sr:
    glm::vec4( 0.64,   0.67,   0.68,   2.08 ), // Y:
    glm::vec4( 0.64,   0.67,   0.68,   1.89 ), // Zr:
    glm::vec4( 0.64,   0.67,   0.68,   1.73 ), // Nb:
    glm::vec4( 0.64,   0.67,   0.68,   1.56 ), // Mo:
    glm::vec4( 0.64,   0.67,   0.68,   1.65 ), // Tc:
    glm::vec4( 0.64,   0.67,   0.68,   1.63 ), // Ru:
    glm::vec4( 0.64,   0.67,   0.68,   1.65 ), // Rh:
    glm::vec4( 0.64,   0.67,   0.68,   1.68 ), // Pd:
    glm::vec4( 0.64,   0.67,   0.68,   1.75 ), // Ag:
    glm::vec4( 0.64,   0.67,   0.68,   1.79 ), // Cd:
    glm::vec4( 0.64,   0.67,   0.68,   1.93 ), // In:
    glm::vec4( 0.64,   0.67,   0.68,   1.71 ), // Sn:
    glm::vec4( 0.64,   0.67,   0.68,   1.75 ), // Sb:
    glm::vec4( 0.64,   0.67,   0.68,   1.73 ), // Te:
    glm::vec4( 0.50,   0.10,   0.50,   1.63 ), // I:
    glm::vec4( 0.64,   0.67,   0.68,   1.70 ), // Xe:
    glm::vec4( 0.64,   0.67,   0.68,   2.95 ), // Cs:
    glm::vec4( 0.64,   0.67,   0.68,   2.47 ), // Ba:
    glm::vec4( 0.64,   0.67,   0.68,   2.17 ), // La:
    glm::vec4( 0.80,   0.80,   0.00,   1.82 ), // Ce:
    glm::vec4( 0.64,   0.67,   0.68,   2.12 ), // Pr:
    glm::vec4( 0.64,   0.67,   0.68,   2.11 ), // Nd:
    glm::vec4( 0.64,   0.67,   0.68,   2.11 ), // Pm:
    glm::vec4( 0.64,   0.67,   0.68,   2.10 ), // Sm:
    glm::vec4( 0.64,   0.67,   0.68,   2.29 ), // Eu:
    glm::vec4( 1.00,   0.84,   0.00,   2.09 ), // Gd:
    glm::vec4( 0.64,   0.67,   0.68,   2.06 ), // Tb:
    glm::vec4( 0.64,   0.67,   0.68,   2.05 ), // Dy:
    glm::vec4( 0.64,   0.67,   0.68,   2.04 ), // Ho:
    glm::vec4( 0.64,   0.67,   0.68,   2.03 ), // Er:
    glm::vec4( 0.64,   0.67,   0.68,   2.03 ), // Tm:
    glm::vec4( 0.64,   0.67,   0.68,   2.24 ), // Yb:
    glm::vec4( 0.64,   0.67,   0.68,   2.02 ), // Lu:
    glm::vec4( 0.64,   0.67,   0.68,   1.86 ), // Hf:
    glm::vec4( 0.64,   0.67,   0.68,   1.73 ), // Ta:
    glm::vec4( 0.64,   0.67,   0.68,   1.58 ), // W:
    glm::vec4( 0.64,   0.67,   0.68,   1.67 ), // Re:
    glm::vec4( 0.64,   0.67,   0.68,   1.64 ), // Os:
    glm::vec4( 0.64,   0.67,   0.68,   1.66 ), // Ir:
    glm::vec4( 0.64,   0.67,   0.68,   1.69 ), // Pt:
    glm::vec4( 0.90,   0.80,   0.00,   1.74 ), // Au:
    glm::vec4( 0.64,   0.67,   0.68,   1.80 ), // Hg:
    glm::vec4( 0.64,   0.67,   0.68,   2.00 ), // Tl:
    glm::vec4( 0.64,   0.67,   0.68,   2.05 ), // Pb:
    glm::vec4( 0.64,   0.67,   0.68,   1.85 ), // Bi:
    glm::vec4( 0.64,   0.67,   0.68,   1.48 ), // Po:
    glm::vec4( 0.80,   0.20,   0.20,   5.43 ), // At:
    glm::vec4( 0.64,   0.67,   0.68,   1.90 ), // Rn:
    glm::vec4( 0.64,   0.67,   0.68,   3.15 ), // Fr:
    glm::vec4( 0.64,   0.67,   0.68,   2.52 ), // Ra:
    glm::vec4( 0.64,   0.67,   0.68,   2.18 ), // Ac:
    glm::vec4( 0.64,   0.67,   0.68,   2.10 ), // Th:
    glm::vec4( 0.64,   0.67,   0.68,   1.91 ), // Pa:
    glm::vec4( 0.64,   0.67,   0.68,   1.69 ), // U:
    glm::vec4( 0.64,   0.67,   0.68,   0.95 ), // Np:
    glm::vec4( 0.64,   0.67,   0.68,   1.81 ), // Pu:
    glm::vec4( 0.64,   0.67,   0.68,   1.61 ), // Am:
    glm::vec4( 0.64,   0.67,   0.68,   0.91 ), // Cm:
    glm::vec4( 0.64,   0.67,   0.68,   0.90 ), // Bk:
    glm::vec4( 0.10,   0.70,   0.30,   1.85 ), // Cf:
    glm::vec4( 0.10,   0.30,   0.70,   2.22 ), // Es:
    glm::vec4( 0.64,   0.67,   0.68,   0.87 ), // Fm:
    glm::vec4( 0.64,   0.67,   0.68,   0.86 ), // Md:
    glm::vec4( 0.64,   0.67,   0.68,   0.85 ), // No:
    glm::vec4( 0.64,   0.67,   0.68,   0.84 ), // Lr:
    glm::vec4( 0.64,   0.67,   0.68,   0.84 ), // Rf:
    glm::vec4( 0.90,   0.80,   0.00,   1.85 ), // Db:
    glm::vec4( 0.64,   0.67,   0.68,   0.84 ), // Sg:
    glm::vec4( 0.64,   0.67,   0.68,   0.84 ), // Bh:
    glm::vec4( 0.64,   0.67,   0.68,   0.84 ), // Hs:
    glm::vec4( 0.64,   0.67,   0.68,   0.84 ), // Mt:
    glm::vec4( 1.00,   0.20,   0.10,   0.40 )  // X:
};

// Atom name
const std::string name[110] =
{
    "Hydrogen",
    "Helium",
    "Lithium",
    "Beryllium",
    "Boron",
    "Carbon",
    "Nitrogen",
    "Oxygen",
    "Fluorine",
    "Neon",
    "Sodium",
    "Magnesium",
    "Aluminium",
    "Silicon",
    "Phosphorus",
    "Sulfur",
    "Chlorine",
    "Argon",
    "Potassium",
    "Calcium",
    "Scandium",
    "Titanium",
    "Vanadium",
    "Chromium",
    "Manganese",
    "Iron",
    "Cobalt",
    "Nickel",
    "Copper",
    "Zinc",
    "Gallium",
    "Germanium",
    "Arsenic",
    "Selenium",
    "Bromine",
    "Krypton",
    "Rubidium",
    "Strontium",
    "Yttrium",
    "Zirconium",
    "Niobium",
    "Molybdenum",
    "Technetium",
    "Ruthenium",
    "Rhodium",
    "Palladium",
    "Silver",
    "Cadmium",
    "Indium",
    "Tin",
    "Antimony",
    "Tellurium",
    "Iodine",
    "Xenon",
    "Caesium",
    "Barium",
    "Lanthanum",
    "Cerium",
    "Praseodymium",
    "Neodymium",
    "Promethium",
    "Samarium",
    "Europium",
    "Gadolinium",
    "Terbium",
    "Dysprosium",
    "Holmium",
    "Erbium",
    "Thulium",
    "Ytterbium",
    "Lutetium",
    "Hafnium",
    "Tantalum",
    "Tungsten",
    "Rhenium",
    "Osmium",
    "Iridium",
    "Platinum",
    "Gold",
    "Mercury",
    "Thallium",
    "Lead",
    "Bismuth",
    "Polonium",
    "Astatine",
    "Radon",
    "Francium",
    "Radium",
    "Actinium",
    "Thorium",
    "Protactinium",
    "Uranium",
    "Neptunium",
    "Plutonium",
    "Americium",
    "Curium",
    "Berkelium",
    "Californium",
    "Einsteinium",
    "Fermium",
    "Mendelevium",
    "Nobelium",
    "Lawrencium",
    "Rutherfordium",
    "Dubnium",
    "Seaborgium",
    "Bohrium",
    "Hassium",
    "Meitnerium",
    "dummy",
};

const glm::vec2 ptable_position[] =
{
    glm::vec2(  7,  30), // H:
    glm::vec2(517,  30), // He:
    glm::vec2(  7,  50), // Li:
    glm::vec2( 37,  50), // Be:
    glm::vec2(367,  50), // B:
    glm::vec2(397,  50), // C:
    glm::vec2(427,  50), // N:
    glm::vec2(457,  50), // O:
    glm::vec2(487,  50), // F:
    glm::vec2(517,  50), // Ne:
    glm::vec2(  7,  70), // Na:
    glm::vec2( 37,  70), // Mg:
    glm::vec2(367,  70), // Al:
    glm::vec2(397,  70), // Si:
    glm::vec2(427,  70), // P:
    glm::vec2(457,  70), // S:
    glm::vec2(487,  70), // Cl:
    glm::vec2(517,  70), // Ar:
    glm::vec2(  7,  90), // K:
    glm::vec2( 37,  90), // Ca:
    glm::vec2( 67,  90), // Sc:
    glm::vec2( 97,  90), // Ti:
    glm::vec2(127,  90), // V:
    glm::vec2(157,  90), // Cr:
    glm::vec2(187,  90), // Mn:
    glm::vec2(217,  90), // Fe:
    glm::vec2(247,  90), // Co:
    glm::vec2(277,  90), // Ni:
    glm::vec2(307,  90), // Cu:
    glm::vec2(337,  90), // Zn:
    glm::vec2(367,  90), // Ga:
    glm::vec2(397,  90), // Ge:
    glm::vec2(427,  90), // As:
    glm::vec2(457,  90), // Se:
    glm::vec2(487,  90), // Br:
    glm::vec2(517,  90), // Kr:
    glm::vec2(  7, 110), // Rb:
    glm::vec2( 37, 110), // Sr:
    glm::vec2( 67, 110), // Y:
    glm::vec2( 97, 110), // Zr:
    glm::vec2(127, 110), // Nb:
    glm::vec2(157, 110), // Mo:
    glm::vec2(187, 110), // Tc:
    glm::vec2(217, 110), // Ru:
    glm::vec2(247, 110), // Rh:
    glm::vec2(277, 110), // Pd:
    glm::vec2(307, 110), // Ag:
    glm::vec2(337, 110), // Cd:
    glm::vec2(367, 110), // In:
    glm::vec2(397, 110), // Sn:
    glm::vec2(427, 110), // Sb:
    glm::vec2(457, 110), // Te:
    glm::vec2(487, 110), // I:
    glm::vec2(517, 110), // Xe:
    glm::vec2(  7, 130), // Cs:
    glm::vec2( 37, 130), // Ba:
    glm::vec2( 67, 130), // La:
    glm::vec2(127, 190), // Ce:
    glm::vec2(157, 190), // Pr:
    glm::vec2(187, 190), // Nd:
    glm::vec2(217, 190), // Pm:
    glm::vec2(247, 190), // Sm:
    glm::vec2(277, 190), // Eu:
    glm::vec2(307, 190), // Gd:
    glm::vec2(337, 190), // Tb:
    glm::vec2(367, 190), // Dy:
    glm::vec2(397, 190), // Ho:
    glm::vec2(427, 190), // Er:
    glm::vec2(457, 190), // Tm:
    glm::vec2(487, 190), // Yb:
    glm::vec2(517, 190), // Lu:
    glm::vec2( 97, 130), // Hf:
    glm::vec2(127, 130), // Ta:
    glm::vec2(157, 130), // W:
    glm::vec2(187, 130), // Re:
    glm::vec2(217, 130), // Os:
    glm::vec2(247, 130), // Ir:
    glm::vec2(277, 130), // Pt:
    glm::vec2(307, 130), // Au:
    glm::vec2(337, 130), // Hg:
    glm::vec2(367, 130), // Tl:
    glm::vec2(397, 130), // Pb:
    glm::vec2(427, 130), // Bi:
    glm::vec2(457, 130), // Po:
    glm::vec2(487, 130), // At:
    glm::vec2(517, 130), // Rn:
    glm::vec2(  7, 150), // Fr:
    glm::vec2( 37, 150), // Ra:
    glm::vec2( 67, 150), // Ac:
    glm::vec2(127, 210), // Th:
    glm::vec2(157, 210), // Pa:
    glm::vec2(187, 210), // U:
    glm::vec2(217, 210), // Np:
    glm::vec2(247, 210), // Pu:
    glm::vec2(277, 210), // Am:
    glm::vec2(307, 210), // Cm:
    glm::vec2(337, 210), // Bk:
    glm::vec2(367, 210), // Cf:
    glm::vec2(397, 210), // Es:
    glm::vec2(427, 210), // Fm:
    glm::vec2(457, 210), // Md:
    glm::vec2(487, 210), // No:
    glm::vec2(517, 210), // Lr:
    glm::vec2( 97, 150), // Rf:
    glm::vec2(127, 150), // Db:
    glm::vec2(157, 150), // Sg:
    glm::vec2(187, 150), // Bh:
    glm::vec2(217, 150), // Hs:
    glm::vec2(247, 150), // Mt:
    glm::vec2( 7,  210), // X:
    glm::vec2(307, 150),
    glm::vec2(337, 150),
    glm::vec2(367, 150),
    glm::vec2(397, 150),
    glm::vec2(427, 150),
    glm::vec2(457, 150),
    glm::vec2(487, 150),
    glm::vec2(517, 150),
    glm::vec2(  7, 150)
};

/*
const glm::vec2 ptable_position[] = {
glm::vec2(13, 5 ), // H
glm::vec2(64, 5 ), // He
glm::vec2(13, 7 ), // Li
glm::vec2(16, 7 ), // Be
glm::vec2(49, 7 ), // B
glm::vec2(52, 7 ), // C
glm::vec2(55, 7 ), // N
glm::vec2(58, 7 ), // O
glm::vec2(61, 7 ), // F
glm::vec2(64, 7 ), // Ne
glm::vec2(13, 9 ), // Na
glm::vec2(16, 9 ), // Mg
glm::vec2(49, 9 ), // Al
glm::vec2(52, 9 ), // Si
glm::vec2(55, 9 ), // P
glm::vec2(58, 9 ), // S
glm::vec2(61, 9 ), // Cl
glm::vec2(64, 9 ), // Ar
glm::vec2(13, 11), // K
glm::vec2(16, 11), // Ca
glm::vec2(19, 11), // Sc
glm::vec2(22, 11), // Ti
glm::vec2(25, 11), // V
glm::vec2(28, 11), // Cr
glm::vec2(31, 11), // Mn
glm::vec2(34, 11), // Fe
glm::vec2(37, 11), // Co
glm::vec2(40, 11), // Ni
glm::vec2(43, 11), // Cu
glm::vec2(46, 11), // Zn
glm::vec2(49, 11), // Ga
glm::vec2(52, 11), // Ge
glm::vec2(55, 11), // As
glm::vec2(58, 11), // Se
glm::vec2(61, 11), // Br
glm::vec2(64, 11), // Kr
glm::vec2(13, 13), // Rb
glm::vec2(16, 13), // Sr
glm::vec2(19, 13), // Y
glm::vec2(22, 13), // Zr
glm::vec2(25, 13), // Nb
glm::vec2(28, 13), // Mo
glm::vec2(31, 13), // Tc
glm::vec2(34, 13), // Ru
glm::vec2(37, 13), // Rh
glm::vec2(40, 13), // Pd
glm::vec2(43, 13), // Ag
glm::vec2(46, 13), // Cd
glm::vec2(49, 13), // In
glm::vec2(52, 13), // Sn
glm::vec2(55, 13), // Sb
glm::vec2(58, 13), // Te
glm::vec2(61, 13), // I
glm::vec2(64, 13), // Xe
glm::vec2(13, 15), // Cs
glm::vec2(16, 15), // Ba
glm::vec2(19, 15), // La
glm::vec2(25, 21), // Ce
glm::vec2(28, 21), // Pr
glm::vec2(31, 21), // Nd
glm::vec2(34, 21), // Pm
glm::vec2(37, 21), // Sm
glm::vec2(40, 21), // Eu
glm::vec2(43, 21), // Gd
glm::vec2(46, 21), // Tb
glm::vec2(49, 21), // Dy
glm::vec2(52, 21), // Ho
glm::vec2(55, 21), // Er
glm::vec2(58, 21), // Tm
glm::vec2(61, 21), // Yb
glm::vec2(64, 21), // Lu
glm::vec2(22, 15), // Hf
glm::vec2(25, 15), // Ta
glm::vec2(28, 15), // W
glm::vec2(31, 15), // Re
glm::vec2(34, 15), // Os
glm::vec2(37, 15), // Ir
glm::vec2(40, 15), // Pt
glm::vec2(43, 15), // Au
glm::vec2(46, 15), // Hg
glm::vec2(49, 15), // Tl
glm::vec2(52, 15), // Pb
glm::vec2(55, 15), // Bi
glm::vec2(58, 15), // Po
glm::vec2(61, 15), // At
glm::vec2(64, 15), // Rn
glm::vec2(13, 17), // Fr
glm::vec2(16, 17), // Ra
glm::vec2(19, 17), // Ac
glm::vec2(25, 23), // Th
glm::vec2(28, 23), // Pa
glm::vec2(31, 23), // U
glm::vec2(34, 23), // Np
glm::vec2(37, 23), // Pu
glm::vec2(40, 23), // Am
glm::vec2(43, 23), // Cm
glm::vec2(46, 23), // Bk
glm::vec2(49, 23), // Cf
glm::vec2(52, 23), // Es
glm::vec2(55, 23), // Fm
glm::vec2(58, 23), // Md
glm::vec2(61, 23), // No
glm::vec2(64, 23), // Lr
glm::vec2(22, 17), // Rf
glm::vec2(25, 17), // Db
glm::vec2(28, 17), // Sg
glm::vec2(31, 17), // Bh
glm::vec2(34, 17), // Hs
glm::vec2(37, 17), // Mt
glm::vec2(40, 17), // Ds
glm::vec2(43, 17), // Rg
glm::vec2(46, 17), // Cn
glm::vec2(49, 17), // ut
glm::vec2(52, 17), // Fl
glm::vec2(55, 17), // up
glm::vec2(58, 17), // Lv
glm::vec2(61, 17), // us
glm::vec2(64, 17)  // uo
};*/

#endif

