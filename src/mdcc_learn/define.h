#ifndef __DEFINE_H__
#define __DEFINE_H__

#include <string>
#include <ctime>
#include <map>
//using namespace std;


//About this program.
#define EXE "mdcc_learn"
#define ABOUT_ME "mdcc_learn ver.0.13.b 02-Feb-2016"
#define PI 3.14159265


typedef unsigned short ifp_u;
enum {
  M_TEST=0,
  M_EM,
  M_VBGMM,
  M_VBFULL,
  M_CAL_GM_DIST,
  M_CAL_GE_DIST,
  M_DUMMY
};

enum {
  DATA_TABLE=0,
  DATA_KKTRJ,
  DATA_DUMMY
};

const std::string sybNew[] = {
  "H",   "CAR",   "C2",   "C3",   "C1",
  "CCA",   "O2",   "O3",   "O2",   "O3",
  "O3",   "N1",   "N2",   "N3",   "NPL",
  "N3+",   "NAM",   "NAR",   "F",   "CL",
  "BR",   "I",   "AL",   "A",   "CA",
  "DU",   "DU",   "H",   "H",   "CL",
  "Q",   "DU",   "K",   "LI",   "LP",
  "NA",   "P3",   "S2",   "S3",   "SO",
  "SO2",   "SI",   "P4",   "D",   "Du",
  "Dum"
};

const std::string aminoacid [] = {
  "ALA", "ASN", "CYS", "GLN", "GLY",
  "ILE", "LEU", "MET", "PHE", "PRO",
  "SER", "THR", "TRP", "TYR", "VAL",
  "ASP", "GLU", "ARG", "HIS", "LYS",
  "DUM"
};

const std::string aaatomname[21][15] = {
  {"N","CA","C","O","CB", "DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //ALA 5  5
  {"N","CA","C","O","CB", "CG", "OD1","ND2","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //ASN 8  13
  {"N","CA","C","O","CB", "SG", "DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //CYS 6  19
  {"N","CA","C","O","CB", "CG", "CD", "OE1","NE2","DUM","DUM","DUM","DUM","DUM","DUM"},  //GLN 9  28
  {"N","CA","C","O","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //GLY 4  32
  {"N","CA","C","O","CB", "CG1","CG2","CD1","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //ILE 7  40
  {"N","CA","C","O","CB", "CG", "CD1","CD2","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //LEU 7  48
  {"N","CA","C","O","CB", "CG", "SD", "CE", "DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //MET 7  56
  {"N","CA","C","O","CB", "CG", "CD1","CD2","CE1","CE2","CZ", "DUM","DUM","DUM","DUM"},  //PHE 11 67
  {"N","CA","C","O","CB", "CG", "CD", "DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //PRO 7  74
  {"N","CA","C","O","CB", "OG", "DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //SER 6  80
  {"N","CA","C","O","CB", "OG1","CG2","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //THR 7  87
  {"N","CA","C","O","CB", "CG", "CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2","DUM"},  //TRP 14 101
  {"N","CA","C","O","CB", "CG", "CD1","CD2","CE1","CE2","CZ", "OH", "DUM","DUM","DUM"},  //TYR 12 113
  {"N","CA","C","O","CB", "CG1","CG2","DUM","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //VAL 7  120
  {"N","CA","C","O","CB", "CG", "OD1","OD2","DUM","DUM","DUM","DUM","DUM","DUM","DUM"},  //ASP 8  128
  {"N","CA","C","O","CB", "CG", "CD", "OE1","OE2","DUM","DUM","DUM","DUM","DUM","DUM"},  //GLU 9  137
  {"N","CA","C","O","CB", "CG", "CD", "NE", "CZ", "NH1","NH2","DUM","DUM","DUM","DUM"},  //ARG 11 148
  {"N","CA","C","O","CB", "CG", "ND1","CD2","CE1","NE2","DUM","DUM","DUM","DUM","DUM"},  //HIS 10 158
  {"N","CA","C","O","CB", "CG", "CD", "CD", "CE", "NZ", "DUM","DUM","DUM","DUM","DUM"},  //LYS 10 168
  {"DUM"}  //DUM
};

const std::string aaatomname_serial[168] = {
  "N","CA","C","O","CB",
  "N","CA","C","O","CB", "CG", "OD1","ND2",
  "N","CA","C","O","CB", "SG", 
  "N","CA","C","O","CB", "CG", "CD", "OE1","NE2",
  "N","CA","C","O",
  "N","CA","C","O","CB", "CG1","CG2","CD1",
  "N","CA","C","O","CB", "CG", "CD1","CD2",
  "N","CA","C","O","CB", "CG", "SD", "CE", 
  "N","CA","C","O","CB", "CG", "CD1","CD2","CE1","CE2","CZ", 
  "N","CA","C","O","CB", "CG", "CD", 
  "N","CA","C","O","CB", "OG", 
  "N","CA","C","O","CB", "OG1","CG2",
  "N","CA","C","O","CB", "CG", "CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2",
  "N","CA","C","O","CB", "CG", "CD1","CD2","CE1","CE2","CZ", "OH", 
  "N","CA","C","O","CB", "CG1","CG2",
  "N","CA","C","O","CB", "CG", "OD1","OD2",
  "N","CA","C","O","CB", "CG", "CD", "OE1","OE2",
  "N","CA","C","O","CB", "CG", "CD", "NE", "CZ", "NH1","NH2",
  "N","CA","C","O","CB", "CG", "ND1","CD2","CE1","NE2",
  "N","CA","C","O","CB", "CG", "CD", "CD", "CE", "NZ",
};

#define DBG 1

#endif
  
