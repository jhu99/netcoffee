/* macro.h
Author: Jialu Hu
Date: 28.06.2012*/

#ifndef MACRO_H_
#define MACRO_H_
const unsigned int MAXIMUM_RECORDS = 500000;//maximum number of records of alignment nodes
const unsigned int NUM_SPECIES=4;
const int LISTSIZE=3;
const unsigned int NUM_INITIAL=1000;// optional; initial step can also be controlled by inner variables.
const unsigned int NUM_ALIGNMENT_NODE=8000;
const int NUM_NODE_REVERSE=3000;
const int NUM_EDGE_REVERSE=10000;
const int NUM_PRO_BARS=95;
const int NUM_PRO_INTERVAL=2;
const int NUM_OFFSET_BARS=90;
const int NUM_GO_SHARED=10;
const float TEMP_MAX=100;
const float TEMP_MIN=10;
const unsigned K_MAX=100;//t_coffee 1000
const int N_MAX=2000;
const float T_CON=0.005;
const float _T_CON=0.003;
const float T_ALPHA=0.5;
const float FACTOR_EDGE=0.10;
const float T_THRESHOLD=0.35;
#endif //MACRO_H_
