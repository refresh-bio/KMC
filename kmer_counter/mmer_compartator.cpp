//
// Created by danflomin on 08/02/2021.
//
#include "mmer_compartator.h"
#include <iostream>

uint32 CMmerNorm::norm[];

CMmerNorm::CMmerNorm(uint32* stats)
{
    std::cout << "before stats" << std::endl;
    Init(stats);
    std::cout << "after stats" << std::endl;
}
