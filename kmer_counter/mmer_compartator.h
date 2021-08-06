//
// Created by danflomin on 08/02/2021.
//

#ifndef KMC_MMER_COMPARTATOR_H
#define KMC_MMER_COMPARTATOR_H

#include "defs.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unistd.h>




class CMmerNorm
{
    static uint32 norm[1<<18];

    class MmerComp
    {
        uint32 stats[1<<18];//uint32* stats;
        uint32 uhs[1<<18];
        uint32 my_norm[1<<18];
        char num_codes[256];
    public:
        bool operator()(uint32 mmer1_idx, uint32 mmer2_idx)
        {
            uint32 mmer1 = my_norm[mmer1_idx];
            uint32 mmer2 = my_norm[mmer2_idx];

            if(uhs[mmer1] == 1 && uhs[mmer2] == 0)
                return true;
            if(uhs[mmer1] == 0 && uhs[mmer2] == 1)
                return false;

            if(is_allowed(mmer1) && !is_allowed(mmer2))
                return true;
            if(!is_allowed(mmer1) && is_allowed(mmer2))
                return false;

            if(stats[mmer1] > stats[mmer2])
                return false;
            if(stats[mmer1] < stats[mmer2])
                return true;

            // frequency is equal
            return mmer1 < mmer2;
        }
        MmerComp(uint32* _signature_occurrences)
        {
            // stats = _signature_occurrences;
            for (int j = 0; j < (1<<18); j++)
            {
                stats[j] = _signature_occurrences[j] + _signature_occurrences[get_reversed(j)];
            }

            for (int i = 0; i < 256; i++)
                num_codes[i] = -1;
            num_codes['A'] = num_codes['a'] = 0;
            num_codes['C'] = num_codes['c'] = 1;
            num_codes['G'] = num_codes['g'] = 2;
            num_codes['T'] = num_codes['t'] = 3;
            std::cout << "before  bitset" << std::endl << std::flush;
            init_uhs_bitset();
            std::cout << "after bitset" << std::endl<< std::flush;

            for (int j = 0; j < (1<<18); j++)
            {
                my_norm[j] = norm[j];
            }
            for (int j = 0; j < (1<<18); j++)
            {
                _signature_occurrences[j] = 0;
                if(uhs[j] == 1)
                    _signature_occurrences[MAX(j,get_reversed(j))] = stats[j];       
                    _signature_occurrences[MIN(j,get_reversed(j))] = 0;       
                            
            }
        }
        void init_uhs_bitset()
        {
            for(int i = 0; i< 1 << 18; i++)
            {
                uhs[i] = 0;
            }


            std::ifstream infile("res_9.txt");
            int str_value;
            std::string line;
            std::cout << "lolz0" << std::endl<< std::flush;
            while (std::getline(infile, line))
            {
                const char *seq = line.c_str();
                str_value = (num_codes[seq[0]] << 16) + (num_codes[seq[1]] << 14) + (num_codes[seq[2]] << 12) + (num_codes[seq[3]] << 10) + (num_codes[seq[4]] << 8) + (num_codes[seq[5]] << 6) + (num_codes[seq[6]] << 4) + (num_codes[seq[7]]<< 2) + (num_codes[seq[8]]);
                uhs[str_value] = 1;
            }
            std::cout << "lolz6" << std::endl<< std::flush;
            infile.close();
            std::cout << "lolz7" << std::endl<< std::flush;
        }
        bool is_allowed(uint32 mmer)
        {
            if ((mmer & 0x3f) == 0x3f)            // TTT suffix
                return false;
            if ((mmer & 0x3f) == 0x3b)            // TGT suffix
                return false;
            if ((mmer & 0x3c) == 0x3c)            // TG* suffix !!!! consider issue #152
                return false;

            for (uint32 j = 0; j < 9 - 3; ++j)
                if ((mmer & 0xf) == 0)                // AA inside
                    return false;
                else
                    mmer >>= 2;

            if (mmer == 0)            // AAA prefix
                return false;
            if (mmer == 0x04)        // ACA prefix
                return false;
            if ((mmer & 0xf) == 0)    // *AA prefix
                return false;

            return true;
        }
        
        uint32 get_reversed(uint32 mmer)
    {
        uint32 rev = 0;
        uint32 shift = 9*2 - 2;
        for(uint32 i = 0 ; i < 9 ; ++i)
        {
            rev += (3 - (mmer & 3)) << shift;
            mmer >>= 2;
            shift -= 2;
        }
        return rev;
    }


    };

    void merge(uint32 arr[], int l, int m, int r, MmerComp* comp)
    {
        int n1 = m - l + 1;
        int n2 = r - m;

        // Create temp arrays
        int L[n1], R[n2];

        // Copy data to temp arrays L[] and R[]
        for (int i = 0; i < n1; i++)
            L[i] = arr[l + i];
        for (int j = 0; j < n2; j++)
            R[j] = arr[m + 1 + j];

        // Merge the temp arrays back into arr[l..r]
        // Initial index of first subarray
        int i = 0;
        // Initial index of second subarray
        int j = 0;
        // Initial index of merged subarray
        int k = l;

        while (i < n1 && j < n2) {
            if (comp->operator()(L[i],R[j]) || norm[L[i]] == norm[R[j]]) {
                arr[k] = L[i];
                i++;
            }
            else {
                arr[k] = R[j];
                j++;
            }
            k++;
        }

        // Copy the remaining elements of
        // L[], if there are any
        while (i < n1) {
            arr[k] = L[i];
            i++;
            k++;
        }

        // Copy the remaining elements of
        // R[], if there are any
        while (j < n2) {
            arr[k] = R[j];
            j++;
            k++;
        }
    }
    void mergeSort(uint32 arr[],int l,int r, MmerComp* comp){
        if(l>=r){
            return;//returns recursively
        }
        int m =l+ (r-l)/2;
        mergeSort(arr,l,m, comp);
        mergeSort(arr,m+1,r, comp);
        merge(arr,l,m,r, comp);
    }


public:
    void Init(uint32* stats)
    {
        std::cout << "in INITTT" << std::endl << std::flush;
        uint32 temp[1 << 18];
        uint32 special = 1 << 18;
        uint32 i;

//        TODO: TAKE CARE OF CANONICAL KMERS. COMBINE stats[i] and stats[rev]. norm[i] == norm[rev]. I think that now they are 1 idx apart.

        for(i = 0 ; i < special ; ++i)
        {
            uint32 rev = get_rev(i);
            norm[i] = MIN(i, rev);
            temp[i] = i;
        }

        MmerComp* comp = new MmerComp(stats);
        mergeSort(temp, 0, special-1, comp);
        std::cout << "after sort" << std::endl;

        for(uint32 i = 0 ; i < special ; ++i)
        {
            temp[i] = norm[temp[i]];
        }
        std::cout << "after seconds for loop" << std::endl;
        for(uint32 i = 0 ; i < special ; ++i)
        {
            norm[i] = temp[i];
        }
        std::cout << "after reverse normalizing for loop" << std::endl;
        for(uint32 i = 0 ; i < special ; ++i)
        {
            uint32 rev = get_rev(i);
            uint32 min_value = MIN(norm[i], norm[rev]);
            norm[i] = min_value;
            norm[rev] = min_value;
        }
        std::cout << "after final for loop" << std::endl;
    }

    CMmerNorm(uint32* stats);

    inline uint32 get_norm_value(uint32 value)
    {
        return norm[value];
    }
    inline uint32 get_rev(uint32 mmer)
    {
        uint32 rev = 0;
        uint32 shift = 9*2 - 2;
        for(uint32 i = 0 ; i < 9 ; ++i)
        {
            rev += (3 - (mmer & 3)) << shift;
            mmer >>= 2;
            shift -= 2;
        }
        return rev;
    }
};

// ***** EOF


#endif //KMC_MMER_COMPARTATOR_H


