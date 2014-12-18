/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.1
  Date   : 2014-12-18
*/

#ifndef _S_MAPPER_H
#define _S_MAPPER_H
#include "defs.h"
#include "mmer.h"
#include "params.h"

#ifdef DEVELOP_MODE
#include "develop.h"
#endif


class CSignatureMapper
{
	uint32 map_size;
	int32* signature_map;
	uint32 signature_len;
	uint32 special_signature;
	CMemoryPool* pmm_stats;
	uint32 n_bins;

	class Comp
	{
		uint32* signature_occurences;
	public:
		Comp(uint32* _signature_occurences) : signature_occurences(_signature_occurences){}
		bool operator()(int i, int j)
		{
			return signature_occurences[i] > signature_occurences[j];
		}
	};
	
public:	
	void Init(uint32* stats)
	{
		uint32 *sorted;
		pmm_stats->reserve(sorted);
		for (uint32 i = 0; i < map_size ; ++i)
			sorted[i] = i;
		sort(sorted, sorted + map_size, Comp(stats));

		list<pair<uint32, uint64>> _stats;
		for (uint32 i = 0; i < map_size ; ++i)
		{
			if (CMmer::is_allowed(sorted[i], signature_len))
				_stats.push_back(make_pair(sorted[i], stats[sorted[i]]));
		}

		list<pair<uint32, uint64>> group;
		uint32 bin_no = 0;
		//counting sum
		double sum = 0.0;
		for (auto &i : _stats)
		{
			i.second += 1000;
			sum += i.second;
		}

		double mean = sum / n_bins;
		double max_bin_size = 1.1 * mean;
		uint32 n = n_bins - 1; //one is needed for disabled signatures
		uint32 max_bins = n_bins - 1;
		while (_stats.size() > n)
		{
			pair<uint32, uint64>& max = _stats.front();

			if (max.second > mean)
			{
				signature_map[max.first] = bin_no++;				
				sum -= max.second;
				mean = sum / (max_bins - bin_no);
				max_bin_size = 1.1 * mean;

				_stats.pop_front();
				--n;
			}
			else
			{
				//heuristic
				group.clear();
				double tmp_sum = 0.0;
				uint32 in_current = 0;
				for (auto it = _stats.begin(); it != _stats.end();)
				{
					if (tmp_sum + it->second < max_bin_size)
					{
						tmp_sum += it->second;
						group.push_back(*it);
						it = _stats.erase(it);
						++in_current;
					}
					else
						++it;
				}

				for (auto i = group.begin(); i != group.end(); ++i)
				{
					signature_map[i->first] = bin_no;
				}
				--n;
				++bin_no;

				sum -= tmp_sum;
				mean = sum / (max_bins - bin_no);
				max_bin_size = 1.1 * mean;
			}
		}
		if (_stats.size() > 0)
		{
			for (auto i = _stats.begin(); i != _stats.end(); ++i)
			{
				signature_map[i->first] = bin_no++;
				//cout << "rest bin: " << i->second << "\n";
			}
		}
		signature_map[special_signature] = bin_no;
		pmm_stats->free(sorted);

#ifdef DEVELOP_MODE
		map_log(signature_len, map_size, signature_map);
#endif

	}
	CSignatureMapper(CMemoryPool* _pmm_stats, uint32 _signature_len, uint32 _n_bins)
	{
		n_bins = _n_bins;
		pmm_stats = _pmm_stats;
		signature_len = _signature_len;
		special_signature = 1 << 2 * signature_len;
		map_size = (1 << 2 * signature_len) + 1;
		signature_map = new int32[map_size];		
		fill_n(signature_map, map_size, -1);
	}
	inline int32 get_bin_id(uint32 signature)
	{
		return signature_map[signature];
	}

	inline int32 get_max_bin_no()
	{
		return signature_map[special_signature];
	}

	~CSignatureMapper()
	{
		delete [] signature_map;
	}

};

#endif 