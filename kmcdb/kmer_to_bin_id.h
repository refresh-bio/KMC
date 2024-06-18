#ifndef KMER_TO_BIN_ID_H_
#define KMER_TO_BIN_ID_H_

#include "mmer.h"
#include "hashers.h"

namespace kmcdb
{
	namespace detail
	{
		//mkokot_TODO: na razie sobie to tak pisze, ale pewnie trzeba bedzie to gdzies przeneisc
		template<unsigned SIZE, typename MMER_T>
		uint64_t get_signature(const CKmer<SIZE>& kmer, uint64_t kmer_len, uint64_t mmer_len)
		{
			MMER_T min_mmer(static_cast<uint32_t>(mmer_len));
			MMER_T cur_mmer(static_cast<uint32_t>(mmer_len));

		
			auto pos = static_cast<uint32_t>(2 * kmer_len - 2);
			for (uint32_t i = 0; i < mmer_len; ++i, pos -= 2)
				cur_mmer.insert(kmer.get_2bits(pos));

			MMER_T min_mmr(cur_mmer);
			for (uint64_t i = mmer_len; i < kmer_len; ++i, pos -=2)
			{
				cur_mmer.insert(kmer.get_2bits(pos));

				if (cur_mmer < min_mmr)
					min_mmr = cur_mmer;
			}
			return min_mmr.get();
		}

		//mkokot_TODO: move to base class
		template<unsigned SIZE>
		uint64_t get_signature(
			SignatureSelectionScheme signature_selection_scheme,
			const CKmer<SIZE>& kmer,
			uint64_t kmer_len,
			uint64_t signature_len
		)
		{
			switch (signature_selection_scheme)
			{
			case SignatureSelectionScheme::MinHash:
				return detail::get_signature<SIZE, MmerMinHash<MurMur64Hash>>(kmer,
					kmer_len, signature_len);
			}
			throw std::runtime_error("signature_selection_scheme not covered by switch statement");
		}

		inline uint64_t get_bin_id(
			SignatureToBinMapping signature_to_bin_mapping,
			uint64_t signature,
			uint64_t num_bins)
		{
			switch (signature_to_bin_mapping)
			{
			case SignatureToBinMapping::Modulo: return signature % num_bins;
			}
			throw std::runtime_error("signature_to_bin_mapping not covered by switch statement");
		}
	}

	template<unsigned SIZE>
	uint64_t get_bin_id(const CKmer<SIZE>& kmer,
		uint64_t kmer_len,
		uint64_t signature_len,
		uint64_t num_bins,
		SignatureSelectionScheme signature_selection_scheme,
		SignatureToBinMapping signature_to_bin_mapping)
	{
		return detail::get_bin_id(signature_to_bin_mapping,
			detail::get_signature(signature_selection_scheme, kmer, kmer_len, signature_len),
			num_bins);
	}
}

#endif // ! KMER_TO_BIN_ID_H_