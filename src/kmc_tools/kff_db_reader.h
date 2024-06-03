#ifndef _KFF_DB_READER_H
#define _KFF_DB_READER_H

#include "kmer.h"
#include "defs.h"
#include "config.h"
#include "bundle.h"
#include "kmer_file_header.h"
#include "queues.h"
#include "kff_info_reader.h"
#include "kff_kmc2_reader_utils.h"
#include <cassert>
#include <memory>
#include <algorithm>

//TODO KFF: consider minimizers sections

//Forward declaration
template<unsigned SIZE> class CKFFDbReaderSorted;
template<unsigned SIZE> class CKFFDbReaderSeq;


template<unsigned SIZE>
class CKFFDbReader : public CInput<SIZE>
{
	CPercentProgress& percent_progress;
	uint32 progress_id;
	std::unique_ptr<CKFFDbReaderSorted<SIZE>> db_reader_sorted;
	std::unique_ptr<CKFFDbReaderSeq<SIZE>> db_reader_sequential;
public:
	CKFFDbReader(const CKmerFileHeader& header, const CInputDesc& desc, CPercentProgress& percent_progress, KmerDBOpenMode open_mode):
		percent_progress(percent_progress)

	{			
		uint64_t total_kmers{};
		for (const auto& scope : header.kff_file_struct.scopes)
			for (const auto& s : scope.data_sections)
				total_kmers += s.nb_blocks;

		progress_id = percent_progress.RegisterItem(total_kmers);

		switch (open_mode)
		{
		case KmerDBOpenMode::sorted:
			db_reader_sorted = std::make_unique <CKFFDbReaderSorted<SIZE>>(header, desc);
			break;
		case KmerDBOpenMode::sequential:			
		case KmerDBOpenMode::counters_only:
			db_reader_sequential = std::make_unique<CKFFDbReaderSeq<SIZE>>(header, desc);
			break;		
		default: //should never be here
			std::cerr << "Error: unknow open mode \n";
			exit(1);
		}
	}

	void NextBundle(CBundle<SIZE>& bundle) override
	{
		db_reader_sorted->NextBundle(bundle, this->finished);
		percent_progress.UpdateItem(progress_id, bundle.Size());
		if (this->finished)
		{
			percent_progress.Complete(progress_id);
		}		
	}

	void IgnoreRest() override
	{
		if (this->finished)
			return;
		db_reader_sorted->IgnoreRest();
		this->finished = true;
	}

	bool NextCounter(uint32& counter)
	{
		if (db_reader_sequential->NextCounter(counter))
		{
			percent_progress.UpdateItem(progress_id);
			return true;
		}
		percent_progress.Complete(progress_id);
		return false;
	}

	bool NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter)
	{
		if (db_reader_sequential->NextKmerSequential(kmer, counter))
		{
			percent_progress.UpdateItem(progress_id);
			return true;
		}
		percent_progress.Complete(progress_id);
		return false;					
	}

	~CKFFDbReader()
	{
		//fclose(file);
	}
};

template<unsigned SIZE> class CSection;

struct CSectionBuff
{
	uchar* buf;
	uint32 size;

	CSectionBuff() :
		buf(nullptr), size(0)
	{
	}

	CSectionBuff(uchar* buf, uint32 size) :buf(buf), size(size)
	{

	}

	CSectionBuff& operator=(CSectionBuff&& rhs) noexcept
	{
		if (this != &rhs)
		{
			buf = rhs.buf;
			size = rhs.size;
			rhs.buf = nullptr;
			rhs.size = 0;
		}
		return *this;
	}

	CSectionBuff(CSectionBuff&& rhs) noexcept
	{
		buf = rhs.buf;
		size = rhs.size;
		rhs.buf = nullptr;
		rhs.size = 0;
	}

	CSectionBuff(const CSectionBuff&) = delete;
	CSectionBuff& operator=(const CSectionBuff&) = delete;
};

template<unsigned SIZE>
class CSectionBuffProvider
{
	std::vector<CSectionBuff> internal_bufs;
	uint32 sections_left_to_read = 0;
	uchar* buf_sections, * buf_internal;

	using desc_t = std::tuple<uint64, uint64, uint32, bool>;//in_file_pos, in_file_pos_end, rec_size, is_empty
	using to_read_t = std::tuple<uint32, uint64, uchar*, uint32>;//id, file_pos, buffer to read, size to read

	std::vector<desc_t> desc;
	std::queue<to_read_t, std::list<to_read_t>> to_read;

	mutable std::mutex mtx;
	std::condition_variable cv_pop;
	std::condition_variable cv_get_next_to_read;

	bool forced_to_finish = false;
public:
	void init(std::vector<CSection<SIZE>>& sections)
	{
		uint64 needed_mem = 0;
		internal_bufs.resize(sections.size());

		for (uint32 i = 0 ; i < sections.size(); ++i)
		{
			const auto& section = sections[i];
			auto rec_size = section.GetRecSize();
			uint32 max_section_bytes = SINGLE_SECTION_BUFF_SIZE_FOR_KFF_READER / rec_size * rec_size;

			auto n_kmers = section.GetTotKmers();
			uint32 mem = (uint32)MIN(n_kmers * rec_size, max_section_bytes);

			auto start_pos_in_file = section.GetStartPosInFile();
			auto end_pos_in_file = start_pos_in_file + n_kmers * rec_size;
			desc.push_back(std::make_tuple(start_pos_in_file, end_pos_in_file, rec_size, true));

			internal_bufs[i] = CSectionBuff(nullptr, mem);
			needed_mem += mem;
		}
		

		sections_left_to_read = (uint32)sections.size();

		buf_sections = new uchar[needed_mem];
		buf_internal = new uchar[needed_mem];

		internal_bufs[0].buf = buf_internal;

		uchar* ptr = buf_sections;
		sections[0].set_section_buff(CSectionBuff(ptr, internal_bufs[0].size));

		for (uint32 i = 1; i < internal_bufs.size(); ++i)
		{
			internal_bufs[i].buf = internal_bufs[i - 1].buf + internal_bufs[i - 1].size;
			ptr += internal_bufs[i - 1].size;
			sections[i].set_section_buff(CSectionBuff(ptr, internal_bufs[i].size));
		}

		for (uint32 section_id = 0; section_id < desc.size(); ++section_id)
		{
			auto rec_size = std::get<2>(desc[section_id]);

			uint64 kmers_left = (std::get<1>(desc[section_id]) - std::get<0>(desc[section_id])) / rec_size;
			if (kmers_left)
			{

				uint32 max_section_bytes = SINGLE_SECTION_BUFF_SIZE_FOR_KFF_READER / rec_size * rec_size;

				uint32 kmers_to_read = (uint32)MIN(kmers_left, max_section_bytes / rec_size);
				internal_bufs[section_id].size = kmers_to_read * rec_size;
				to_read.push(std::make_tuple(section_id, std::get<0>(desc[section_id]), internal_bufs[section_id].buf, internal_bufs[section_id].size));
				std::get<0>(desc[section_id]) += kmers_to_read * rec_size;
			}
			else
			{
				--sections_left_to_read;
			}
		}
	}

	void pop(uint32 id, CSectionBuff& section_buf)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this, id] {return !std::get<3>(desc[id]); });

		std::swap(section_buf, internal_bufs[id]);
		std::get<3>(desc[id]) = true;

		uint32 rec_size = std::get<2>(desc[id]);
		uint64 kmers_left = (std::get<1>(desc[id]) - std::get<0>(desc[id])) / rec_size;
		

		if (kmers_left)
		{
			uint32 max_section_bytes = SINGLE_SECTION_BUFF_SIZE_FOR_KFF_READER / rec_size * rec_size;

			uint32 kmers_to_read = (uint32)MIN(kmers_left, max_section_bytes / rec_size);
			internal_bufs[id].size = kmers_to_read * rec_size;
			bool was_empty = to_read.empty();
			to_read.push(std::make_tuple(id, std::get<0>(desc[id]), internal_bufs[id].buf, internal_bufs[id].size));
			std::get<0>(desc[id]) += kmers_to_read * rec_size;
			if (was_empty)
				cv_get_next_to_read.notify_all();
		}
		else
		{
			--sections_left_to_read;
			if (!sections_left_to_read)
				cv_get_next_to_read.notify_all();
		}
	}

	void notify_section_filled(uint32 id)
	{
		std::lock_guard<std::mutex> lck(mtx);
		std::get<3>(desc[id]) = false;
		cv_pop.notify_all();
	}
	bool get_next_to_read(uint32& id, uint64& file_pos, uchar*& buf, uint32& size)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_get_next_to_read.wait(lck, [this] {return !to_read.empty() || !sections_left_to_read || forced_to_finish; });
		if (forced_to_finish || (to_read.empty() && !sections_left_to_read))
			return false;

		std::tie(id, file_pos, buf, size) = to_read.front();
		to_read.pop();
		return true;
	}

	void force_to_finish()
	{
		std::lock_guard<std::mutex> lck(mtx);
		forced_to_finish = true;
		cv_get_next_to_read.notify_all();
	}

	~CSectionBuffProvider()
	{
		delete[] buf_sections;
		delete[] buf_internal;
	}
};


template<unsigned SIZE> class CSectionKmerPQ;
template<unsigned SIZE>
class CSection
{
	uint32 id;
	CSectionBuff section_buf;
	uint32 pos = 0;
	CSectionBuffProvider<SIZE>& section_buf_provider;
	uint64 tot_kmers;
	uint32 kmer_bytes, counter_bytes;
	uint64 start_pos_in_file;
	uint64 kmers_left;
	uint32 cutoff_range;
	uint32 cutoff_min;

	friend class CSectionKmerPQ<SIZE>;
public:
	uint32 GetKmerBytes() const
	{
		return kmer_bytes;
	}

	uint32 GetRecSize() const
	{
		return kmer_bytes + counter_bytes;
	}
	uint64 GetKmersLeft() const
	{
		return kmers_left;
	}
	void DecKmersLeft()
	{
		--kmers_left;
	}
	uint64 GetStartPosInFile() const
	{
		return start_pos_in_file;
	}

	uint64 GetTotKmers() const
	{
		return tot_kmers;
	}

	void ReloadIfNeeded()
	{
		if (pos >= section_buf.size)
			reload_section_buf();
	}

	uchar* GetNextRecord()
	{
		uchar* res = section_buf.buf + pos;
		pos += GetRecSize();
		return res;
	}

	CSection(
		uint32 id, 
		CSectionBuffProvider<SIZE>& section_buf_provider, 
		uint64 tot_kmers, 
		uint32 kmer_bytes, 
		uint32 counter_bytes, 
		uint64 start_pos_in_file, 
		uint32 cutoff_range, 
		uint32 cutoff_min
	) :
		id(id),
		section_buf_provider(section_buf_provider),		
		tot_kmers(tot_kmers),		
		kmer_bytes(kmer_bytes),
		counter_bytes(counter_bytes),
		start_pos_in_file(start_pos_in_file),
		kmers_left(tot_kmers),
		cutoff_range(cutoff_range),
		cutoff_min(cutoff_min)
	{

	}
	void reload_section_buf()
	{
		section_buf_provider.pop(id, section_buf);
		pos = 0;
	}
	bool NextKmer(CKmer<SIZE>& kmer, uint32& counter)
	{
		while (kmers_left)
		{
			ReloadIfNeeded();
			
			uchar* record = GetNextRecord();
			
			kmer.load(record, kmer_bytes); //TODO KFF: consider load_fast
			
			if (counter_bytes == 0)
				counter = 1;
			else
			{
				//TODO KFF: consider faster load like in kmc1 db reader
				counter = 0;
				uchar* counter_ptr = record; //load shifts buffer
				for (uint64_t i = 0; i < counter_bytes; ++i)
					counter += counter_ptr[i] << (8 * (counter_bytes - i - 1));
			}
			
			--kmers_left;		
			if (counter - cutoff_min <= cutoff_range)
				return true;
			
		}
		return false;
	}

	void set_section_buff(CSectionBuff&& _section_buf)
	{
		section_buf = std::move(_section_buf);
		pos = section_buf.size; //force reload
	}
};


//************************************************************************************************************
// CSectionMergerParent - Merger of k-mers produced by CMergerChilds
//************************************************************************************************************
template<unsigned SIZE> class CSectionMergerParent
{
public:
	CSectionMergerParent(std::vector<CCircularQueue<SIZE>*>& input_queues, CCircularQueue<SIZE>& output_queue, uint32 n_subthreads) :
		input_queues(input_queues),
		output_queue(output_queue),
		n_subthreads(n_subthreads)
	{
		input_bundles.resize(input_queues.size());
	}

	void operator()()
	{
		if (n_subthreads > 1)
		{
			ProcessWithSubthreads();
		}

		else if (input_queues.size() == 2)
		{
			Process2Inputs();
		}
		else
		{
			ProcessMoreInputs();
		}
	}

private:
	void Process2Inputs()
	{
		CBundleData<SIZE> b1, b2;
		CCircularQueue<SIZE>* q1, * q2;
		q1 = input_queues[0];
		q2 = input_queues[1];
		bool q1_empty = !q1->pop(b1);
		bool q2_empty = !q2->pop(b2);

		if (q1_empty && q2_empty)
		{
			output_queue.mark_completed();
			return;
		}
		if (q1_empty || q2_empty)
		{
			CCircularQueue<SIZE>* q = q1_empty ? q2 : q1;
			CBundleData<SIZE>& b = q1_empty ? b2 : b1;
			while (true)
			{
				if (!output_queue.push(b))
					break;
				if (!q->pop(b))
					break;
			}
			output_queue.mark_completed();
			return;
		}

		uint32 get1 = 0;
		uint32 get2 = 0;

		CKmer<SIZE> kmer2 = b2.kmers_with_counters[get2].kmer;
		uint32 counter2 = b2.kmers_with_counters[get2].counter;
		CKmer<SIZE> kmer1 = b1.kmers_with_counters[get1].kmer;
		uint32 counter1 = b1.kmers_with_counters[get1].counter;

		uint32 left1 = b1.NRecLeft();
		uint32 left2 = b2.NRecLeft();

		uint32 out_insert_pos = 0;
		uint32 out_size = output_bundle.size;

		while (true)
		{
			if (kmer1 < kmer2)
			{
				output_bundle.kmers_with_counters[out_insert_pos].kmer = kmer1;
				output_bundle.kmers_with_counters[out_insert_pos++].counter = counter1;
				if (out_insert_pos == out_size)
				{
					output_bundle.insert_pos = out_insert_pos;

					if (!output_queue.push(output_bundle))
						break;
					out_insert_pos = 0;
					out_size = output_bundle.size;
				}

				++get1;
				if (--left1)
				{
					kmer1 = b1.kmers_with_counters[get1].kmer;
					counter1 = b1.kmers_with_counters[get1].counter;
				}
				else
				{
					b1.get_pos = get1;
					if (q1->pop(b1))
					{
						get1 = 0;
						kmer1 = b1.kmers_with_counters[get1].kmer;
						counter1 = b1.kmers_with_counters[get1].counter;
						left1 = b1.NRecLeft();
					}
					else
						break;

				}
			}
			else
			{
				output_bundle.kmers_with_counters[out_insert_pos].kmer = kmer2;
				output_bundle.kmers_with_counters[out_insert_pos++].counter = counter2;
				if (out_insert_pos == out_size)
				{
					output_bundle.insert_pos = out_insert_pos;
					if (!output_queue.push(output_bundle))
						break;
					out_insert_pos = 0;
					out_size = output_bundle.size;
				}

				++get2;
				if (--left2)
				{
					kmer2 = b2.kmers_with_counters[get2].kmer;
					counter2 = b2.kmers_with_counters[get2].counter;
				}
				else
				{
					b2.get_pos = get2;
					if (q2->pop(b2))
					{
						get2 = 0;
						kmer2 = b2.kmers_with_counters[get2].kmer;
						counter2 = b2.kmers_with_counters[get2].counter;
						left2 = b2.NRecLeft();
					}
					else
						break;
				}
			}
		}

		b1.get_pos = get1;
		b2.get_pos = get2;
		output_bundle.insert_pos = out_insert_pos;

		if (b1.Empty())
			q1->pop(b1);
		if (!b1.Empty())
		{
			while (true)
			{
				if (b1.Empty())
				{
					if (!q1->pop(b1))
						break;
				}
				output_bundle.Insert(b1.TopKmer(), b1.TopCounter());
				b1.Pop();
				if (output_bundle.Full())
				{
					if (!output_queue.push(output_bundle))
						break;
				}
			}
		}

		if (b2.Empty())
			q2->pop(b2);
		if (!b2.Empty())
		{
			while (true)
			{
				if (b2.Empty())
				{
					if (!q2->pop(b2))
						break;
				}
				output_bundle.Insert(b2.TopKmer(), b2.TopCounter());
				b2.Pop();
				if (output_bundle.Full())
				{
					if (!output_queue.push(output_bundle))
						break;
				}
			}
		}

		if (!output_bundle.Empty())
			output_queue.push(output_bundle);
		output_queue.mark_completed();
	}
	void ProcessMoreInputs()
	{
		//init
		auto q_iter = input_queues.begin();
		auto b_iter = input_bundles.begin();
		for (; q_iter != input_queues.end();)
		{
			if (!(*q_iter)->pop(*b_iter))
			{
				q_iter = input_queues.erase(q_iter);
				b_iter = input_bundles.erase(b_iter);
			}
			else
				++q_iter, ++b_iter;
		}
		uint32 index_of_min = 0;
		CKmer<SIZE> min_kmer;

		uint32 output_bundle_insert_pos = output_bundle.insert_pos;
		uint32 output_bundle_size = output_bundle.size;
		decltype(output_bundle.kmers_with_counters) kmers_counters = output_bundle.kmers_with_counters;

		while (input_bundles.size())
		{
			index_of_min = 0;
			min_kmer = input_bundles[index_of_min].TopKmer();

			if (input_bundles.size() == 4)
			{
				uint32 tmp_min = 2;
				CKmer<SIZE> tmp_kmer = input_bundles[tmp_min].TopKmer();

				if (input_bundles[1].TopKmer() < min_kmer)
				{
					min_kmer = input_bundles[1].TopKmer();
					index_of_min = 1;
				}
				if (input_bundles[3].TopKmer() < tmp_kmer)
				{
					tmp_kmer = input_bundles[3].TopKmer();
					tmp_min = 3;
				}
				if (tmp_kmer < min_kmer)
				{
					index_of_min = tmp_min;
				}
			}
			else if (input_bundles.size() == 3)
			{
				if (input_bundles[1].TopKmer() < min_kmer)
				{
					min_kmer = input_bundles[1].TopKmer();
					index_of_min = 1;
				}
				if (input_bundles[2].TopKmer() < min_kmer)
				{
					min_kmer = input_bundles[2].TopKmer();
					index_of_min = 2;
				}
			}
			else if (input_bundles.size() == 2)
			{
				if (input_bundles[1].TopKmer() < min_kmer)
				{
					min_kmer = input_bundles[1].TopKmer();
					index_of_min = 1;
				}
			}
			else if (input_bundles.size() == 1)
			{
			}
			else
			{
				for (uint32 i = 1; i < input_bundles.size(); ++i)
				{
					if (input_bundles[i].TopKmer() < min_kmer)
					{
						index_of_min = i;
						min_kmer = input_bundles[index_of_min].TopKmer();
					}
				}
			}

			kmers_counters[output_bundle_insert_pos].kmer = input_bundles[index_of_min].TopKmer();
			kmers_counters[output_bundle_insert_pos++].counter = input_bundles[index_of_min].TopCounter();
			//output_bundle.Insert(input_bundles[index_of_min].TopKmer(), input_bundles[index_of_min].TopCounter());

			input_bundles[index_of_min].Pop();
			if (input_bundles[index_of_min].Empty())
			{

				if (!input_queues[index_of_min]->pop(input_bundles[index_of_min]))
				{
					input_queues.erase(input_queues.begin() + index_of_min);
					input_bundles.erase(input_bundles.begin() + index_of_min);
				}
			}


			if (output_bundle_insert_pos == output_bundle_size)
				//if (output_bundle.Full())
			{
				output_bundle.insert_pos = output_bundle_insert_pos;
				if (!output_queue.push(output_bundle))
				{
					output_bundle_insert_pos = output_bundle.insert_pos; //0
					kmers_counters = output_bundle.kmers_with_counters;
					break;
				}
				output_bundle_insert_pos = output_bundle.insert_pos; //0
				kmers_counters = output_bundle.kmers_with_counters;
			}
		}
		output_bundle.insert_pos = output_bundle_insert_pos;
		if (!output_bundle.Empty())
			output_queue.push(output_bundle);
		output_queue.mark_completed();
	}

	uint32 LastLowerPlus1(CBundleData<SIZE>& b, uint32 start, uint32 end, CKmer<SIZE> kmer)
	{
		auto ub = std::upper_bound(b.kmers_with_counters + start, b.kmers_with_counters + end, kmer,
			[](const CKmer<SIZE>& kmer, typename CBundleData<SIZE>::CKmerWithCounter& kc)
		{
			return kmer < kc.kmer;
		});
		return ub - b.kmers_with_counters;
	}
	void ProcessWithSubthreads()
	{
		std::vector<CBundleData<SIZE>> input_bundles(input_queues.size());

		CBundleData<SIZE> output_bundle(KFF_DB_READER_BUNDLE_CAPACITY);
		uint32 curr_end = output_bundle.insert_pos;

		auto q_iter = input_queues.begin();
		auto b_iter = input_bundles.begin();
		for (; q_iter != input_queues.end();)
		{
			if (!(*q_iter)->pop(*b_iter))
			{
				q_iter = input_queues.erase(q_iter);
				b_iter = input_bundles.erase(b_iter);
			}
			else
				++q_iter, ++b_iter;
		}

		std::vector<CParentSubthreadDescQueue<SIZE>> task_descs(n_subthreads);
		std::vector<CParentSubthreadPartDesc> descs(input_queues.size());
		for (uint32 i = 0; i < input_queues.size(); ++i)
		{
			descs[i].start = descs[i].part_end = input_bundles[i].get_pos;
			descs[i].end = input_bundles[i].insert_pos;
		}

		std::vector<CMergerParentSubthread<SIZE>> subtasks;
		std::vector<std::thread> subthreads;

		subtasks.reserve(n_subthreads);
		subthreads.reserve(n_subthreads);

		CParentSubthreadSynchronizer syncer;
		for (uint32 i = 0; i < n_subthreads; ++i)
		{
			subtasks.emplace_back(task_descs[i], syncer);
			subthreads.push_back(std::thread(std::ref(subtasks.back())));
		}

		while (input_bundles.size())
		{
			//prepare threads					
			for (uint32 th = 0; th < n_subthreads; ++th)
			{
				bool any_empty = false;
				for (uint32 i = 0; i < descs.size(); ++i)
				{
					descs[i].part_end = descs[i].start + MIN((output_bundle.size - curr_end) / input_bundles.size() / (n_subthreads - th), descs[i].left());

					//if any is empty it must be refilled or removed
					if (descs[i].part_end == descs[i].start)
					{
						any_empty = true;
						break;
					}
				}
				if (any_empty)
					break;

				//find min kmer
				uint32 min_kmer_i = 0;
				CKmer<SIZE> min_kmer = input_bundles[0].kmers_with_counters[descs[0].part_end - 1].kmer;

				for (uint32 i = min_kmer_i + 1; i < input_bundles.size(); ++i)
				{
					if (input_bundles[i].kmers_with_counters[descs[i].part_end - 1].kmer < min_kmer)
					{
						min_kmer = input_bundles[i].kmers_with_counters[descs[i].part_end - 1].kmer;
						min_kmer_i = i;
					}
				}

				uint32 prev_end = curr_end;

				//correct part_end according to min kmer
				for (uint32 i = 0; i < descs.size(); ++i)
				{
					if (i != min_kmer_i)
					{
						descs[i].part_end = LastLowerPlus1(input_bundles[i], descs[i].start, descs[i].part_end, min_kmer);
					}
					curr_end += descs[i].part_end - descs[i].start;
				}

				task_descs[th].desc.desc = descs;
				task_descs[th].desc.inputs = &input_bundles;
				task_descs[th].desc.out = &output_bundle;
				task_descs[th].desc.o_start = prev_end;
				syncer.increment();
				task_descs[th].start();

				for (uint32 i = 0; i < descs.size(); ++i)
				{
					descs[i].start = descs[i].part_end;
				}
			}

			syncer.wait(); //BARIER

			//send output
			output_bundle.insert_pos = curr_end;
			if (!output_bundle.Empty())
			{
				if (!output_queue.push(output_bundle))
				{
					break;
				}
			}

			curr_end = output_bundle.insert_pos;

			auto q_iter = input_queues.begin();
			auto b_iter = input_bundles.begin();
			auto d_iter = descs.begin();

			for (; b_iter != input_bundles.end();)
			{
				b_iter->get_pos = d_iter->start;
				if ((*b_iter).Empty())
				{
					if ((*q_iter)->pop(*b_iter))
					{
						d_iter->start = d_iter->part_end = b_iter->get_pos;
						d_iter->end = b_iter->insert_pos;
					}
					else
					{
						d_iter = descs.erase(d_iter);
						b_iter = input_bundles.erase(b_iter);
						q_iter = input_queues.erase(q_iter);
						continue;
					}
				}
				++q_iter, ++b_iter, ++d_iter;
			}

		}

		for (auto& t : task_descs)
			t.mark_completed();

		for (auto& t : subthreads)
			t.join();

		output_queue.mark_completed();
	}

	std::vector<CBundleData<SIZE>> input_bundles;
	std::vector<CCircularQueue<SIZE>*>& input_queues;

	CBundleData<SIZE> output_bundle;
	CCircularQueue<SIZE>& output_queue;
	uint32 n_subthreads;
};

template<unsigned SIZE>
class CSectionBlockReader
{
	FILE* file;	
	CSectionBuffProvider<SIZE>& section_buf_provider;
public:
	CSectionBlockReader(FILE* file, CSectionBuffProvider<SIZE>& section_buf_provider) :
		file(file),
		section_buf_provider(section_buf_provider)
	{

	}

	void operator()()
	{
		uint32 id;
		uint64 file_pos;
		uchar* buf;
		uint32 size;

		while (section_buf_provider.get_next_to_read(id, file_pos, buf, size))
		{
			my_fseek(file, file_pos, SEEK_SET);
			if (fread(buf, 1, size, file) != size)
			{
				std::cerr << "Error while reading suffix file\n";
				exit(1);
			}
			section_buf_provider.notify_section_filled(id);
		}
	}
};

//************************************************************************************************************
// CSectionKmerPQ - Priority Queue of k-mers - binary heap. K-mers from sections are processed by this priority queue
//************************************************************************************************************
template<unsigned SIZE> class CSectionKmerPQ
{
public:
	CSectionKmerPQ(uint32 _no_of_sections)
	{
		elems.resize(_no_of_sections + 1);
		descs.resize(_no_of_sections + 1);
		pos = 1;
		desc_pos = 0;
	}
	inline void init_add(CSection<SIZE>* section)
	{
		CKmer<SIZE> kmer;
		uint32 counter;
		if (section->NextKmer(kmer, counter))
		{
			descs[desc_pos] = std::make_pair(section, counter);
			elems[pos] = std::make_pair(kmer, desc_pos);
			uint32 child_pos = pos++;

			while (child_pos > 1 && elems[child_pos].first < elems[child_pos / 2].first)
			{
				swap(elems[child_pos], elems[child_pos / 2]);
				child_pos /= 2;
			}

			++desc_pos;
		}
	}

	void Process(CCircularQueue<SIZE>& output_queue)
	{
		if (pos <= 1)
		{
			output_queue.mark_completed();
			return;
		}
		CBundleData<SIZE> bundle_data;
		uint32 desc_id = 0;
		CSection<SIZE>* section = descs[elems[1].second].first;
		CKmer<SIZE> kmer;
		uint32 counter;
		uchar* record = nullptr;

		
		//uint32 suffix_bytes = section->suffix_bytes;
		uint32 counter_size = section->counter_bytes;
		//uint32 counter_mask = section->counter_mask;
		//uint32 record_size = section->record_size;
		uint32 cutoff_min = section->cutoff_min;
		uint32 cutoff_range = section->cutoff_range;

		//bool endian = CConfig::GetInstance().IsLittleEndian();

		while (true)
		{
			if (pos <= 1)
				break;
			bundle_data.Insert(elems[1].first, descs[elems[1].second].second);


			//UPDATE HEAP!
			desc_id = elems[1].second;
			section = descs[desc_id].first;

			bool exists = false;


			while (section->GetKmersLeft())
			{
				section->ReloadIfNeeded();
				
				record = section->GetNextRecord();
				uint32 kmer_bytes = section->GetKmerBytes();

				kmer.load(record, kmer_bytes); //TODO KFF: consider load_fast

				if (counter_size == 0)
					counter = 1;
				else
				{
					//TODO KFF: consider faster load like in kmc1/kmc2 db reader
					counter = 0;
					uchar* counter_ptr = record;
					for (uint64_t i = 0; i < counter_size; ++i)
						counter += counter_ptr[i] << (8 * (counter_size - i - 1));
				}

//				kmer.load_fast(record, suffix_bytes, endian);
//				kmer.set_prefix(bin->prefix, suffix_bytes);

				//CCounterBuilder::build_counter(counter, record, counter_size, counter_mask, endian);
				//bin->pos += record_size;

				section->DecKmersLeft();
				
				if (counter - cutoff_min <= cutoff_range)
				{
					exists = true;
					break;
				}
			}

			if (!exists)
			{
				kmer.set(elems[--pos].first);
				desc_id = elems[pos].second;
			}
			else
				descs[desc_id].second = counter;

			uint32 parent, less;
			parent = less = 1;
			while (true)
			{
				if (parent * 2 >= pos)
					break;
				if (parent * 2 + 1 >= pos)
					less = parent * 2;
				else if (elems[parent * 2].first < elems[parent * 2 + 1].first)
					less = parent * 2;
				else
					less = parent * 2 + 1;
				if (elems[less].first < kmer)
				{
					elems[parent] = elems[less];
					parent = less;
				}
				else
					break;
			}
			elems[parent] = std::make_pair(kmer, desc_id);

			if (bundle_data.Full())
			{
				if (!output_queue.push(bundle_data))
					break;
			}
		}
		if (!bundle_data.Empty())
			output_queue.push(bundle_data);
		output_queue.mark_completed();
	}

	inline void reset()
	{
		pos = 1;
		desc_pos = 0;
	}

private:
	using elem_t = std::pair<CKmer<SIZE>, uint32>;//kmer, desc_id
	using desc_t = std::pair<CSection<SIZE>*, uint32>;//bin, counter
	std::vector<elem_t> elems;
	std::vector<desc_t> descs;
	uint32 pos, desc_pos;
};

//************************************************************************************************************
// CSectionMergerChild - Merger of k-mers from sections
//************************************************************************************************************
template<unsigned SIZE> class CSectionMergerChild
{
	using section_iter = typename std::vector<CSection<SIZE>>::iterator;
public:
	CSectionMergerChild(section_iter begin, section_iter end, CCircularQueue<SIZE>& output_queue) :
		sections(begin, end),
		output_queue(output_queue)
	{

	}
	void operator()()
	{
		CSectionKmerPQ<SIZE> kmers_pq(static_cast<uint32>(sections.size()));
		for (uint32 i = 0; i < sections.size(); ++i)
			kmers_pq.init_add(&sections[i].get());

		kmers_pq.Process(output_queue);
	}

private:
	std::vector<std::reference_wrapper<CSection<SIZE>>> sections;
	CCircularQueue<SIZE>& output_queue;
};

//************************************************************************************************************
// CKFFDbReaderSorted - Produce k-mers in sorted order from KFF database
//************************************************************************************************************
template<unsigned SIZE> class CKFFDbReaderSorted
{
	FILE* file;
	std::vector<CSection<SIZE>> sections;
	CSectionBuffProvider<SIZE> section_buf_provider;
	std::unique_ptr<CSectionBlockReader<SIZE>> section_block_reader;
	std::thread section_block_reader_th;

	uint32 n_child_threads;
	uint32 n_parent_threads;

	CCircularQueue<SIZE>* output_queue;
	std::vector<CCircularQueue<SIZE>*> childs_parent_queues;

	CSectionMergerParent<SIZE>* parent = nullptr;
	std::thread parent_thread;

	std::vector<CSectionMergerChild<SIZE>*> childs;
	std::vector<std::thread> childs_threads;

	friend class CSection<SIZE>;
public:
	CKFFDbReaderSorted(const CKmerFileHeader& header, const CInputDesc& desc);

	void NextBundle(CBundle<SIZE>& bundle, bool& finished) 
	{
		if (output_queue->pop(bundle.Data()))
		{
			return;
		}

		for (auto& child_thread : childs_threads)
			child_thread.join();
		for (auto& child : childs)
			delete child;

		if (parent_thread.joinable()) //for small number of threads there is only one child thread and parent threads is not needed
			parent_thread.join();
		if (parent) //as above
			delete parent;

		delete output_queue;

		for (auto& q : childs_parent_queues)
			delete q;

		section_block_reader_th.join();

		
		finished = true;
	}

	void IgnoreRest() 
	{
		output_queue->force_finish();

		for (auto& q : childs_parent_queues)
			q->force_finish();

		for (auto& child_thread : childs_threads)
			child_thread.join();
		for (auto& child : childs)
			delete child;

		if (parent_thread.joinable()) //for small number of threads there is only one child thread and parent threads is not needed
			parent_thread.join();

		if (parent) //as above
			delete parent;

		delete output_queue;

		for (auto& q : childs_parent_queues)
			delete q;

		
		section_buf_provider.force_to_finish();
		section_block_reader_th.join();
	}

	~CKFFDbReaderSorted() 
	{
		fclose(file);
	}
};

template<unsigned SIZE>
CKFFDbReaderSorted<SIZE>::CKFFDbReaderSorted(const CKmerFileHeader& header, const CInputDesc& desc)
{
	file = fopen(desc.file_src.c_str(), "rb");
	if (!file)
		file = fopen((desc.file_src + ".kff").c_str(), "rb");

	if (!file)
	{
		std::cerr << "Error: cannot open file " << desc.file_src << "\n";
		exit(1);
	}
	setvbuf(file, NULL, _IONBF, 0);

	uint32 section_id{};
	uint32 cutoff_range = desc.cutoff_max - desc.cutoff_min;
	uint32 cutoff_min = desc.cutoff_min;

	for (const auto& scope : header.kff_file_struct.scopes)
	{
		for (const auto& section : scope.data_sections)
		{			
			uint32 kmer_bytes = (scope.kmer_size + 3) / 4;
			sections.emplace_back(section_id++, section_buf_provider, section.nb_blocks, kmer_bytes, scope.data_size, section.data_start_pos, cutoff_range, cutoff_min);
		}
	}

	section_buf_provider.init(sections);

	section_block_reader = std::make_unique<CSectionBlockReader<SIZE>>(file, section_buf_provider);
	section_block_reader_th = std::thread(std::ref(*section_block_reader.get()));

	uint32 n_threads = desc.threads;

	if (n_threads < 3)
	{
		n_threads = n_child_threads = 1;
		output_queue = new CCircularQueue<SIZE>(DEFAULT_CIRCULAL_QUEUE_CAPACITY);
		childs.push_back(new CSectionMergerChild<SIZE>(sections.begin(), sections.end(), *output_queue));
		childs_threads.push_back(std::thread(std::ref(*childs.front())));
		return;
	}

	else if (n_threads == 3)
	{
		n_child_threads = 2;
		n_parent_threads = 1;
	}
	//based on experiment on 24 core machine
	else if (n_threads < 6)
	{
		n_child_threads = 3;
		n_parent_threads = n_threads - n_child_threads;
	}
	else if (n_threads < 9)
	{
		n_child_threads = 4;
		n_parent_threads = n_threads - n_child_threads;
	}
	else if (n_threads < 11)
	{
		n_child_threads = 5;
		n_parent_threads = n_threads - n_child_threads;
	}
	else if (n_threads < 14)
	{
		n_child_threads = 6;
		n_parent_threads = n_threads - n_child_threads;
	}
	else if (n_threads < 17)
	{
		n_child_threads = 7;
		n_parent_threads = n_threads - n_child_threads;
	}
	else
	{
		n_child_threads = (n_threads - 17) / 5 + 8;
		n_parent_threads = n_threads - n_child_threads;
	}
	childs_parent_queues.reserve(n_child_threads);
	childs.reserve(n_child_threads);

	uint32 bundle_size = KFF_DB_READER_BUNDLE_CAPACITY;

	if (n_parent_threads < 2)
	{
		bundle_size = BUNDLE_CAPACITY;
	}

	for (uint32 i = 0; i < n_child_threads; ++i)
		childs_parent_queues.push_back(new CCircularQueue<SIZE>(DEFAULT_CIRCULAL_QUEUE_CAPACITY, bundle_size));


	uint32 sections_per_thread = sections.size() / n_child_threads;

	for (uint32 i = 0; i < n_child_threads - 1; ++i)
	{
		childs.push_back(new CSectionMergerChild<SIZE>(sections.begin() + i * sections_per_thread, sections.begin() + (i + 1) * sections_per_thread, *childs_parent_queues[i]));
		childs_threads.push_back(std::thread(std::ref(*childs.back())));
	}

	//last one
	childs.push_back(new CSectionMergerChild<SIZE>(sections.begin() + (n_child_threads - 1) * sections_per_thread, sections.end(), *childs_parent_queues.back()));
	childs_threads.push_back(std::thread(std::ref(*childs.back())));

	output_queue = new CCircularQueue<SIZE>(DEFAULT_CIRCULAL_QUEUE_CAPACITY, bundle_size);

	parent = new CSectionMergerParent<SIZE>(childs_parent_queues, *output_queue, n_parent_threads);
	parent_thread = std::thread(std::ref(*parent));

}


//************************************************************************************************************
// CKFFDbReaderSeqCounter_Base - Base class for classes to access k-mers one by one (not sorted) or 
// for counters only from KFF database
//************************************************************************************************************
template <unsigned SIZE>
class CKFFDbReaderSeq
{

	const CKmerFileHeader& header;
	const CInputDesc& desc;

	uint64 buff_bytes{};	
	uint64 buff_pos{}; //in_recs
	uchar* buff = nullptr;

	FILE* file{};

	struct CSectionInternal
	{	
		uint64_t kmers_left;
		uint64_t kmer_bytes;
		uint64_t counter_bytes;
		uint64_t start_in_file;
		bool started = false;
		
	public:
		CSectionInternal(uint64_t n_kmers, uint64_t kmer_bytes, uint64_t counter_bytes, uint64_t start_in_file) :
			kmers_left(n_kmers),
			kmer_bytes(kmer_bytes),
			counter_bytes(counter_bytes),
			start_in_file(start_in_file)
		{

		}
	};
	std::vector<CSectionInternal> internal_sections;
	uint64_t cur_section_id{};

	bool reload_buff()
	{
		while (cur_section_id < internal_sections.size())
		{
			auto& cur_section = internal_sections[cur_section_id];
			if (cur_section.kmers_left == 0)
			{
				++cur_section_id;
				continue;
			}
			if (!cur_section.started)
			{
				fseek(file, cur_section.start_in_file, SEEK_SET);
				cur_section.started = true;
			}

			auto rec_size_bytes = cur_section.kmer_bytes + cur_section.counter_bytes;
			uint64_t max_rec_to_read = KFF_DB_READER_BUFF_BYTES / rec_size_bytes;
			uint64_t recs_left = cur_section.kmers_left;

			auto buff_recs = MIN(max_rec_to_read, recs_left);

			buff_bytes = buff_recs * rec_size_bytes;
			fread(buff, 1, buff_bytes, file);
			cur_section.kmers_left -= buff_recs;
			buff_pos = 0;
			return true;
		}
		return false;
	}


	bool buff_empty()
	{
		return buff_pos == buff_bytes;
	}

	void read_kmer_from_buff(CKmer<SIZE>& kmer, uint32& counter)
	{
		auto kmer_bytes = internal_sections[cur_section_id].kmer_bytes;
		auto counter_bytes = internal_sections[cur_section_id].counter_bytes;
		
		auto record = buff + buff_pos;
		kmer.load(record, kmer_bytes);

		if (counter_bytes == 0)
			counter = 1;
		else
		{
			//TODO KFF: consider faster load like in kmc1/2 db reader
			counter = 0;
			uchar* counter_ptr = record; 
			for (uint64_t i = 0; i < counter_bytes; ++i)
				counter += counter_ptr[i] << (8 * (counter_bytes - i - 1));
		}

		buff_pos += kmer_bytes + counter_bytes;
	}

	void read_counter_from_buff(uint32& counter)
	{
		auto kmer_bytes = internal_sections[cur_section_id].kmer_bytes;
		auto counter_bytes = internal_sections[cur_section_id].kmer_bytes;

		auto record = buff + buff_pos + kmer_bytes;
		
		if (counter_bytes == 0)
			counter = 1;
		else
		{			
			//TODO KFF: consider faster load like in kmc1/2 db reader
			counter = 0;
			uchar* counter_ptr = record; //load shifts buffer
			for (uint64_t i = 0; i < counter_bytes; ++i)
				counter += counter_ptr[i] << (8 * (counter_bytes - i - 1));
		}

		buff_pos += kmer_bytes + counter_bytes;
	}
public:
	CKFFDbReaderSeq(const CKmerFileHeader& header, const CInputDesc& desc):
		header(header),
		desc(desc)
	{
		file = fopen(desc.file_src.c_str(), "rb");
		if (!file)
			file = fopen((desc.file_src + ".kff").c_str(), "rb");

		if (!file)
		{
			std::cerr << "Error: cannot open file " << desc.file_src << "\n";
			exit(1);
		}
		setvbuf(file, NULL, _IONBF, 0);
		buff = new uchar[KFF_DB_READER_BUFF_BYTES];	
				
		for (const auto& scope : header.kff_file_struct.scopes)
		{			
			for (const auto& s : scope.data_sections)
			{		
				uint64 kmer_bytes = (scope.kmer_size + 3) / 4;				
				internal_sections.emplace_back(s.nb_blocks, kmer_bytes, scope.data_size, s.data_start_pos);
			}
		}		
	}
	~CKFFDbReaderSeq()
	{		
		fclose(file);
		delete[] buff;
	}
	
	bool NextKmerSequential(CKmer<SIZE>& kmer, uint32& counter)
	{
		while (true)
		{
			if (buff_empty())
				if (!reload_buff())
					return false;

			read_kmer_from_buff(kmer, counter);			
			if (counter >= this->desc.cutoff_min && counter <= this->desc.cutoff_max)
				return true;
		}	
	}

	bool NextCounter(uint32& counter)
	{
		while (true)
		{
			if (buff_empty())
				if (!reload_buff())
					return false;

			read_counter_from_buff(counter);
			if (counter >= this->desc.cutoff_min && counter <= this->desc.cutoff_max)
				return true;
		}
	}
};

#endif