#ifndef WRITERS_H_
#define WRITERS_H_
#include "bin_writers.h"
#include "capture_ostream.h"
#include "capture_sys_hw_cmd.h"
#include <iostream>

//mkokot_TODO: poniewaz generalnie te klasy nie wszystko kontroluja
//moze warto by bylo dodac taki wariant ktory dziala dla dlugosci k-mera = 0
//wtedy nie przechowujemy k-merow tylko cokolwiek chcemy
//np jakies sumaryczne dan
//albo jakbysmy na przyklad mieli jedna baze z k-merami i danymi, a potem dla tych samych k-merow chcieli jakies inne dane
//to oczywisce zrzuca na uzytkownika koniecznosc zachowania odpowiedniej kolejnosci
//no ale tak samo ma jak np. robi dwa wektory w ktorych odpowiadajace w jakis sposob elemetny maja byc na tych samych pozycjach
//w wariantach odczytu random access moze byc warto dodac metody GetIth, ktore zwroca i-ty k-mer w binie
//w wariantach odczytu listing moze byc warto dodac metode Skip, albo ogolniejsza, SetToIth
//ktora by ustawiala na i-ty k-mer i iterowanie by lecialo od tego wlasnie ustawienia
//to by mozna bylo zastosowac jakby sie chcialo np. losowo samplowac k-mery zeby nie bylo koniecznosci iterowania po wszystkich
namespace kmcdb
{
	namespace detail
	{
		class HistoryWriter
		{
			archive_output_t* archive;
			int stream_id;
			HistoryItem cur_item;

			std::unique_ptr<detail::CaptureOstream> capture_std_cout;
			std::unique_ptr<detail::CaptureOstream> capture_std_cerr;

			bool is_finished = false;
		private:
			static uint64_t get_cur_time()
			{
				return std::chrono::duration_cast<std::chrono::milliseconds>(
					std::chrono::system_clock::now().time_since_epoch()).count();
			}

			//returns true on success, i.e. the state of capturing was changed
			bool capture_ostream_impl(
				bool do_capture,
				std::unique_ptr<detail::CaptureOstream>& inst,
				std::ostream& stream,
				std::string& capture_to)
			{
				if (do_capture)
				{
					//if is already capturing do noting
					if (inst)
						return false;
					inst = std::make_unique<detail::CaptureOstream>(stream, capture_to);
					return true;
				}

				//if was not capturing do nothing
				if (!inst)
					return false;

				//by releasing the pointer we stop capturing
				inst.reset();
				return true;
			}

			void capture_command_line_system_and_hw_info()
			{
				const auto v_cmd = detail::get_command_line_args();
				std::string cmd;
				for (const auto& item : v_cmd)
					cmd += item + " ";

				//remove last ' '
				if (!cmd.empty())
					cmd.pop_back();

				//mkokot_TODO: probably better store as vector<string>
				//this is because we may have parameters with whitespaces
				//and this way we lost this info
				cur_item.command_line = std::move(cmd);

				OSInfo os_info;
				CPUInfo cpu_info;

				//os_info.
				std::string system_info;
				system_info += "{\n";

				system_info += "\t\"OS\":\n";
				system_info += "\t{\n";
				system_info += "\t\t\"Name\": \"" + os_info.GetName() + "\"\n";
				system_info += "\t\t\"Kernel\": \"" + os_info.GetKernel() + "\"\n";
				system_info += "\t\t\"Version\": \"" + os_info.GetVersion() + "\"\n";
				system_info += "\t\t\"LittleEndian\": " + std::to_string(os_info.IsLittleEndian()) + "\n";

				system_info += "\t}\n"; //OS

				system_info += "\t\"CPU\":\n";
				system_info += "\t{\n";
				system_info += "\t\t\"ModelName\": \"" + cpu_info.GetModelName() + "\"\n";
				system_info += "\t\t\"Vendor\": \"" + cpu_info.GetVendor() + "\"\n";
				system_info += "\t\t\"Architecture\": \"" + cpu_info.GetArchitecture() + "\"\n";
				system_info += "\t\t\"NumSockets\": " + std::to_string(cpu_info.GetNumSockets()) + "\n";
				system_info += "\t\t\"NumPhysicalCores\": " + std::to_string(cpu_info.GetNumPhysicalCores()) + "\n";
				system_info += "\t\t\"NumLogicalCores\": " + std::to_string(cpu_info.GetNumLogicalCores()) + "\n";

				system_info += "\t}\n"; //CPU

				system_info += "\t\"RAM\":\n";
				system_info += "\t{\n";
				system_info += "\t\t\"Total[B]\": " + std::to_string(detail::getRAMPhysicalTotal()) + "\n";

				system_info += "\t}\n"; //RAM

				system_info += "}\n";

				cur_item.system_info = std::move(system_info);
			}

		public:
			//org_path - optional path to a database from which the current one is being created
			HistoryWriter(archive_output_t* archive, const std::string& org_path):
				archive(archive),
				stream_id(archive->register_stream(stream_names::HISTORY))
			{
				if (stream_id == -1)
					throw std::runtime_error("Cannot create stream " + stream_names::HISTORY);

				cur_item.open_time = get_cur_time();

				//mkokot_TODO: by default false, because there are issues if for example multiple threads creates multiple kmcdbs at once and each is capturing
				//maybe we can fix this somehow, but for now lets keep it turned off
				CaptureStdCout(false);
				CaptureStdCerr(false);

				//if there were some history, rewrite it
				if (!org_path.empty())
				{
					MetadataReader org_metadata_reader(org_path, false);
					const auto org_archive = get_archive_from_metadata_reader(org_metadata_reader);

					const auto org_history_stream_id = org_archive->get_stream_id(stream_names::HISTORY);
					if (org_history_stream_id == -1)
						throw std::runtime_error("Cannot find stream " + stream_names::HISTORY + " in " + org_path);

					//mkokot_TODO: I assume a single part is a whole serialized history item... not the nicest solution but lets use it
					std::vector<uint8_t> serialized;
					uint64_t meta;
					while (org_archive->get_part(org_history_stream_id, serialized, meta))
						if (!archive->add_part(stream_id, serialized, meta))
							throw std::runtime_error("Error rewriting history from " + org_path);
				}
			}

			void Finish()
			{
				assert(!is_finished);
				is_finished = true;

				capture_std_cout.reset();
				capture_std_cerr.reset();

				capture_command_line_system_and_hw_info();
				cur_item.mem_peak_bytes = getPeakRSS();
				cur_item.close_time = get_cur_time();

				std::vector<uint8_t> serialized;
				//previous items were serialized/copied from prev archive in ctor
				cur_item.serialize(serialized);

				if (!archive->add_part(stream_id, serialized))
					throw std::runtime_error("Error saving history item");
			}

			bool CaptureStdCout(bool do_capture)
			{
				return capture_ostream_impl(do_capture, capture_std_cout, std::cout,  cur_item.std_cout);
			}
			bool CaptureStdCerr(bool do_capture)
			{
				return capture_ostream_impl(do_capture, capture_std_cerr, std::cerr, cur_item.std_cerr);
			}

			HistoryItem& GetHistoryItem()
			{
				return cur_item;
			}
			~HistoryWriter()
			{
				if (is_finished)
					return;
				try
				{
					Finish();
				}
				catch (...)
				{
				}
			}
		};
	}

	template<typename VALUE_T, typename BIN_T>
	class WriterBase
	{
	private:
		const std::string& path;
		bool is_closed = false;
		
	protected:
		std::vector<std::unique_ptr<BIN_T>> bins{};

		archive_output_t archive;
		detail::Metadata metadata;

	private:
		detail::HistoryWriter history_writer;

		archive_output_t* open_archive()
		{
			archive_output_t::params_t archive_params;
			archive_params.archive_version = 2;
			archive_params.parts_metadata_empty = true;
			archive_params.parts_metadata_fixed_size = true;
			//mkokot_TODO: rozwazyc czy wyeksponowac mozliwosc wyboru:
			//  -  czy buforowane czy niebuforowane
			//  -  czy z reopen czy bez
			if (!archive.open_file_unbuffered(path, false, archive_params))
				throw std::runtime_error("Cannot open file " + path);

			return &archive;
		}
	protected:

		//path_created_from - optional path to a database from which the current one is being created
		//sample_names - not obligatory but possible, if given must be of appropriate size
		template<typename REPRESENTATION_CONFIG_T>
		WriterBase(const Config& config,
			const std::string& path,
			detail::KmersRepresentation kmers_representation,
			const REPRESENTATION_CONFIG_T& representation_config,
			const std::string& path_created_from,
			const std::vector<std::string>& sample_names
			):
			path(path),
			history_writer(open_archive(), path_created_from)
		{
			detail::IterateValues(VALUE_T{}, [&]<typename T>(auto idx, T& /*val*/)
			{
				if constexpr (!std::is_integral_v<T>)
				{
					if (sizeof(T) != config.num_bytes_single_value[idx])
						throw std::runtime_error("For float/double VALUE_T num_bytes_single_value must be sizeof(VALUE_T)");
				}
			});

			metadata.config = config;
			metadata.representation_config = representation_config;

			metadata.value_types = detail::ToValueTypes<VALUE_T>();
			metadata.kmers_representation = kmers_representation;

			//store sample names
			if (!sample_names.empty())
			{
				if (config.num_samples != sample_names.size())
					throw std::runtime_error("Error: wrong number of sample names, should be 0 or " + std::to_string(config.num_samples) + ", but is " + std::to_string(sample_names.size()));

				const auto sample_names_stream_id = archive.register_stream(stream_names::SAMPLE_NAMES);

				if (sample_names_stream_id == -1)
					throw std::runtime_error("Cannot create stream " + stream_names::SAMPLE_NAMES);

				std::vector<uint8_t> serialized;
				for (const auto& sample_name : sample_names)
					refresh::serialization::serialize_string_with_len(sample_name, serialized);

				if (!archive.add_part(sample_names_stream_id, serialized))
					throw std::runtime_error("Error saving stream names");
			}

		}

	public:

		bool CaptureStdCout(bool do_capture)
		{
			return history_writer.CaptureStdCout(do_capture);
		}

		bool CaptureStdCerr(bool do_capture)
		{
			return history_writer.CaptureStdCerr(do_capture);
		}

		void AppendAdditionalInfo(const std::string& info)
		{
			history_writer.GetHistoryItem().info += info;
		}

		BIN_T* GetBin(uint32_t id)
		{
			return bins[id].get();
		}

		//may be called only once, and before close
		//sample_names must be of a valid size
		void SetSampleNames(const std::vector<std::string>& sample_names)
		{
			assert(sample_names.size() == metadata.config.num_samples);
			assert(!is_closed);

		}


		//dtor will also close, but the possible exceptions will be silenced
		//having this method the caller may catch the exception and handle it
		void Close()
		{
			if (is_closed)
				return;
			is_closed = true;

			history_writer.Finish();

			for (auto& bin : bins)
				bin->Close();

			std::vector<uint8_t> metadata_serialized;
			metadata.serialize(metadata_serialized);
			auto res = archive.add_part(archive.register_stream(stream_names::METADATA), metadata_serialized, 0);
			if (!res)
				throw std::runtime_error("Cannot store metadata");

			archive.close();
		}
		~WriterBase() noexcept
		{
			try
			{
				Close();
			}
			catch (...)
			{
			}
		}
	};

	//just k-mer, no lut
	template<typename VALUE_T>
	class WriterSortedPlain : public WriterBase<VALUE_T, BinWriterSortedPlain<VALUE_T>>
	{
		using BIN_T = BinWriterSortedPlain<VALUE_T>;
		using WriterBase<VALUE_T, BIN_T>::archive;
		using WriterBase<VALUE_T, BIN_T>::bins;
	public:
		WriterSortedPlain(
			const Config& config,
			const ConfigSortedPlain& representation_config,
			const std::string& path,
			const std::string& path_created_from = "",
			const std::vector<std::string>& sample_names = {},
			size_t max_part_size = 1ull << 20) :

			WriterBase<VALUE_T, BIN_T>(
				config,
				path,
				detail::KmersRepresentation::SortedPlain,
				representation_config,
				path_created_from,
				sample_names)
		{
			bins.reserve(config.num_bins);
			for (uint64_t i = 0; i < config.num_bins; ++i)
				bins.emplace_back(std::make_unique<BinWriterSortedPlain<VALUE_T>>(
					i,
					&archive,
					config.kmer_len,
					config.num_samples,
					config.num_bytes_single_value,
					max_part_size));

		}
	};

	template<typename VALUE_T>
	class WriterSortedWithLUTRaw : public WriterBase<VALUE_T, BinWriterSortedWithLUTRaw<VALUE_T>>
	{
		using BIN_T = BinWriterSortedWithLUTRaw<VALUE_T>;
		using WriterBase<VALUE_T, BIN_T>::archive;
		using WriterBase<VALUE_T, BIN_T>::bins;
		using WriterBase<VALUE_T, BIN_T>::metadata;
	public:
		WriterSortedWithLUTRaw(
			const Config& config,
			const ConfigSortedWithLUT& representation_config,
			const std::string& path,
			const std::string& path_created_from = "",
			const std::vector<std::string>& sample_names = {}
			) :
			WriterBase<VALUE_T, BIN_T>(
				config,
				path,
				detail::KmersRepresentation::SortedWithLUT,
				representation_config,
				path_created_from,
				sample_names)
		{
			bins.reserve(config.num_bins);
			for (uint64_t i = 0; i < config.num_bins; ++i)
				bins.emplace_back(std::make_unique<BinWriterSortedWithLUTRaw<VALUE_T>>(
					i,
					&archive,
					config.kmer_len,
					config.num_samples,
					representation_config.lut_prefix_len,
					config.num_bytes_single_value));
		}

		void ChangeLutPrefixLen(uint64_t new_lut_prefix_len)
		{
			std::get<ConfigSortedWithLUT>(metadata.representation_config).lut_prefix_len = new_lut_prefix_len;
			for (auto& bin : bins)
				bin->ChangeLutPrefixLen(new_lut_prefix_len);
		}
	};

	template<typename VALUE_T>
	class WriterSortedWithLUT : public WriterBase<VALUE_T, BinWriterSortedWithLUT<VALUE_T>>
	{
		using BIN_T = BinWriterSortedWithLUT<VALUE_T>;
		using WriterBase<VALUE_T, BIN_T>::archive;
		using WriterBase<VALUE_T, BIN_T>::bins;
	public:
		WriterSortedWithLUT(
			const Config& config,
			const ConfigSortedWithLUT& representation_config,
			const std::string& path,
			const std::string& path_created_from,
			const std::vector<std::string>& sample_names = {},
			size_t max_part_size = 1ull << 23) :
			WriterBase<VALUE_T, BIN_T>(
				config,
				path,
				detail::KmersRepresentation::SortedWithLUT,
				representation_config,
				path_created_from,
				sample_names)
		{
			bins.reserve(config.num_bins);
			for (uint64_t i = 0; i < config.num_bins; ++i)
				bins.emplace_back(std::make_unique<BinWriterSortedWithLUT<VALUE_T>>(
					i,
					&archive,
					config.kmer_len,
					config.num_samples,
					representation_config.lut_prefix_len,
					config.num_bytes_single_value, max_part_size));
		}
	};
}
#endif // ! WRITERS_H_
