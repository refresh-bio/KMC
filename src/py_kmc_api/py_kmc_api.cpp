#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <kmc_file.h>

namespace py = pybind11;

struct Count
{
	uint64 value;	
};

struct CountVec
{
	std::vector<uint32> value;
};

struct LongKmerRepresentation
{
	std::vector<uint64> value;
};

PYBIND11_MODULE(py_kmc_api, m) {
	m.doc() = "Python wrapper for KMC_API."; // optional module docstring
	
	py::class_<LongKmerRepresentation>(m, "LongKmerRepresentation")
		.def(py::init<>())
		.def_readwrite("value", &LongKmerRepresentation::value);

	py::class_<CountVec>(m, "CountVec")
		.def(py::init<>())
		.def_readwrite("value", &CountVec::value);

	py::class_<Count>(m, "Count")
		.def(py::init<>())		
		.def_readwrite("value", &Count::value);

	py::class_<CKMCFileInfo>(m, "KMCFileInfo")
		.def(py::init<>())
		.def_readwrite("kmer_length", &CKMCFileInfo::kmer_length)
		.def_readwrite("mode", &CKMCFileInfo::mode)
		.def_readwrite("counter_size", &CKMCFileInfo::counter_size)
		.def_readwrite("counter_size", &CKMCFileInfo::counter_size)
		.def_readwrite("lut_prefix_length", &CKMCFileInfo::lut_prefix_length)
		.def_readwrite("signature_len", &CKMCFileInfo::signature_len)
		.def_readwrite("min_count", &CKMCFileInfo::min_count)
		.def_readwrite("max_count", &CKMCFileInfo::max_count)
		.def_readwrite("both_strands", &CKMCFileInfo::both_strands)
		.def_readwrite("total_kmers", &CKMCFileInfo::total_kmers);


	py::class_<CKmerAPI>(m, "KmerAPI")
		.def(py::init<uint32>(), py::arg("length") = 1)
		.def(py::init<const CKmerAPI&>())
		.def("assign", &CKmerAPI::operator=)
		.def(py::self == py::self)
		.def(py::self < py::self)
		.def("get_asci_symbol", &CKmerAPI::get_asci_symbol)
		.def("get_num_symbol", &CKmerAPI::get_num_symbol)
		.def("to_string", [](CKmerAPI& ptr) {return ptr.to_string(); })
		.def("to_string", [](CKmerAPI& ptr, char* str) { ptr.to_string(str); })
		.def("__str__", [](CKmerAPI& ptr) {return ptr.to_string(); })
		.def("to_string", [](CKmerAPI& ptr, std::string& str) { ptr.to_string(str); })
		.def("to_long", [](CKmerAPI& ptr, LongKmerRepresentation& res) {ptr.to_long(res.value); })
		.def("reverse", &CKmerAPI::reverse)
		.def("get_signature", &CKmerAPI::get_signature)
		.def("from_string", [](CKmerAPI& ptr, const char* str) { return ptr.from_string(str); })
		.def("from_string", [](CKmerAPI& ptr, const std::string& str) { return ptr.from_string(str); });
		


	py::class_<CKMCFile>(m, "KMCFile")
		.def(py::init<>())
		.def("OpenForRA", &CKMCFile::OpenForRA)
		.def("OpenForListing", &CKMCFile::OpenForListing)
		.def("ReadNextKmer", [](CKMCFile& ptr, CKmerAPI& kmer, Count& count) {return ptr.ReadNextKmer(kmer, count.value); })
		.def("Close", &CKMCFile::Close)
		.def("SetMinCount", &CKMCFile::SetMinCount)
		.def("GetMinCount", &CKMCFile::GetMinCount)
		.def("SetMaxCount", &CKMCFile::SetMaxCount)
		.def("GetMaxCount", &CKMCFile::GetMaxCount)
		.def("GetBothStrands", &CKMCFile::GetBothStrands)
		.def("KmerCount", &CKMCFile::KmerCount)
		.def("KmerLength", &CKMCFile::KmerLength)
		.def("RestartListing", &CKMCFile::RestartListing)
		.def("Eof", &CKMCFile::Eof)
		.def("CheckKmer", [](CKMCFile& ptr, CKmerAPI& kmer, Count& count) { return ptr.CheckKmer(kmer, count.value); })
		.def("IsKmer", &CKMCFile::IsKmer)
		.def("ResetMinMaxCounts", &CKMCFile::ResetMinMaxCounts)
		.def("Info", [](CKMCFile& ptr, CKMCFileInfo& info) {return ptr.Info(info); })
		.def("Info", [](CKMCFile& ptr) { CKMCFileInfo info; ptr.Info(info); return info; })
		.def("GetCountersForRead", [](CKMCFile& ptr, const std::string& read, CountVec& counters) {
			return ptr.GetCountersForRead(read, counters.value);
		})
		;
		
}
