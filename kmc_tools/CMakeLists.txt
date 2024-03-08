# -*- cmake -*-

project(KMC_TOOLS
        LANGUAGES CXX)

add_executable(kmc_tools
  kmer_file_header.cpp kmc_tools.cpp nc_utils.cpp 
  parameters_parser.cpp parser.cpp tokenizer.cpp 
  fastq_filter.cpp fastq_reader.cpp fastq_writer.cpp 
  percent_progress.cpp kff_info_reader.cpp)
target_link_libraries(kmc_tools kmc_api kff_core zlib)

install(TARGETS kmc_tools
        DESTINATION bin
        COMPONENT kmc)
      
