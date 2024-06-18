//Original preamble
/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

 /*
  * Author of changes: Marek Kokot
  * Info: I first find this solution on SO: https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c - but is seems the links there are not working
  * To the best of my understanding I may use this file
  * I have added `inline` to these functions to may it use directly from header file
  * I have also added function getRAMPhysicalTotal based on https://github.com/DigitalInBlue/Celero/blob/master/src/Memory.cpp
  * A lot of code here is taken directly from https://github.com/lfreist/hwinfo
  *	For this project I encapsulated everything in appropriate namespace
  */

#ifndef CAPTURE_SYS_HW_CMD_H
#define CAPTURE_SYS_HW_CMD_H


#include <sstream>
#include <fstream>
#include <cstring>
#include <vector>
#include <set>

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#include <WbemIdl.h>
#include <comdef.h>
#pragma comment(lib, "wbemuuid.lib")

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>
#include <sys/utsname.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#include <sys/sysctl.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#include <sys/sysinfo.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif


#ifdef __APPLE__
#include <crt_externs.h>
#endif 

namespace kmcdb::detail
{
    inline std::vector<std::string> split(const std::string& input, const char delimiter)
    {
        std::vector<std::string> result;
        size_t shift = 0;
        while (true)
        {
            size_t match = input.find(delimiter, shift);
            result.emplace_back(input.substr(shift, match - shift));
            if (match == std::string::npos)
                break;

            shift = match + 1;
        }
        return result;
    }

    inline std::vector<std::string> split(const std::string& input, const std::string& delimiter)
    {
        std::vector<std::string> result;
        size_t shift = 0;
        while (true)
        {
            size_t match = input.find(delimiter, shift);
            result.emplace_back(input.substr(shift, match - shift));
            if (match == std::string::npos)
                break;

            shift = match + delimiter.size();
        }
        return result;
    }

    inline void strip(std::string& input)
    {
        if (input.empty())
            return;

        // optimization for input size == 1
        if (input.size() == 1)
        {
            if (input[0] == ' ' || input[0] == '\t' || input[0] == '\n')\
            {
                input = "";
                return;
            }
            else {
                return;
            }
        }
        size_t start_index = 0;
        while (true)
        {
            char c = input[start_index];
            if (c != ' ' && c != '\t' && c != '\n')
            {
                break;
            }
            start_index++;
        }
        size_t end_index = input.size() - 1;
        while (true)
        {
            char c = input[end_index];
            if (c != ' ' && c != '\t' && c != '\n')
                break;

            end_index--;
        }

        if (end_index < start_index)
        {
            input.assign("");
            return;
        }
        input.assign(input.begin() + start_index, input.begin() + end_index + 1);
    }

    inline std::string wstring_to_std_string(const std::wstring& ws)
    {
        std::string str_locale = setlocale(LC_ALL, "");
        const wchar_t* wch_src = ws.c_str();

#ifdef _MSC_VER
        size_t n_dest_size;
        wcstombs_s(&n_dest_size, nullptr, 0, wch_src, 0);
        n_dest_size++;  // Increase by one for null terminator

        char* ch_dest = new char[n_dest_size];
        memset(ch_dest, 0, n_dest_size);

        size_t n_convert_size;
        wcstombs_s(&n_convert_size, ch_dest, n_dest_size, wch_src,
            n_dest_size - 1);  // subtract one to ignore null terminator

        std::string result_text = ch_dest;
        delete[] ch_dest;
#else
        size_t n_dest_size = wcstombs(NULL, wch_src, 0) + 1;
        char* ch_dest = new char[n_dest_size];
        memset(ch_dest, 0, n_dest_size);
        wcstombs(ch_dest, wch_src, n_dest_size);
        std::string result_text = ch_dest;
        delete[] ch_dest;
#endif

        setlocale(LC_ALL, str_locale.c_str());
        return result_text;
    }


    /**
      * Returns the peak (maximum so far) resident set size (physical
      * memory use) measured in bytes, or zero if the value cannot be
      * determined on this OS.
      */
    inline size_t getPeakRSS()
    {
#if defined(_WIN32)
        /* Windows -------------------------------------------------- */
        PROCESS_MEMORY_COUNTERS info;
        GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
        return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
        /* AIX and Solaris ------------------------------------------ */
        struct psinfo psinfo;
        int fd = -1;
        if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
            return (size_t)0L;      /* Can't open? */
        if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
        {
            close(fd);
            return (size_t)0L;      /* Can't read? */
        }
        close(fd);
        return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
        /* BSD, Linux, and OSX -------------------------------------- */
        struct rusage rusage;
        getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
        return (size_t)rusage.ru_maxrss;
#else
        return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
        /* Unknown OS ----------------------------------------------- */
        return (size_t)0L;          /* Unsupported. */
#endif
    }

    /**
     * Returns the current resident set size (physical memory use) measured
     * in bytes, or zero if the value cannot be determined on this OS.
     */
    inline size_t getCurrentRSS()
    {
#if defined(_WIN32)
        /* Windows -------------------------------------------------- */
        PROCESS_MEMORY_COUNTERS info;
        GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
        return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
        /* OSX ------------------------------------------------------ */
        struct mach_task_basic_info info;
        mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
        if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
            (task_info_t)&info, &infoCount) != KERN_SUCCESS)
            return (size_t)0L;      /* Can't access? */
        return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
        /* Linux ---------------------------------------------------- */
        long rss = 0L;
        FILE* fp = NULL;
        if ((fp = fopen("/proc/self/statm", "r")) == NULL)
            return (size_t)0L;      /* Can't open? */
        if (fscanf(fp, "%*s%ld", &rss) != 1)
        {
            fclose(fp);
            return (size_t)0L;      /* Can't read? */
        }
        fclose(fp);
        return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
        /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
        return (size_t)0L;          /* Unsupported. */
#endif
    }

    //based on https://github.com/DigitalInBlue/Celero/blob/master/src/Memory.cpp
    inline int64_t getRAMPhysicalTotal()
    {
#ifdef _WIN32
        MEMORYSTATUSEX memInfo;
        memInfo.dwLength = sizeof(MEMORYSTATUSEX);
        GlobalMemoryStatusEx(&memInfo);
        return static_cast<int64_t>(memInfo.ullTotalPhys);
#elif defined(__APPLE__)
        int64_t physicalMemory;
        size_t length = sizeof(physicalMemory);
        if (sysctlbyname("hw.memsize", &physicalMemory, &length, nullptr, 0) == 0)
        {
            return static_cast<size_t>(physicalMemory);
        }
        return 0;
#else
        struct sysinfo memInfo;
        sysinfo(&memInfo);
        return memInfo.totalram * memInfo.mem_unit;
#endif
    }


    class OSInfo
    {

        std::string name = "<unknown name>";
        std::string version = "<unknown version>";
        std::string kernel = "<unknown kernel>";
        bool littleEndian = false;

    public:
        OSInfo()
        {
            // Get endian. This is platform independent...
            char16_t dummy = 0x0102;
            littleEndian = ((char*)&dummy)[0] == 0x02;
#ifdef _WIN32

            struct locator_and_service_t
            {
                IWbemLocator* locator = nullptr;
                IWbemServices* service = nullptr;
                ~locator_and_service_t()
                {
                    if (locator) locator->Release();
                    if (service) service->Release();
                    CoUninitialize();
                }

            };
            locator_and_service_t locator_and_service;

            auto res = CoInitializeSecurity(nullptr, -1, nullptr, nullptr, RPC_C_AUTHN_LEVEL_DEFAULT, RPC_C_IMP_LEVEL_IMPERSONATE,
                nullptr, EOAC_NONE, nullptr);
            res &= CoInitializeEx(nullptr, COINIT_MULTITHREADED);
            res &= CoCreateInstance(CLSID_WbemLocator, nullptr, CLSCTX_INPROC_SERVER, IID_IWbemLocator, (LPVOID*)&locator_and_service.locator);
            if (locator_and_service.locator)
            {
                res &= locator_and_service.locator->ConnectServer(
                    _bstr_t("ROOT\\CIMV2"), nullptr, nullptr, nullptr, 0,
                    nullptr, nullptr, &locator_and_service.service);

                if (locator_and_service.service)
                    res &= CoSetProxyBlanket(locator_and_service.service, RPC_C_AUTHN_WINNT, RPC_C_AUTHZ_NONE, nullptr, RPC_C_AUTHN_LEVEL_CALL,
                        RPC_C_IMP_LEVEL_IMPERSONATE, nullptr, EOAC_NONE);
            }
            if (!SUCCEEDED(res)) {
                return;
            }

            const std::wstring query_string(L"SELECT Caption, BuildNumber, Version FROM Win32_OperatingSystem");

            if (locator_and_service.service == nullptr)
                return;

            IEnumWbemClassObject* enumerator = nullptr;

            res = locator_and_service.service->ExecQuery(bstr_t(L"WQL"), bstr_t(std::wstring(query_string.begin(), query_string.end()).c_str()),
                WBEM_FLAG_FORWARD_ONLY | WBEM_FLAG_RETURN_IMMEDIATELY, nullptr, &enumerator);

            if (!SUCCEEDED(res))
                return;

            ULONG u_return = 0;
            IWbemClassObject* obj = nullptr;
            enumerator->Next(WBEM_INFINITE, 1, &obj, &u_return);
            if (!u_return) {
                return;
            }
            VARIANT vt_prop;
            HRESULT hr;
            hr = obj->Get(L"Caption", 0, &vt_prop, nullptr, nullptr);
            if (SUCCEEDED(hr) && (V_VT(&vt_prop) == VT_BSTR)) {
                name = detail::wstring_to_std_string(vt_prop.bstrVal);
            }

            hr = obj->Get(L"BuildNumber", 0, &vt_prop, nullptr, nullptr);
            if (SUCCEEDED(hr) && (V_VT(&vt_prop) == VT_BSTR)) {
                version = detail::wstring_to_std_string(vt_prop.bstrVal);
            }
            hr = obj->Get(L"Version", 0, &vt_prop, nullptr, nullptr);
            if (SUCCEEDED(hr) && (V_VT(&vt_prop) == VT_BSTR)) {
                kernel = detail::wstring_to_std_string(vt_prop.bstrVal);
            }
            VariantClear(&vt_prop);
            obj->Release();

#elif defined(__APPLE__)
            size_t size = 1024;

            name.resize(size);
            if (sysctlbyname("kern.ostype", (void*)(name.data()), &size, nullptr, 0) == 0)
            {
                name.resize(size);  // trim the string to the actual size
            }

            size = 1024;

            version.resize(size);
            if (sysctlbyname("kern.osproductversion", (void*)(version.data()), &size, nullptr, 0) == 0)
            {
                version.resize(size);  // trim the string to the actual size
            }
#else
            std::string line;
            std::ifstream stream("/etc/os-release");
            if (!stream) {
                name = "Linux";
            }
            while (std::getline(stream, line)) {
                if (line.starts_with("PRETTY_NAME")) {
                    line = line.substr(line.find('=') + 1, line.length());
                    // remove \" at begin and end of the substring result
                    name = { line.begin() + 1, line.end() - 1 };
                }
                if (line.starts_with("VERSION=")) {
                    line = line.substr(line.find('=') + 1, line.length());
                    // remove \" at begin and end of the substring result
                    version = { line.begin() + 1, line.end() - 1 };
                }
            }
            stream.close();

            static utsname info;
            if (uname(&info) == 0) {
                kernel = info.release;
            }
#endif
        }
        const std::string& GetName() const
        {
            return name;
        }
        const std::string& GetVersion() const
        {
            return version;
        }
        const std::string& GetKernel() const
        {
            return kernel;
        }
        bool IsLittleEndian() const
        {
            return littleEndian;
        }
    };


    class CPUInfo
    {
        std::string architecture = "<unknown>";
        std::string modelName = "<unknown>";
        std::string vendor = "<unknown>";
        int numSockets{ -1 };
        int numPhysicalCores{ -1 };
        int numLogicalCores{ -1 };
    public:
        //mkokot_TODO: move to ctor, just use a single cpu, the list does not make sense i think... also it seems the original hwinfo works strange on multi cpu platform
        CPUInfo()
        {
#ifdef _WIN32
            struct locator_and_service_t
            {
                IWbemLocator* locator = nullptr;
                IWbemServices* service = nullptr;
                ~locator_and_service_t()
                {
                    if (locator) locator->Release();
                    if (service) service->Release();
                    CoUninitialize();
                }

            };
            locator_and_service_t locator_and_service;

            auto res = CoInitializeSecurity(nullptr, -1, nullptr, nullptr, RPC_C_AUTHN_LEVEL_DEFAULT, RPC_C_IMP_LEVEL_IMPERSONATE,
                nullptr, EOAC_NONE, nullptr);
            res &= CoInitializeEx(nullptr, COINIT_MULTITHREADED);
            res &= CoCreateInstance(CLSID_WbemLocator, nullptr, CLSCTX_INPROC_SERVER, IID_IWbemLocator, (LPVOID*)&locator_and_service.locator);
            if (locator_and_service.locator)
            {
                res &= locator_and_service.locator->ConnectServer(
                    _bstr_t("ROOT\\CIMV2"), nullptr, nullptr, nullptr, 0,
                    nullptr, nullptr, &locator_and_service.service);

                if (locator_and_service.service)
                    res &= CoSetProxyBlanket(locator_and_service.service, RPC_C_AUTHN_WINNT, RPC_C_AUTHZ_NONE, nullptr, RPC_C_AUTHN_LEVEL_CALL,
                        RPC_C_IMP_LEVEL_IMPERSONATE, nullptr, EOAC_NONE);
            }
            if (!SUCCEEDED(res)) {
                return;
            }

            const std::wstring query_string(
                L"SELECT Name, Architecture, Manufacturer, NumberOfCores, NumberOfLogicalProcessors "
                L"FROM Win32_Processor");
            if (locator_and_service.service == nullptr)
                return;

            IEnumWbemClassObject* enumerator = nullptr;

            res = locator_and_service.service->ExecQuery(bstr_t(L"WQL"), bstr_t(std::wstring(query_string.begin(), query_string.end()).c_str()),
                WBEM_FLAG_FORWARD_ONLY | WBEM_FLAG_RETURN_IMMEDIATELY, nullptr, &enumerator);

            if (!SUCCEEDED(res))
                return;

            ULONG u_return = 0;
            IWbemClassObject* obj = nullptr;

            if (!enumerator)
                return;

            //process first
            int iterations = 0;
            while (true)
            {
                enumerator->Next(WBEM_INFINITE, 1, &obj, &u_return);
                if (!u_return)
                    break;

                ++iterations;

                VARIANT vt_prop;
                HRESULT hr;
                hr = obj->Get(L"Name", 0, &vt_prop, nullptr, nullptr);
                if (SUCCEEDED(hr) && (V_VT(&vt_prop) == VT_BSTR)) {
                    if (iterations == 1) {
                        modelName = detail::wstring_to_std_string(vt_prop.bstrVal);
                        detail::strip(modelName);
                    }
                    else if (modelName != detail::wstring_to_std_string(vt_prop.bstrVal))
                        modelName = "<different model names>";
                }
                hr = obj->Get(L"Manufacturer", 0, &vt_prop, nullptr, nullptr);
                if (SUCCEEDED(hr) && (V_VT(&vt_prop) == VT_BSTR)) {
                    if (iterations == 1)
                        vendor = detail::wstring_to_std_string(vt_prop.bstrVal);
                    else if (vendor != detail::wstring_to_std_string(vt_prop.bstrVal))
                        modelName = "<different vendors>";
                }
                hr = obj->Get(L"NumberOfCores", 0, &vt_prop, nullptr, nullptr);
                if (SUCCEEDED(hr) && (V_VT(&vt_prop) == VT_I4)) {
                    if (iterations == 1)
                        numPhysicalCores = vt_prop.intVal;
                    else if (numPhysicalCores != vt_prop.intVal)
                        numPhysicalCores = -1;
                }
                hr = obj->Get(L"Architecture", 0, &vt_prop, nullptr, nullptr);
                if (SUCCEEDED(hr) && (V_VT(&vt_prop) == VT_I4)) {
                    //https://learn.microsoft.com/en-us/windows/win32/cimwin32prov/win32-processor
                    std::string arch;
                    switch (vt_prop.intVal)
                    {
                    case 0: arch = "x86"; break;
                    case 1: arch = "MIPS"; break;
                    case 2: arch = "Alpha"; break;
                    case 3: arch = "PowerPC"; break;
                    case 5: arch = "ARM"; break;
                    case 6: arch = "ia64"; break;
                    case 9: arch = "x64"; break;
                    case 12: arch = "ARM64 "; break;
                    }
                    if (iterations == 1)
                        architecture = arch;
                    else if (architecture != arch)
                        architecture = "<different architectures>";

                }

                hr = obj->Get(L"NumberOfLogicalProcessors", 0, &vt_prop, nullptr, nullptr);
                if (SUCCEEDED(hr) && (V_VT(&vt_prop) == VT_I4)) {
                    if (iterations == 1)
                        numLogicalCores = vt_prop.intVal;
                    else if (numLogicalCores != vt_prop.intVal)
                        numLogicalCores = -1;
                }

                VariantClear(&vt_prop);
                obj->Release();
            }
            numSockets = iterations;

#elif defined(__APPLE__)
            size_t size = 1024;
            std::string tmp;
            int num;


            //model
            size = 1024;
            tmp.resize(size);
            if (sysctlbyname("machdep.cpu.brand_string", tmp.data(), &size, NULL, 0) == 0)
            {
                tmp.resize(size);
                modelName = tmp;
            }

            //vendor
            size = 1024;
            tmp.resize(size);
            if (sysctlbyname("machdep.cpu.vendor", tmp.data(), &size, NULL, 0) == 0)
            {
                tmp.resize(size);
                vendor = tmp;
            }

            //architecture
            size = 1024;
            tmp.resize(size);
            if (sysctlbyname("hw.machine", tmp.data(), &size, nullptr, 0) == 0)
            {
                tmp.resize(size);
                architecture = tmp;
            }

            //numPhysicalCores
            size = sizeof(num);
            if (sysctlbyname("hw.physicalcpu", &num, &size, nullptr, 0) == 0)
                numPhysicalCores = num;

            //numLogicalCores
            size = sizeof(num);
            if (sysctlbyname("hw.logicalcpu", &num, &size, nullptr, 0) == 0)
                numLogicalCores = num;

            //numSockets
            size = sizeof(num);
            if (sysctlbyname("hw.packages", &num, &size, nullptr, 0) == 0)
                numSockets = num;

#else
            //try using lscpu
            auto exec_cmd = [](const char* cmd) -> std::string
                {
                    constexpr int buff_size = 128;
                    char buff[buff_size];
                    std::string result;
                    FILE* pipe = popen(cmd, "r");
                    if (!pipe)
                        return "";

                    while (fgets(buff, buff_size, pipe) != nullptr)
                        result += buff;

                    pclose(pipe);
                    return result;
                };

            auto lscpu = exec_cmd("lscpu");

            if (!lscpu.empty())
            {
                std::istringstream iss(lscpu);
                std::string line;

                int threads_per_socket{};
                while (std::getline(iss, line))
                {

                    auto parts = detail::split(line, ':');
                    if (parts.size() < 2)
                        continue;
                    auto key = std::move(parts[0]);
                    auto value = std::move(parts[1]);
                    detail::strip(key);
                    detail::strip(value);

                    if (key == "Architecture")
                        architecture = value;
                    else if (key == "Vendor ID")
                        vendor = value;
                    else if (key == "Model name")
                        modelName = value;
                    else if (key == "Socket(s)")
                        numSockets = std::stoi(value);
                    else if (key == "Core(s) per socket")
                        numPhysicalCores = std::stoi(value);
                    else if (key == "Thread(s) per core")
                        threads_per_socket = std::stoi(value);
                }
                numLogicalCores = threads_per_socket * numPhysicalCores;
                return;
            }
            else
            {
                //try using the same as in hwinfo
                std::ifstream cpuinfo("/proc/cpuinfo");
                if (!cpuinfo.is_open())
                    return;

                std::string file((std::istreambuf_iterator<char>(cpuinfo)), (std::istreambuf_iterator<char>()));
                cpuinfo.close();
                auto cpu_blocks_string = detail::split(file, "\n\n");

                if (cpu_blocks_string.empty())
                    return;

                const auto& first_block = cpu_blocks_string[0];
                std::set<int> physical_ids;
                std::istringstream iss(first_block);
                std::string line;
                while (std::getline(iss, line))
                {
                    auto parts = detail::split(line, ':');
                    if (parts.size() < 2)
                        continue;

                    auto key = std::move(parts[0]);
                    auto value = std::move(parts[1]);
                    detail::strip(key);
                    detail::strip(value);

                    if (key == "vendor_id")
                        vendor = value;
                    else if (key == "model name")
                        modelName = value;
                    else if (key == "siblings")
                        numLogicalCores = std::stoi(value);
                    else if (key == "cpu cores")
                        numPhysicalCores = std::stoi(value);
                    else if (key == "physical id")
                        physical_ids.insert(std::stoi(value));
                }

                //walk through remaining to get the number of sockets which is I guess the number of unique physical ids...
                for (int i = 1; i < std::ssize(cpu_blocks_string); ++i)
                {
                    const auto& cur_block = cpu_blocks_string[i];
                    std::istringstream iss(cur_block);
                    std::string line;
                    while (std::getline(iss, line))
                    {
                        auto parts = detail::split(line, ':');
                        if (parts.size() < 2)
                            continue;

                        auto key = std::move(parts[0]);
                        auto value = std::move(parts[1]);
                        detail::strip(key);
                        detail::strip(value);


                        if (key == "vendor_id")
                        {
                            if (vendor != value) //I think this should not happen
                                vendor = "<different vendors>";
                        }
                        else if (key == "model name")
                        {
                            if (modelName != value)
                                modelName = "<different model names>";
                        }
                        else if (key == "siblings")
                        {
                            if (numLogicalCores != std::stoi(value))
                                numLogicalCores = -1;
                        }
                        else if (key == "cpu cores")
                        {
                            if (numPhysicalCores != std::stoi(value))
                                numPhysicalCores = -1;
						}
                        else if (key == "physical id")
                        {
                            physical_ids.insert(std::stoi(value));
                        }
                    }
                }
                numSockets = physical_ids.size();
            }
#endif
        }

        const std::string& GetModelName() const
        {
            return modelName;
        }

        const std::string& GetArchitecture() const
        {
            return architecture;
        }

        const std::string& GetVendor() const
        {
            return vendor;
        }

        int GetNumSockets() const
        {
            return numSockets;
        }

        int GetNumPhysicalCores() const
        {
            return numPhysicalCores;
        }

        int GetNumLogicalCores() const
        {
            return numLogicalCores;
        }
    };

#ifdef _WIN32
    inline std::vector<std::string> get_command_line_args()
    {
        std::vector<std::string> args;
        LPWSTR cmdline = GetCommandLineW();

        int argc;
        LPWSTR* argv = CommandLineToArgvW(cmdline, &argc);

        if (argv)
        {
            for (int i = 0; i < argc; ++i)
                args.push_back(detail::wstring_to_std_string(argv[i]));
            LocalFree(argv);
        }

        return args;
    }
#elif defined(__APPLE__)

	inline std::vector<std::string> get_command_line_args() {
	    std::vector<std::string> args;
	    char** argv = *_NSGetArgv();
	    int argc = *_NSGetArgc();

	    for (int i = 0; i < argc; ++i) {
	        args.push_back(argv[i]);
	    }

	    return args;
}
#else

	inline std::vector<std::string> get_command_line_args()
	{
	    std::ifstream cmdline_file("/proc/self/cmdline");
	    //std::string cmdline;

	    std::vector<std::string> args;
	    std::string arg;
	    while (std::getline(cmdline_file, arg, '\0'))
	        args.push_back(arg);

	    return args;
	}
#endif
}

#endif // ! CAPTURE_SYS_HW_CMD_H
