#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <bit>

namespace refresh 
{
    namespace hash 
    {
        template <typename T>
        struct MurMur;

        // **********************************************************************************
        template <>
        struct MurMur<uint32_t>
        {
            std::size_t operator()(uint32_t h) const noexcept
            {
                h ^= h >> 16;
                h *= 0x85ebca6b;
                h ^= h >> 13;
                h *= 0xc2b2ae35;
                h ^= h >> 16;

                return static_cast<std::size_t>(h);
            }
        };

        // **********************************************************************************
        template <>
        struct MurMur<uint64_t>
        {
            std::size_t operator()(uint64_t h) const noexcept
            {
                h ^= h >> 33;
                h *= 0xff51afd7ed558ccdULL;
                h ^= h >> 33;
                h *= 0xc4ceb9fe1a85ec53ULL;
                h ^= h >> 33;

                return static_cast<std::size_t>(h);
            }
        };

        // **********************************************************************************
        template <typename T1, typename T2>
        struct MurMur<std::pair<T1, T2>>
        {
            std::size_t operator()(const std::pair<T1, T2>& x) const noexcept
            {
                size_t h1 = MurMur<T1>{}(x.first);
                size_t h2 = MurMur<T2>{}(x.second);

                return (h1 * static_cast<std::size_t>(0x87c37b91114253d5ULL)) ^ h2;
            }
        };

        // **********************************************************************************
        template <>
        struct MurMur<std::string_view>
        {
        private:
            static uint64_t load64(const char*& p)
            {
                uint64_t x = (uint64_t)(*p++);
                x <<= 8; x += (uint64_t)(*p++);
                x <<= 8; x += (uint64_t)(*p++);
                x <<= 8; x += (uint64_t)(*p++);
                x <<= 8; x += (uint64_t)(*p++);
                x <<= 8; x += (uint64_t)(*p++);
                x <<= 8; x += (uint64_t)(*p++);
                x <<= 8; x += (uint64_t)(*p++);
                return x;
            }

        public:
            std::size_t operator()(const std::string_view& s) const
            {
                uint64_t h1 = 0;
                uint64_t h2 = 0;

                const uint64_t c1 = 0x87c37b91114253d5ULL;
                const uint64_t c2 = 0x4cf5ad432745937fULL;

                const char* data = s.data();

                for (std::size_t i = 0; i < s.size() / 16; i++)
                {
                    uint64_t k1 = load64(data);
                    uint64_t k2 = load64(data);

                    k1 *= c1; k1 = std::rotl<uint64_t>(k1, 31); k1 *= c2; h1 ^= k1;
                    h1 = std::rotl<uint64_t>(h1, 27); h1 += h2; h1 = h1 * 5 + 0x52dce729;

                    k2 *= c2; k2 = std::rotl<uint64_t>(k2, 33); k2 *= c1; h2 ^= k2;
                    h2 = std::rotl<uint64_t>(h2, 31); h2 += h1; h2 = h2 * 5 + 0x38495ab5;
                }

                std::size_t tail = s.size() % 16;
                uint64_t k1 = 0;
                uint64_t k2 = 0;

                switch (tail & 15)
                {
                case 15: k2 ^= ((uint64_t)data[14]) << 48; [[fallthrough]];
                case 14: k2 ^= ((uint64_t)data[13]) << 40; [[fallthrough]];
                case 13: k2 ^= ((uint64_t)data[12]) << 32; [[fallthrough]];
                case 12: k2 ^= ((uint64_t)data[11]) << 24; [[fallthrough]];
                case 11: k2 ^= ((uint64_t)data[10]) << 16; [[fallthrough]];
                case 10: k2 ^= ((uint64_t)data[9]) << 8; [[fallthrough]];
                case  9: k2 ^= ((uint64_t)data[8]) << 0;
                    k2 *= c2; k2 = std::rotl<uint64_t>(k2, 33); k2 *= c1; h2 ^= k2;
                    [[fallthrough]];
                case  8: k1 ^= ((uint64_t)data[7]) << 56; [[fallthrough]];
                case  7: k1 ^= ((uint64_t)data[6]) << 48; [[fallthrough]];
                case  6: k1 ^= ((uint64_t)data[5]) << 40; [[fallthrough]];
                case  5: k1 ^= ((uint64_t)data[4]) << 32; [[fallthrough]];
                case  4: k1 ^= ((uint64_t)data[3]) << 24; [[fallthrough]];
                case  3: k1 ^= ((uint64_t)data[2]) << 16; [[fallthrough]];
                case  2: k1 ^= ((uint64_t)data[1]) << 8; [[fallthrough]];
                case  1: k1 ^= ((uint64_t)data[0]) << 0;
                    k1 *= c1; k1 = std::rotl<uint64_t>(k1, 31); k1 *= c2; h1 ^= k1;
                }

                h1 ^= (uint64_t)s.size(); 
                h2 ^= (uint64_t)s.size();
                h1 += h2; 
                h2 += h1;
                h1 = (uint64_t) MurMur<uint64_t>{}(h1);
                h2 = (uint64_t) MurMur<uint64_t>{}(h2);
                h1 += h2; 
                h2 += h1;

                return std::size_t(h1 ^ h2);
            }
        };

        template <>
        struct MurMur<std::string>
        {
            std::size_t operator()(const std::string& s) const
            {
                return MurMur<std::string_view>{}(s);
            }
        };

        // **********************************************************************************
        //
        // **********************************************************************************

        template <typename T>
        struct MurMurSeeded;

        // **********************************************************************************
        template <>
        struct MurMurSeeded<uint32_t>
        {
            std::size_t operator()(uint32_t h) const noexcept
            {
                h ^= 0x85ebca6bUL; 
                h ^= h >> 16;
                h *= 0x85ebca6b;
                h ^= h >> 13;
                h *= 0xc2b2ae35;
                h ^= h >> 16;

                return std::size_t(h);
            }
        };

        // **********************************************************************************
        template <>
        struct MurMurSeeded<uint64_t>
        {
            std::size_t operator()(uint64_t h) const noexcept
            {
                h ^= 0x87c37b91114253d5ULL;
                h ^= h >> 33;
                h *= 0xff51afd7ed558ccdULL;
                h ^= h >> 33;
                h *= 0xc4ceb9fe1a85ec53ULL;
                h ^= h >> 33;

                return std::size_t(h);
            }
        };

        // **********************************************************************************
        template <typename T1, typename T2>
        struct MurMurSeeded<std::pair<T1, T2>>
        {
            std::size_t operator()(const std::pair<T1, T2>&x) const noexcept
            {
                size_t h1 = MurMurSeeded<T1>{}(x.first);
                size_t h2 = MurMurSeeded<T2>{}(x.second);

                return (h1 * std::size_t(0x87c37b91114253d5ULL)) ^ h2;
            }
        };

        // **********************************************************************************
        template <>
        struct MurMurSeeded<std::string>
        {
            std::size_t operator()(const std::string& s) const
            {
                return MurMur<std::string>{}(s);
            }
        };

        // **********************************************************************************
        template <>
        struct MurMurSeeded<std::string_view>
        {
            std::size_t operator()(const std::string_view& s) const
            {
                return MurMur<std::string_view>{}(s);
            }
        };

        // **********************************************************************************
        // 
        // **********************************************************************************
        template <typename T>
        struct MurMurSimplified;

        // **********************************************************************************
        template <>
        struct MurMurSimplified<uint32_t>
        {
            std::size_t operator()(uint32_t h) const noexcept
            {
                h ^= h >> 16;
                h *= 0x85ebca6b;
                h ^= h >> 13;

                return std::size_t(h);
            }
        };

        // **********************************************************************************
        template <>
        struct MurMurSimplified<uint64_t>
        {
            std::size_t operator()(uint64_t h) const noexcept
            {
                h ^= h >> 33;
                h *= 0xff51afd7ed558ccdULL;
                h ^= h >> 33;

                return std::size_t(h);
            }
        };

    } // namespace hash
} // namespace refresh
