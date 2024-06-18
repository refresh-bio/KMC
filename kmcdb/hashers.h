#ifndef HASHERS_H
#define HASHERS_H

namespace kmcdb
{
	struct MurMur64Hash
	{
		size_t operator()(size_t h) const noexcept
		{
			h ^= h >> 33;
			h *= 0xff51afd7ed558ccdL;
			h ^= h >> 33;
			h *= 0xc4ceb9fe1a85ec53L;
			h ^= h >> 33;

			return h;
		}
	};
}

#endif // ! HASHERS_H
