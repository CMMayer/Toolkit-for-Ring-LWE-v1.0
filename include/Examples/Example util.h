#pragma once

#ifndef EU_H
#define EU_H

#include <bitset>
#include "RLWE/RingLweCryptographyField.h"

// Convert char array to bit-vector
coordinate_vec toBitVec(char* text, int n)
{
	coordinate_vec coords = coordinate_vec(n);
	for (int i = 0; i < 9; i++)
	{
		std::bitset<8> bs = std::bitset<8>(text[i]);
		for (int j = 0; j < 8; j++)
		{
			coords[i * 8 + j] = bs[j];
		}
	}
	return coords;
}

// Convert bit-vector to char array
char * toCharArray(coordinate_vec& coords)
{
	int n = coords.size();
	char* text = new char[(n / 8) + 1];
	text[(n / 8)] = '\0';
	for (int i = 0; i < n / 8; i++)
	{
		std::bitset<8> bs = std::bitset<8>();
		for (int j = 0; j < 8; j++)
		{
			bs[j] = coords[i * 8 + j];
		}
		text[i] = bs.to_ulong();
	}
	return text;
}

#endif // !EU_H