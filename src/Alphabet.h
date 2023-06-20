#pragma once
#include <array>
#include <vector>

struct Alphabet {
	std::vector<std::array<unsigned int, 3>> wordbook;
	std::array<int, 3> order;
	Alphabet(const std::vector<std::array<unsigned int, 3>>& _wordbook, const std::array<int, 3>& _order) : wordbook(_wordbook), order(_order) {}
	Alphabet(unsigned int alphabetSize, int type);
	std::vector<int> decode(const std::vector<int>& sequence);
};
