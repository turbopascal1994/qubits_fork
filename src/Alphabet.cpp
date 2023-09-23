#include "Alphabet.h"
#include <cassert>
#include <algorithm>

void generateBipolarAlphabet(unsigned int alphabetSize, std::vector<std::array<unsigned int, 3>>& wordbook, std::array<int, 3>& order) {
	order = { 1, 0 };
	assert(0);
}

void generateUnipolarAlphabet(unsigned int alphabetSize, std::vector<std::array<unsigned int, 3>>& wordbook, std::array<int, 3>& order) {
	order = { 1, 0, -1 };
	for (unsigned int i = 0; i < 100; ++i) {
		for (unsigned int j = 0; j < 100; ++j) {
			for (unsigned int k = 0; k < 4; ++k) {
				if (i + j + k >= 3 && i + j + k <= 7 && i + j > 0) {
					wordbook.push_back({ i, j, k });
				}
			}
		}
	}
	std::sort(wordbook.begin(), wordbook.end(), [](std::array<unsigned int, 3> &lb, std::array<unsigned int, 3> & rb) {
		int lb_sum = lb[0] + lb[1] + lb[2];
		int rb_sum = rb[0] + rb[1] + rb[2];
		if (lb_sum == rb_sum) {
			if (lb[0] == rb[0]) {
				if (lb[1] == rb[1]) {
					return lb[2] > rb[2];
				}
				return lb[1] > rb[1];
			}
			return lb[0] > rb[0];
		}
		return lb_sum < rb_sum;
	});
	wordbook.resize(alphabetSize);
}

Alphabet::Alphabet(unsigned int alphabetSize, int type) {
	// type == 2 -> bipolar
	// type == 3 -> unipolar
	assert(type == 2 || type == 3);
	if (type == 2) {
		generateBipolarAlphabet(alphabetSize, wordbook, order);
	}
	else if (type == 3) {
		generateUnipolarAlphabet(alphabetSize, wordbook, order);
	}
}

std::vector<int> Alphabet::decode(const std::vector<int>& sequence) {
	std::vector<int> decoded;
	for (auto& i : sequence) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < wordbook[i][j]; ++k) {
				decoded.push_back(order[j]);
			}
		}
	}
	return decoded;
}
