#pragma once
#include<vector>

class BaseIndividual {
public:
	std::vector<int> sequence;
	double fidelity;
	double leak;
	int numberOfCycles;
	int prefix_count;
	BaseIndividual() {}
	BaseIndividual(
		const std::vector<int>& _sequence,
		double _fidelity,
		double _leak,
		int _numberOfCycles,
		int _prefix_count
	) : sequence(_sequence), fidelity(_fidelity),
		leak(_leak), numberOfCycles(_numberOfCycles), prefix_count(_prefix_count) {}
	bool operator<(const BaseIndividual& lb) const {
		return fidelity > lb.fidelity;
	}
};
