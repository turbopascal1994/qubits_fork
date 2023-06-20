#pragma once
#include<vector>

class BaseIndividual {
public:
	std::vector<int> sequence;
	double fidelity;
	double leak;
	int numberOfCycles;
	BaseIndividual() {}
	BaseIndividual(
		const std::vector<int>& _sequence,
		double _fidelity,
		double _leak,
		int _numberOfCycles
	) : sequence(_sequence), fidelity(_fidelity),
		leak(_leak), numberOfCycles(_numberOfCycles) {}
	bool operator<(const BaseIndividual& lb) const {
		return fidelity > lb.fidelity;
	}
};
