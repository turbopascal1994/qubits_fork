#include "ArgsPreprocessor.h"

std::map<std::string, double> ArgsPreprocessor::run(int argc, char** argv) {
	std::map<std::string, double> mp = {
		{"regular_len", 30},
		{"len", 20},
		{"max_len", 120},
		{"w01", 5},
		{"w12", 0.25},
		{"wt", 25},
		{"angle", 0.024},
		{"module", 0.0001},
		{"mp", 0.5},
		{"cp", 0.9},
		{"iter", 500},
		{"type", 3}
	};
	for (int i = 1; i < argc; i += 2) {
		std::string name = argv[i];
		name = name.substr(2, std::string::npos);
		double value;
		if (name == "type") {
			std::string val = argv[i + 1];
			value = val == "bipolar" ? 3 : 2;
		}
		else {
			value = atof(argv[i + 1]);
		}
		mp[name] = value;
	}
	return mp;
}
