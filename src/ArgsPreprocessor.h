#pragma once
#include <map>
#include <string>

class ArgsPreprocessor {
public:
	static std::map<std::string, double> run(int argc, char** argv);
};

