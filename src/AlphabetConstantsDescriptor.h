#pragma once
#include "Alphabet.h"
#include "ConstantsDescriptor.h"

struct AlphabetConstantsDescriptor : public ConstantsDescriptor {
public:
	Alphabet alp;
	AlphabetConstantsDescriptor(
		const Alphabet& _alp,
		double _w01, double _w12, double _wt, double _w,
		double _T, double _Tstep, double _Theta,
		double _neededAngle, int _numberOfCycles,
		double _angleUpperBound, int _type
	) : ConstantsDescriptor(_w01, _w12, _wt, _w, _T, _Tstep, _Theta, _neededAngle, _numberOfCycles, _angleUpperBound, _type),
		alp(_alp) {}
};
