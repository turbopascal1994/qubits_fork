#pragma once

struct ConstantsDescriptor {
public:
	double w01, w12, wt, w, T, Tstep, Theta;
	int numberOfCycles;
	double neededAngle;
	double angleUpperBound;
	int type;
	ConstantsDescriptor(
		double _w01, double _w12, double _wt, double _w,
		double _T, double _Tstep, double _Theta,
		double _neededAngle, int _numberOfCycles,
		double _angleUpperBound, int _type) :
		w01(_w01), w12(_w12), wt(_wt), w(_w),
		T(_T), Tstep(_Tstep), Theta(_Theta),
		neededAngle(_neededAngle), numberOfCycles(_numberOfCycles),
		angleUpperBound(_angleUpperBound), type(_type) {}
};
