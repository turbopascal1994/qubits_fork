#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <complex>
#include <algorithm>
#include <omp.h>
#include "linalg.h"
#include "compareDouble.h"
#include <map>

using namespace std;

const double PI = acos(-1.0);

vector<int> CreateStartSCALLOP(double w01, double wt, double w, int M, double AmpThreshold){
	
	double Tt = PI * 2 / wt;
	double LenImp = Tt * M;

	double a = 0;
	double b = w;

	vector<int> NewSequence(M);
	int i = 0;
	int j = 0;

	while (Less_equal(b, LenImp)){
		int iSeq = 0;
		if (More_equal(sin(w01 * a), AmpThreshold)){
			NewSequence[j] = 1;
			iSeq = 1;
		}
		if (Less_equal(sin(w01 * a), -AmpThreshold)){
			NewSequence[j] = -1;
			iSeq = 1;
		}
		if (iSeq == 0){
			NewSequence[j] = 0;
		}
		a = Tt * (i + 1) - w / 2;
		b = Tt * (i + 1) + w / 2;
		i++;
		j++;
	}
	return NewSequence;
}

vector<double> CreateAmpThresholds(const double w01, const double wt, const double w, const int M){
	double Tt = 2.0 * PI / wt;
	double LenImp = M * Tt + w;
	double StartA = 0.1;
	vector<double> AmpsArrayRAW;
	int PrevN = -1;
	for(double AmpThreshold = 0.1;Less_equal(AmpThreshold, 0.9, 1e-4); AmpThreshold += 0.01){
		int PulsesN = 0;
		if (Equal(AmpThreshold, StartA)){
			PrevN = 0;
		}
		double a = 0;
		double b = w;
		int i = 1;
		int j = 1;
		while(Less_equal(b, LenImp)){
			if (More_equal(sin(w01 * a), AmpThreshold)){
				PulsesN++;
			}
			if (Less_equal(sin(w01 * a), -AmpThreshold)){
				PulsesN++;
			}
			a = i * Tt - w / 2;
			b = i * Tt + w / 2;
			++i;
			++j;
		}
		if (PulsesN != PrevN){
			AmpsArrayRAW.push_back(AmpThreshold);
			PrevN = PulsesN;
		}
	}

	int TrueCount = 0;
	for(int i = 0;i < AmpsArrayRAW.size();i++){
		if (Not_equal(AmpsArrayRAW[i], 0.0)){
			TrueCount++;
		}
	}
	vector<double> AmpsArray(TrueCount);
	for(int i = 0;i < TrueCount;i++){
		AmpsArray[i] = AmpsArrayRAW[i];
	}
	return AmpsArray;
}

const int DT_SIZE = 106480;
// const int DT_SIZE = 96000;

double* dt;
vector<complex<double>> Id, H0, EigVec(9), EigVal(3), WF1, WF3, Hmatrix, HrPlus(9), HrMinus(9), HrZero(9), InitStates;
const double h = 1.054e-34;
const double C1 = 1e-12;
const double F0 = 2.06e-15;

double calcProbability(int N,
						const int* imp,
						double Amp,
						vector<complex<double>>& UPlus,
						vector<complex<double>>& UMinus,
						vector<complex<double>>& UZero,
						vector<complex<double>>& WF,
						vector<complex<double>>& WF1) {
	vector<complex<double>> newWF(N);
	vector<complex<double>> res(1);
	vector<pair<int, int>> counter;
	int now_ind = -2;
	int cnt = 0;

	auto make_divisible = [](int a) {
		for (int i = a - 2; i <= a + 2; i++) {
			if (i != 0 && i % 80 == 0) return i;
		}
		exit(1);
	};

	for (int j = 0; j < DT_SIZE; j++) {
		int cur_ind;
		if (imp[j] == 1) {
			cur_ind = 1;
		}
		else if (imp[j] == -1) {
			cur_ind = -1;
		}
		else if (imp[j] == 0) {
			cur_ind = 0;
		}
		else exit(1);

		if (now_ind != cur_ind) {
			if (cnt) counter.push_back({ now_ind, make_divisible(cnt) });
			now_ind = cur_ind;
			cnt = 1;
		}
		else {
			cnt++;
		}
	}
	if (cnt) counter.push_back({ now_ind, make_divisible(cnt) });

	auto UPlus80 = pow_matrix(UPlus, 80, N);
	auto UMinus80 = pow_matrix(UMinus, 80, N);
	auto UZero80 = pow_matrix(UZero, 80, N);
	auto UZero720 = pow_matrix(UZero, 720, N);

	auto UPlusZero = mult(UPlus80, UZero720, N, N, N, N, N, N);
	auto UMinusZero = mult(UMinus80, UZero720, N, N, N, N, N, N);
	auto UZeroZero = mult(UZero80, UZero720, N, N, N, N, N, N);

	vector<pair<int, int>> new_counter;
	int last = 720;
	for (int i = 0; i < counter.size(); i++) {
		while (counter[i].second) {
			if (last == 720) {
				last = 80;
			}
			else {
				last = 720;
			}
			new_counter.push_back({ counter[i].first, last });
			counter[i].second -= last;
		}
	}

	vector<complex<double>> resU(N * N);
	for (int i = 0; i < N; i++) resU[i * N + i] = { 1, 0 };
	new_counter.resize(240);
	
	for (size_t i = 0; i < new_counter.size(); ++i) {
		vector<complex<double>> U;
		if (i + 1 < new_counter.size()) {
			if (new_counter[i].first == -1) {
				if (new_counter[i + 1].first != 0) exit(1);
				U = UMinusZero;
			}
			else if (new_counter[i].first == 1) {
				if (new_counter[i + 1].first != 0) exit(1);
				U = UPlusZero;
			}
			else {
				if (new_counter[i + 1].first != 0) exit(1);
				U = UZeroZero;
			}
			++i;
		}
		else {
			if (new_counter[i].first != 0) exit(1);
			U = UZero80;
		}
		// resU = mult(resU, U, N, N, N, N, N, N);
		WF = mult(U, WF, N, 1, N, N, 1, 1);
	}


	res[0] = { 0, 0 };
	for (int i = 0; i < N; i++) {
		res[0] += WF[i] * conj(WF1[i]);
	}
	return norm(res[0]);
}

void calcImpSequence(double w, double Amp, double T, int NumberOfCycles, vector<int>& SignalString, int* y) {
	double a, b;
	int CellsNumber = SignalString.size();
	int k = 1, kk = 1;
	fill(y, y + DT_SIZE, 0);
	while (kk <= NumberOfCycles) {
		while (k <= CellsNumber) {
			a = T * (CellsNumber * (kk - 1) + k - 1);
			b = w + a;
			for (size_t i = 0; i < DT_SIZE; ++i) {
				if (More_equal(dt[i], a) && Less_equal(dt[i], b)) {
					y[i] += SignalString[k - 1];
				}
			}
			k++;
		}
		kk++;
		k = 1;
	}
}

void precalcFidelityAndRotateCheck(double tstep, double w01, double w12){
	Id = { {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0} };
	dt = new double[DT_SIZE];
	double iter = 0;
	for (size_t i = 0; i < DT_SIZE; ++i) {
		dt[i] = iter;
		iter += tstep;
	}

	H0 = { {0, 0}, {0, 0}, {0, 0}, {0, 0}, {h * w01, 0}, {0, 0}, {0, 0}, {0, 0}, {h * (w01 + w12), 0} };
	eig(H0, EigVec, EigVal, 3);

	vector<int> indices = { 0, 1, 2 };
	sort(indices.begin(), indices.end(), [&](int a, int b) {
		return Less(EigVal[a].real(), EigVal[b].real());
	});
	WF1 = { EigVec[indices[0]], EigVec[indices[0] + 3], EigVec[indices[0] + 6] };
	WF3 = { EigVec[indices[2]], EigVec[indices[2] + 3], EigVec[indices[2] + 6] };

	Hmatrix = { {0, 0}, {-1, 0}, {0, 0}, {1, 0}, {0, 0}, {-sqrt(2), 0}, {0, 0}, {sqrt(2), 0}, {0, 0} };
	InitStates = {
		{1, 0}, {0, 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0},
		{0, 0}, {1, 0}, {1.0 / sqrt(2), 0}, {-1.0 / sqrt(2), 0}, {0, 1.0 / sqrt(2)}, {0, -1.0 / sqrt(2)},
		{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
	};
}

double Fidelity(double w01, double w12, double w, double T, double Len, double tstep, vector<int>& SignalString, int NumberOfCycles, double Theta) {
	int CellsNumber = SignalString.size();
	
	double V = F0 / w;
	double Cc = Theta / (F0 * sqrt(2.0 * w01 / (h * C1)));
	double Amp = Cc * V * sqrt(h * w01 / (2.0 * C1));
	int* imp = new int[DT_SIZE];
	calcImpSequence(w, Amp, T, NumberOfCycles, SignalString, imp);
	double sum = 0;
	for (size_t i = 0; i < H0.size(); ++i) {
		HrPlus[i] = { H0[i].real(), Amp * Hmatrix[i].real() };
		HrMinus[i] = { H0[i].real(), -Amp * Hmatrix[i].real() };
		HrZero[i] = H0[i];
	}

	auto UPlus = getUMatrix(Id, HrPlus, tstep, h, 3);
	auto UMinus = getUMatrix(Id, HrMinus, tstep, h, 3);
	auto UZero = getUMatrix(Id, HrZero, tstep, h, 3);

	for(size_t IS = 0;IS < 6;++IS){
		vector<complex<double>> WF = { InitStates[IS], InitStates[IS + 6], InitStates[IS + 12] };
		sum += calcProbability(3, imp, Amp, UPlus, UMinus, UZero, WF, WF3);
	}
	return sum / 6;
}

double RotateCheck(double w01, double w12, double w, double T, double Len, double tstep, int NumberOfCycles, vector<int>& InputSequence, double Theta){
	double V = F0 / w;
	double Cc = Theta / (F0 * sqrt(2.0 * w01 / (h * C1)));
	double Amp = Cc * V * sqrt(h * w01 / (2.0 * C1));
	int CellsNumber = InputSequence.size();
	int* imp = new int[DT_SIZE];
	calcImpSequence(w, Amp, T, NumberOfCycles, InputSequence, imp);

	vector<complex<double>> WF(WF1);
	for(int i = 0;i < H0.size();i++){
		HrPlus[i] = { H0[i].real(), Amp * Hmatrix[i].real() };
		HrMinus[i] = { H0[i].real(), -Amp * Hmatrix[i].real() };
		HrZero[i] = H0[i];
	}
	
	auto UPlus = getUMatrix(Id, HrPlus, tstep, h, 3);
	auto UMinus = getUMatrix(Id, HrMinus, tstep, h, 3);
	auto UZero = getUMatrix(Id, HrZero, tstep, h, 3);
	
	return calcProbability(3, imp, Amp, UPlus, UMinus, UZero, WF, WF1);
}

double NewThetaOptimizer(double w01, double w12, double w, double T, double Len, double tstep, vector<int> & InputSequence, int NumberOfCycles, double UnOptTheta){
	double SmolStep = 0.0005;
	double BigStep = 0.005;

	double PCrit = 0.01;
	double StepCrit = 0.15;
	double P = 0;
	double Theta = UnOptTheta;

	while(More(abs(P - 0.5), PCrit)){
		P = RotateCheck(w01, w12, w, T, Len, tstep, NumberOfCycles, InputSequence, Theta);
		double ThetaStep = 0;
		if (More(abs(P - 0.5), PCrit)){
			if (More_equal(abs(P - 0.5), StepCrit)){
				if (More(P, 0.5)){
					ThetaStep = BigStep;
				} else {
					ThetaStep = -BigStep;
				}
			} else {
				if (More(P, 0.5)){
					ThetaStep = SmolStep;
				} else {
					ThetaStep = -SmolStep;
				}
			}
		}
		Theta += ThetaStep;
	}
	return Theta;
}

void WriteSequence(vector<int> & InputString){
	for(auto &i : InputString){
		cout << i;
	}
	cout << '\n';
}

vector<int> ChangingOneElement(vector<int> &InputSequence, int index, int type) {
	auto Sequence = InputSequence;

	if (InputSequence[index] == -1){
		if (type == 0){
			Sequence[index] = 1;
		} else {
			Sequence[index] = 0;
		}
	} else if (InputSequence[index] == 1){
		if (type == 0) {
			Sequence[index] = -1;
		}
		else {
			Sequence[index] = 0;
		}
	} else {
		if (type == 0) {
			Sequence[index] = 1;
		}
		else {
			Sequence[index] = -1;
		}
	}
	return Sequence;
}

void GenSearch(double w01, double w12, double w, double T, double Len,
			   double tstep, vector<int> InputString, int NumberOfCycles,
			   double StartTheta, vector<int> & BestSequenceOverall, double & BestLeakOverall, double & BestAngleOverall){
	int Gen = 1;
	double NewLeak = 0;
	double GoldLeak = 1;
	vector<int> BestSequence;
	double Theta = StartTheta;
	cout << "Стартовая последовательность\n";
	WriteSequence(InputString);
	int SequencesNumber = 2 * InputString.size();
	vector<double> FArray(SequencesNumber, 0);
	vector<double> AnglesArray(SequencesNumber, 0);

	while(Less(NewLeak, GoldLeak)){
		cout << "Поколение #" << Gen << '\n';
		if (Gen > 1){
			GoldLeak = NewLeak;
			InputString = BestSequence;
		}
#pragma omp parallel for
		for(size_t j = 0;j < SequencesNumber;++j){
			auto SignalString = ChangingOneElement(InputString, j / 2, j % 2);
			double OptimizedTheta = NewThetaOptimizer(w01, w12, w, T, Len, tstep, SignalString, NumberOfCycles, Theta);
			AnglesArray[j] = OptimizedTheta;
			FArray[j] = Fidelity(w01, w12, w, T, Len, tstep, SignalString, NumberOfCycles, OptimizedTheta);
		}
		double MinLeak = FArray[0];
		int MinLeakN = 0;
		for(size_t j = 1;j < FArray.size();++j){
			if (Less(FArray[j], MinLeak)){
				MinLeak = FArray[j];
				MinLeakN = j;
			}
		}
		BestSequence = ChangingOneElement(InputString, MinLeakN / 2, MinLeakN % 2);
		auto BestAngle = AnglesArray[MinLeakN];
		cout << "Минимальная утечка у последовательности ";
		for (auto& i : BestSequence) cout << i;
		cout << '\n';
		cout << "F = " << MinLeak << '\n';
		cout << "Угол Th = " << BestAngle << '\n';
		NewLeak = MinLeak;
		if (Less(NewLeak, GoldLeak)){
			BestLeakOverall = MinLeak;
			BestSequenceOverall = BestSequence;
			BestAngleOverall = BestAngle;
		}
		++Gen;
	}
	cout << "Поиск занял " << Gen - 1 << " поколений\n";
	cout << "Самая лучшая последовательность\n";
	WriteSequence(BestSequenceOverall);
	cout << "F = " << BestLeakOverall << '\n';
	cout << "Угол Th = " << BestAngleOverall << '\n';
}

int main(){
	setlocale(LC_ALL, "Rus");
	cout.precision(20);

	auto program_start = chrono::high_resolution_clock::now();
	const double tstep = 5e-14;
	const double w01 = PI * 8e9;
	const double w12 = w01 - PI * 5e8;
	const double wt = PI * 5e10;
	const double w = 4e-12;
	const double T = PI / wt * 2;
	double Theta = 0.001;
	int CellsNumber = 120;
	const int NumberOfCycles = 1;

	const double Len = (double)CellsNumber * T * NumberOfCycles;
	const double NeededAngle = 0.024;
	const double AngleThreshold = 0.0005;

	precalcFidelityAndRotateCheck(tstep, w01, w12);

	auto Amps = CreateAmpThresholds(w01, wt, w, CellsNumber);

	double A1 = Amps[0];
	double A2 = Amps[Amps.size() / 2 - 1];
	
	cout << "Берём наименьшую амплитуду А1.\n";
	vector<int> SignalString = CreateStartSCALLOP(w01, wt, w, CellsNumber, A1);
	double NewTheta = NewThetaOptimizer(w01, w12, w, T, Len, tstep, SignalString, 1, Theta);
	cout << NewTheta << '\n';
	cout << "A1 = " << A1 << '\n';

	vector<int> A1Sequence;
	double A1F, A1Angle;
	auto t_start = chrono::high_resolution_clock::now();
	GenSearch(w01, w12, w, T, Len, tstep, SignalString, NumberOfCycles, NewTheta, A1Sequence, A1F, A1Angle);
	auto t_end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<std::chrono::seconds>(t_end - t_start).count();
	cout << "Время вычисления алгоритма = " << duration << " сек.\n";
	
	Theta = 0.001;
	cout << "Берём наибольшую амплитуду A2.\n";
	SignalString = CreateStartSCALLOP(w01, wt, w, CellsNumber, A2);
	NewTheta = NewThetaOptimizer(w01, w12, w, T, Len, tstep, SignalString, 1, Theta);
	cout << "A2 = " << A2 << '\n';
	vector<int> A2Sequence;
	double A2F, A2Angle;
	GenSearch(w01, w12, w, T, Len, tstep, SignalString, NumberOfCycles, NewTheta, A2Sequence, A2F, A2Angle);


	// Interpolation
	cout << "Интерполируем амплитуду до нужной, А3.\n";
	auto b_int = (A2Angle - A1Angle) / (A2 - A1);
	auto a_int = A1Angle - b_int * A1;
	auto NewAmpRough = (NeededAngle - a_int) / b_int;

	vector<double> AmpChecks(Amps.size());
	for(size_t i = 0;i < AmpChecks.size();++i){
		AmpChecks[i] = abs(Amps[i] - NewAmpRough);
	}
	int MinAmpCheckN = min_element(AmpChecks.begin(), AmpChecks.end()) - AmpChecks.begin();
	double NewAmp = Amps[MinAmpCheckN];
	Theta = 0.001;

	SignalString = CreateStartSCALLOP(w01, wt, w, CellsNumber, NewAmp);
	NewTheta = NewThetaOptimizer(w01, w12, w, T, Len, tstep, SignalString, NumberOfCycles, NewTheta);
	vector<int> A3Sequence;
	double A3F, A3Angle;
	GenSearch(w01, w12, w, T, Len, tstep, SignalString, NumberOfCycles, NewTheta, A3Sequence, A3F, A3Angle);
	cout << "Для новой амлитуды А3 = " << NewAmp << "; угол Th3 = " << A3Angle << '\n';
	cout << "Желаемый угол Th0 = " << NeededAngle << '\n';
	
	vector<int> NewSequence;
	double NewF = A3F, NewAngle = A3Angle;
	if (Equal(NeededAngle, NewAngle)){
		cout << "Изменение длины последовательности ничего не даст\n";
		NewSequence = A3Sequence;
	}
	else {
		cout << "Теперь меняем длину последовательности М.\n";
		if (Less(NewAngle, NeededAngle + AngleThreshold)){
			while(Not_equal(NewAngle, NeededAngle + AngleThreshold)){
				CellsNumber--;
				cout << "M = " << CellsNumber << '\n';
				SignalString = CreateStartSCALLOP(w01, wt, w, CellsNumber, NewAmp);
				GenSearch(w01, w12, w, T, Len, tstep, SignalString, NumberOfCycles, NewTheta, NewSequence, NewF, NewAngle);
				if (Equal(NewAngle, NeededAngle) && Less(NewF, 1e-4)){
					break;
				}
			}
		} else if (More(NewAngle, NeededAngle - AngleThreshold)){
			while (Not_equal(NewAngle, NeededAngle - AngleThreshold)) {
				CellsNumber++;
				cout << "M = " << CellsNumber << '\n';
				SignalString = CreateStartSCALLOP(w01, wt, w, CellsNumber, NewAmp);
				GenSearch(w01, w12, w, T, Len, tstep, SignalString, NumberOfCycles, NewTheta, NewSequence, NewF, NewAngle);
				if (Equal(NewAngle, NeededAngle) && Less(NewF, 1e-4)) {
					break;
				}
			}
		}
	}

	cout << "Желаемый угол Th0 = " << NeededAngle << '\n';
	cout << "Новый угол ThN = " << NewAngle << '\n';
	cout << "Конечная последовательность\n";
	WriteSequence(NewSequence);
	cout << "Утечка F = " << NewF << '\n';
	auto program_end = chrono::high_resolution_clock::now();
	auto program_duration = chrono::duration_cast<std::chrono::seconds>(program_end - program_start).count();
	cout << "Время работы программы = " << program_duration << " сек.\n";
	return 0;
}