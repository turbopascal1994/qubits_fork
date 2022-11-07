#include "kernel.h"
#include<algorithm>


void Kernel::precalcFidelityAndRotateCheck(double tstep, double w01, double w12) {

	Id = { {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0} };

	H0 = { {0, 0}, {0, 0}, {0, 0}, {0, 0}, {h * w01, 0}, {0, 0}, {0, 0}, {0, 0}, {h * (w01 + w12), 0} };
	EigVec.resize(9);
	EigVal.resize(3);
	eig(H0, EigVec, EigVal, 3);

	vector<int> indices = { 0, 1, 2 };
	sort(indices.begin(), indices.end(), [&](int a, int b) {
		return EigVal[a].real() < EigVal[b].real();
	});
	WF1 = { EigVec[indices[0]], EigVec[indices[0] + 3], EigVec[indices[0] + 6] };
	WF3 = { EigVec[indices[2]], EigVec[indices[2] + 3], EigVec[indices[2] + 6] };

	Hmatrix = { {0, 0}, {-1, 0}, {0, 0}, {1, 0}, {0, 0}, {-sqrt(2), 0}, {0, 0}, {sqrt(2), 0}, {0, 0} };
	InitStates = {
		{1, 0}, {0, 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0},
		{0, 0}, {1, 0}, {1.0 / sqrt(2), 0}, {-1.0 / sqrt(2), 0}, {0, 1.0 / sqrt(2)}, {0, -1.0 / sqrt(2)},
		{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
	};
};

double Kernel::calcProbability(int N, const vector<int>& InputSequence, const int CyclePlusMinusSteps, const int CycleZeroSteps, const vector<complex<double>>& UPlus, const vector<complex<double>>& UMinus, const vector<complex<double>>& UZero, vector<complex<double>>& WF, const vector<complex<double>>& WF1) {
	auto UPlus80 = pow_matrix(UPlus, CyclePlusMinusSteps, N);
	auto UMinus80 = pow_matrix(UMinus, CyclePlusMinusSteps, N);
	auto UZero80 = pow_matrix(UZero, CyclePlusMinusSteps, N);
	auto UZero720 = pow_matrix(UZero, CycleZeroSteps, N);

	auto UPlusZero = mult(UPlus80, UZero720, N, N, N, N, N, N);
	auto UMinusZero = mult(UMinus80, UZero720, N, N, N, N, N, N);
	auto UZeroZero = mult(UZero80, UZero720, N, N, N, N, N, N);

	vector<complex<double>> resU(N * N);
	for (int i = 0; i < N; i++) resU[i * N + i] = { 1, 0 };
	for (int i = 0; i < InputSequence.size(); i++) {
		vector<complex<double>> U = UZeroZero;
		if (InputSequence[i] == -1) {
			U = UMinusZero;
		}
		else if (InputSequence[i] == 1) {
			U = UPlusZero;
		}
		WF = mult(U, WF, N, 1, N, N, 1, 1);
	}

	complex<double> res = { 0, 0 };
	for (int i = 0; i < N; i++) {
		res += WF[i] * conj(WF1[i]);
	}
	return norm(res);
};

double Kernel::RotateCheck(const vector<int>& InputSequence, double Theta) {
	double V = F0 / w;
	double Cc = Theta / (F0 * sqrt(2.0 * w01 / (h * C1)));
	double Amp = Cc * V * sqrt(h * w01 / (2.0 * C1));
	int CellsNumber = InputSequence.size();

	//	Число тактов с импульсом и без него на одном периоде тактовой частоты генератора
	int CycleSteps = floor(T / tstep + 0.5);
	int CyclePlusMinusSteps = floor(CycleSteps * w / T + 0.5);
	int CycleZeroSteps = CycleSteps - CyclePlusMinusSteps;

	vector<complex<double>> WF(WF1);
	vector<complex<double>> HrPlus(H0.size());
	vector<complex<double>> HrMinus(H0.size());
	vector<complex<double>> HrZero(H0.size());

	for (int i = 0; i < H0.size(); i++) {
		HrPlus[i] = { H0[i].real(), Amp * Hmatrix[i].real() };
		HrMinus[i] = { H0[i].real(), -Amp * Hmatrix[i].real() };
		HrZero[i] = H0[i];
	}

	auto UPlus = getUMatrix(Id, HrPlus, tstep, h, 3);
	auto UMinus = getUMatrix(Id, HrMinus, tstep, h, 3);
	auto UZero = getUMatrix(Id, HrZero, tstep, h, 3);

	return calcProbability(3, InputSequence, CyclePlusMinusSteps, CycleZeroSteps, UPlus, UMinus, UZero, WF, WF1);
}
vector<int> Kernel::CreateStartSCALLOP(int M, double AmpThreshold) {
	double Tt = PI * 2 / wt;

	vector<int> NewSequence(M);
	for (int A = 0; A < M; A++) {
		double sina = sin(w01 * (A * Tt + w / 2.0));
		if (sina >= AmpThreshold) {
			NewSequence[A] = 1;
		}
		else if (sina <= -AmpThreshold) {
			NewSequence[A] = (TYPE == 3 ? -1 : 1);
		}
		else {
			NewSequence[A] = 0;
		}
	}

	return NewSequence;
}
vector<double> Kernel::CreateAmpThresholds(const int M) {
	double Tt = 2.0 * PI / wt;
	int PrevN = 0;
	vector<double> AmpsArray;
	for (int i = 0; i <= 80; i++) {
		double AmpThreshold = 0.1 + 0.01 * i;
		int PulsesN = 0;

		for (int A = 0; A < M; A++) {
			double sina = sin(w01 * (A * Tt + w / 2.0));
			if (sina >= AmpThreshold) {
				PulsesN++;
			}
			if (sina <= -AmpThreshold) {
				PulsesN++;
			}
		}

		if (PulsesN != PrevN) {
			AmpsArray.push_back(AmpThreshold);
			PrevN = PulsesN;
		}
	}
	return AmpsArray;
}
double Kernel::Fidelity(const vector<int>& SignalString, double Theta) {
	int CellsNumber = SignalString.size();

	double V = F0 / w;
	double Cc = Theta / (F0 * sqrt(2.0 * w01 / (h * C1)));
	double Amp = Cc * V * sqrt(h * w01 / (2.0 * C1));

	//	Число тактов с импульсом и без него на одном периоде тактовой частоты генератора
	int CycleSteps = floor(T / tstep + 0.5);
	int CyclePlusMinusSteps = floor(CycleSteps * w / T + 0.5);
	int CycleZeroSteps = CycleSteps - CyclePlusMinusSteps;

	double sum = 0;
	vector<complex<double>> HrPlus(H0.size());
	vector<complex<double>> HrMinus(H0.size());
	vector<complex<double>> HrZero(H0.size());
	for (size_t i = 0; i < H0.size(); ++i) {
		HrPlus[i] = { H0[i].real(), Amp * Hmatrix[i].real() };
		HrMinus[i] = { H0[i].real(), -Amp * Hmatrix[i].real() };
		HrZero[i] = H0[i];
	}

	auto UPlus = getUMatrix(Id, HrPlus, tstep, h, 3);
	auto UMinus = getUMatrix(Id, HrMinus, tstep, h, 3);
	auto UZero = getUMatrix(Id, HrZero, tstep, h, 3);

	for (size_t IS = 0; IS < 6; ++IS) {
		vector<complex<double>> WF = { InitStates[IS], InitStates[IS + 6], InitStates[IS + 12] };
		sum += calcProbability(3, SignalString, CyclePlusMinusSteps, CycleZeroSteps, UPlus, UMinus, UZero, WF, WF3);
	}
	return sum / 6;
}
double Kernel::NewThetaOptimizer(const vector<int>& InputSequence, double UnOptTheta) {
	double CurrentStep = 0.01;
	double MinStep = 1e-9;
	double PPrecision = 1e-9;
	double CurrentTheta = UnOptTheta;
	double PreviousTheta = CurrentTheta;

	double P = RotateCheck(InputSequence, CurrentTheta);
	for (int iteration = 0; iteration < 1000 && abs(P - 0.5) > PPrecision && CurrentStep >= MinStep; ++iteration) {
		CurrentTheta = PreviousTheta + CurrentStep;
		P = RotateCheck(InputSequence, CurrentTheta);
		if (P < 0.5) {
			CurrentStep /= 2.0;
		}
		else {
			PreviousTheta = CurrentTheta;
		}
	}
	return CurrentTheta;
}
void Kernel::WriteSequence(ostream& fout, const vector<int>& InputString, char End) {
	for (auto& i : InputString) {
		fout << i;
	}
	fout << End;
}
vector<int> Kernel::ChangingOneElement(vector<int> InputSequence, int index, int type) {
	if (TYPE != 3) {
		InputSequence[index] ^= 1;
		return InputSequence;
	}
	if (InputSequence[index] == -1) {
		if (type == 0) {
			InputSequence[index] = 1;
		}
		else {
			InputSequence[index] = 0;
		}
	}
	else if (InputSequence[index] == 1) {
		if (type == 0) {
			InputSequence[index] = -1;
		}
		else {
			InputSequence[index] = 0;
		}
	}
	else {
		if (type == 0) {
			InputSequence[index] = 1;
		}
		else {
			InputSequence[index] = -1;
		}
	}
	return InputSequence;
}
void Kernel::GenSearch(vector<int> InputString, double StartTheta, vector<int>& BestSequenceOverall, double& BestLeakOverall, double& BestAngleOverall) {
	int Gen = 1;
	double NewLeak = 0;
	double GoldLeak = 1;
	vector<int> BestSequence;
	double Theta = StartTheta;
	cout << "Стартовая последовательность\n";
	WriteSequence(cout, InputString);
	int SequencesNumber = InputString.size() * (TYPE == 3 ? 2 : 1);
	vector<double> FArray(SequencesNumber, 0);
	vector<double> AnglesArray(SequencesNumber, 0);

	while (NewLeak < GoldLeak) {
		if (Gen > 1) {
			GoldLeak = NewLeak;
			InputString = BestSequence;
		}
#pragma omp parallel shared(InputString, Theta, AnglesArray, FArray)
		{
			int j;
#pragma omp for private(j)
			for (j = 0; j < SequencesNumber; ++j) {
				auto SignalString = ChangingOneElement(InputString, j / 2, j % 2);
				AnglesArray[j] = NewThetaOptimizer(SignalString, Theta);
				FArray[j] = Fidelity(SignalString, AnglesArray[j]);
			}
		}

		double MinLeak = FArray[0];
		int MinLeakN = 0;
		for (size_t j = 1; j < FArray.size(); ++j) {
			if (FArray[j] < MinLeak) {
				MinLeak = FArray[j];
				MinLeakN = j;
			}
		}
		BestSequence = ChangingOneElement(InputString, MinLeakN / 2, MinLeakN % 2);
		auto BestAngle = AnglesArray[MinLeakN];

		/*cout << "Поколение #" << Gen << '\n';
		cout << "Минимальная утечка у последовательности ";
		for (auto& i : BestSequence) cout << i;
		cout << '\n';
		cout << "F = " << MinLeak << '\n';
		cout << "Угол Th = " << BestAngle << '\n';*/
		NewLeak = MinLeak;
		if (NewLeak < GoldLeak) {
			BestLeakOverall = MinLeak;
			BestSequenceOverall = BestSequence;
			BestAngleOverall = BestAngle;
		}
		++Gen;
	}
	cout << "Поиск занял " << Gen - 1 << " поколений\n";
	cout << "Самая лучшая последовательность\n";
	WriteSequence(cout, BestSequenceOverall);
	cout << "F = " << BestLeakOverall << '\n';
	cout << "Угол Th = " << BestAngleOverall << '\n';
}
;
