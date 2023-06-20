#include "kernel.h"


void Kernel::precalcFidelityAndRotateCheck(double tstep, double w01, double w12) {

	Id = { {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0} };

	H0 = { {0, 0}, {0, 0}, {0, 0}, {0, 0}, {h * w01, 0}, {0, 0}, {0, 0}, {0, 0}, {h * (w01 + w12), 0} };
	EigVec.resize(9);
	EigVecRight.resize(9);
	EigVal.resize(3);
	linalg::eig(H0, EigVec, EigVecRight, EigVal, 3);

	vector<int> indices = { 0, 1, 2 };
	std::sort(indices.begin(), indices.end(), [&](int a, int b) {
		return EigVal[a].real() < EigVal[b].real();
	});
	WF1 = { EigVec[indices[0]], EigVec[indices[0] + 3], EigVec[indices[0] + 6] };
	WF2 = { EigVec[indices[1]], EigVec[indices[1] + 3], EigVec[indices[1] + 6] };
	WF3 = { EigVec[indices[2]], EigVec[indices[2] + 3], EigVec[indices[2] + 6] };

	Hmatrix = { {0, 0}, {-1, 0}, {0, 0}, {1, 0}, {0, 0}, {-sqrt(2), 0}, {0, 0}, {sqrt(2), 0}, {0, 0} };
	InitStates = {
		{1, 0}, {0, 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0},
		{0, 0}, {1, 0}, {1.0 / sqrt(2), 0}, {-1.0 / sqrt(2), 0}, {0, 1.0 / sqrt(2)}, {0, -1.0 / sqrt(2)},
		{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
	};
	IdealStates = {
		{1.0 / sqrt(2), 0}, {-1.0 / sqrt(2), 0}, {0, 0}, {1, 0}, {0.5, -0.5}, {0.5, 0.5},
		{1.0 / sqrt(2), 0}, {1.0 / sqrt(2), 0}, {1, 0}, {0, 0}, {0.5, 0.5}, {0.5, -0.5},
		{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
	};
};

void Kernel::printState(ofstream& out, const vector<complex<double>>& state, int IS) {
	for (int j = 0; j < state.size(); ++j) {
		out << state[j].real() << ' ' << state[j].imag() << ' ';
	}
	out << IS << '\n';
}

Kernel::IntegrationResult Kernel::integration(
	int N,
	const vector<int>& InputSequence, 
	const vector<complex<double>>& UTZero,
	const vector<complex<double>>& UTMinus,
	const vector<complex<double>>& UTPlus,
	vector<complex<double>>& WF,
	int IS,
	ofstream* out_stream
) {
	if (out_stream) {
		printState(*out_stream, WF, IS);
	}
	for (int i = 0; i < InputSequence.size(); i++) {
		vector<complex<double>> U = UTZero;
		if (InputSequence[i] == -1) {
			U = UTMinus;
		}
		else if (InputSequence[i] == 1) {
			U = UTPlus;
		}
		WF = linalg::matmul(U, WF, N, 1, N, N, 1, 1);
		if (out_stream) {
			printState(*out_stream, WF, IS);
		}
	}

	complex<double> l0(0), l1(0), l2(0);
	for (int i = 0; i < N; i++) {
		l0 += WF[i] * conj(WF1[i]);
		l1 += WF[i] * conj(WF2[i]);
		l2 += WF[i] * conj(WF3[i]);
	}
	return {norm(l0), norm(l1), norm(l2), WF};
};

void Kernel::prepareUMatrices(double Theta, vector<complex<double>>& UTZero, vector<complex<double>>& UTMinus, vector<complex<double>>& UTPlus) {
	double V = F0 / w;
	double Cc = Theta / (F0 * sqrt(2.0 * w01 / (h * C1)));
	double Amp = Cc * V * sqrt(h * w01 / (2.0 * C1));

	//	����� ������ � ��������� � ��� ���� �� ����� ������� �������� ������� ����������
	int CycleSteps = floor(T / tstep + 0.5);
	int CyclePlusMinusSteps = floor(CycleSteps * w / T + 0.5);
	int CycleZeroSteps = CycleSteps - CyclePlusMinusSteps;

	vector<complex<double>> HrPlus(H0.size());
	vector<complex<double>> HrMinus(H0.size());
	vector<complex<double>> HrZero(H0.size());
	for (size_t i = 0; i < H0.size(); ++i) {
		HrPlus[i] = { H0[i].real(), Amp * Hmatrix[i].real() };
		HrMinus[i] = { H0[i].real(), -Amp * Hmatrix[i].real() };
		HrZero[i] = H0[i];
	}

	auto UZeroStep = linalg::getUMatrix(Id, HrZero, tstep, h, 3);
	auto UMinusStep = linalg::getUMatrix(Id, HrMinus, tstep, h, 3);
	auto UPlusStep = linalg::getUMatrix(Id, HrPlus, tstep, h, 3);

	auto UZero = linalg::matpow(UZeroStep, CyclePlusMinusSteps, 3);
	auto UMinus = linalg::matpow(UMinusStep, CyclePlusMinusSteps, 3);
	auto UPlus = linalg::matpow(UPlusStep, CyclePlusMinusSteps, 3);

	auto UT = linalg::matpow(UZeroStep, CycleZeroSteps, 3);

	UTZero = linalg::matmul(UT, UZero, 3, 3, 3, 3, 3, 3);
	UTMinus = linalg::matmul(UT, UMinus, 3, 3, 3, 3, 3, 3);
	UTPlus = linalg::matmul(UT, UPlus, 3, 3, 3, 3, 3, 3);
}
double Kernel::RotateCheck(const vector<int>& InputSequence, double Theta) {
	vector<complex<double>> WF(WF1);
	vector<complex<double>> UTZero, UTMinus, UTPlus;
	prepareUMatrices(Theta, UTZero, UTMinus, UTPlus);
	return integration(3, InputSequence, UTZero, UTMinus, UTPlus, WF, 1).level0;
}
vector<int> Kernel::CreateStartSCALLOP(int M, double AmpThreshold, int TYPE) {
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
Kernel::FidelityResult Kernel::Fidelity(const vector<int>& SignalString, double Theta, ofstream* out_stream) {
	double leak = 0, fidelity = 0;
	vector<complex<double>> UTZero, UTMinus, UTPlus;
	prepareUMatrices(Theta, UTZero, UTMinus, UTPlus);
	for (size_t IS = 0; IS < 6; ++IS) {
		vector<complex<double>> WF = { InitStates[IS], InitStates[IS + 6], InitStates[IS + 12] };
		vector<complex<double>> ideal = { IdealStates[IS], IdealStates[IS + 6], IdealStates[IS + 12] };
		auto res = integration(3, SignalString, UTZero, UTMinus, UTPlus, WF, IS, out_stream);
		leak += res.level2;
		complex<double> dot_product(0);
		for (int i = 0; i < 3; ++i) {
			dot_product += res.state[i] * conj(ideal[i]);
		}
		fidelity += norm(dot_product);
	}
	leak /= 6;
	fidelity /= 6;
	return Kernel::FidelityResult( fidelity, leak );
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

vector<int> Kernel::ChangingOneElement(vector<int> InputSequence, int index, int type, int TYPE) {
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