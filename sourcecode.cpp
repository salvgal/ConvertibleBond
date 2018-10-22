#include<iostream>
#include<algorithm>
#include<fstream>
#include<vector>
#include <iomanip>
#include<cmath> 
using namespace std;

class Database {
public:

	Database() { r = 0; q = 0; }
	~Database() {}

	double r;
	double q;
};
class TreeNode {
public:

	TreeNode() { S = 0; V = 0; t = 0; }
	~TreeNode() {}

	double S;
	double V;
	double t;
};
class Derivative {
public:
	virtual ~Derivative() {}

	virtual double TerminalPayoff(double S) const { return 0; }
	virtual int ValuationTests(TreeNode &node) const { return 0; }

	double T;
protected:
	Derivative() { T = 0; }
};

class BinomialTree {
public:
	BinomialTree(int n);
	~BinomialTree();

	int FairValue(int n, const Derivative * p_derivative,
		const Database * p_db, double S, double sigma,
		double t0, double &FV);
	int ImpliedVolatility(int n, const Derivative * p_derivative, const Database * p_db,
		double S, double t0, double target, double & implied_vol, int & num_iter);
private:

	void Clear();
	int Allocate(int n);


	int n_tree;
	TreeNode **tree_nodes;
};

BinomialTree::BinomialTree(int n) {
	n_tree = 0;
	tree_nodes = NULL;
	Allocate(n);
}
BinomialTree :: ~BinomialTree() {
	Clear();
}
void BinomialTree::Clear() {
	if (tree_nodes == NULL) return;
	for (int i = 0; i < n_tree; i++) {
		delete[] tree_nodes[i];
	}
	delete[] tree_nodes;
	tree_nodes = NULL;
}

int BinomialTree::Allocate(int n) {
	if (n < n_tree) return 0;

	Clear();

	n_tree = n;
	tree_nodes = new TreeNode*[n + 1];
	for (int i = 0; i <= n_tree; i++) {
		tree_nodes[i] = new TreeNode[n + 1];
		for (int j = 0; j <= n_tree; j++) {
			TreeNode temp_node;
			tree_nodes[i][j] = temp_node;
		}
	}
	return 0;
}
int BinomialTree::FairValue(int n, const Derivative * p_derivative,
	const Database * p_db, double S, double sigma,
	double t0, double &FV) {
	FV = 0;

	if (n < 1 ||
		S <= 0 ||
		p_derivative == NULL ||
		p_db == NULL ||
		p_derivative->T <= t0 ||
		sigma <= 0) return 1;

	double dt = (p_derivative->T - t0) / double(n);
	double df = exp(-1 * p_db->r * dt);
	double growth = exp((p_db->r - p_db->q) * dt);
	double u = exp(sigma * sqrt(dt));
	double d = 1.0 / u;

	double p_prob = (growth - d) / (u - d);
	double q_prob = 1.0 - p_prob;

	if (p_prob < 0.0 || p_prob > 1.0) {
		return 1;
	}


	Allocate(n);

	TreeNode * node_tmp = tree_nodes[0];
	node_tmp[0].S = S;
	node_tmp[0].t = t0;

	for (int i = 1; i <= n; i++) {
		double t = t0 + i * dt;
		TreeNode * prev = tree_nodes[i - 1];
		node_tmp = tree_nodes[i];
		node_tmp[0].S = prev[0].S * d;
		node_tmp[0].t = t;
		for (int j = 1; j <= n; j++) {
			node_tmp[j].S = node_tmp[j - 1].S * u * u;
			node_tmp[j].t = t;
		}
	}


	int i = n;
	node_tmp = tree_nodes[i];
	for (int j = 0; j <= n; j++) {
		node_tmp[j].V = p_derivative->TerminalPayoff(node_tmp[j].S);
	}

	for (int i = n - 1; i >= 0; i--) {
		node_tmp = tree_nodes[i];
		TreeNode * node_next = tree_nodes[i + 1];
		for (int j = 0; j <= i; j++) {
			node_tmp[j].V = df * (p_prob * node_next[j + 1].V + q_prob * node_next[j].V);

			p_derivative->ValuationTests(node_tmp[j]);
		}
	}

	node_tmp = tree_nodes[0];
	FV = node_tmp[0].V;

	return 0;
}

int BinomialTree::ImpliedVolatility(int n,
	const Derivative * p_derivative,
	const Database * p_db,
	double S, double t0, double target,
	double & implied_vol, int & num_iter) {

	const double tol = 1.0e-6;
	const int max_iter = 100;


	implied_vol = 0;
	num_iter = 0;

	double sigma_low = 0.01;
	double sigma_high = 2.0;
	double sigma = 0.0;
	double FV_low = 0.0;
	double FV_high = 0.0;
	double FV = 0.0;

	sigma = sigma_low;
	FairValue(n, p_derivative, p_db, S, sigma_low, t0, FV_low);
	double diff_FV_low = FV_low - target;
	if (abs(diff_FV_low) <= tol) {
		implied_vol = sigma;
		return 0;
	}
	sigma = sigma_high;
	FairValue(n, p_derivative, p_db, S, sigma_high, t0, FV_high);
	double diff_FV_high = FV_high - target;
	if (abs(diff_FV_high) <= tol) {
		implied_vol = sigma;
		return 0;
	}

	if (diff_FV_low * diff_FV_high > 0) {
		implied_vol = 0;
		return 1;
	}

	for (num_iter = 1; num_iter <max_iter; num_iter++) {
		sigma = 0.5 * (sigma_low + sigma_high);
		FairValue(n, p_derivative, p_db, S, sigma, t0, FV);
		double diff_FV = FV - target;
		if (abs(diff_FV) < tol) {
			implied_vol = sigma;
			return 0;
		}
		else {
			if ((diff_FV_low * diff_FV) > 0) {
				sigma_low = sigma;
			}
			else {
				sigma_high = sigma;
			}
		}

		if (abs(sigma_high - sigma_low) < tol) {
			implied_vol = sigma;
			return 0;
		}
	}
	implied_vol = 0;
	num_iter = max_iter;
	return 1;
}
class ConvertibleBond : public Derivative {
public:
	ConvertibleBond() { K = 0; isCall = false; isAmerican = false; B = 0; }
	virtual ~ConvertibleBond() {}
	virtual double TerminalPayoff(double S) const;
	virtual int ValuationTests(TreeNode &node) const;

	double K;
	double B; 
	bool isCall;
	bool isAmerican;
};
double ConvertibleBond::TerminalPayoff(double S) const {
	if (isCall && S >= K) return S;
	if (isCall && S < K) return K;
	return 0;
}

int ConvertibleBond::ValuationTests(TreeNode &node) const {
	if (node.S >= B) {
		node.V = node.S;
		return 0;
	}
	if (isAmerican) {
		double intrinsic = 0;
		if (isCall) {
			if (node.S >= K) intrinsic = node.S;
		}
		else {
			intrinsic = 0;
		}
		if (intrinsic > node.V) {
			node.V = intrinsic;
		}
	}
	return 0;
}

int test_plot_data() {
	ofstream outfile_conver_bond("convertible_bond_plot_data.txt");
	int rc = 0;
	double U = 0;
	double M = 0;

	double r = 0.05;
	double q = 0;
	double T = 5.0;
	double t0 = 0;

	Database db;
	db.r = r;
	db.q = q;

	double S = 0;
	double K = 100;
	double sigma = 0.5;

	double B = 130;
	int n = 1000;


	ConvertibleBond callableAmericanConvertibleBond;
	callableAmericanConvertibleBond.K = K;
	callableAmericanConvertibleBond.T = T;
	callableAmericanConvertibleBond.B = B;
	callableAmericanConvertibleBond.isCall = true;
	callableAmericanConvertibleBond.isAmerican = true;

	double FV_convertible_bond = 0;
	BinomialTree binom(n);
	cout << endl << endl;
	cout << "sigma = 0.5" << endl << endl;
	for (int i = 1; i <= 150; i++) {
		S = i * 1.00;
		rc = binom.FairValue(n, &callableAmericanConvertibleBond, &db, S, sigma, t0, FV_convertible_bond);
		outfile_conver_bond << std::setw(6) << S << " ";
		outfile_conver_bond << std::setw(16) << FV_convertible_bond << " " << endl;
	}
	outfile_conver_bond.close();
	return 0;
}

int test_gamma_trading() {
	cout << endl << endl;
	cout << "Gamma Trading:" << endl << endl;
	int rc = 0;
	double U = 0;
	double M = 0;

	double r = 0.05;
	double q = 0;
	double T = 5.0;
	double t0 = 0;

	Database db;
	db.r = r;
	db.q = q;

	double S = 0;
	double K = 100;
	double sigma = 0.5;

	double B = 130;
	int n = 1000;


	double delta_S1 = 2358.0 / 10000.0;  
	double delta_S2 = -7201.0 / 10000.0;

	//Day 0
	t0 = 0;
	double S0 = 60;
	double M0 = 90; 
	double target = M0;
	ConvertibleBond convertibleBondDay0;
	convertibleBondDay0.K = K;
	convertibleBondDay0.T = T;
	convertibleBondDay0.B = B;
	convertibleBondDay0.isCall = true;
	convertibleBondDay0.isAmerican = true;

	double implied_volatility = 0;
	int num_iter = 0;
	BinomialTree binomDay0(n);

	rc = binomDay0.ImpliedVolatility(n, &convertibleBondDay0, &db, S0, t0, target, implied_volatility, num_iter);
	cout << "implied volatility (Sigma_0) = " << setprecision(4) << implied_volatility << endl;

	double U_S0Plus1 = 0;
	double U_S0Minus1 = 0;
	binomDay0.FairValue(n, &convertibleBondDay0, &db, S0 + 1, implied_volatility, t0, U_S0Plus1);
	binomDay0.FairValue(n, &convertibleBondDay0, &db, S0 - 1, implied_volatility, t0, U_S0Minus1);
	double delta0_con_bond = (U_S0Plus1 - U_S0Minus1) / 2.0;
	cout << "delta0_con_bond =  " << delta0_con_bond << endl;

	double money0 = delta0_con_bond * S0 - M0;
	cout << "Money0 on day 0 is " << setprecision(4) << money0 << endl << endl; 

	//Day1
	t0 = 0.01;
	double M1 = 90.2;
	double S1 = S0 + delta_S1;
	target = M1;
	ConvertibleBond convertibleBondDay1;
	convertibleBondDay1.K = K;
	convertibleBondDay1.T = T;
	convertibleBondDay1.B = B;
	convertibleBondDay1.isCall = true;
	convertibleBondDay1.isAmerican = true;

	implied_volatility = 0;
	num_iter = 0;
	BinomialTree binomDay1(n);
	rc = binomDay1.ImpliedVolatility(n, &convertibleBondDay1, &db, S1, t0, target, implied_volatility, num_iter);
	cout << "implied volatility (Sigma_1) = " << setprecision(4) << implied_volatility << endl;

	double U_S1Plus1 = 0;
	double U_S1Minus1 = 0;
	binomDay1.FairValue(n, &convertibleBondDay1, &db, S1 + 1, implied_volatility, t0, U_S1Plus1);
	binomDay1.FairValue(n, &convertibleBondDay1, &db, S1 - 1, implied_volatility, t0, U_S1Minus1);
	double delta1_con_bond = (U_S1Plus1 - U_S1Minus1) / 2.0;
	cout << "delta1_con_bond =  " << delta1_con_bond << endl;

	double money1 = money0 + (delta1_con_bond - delta0_con_bond) * S1;
	cout << "Money1 on day 1 is " << setprecision(4) << money1 << endl << endl;;

	//Day2
	t0 = 0.02;
	double M2 = 90.15;
	double S2 = S1 + delta_S2;
	target = M2;
	ConvertibleBond convertibleBondDay2;
	convertibleBondDay2.K = K;
	convertibleBondDay2.T = T;
	convertibleBondDay2.B = B;
	convertibleBondDay2.isCall = true;
	convertibleBondDay2.isAmerican = true;

	implied_volatility = 0;
	num_iter = 0;
	BinomialTree binomDay2(n);

	rc = binomDay2.ImpliedVolatility(n, &convertibleBondDay2, &db, S2, t0, target, implied_volatility, num_iter);
	cout << "implied volatility (Sigma_2) = " << setprecision(4) << implied_volatility << endl;

	double U_S2Plus1 = 0;
	double U_S2Minus1 = 0;
	binomDay2.FairValue(n, &convertibleBondDay1, &db, S2 + 1, implied_volatility, t0, U_S2Plus1);
	binomDay2.FairValue(n, &convertibleBondDay1, &db, S2 - 1, implied_volatility, t0, U_S2Minus1);
	double delta2_con_bond = (U_S2Plus1 - U_S2Minus1) / 2.0;
	cout << "delta2_con_bond =  " << delta2_con_bond << endl;

	double money2 = money1 + (delta2_con_bond - delta1_con_bond)* S2;
	cout << "Money2 on day 2 is " << setprecision(4) << money2 << endl;

	double Profit = money2 + M2 - delta2_con_bond * S2;
	cout << endl << endl;
	cout << "Total profit is " << setprecision(4) << Profit << endl << endl;

	return 0;
}

int main() {
	test_plot_data();
	test_gamma_trading();
	system("PAUSE");
	return 0;
}
