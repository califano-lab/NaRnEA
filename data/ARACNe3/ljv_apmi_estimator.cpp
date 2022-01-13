// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/chi_squared.hpp>
using namespace Rcpp;
using namespace std;

// takes in a vector of type T and a vector if indices (ints)
// returns sub-vector w/ specified indices
// [[Rcpp::export]]
vector<double> getElements(vector<double> val_vec, vector<int> ind_vec) {
	vector<double> ret_vec;
	for (size_t index = 0; index < ind_vec.size(); ++index) {
		ret_vec.push_back(val_vec[ind_vec[index]]);
	}
	return ret_vec;
}

// takes in a vector of type T
// returns a vector of ints w/ indices of elements less than threshold
// [[Rcpp::export]]
vector<int> lowerInds(vector<double> samp_vec, double thresh_val) {
	vector<int> low_inds;
	for (size_t index = 0; index < samp_vec.size(); ++index) {
		if (samp_vec[index] < thresh_val) {
			low_inds.push_back(index);
		}
	}
	return low_inds;
}

// takes in two vectors of doubles - one for X and Y - and an additional vector of doubles for the range of X
// returns a 3D vector; one 2D vector for each quadrant w/ 3 elements; X, Y, x.bounds
// [[Rcpp::export]]
vector<vector<vector<double>>> vectorSplit(vector<vector<double>> quadrant_obj) {
	// make sub-lists
	vector<double> vec_x, vec_y, bounds_x, bounds_y;
	vec_x = quadrant_obj[0]; vec_y = quadrant_obj[1]; bounds_x = quadrant_obj[2]; bounds_y = quadrant_obj[3];
	// get new thresholds
	double thresh_x = (quadrant_obj[2][1] + quadrant_obj[2][0]) / 2;
	double thresh_y = (quadrant_obj[3][1] + quadrant_obj[3][0]) / 2;
	// get lower indices
	vector<int> lInds_x = lowerInds(quadrant_obj[0], thresh_x);
	vector<int> lInds_y = lowerInds(quadrant_obj[1], thresh_y);
	// vector of all indices
	vector<int> all_inds(quadrant_obj[0].size());
	iota(begin(all_inds), end(all_inds), 0);
	// union vector
	vector<int> lInd_union;
	set_union(lInds_x.begin(), lInds_x.end(), lInds_y.begin(), lInds_y.end(),
		back_inserter(lInd_union));

	// 1st quadrant (top right)
	vector<int> q1_inds;
	set_difference(all_inds.begin(), all_inds.end(), lInd_union.begin(), lInd_union.end(),
		back_inserter(q1_inds));
	vector<vector<double>> q1;
	q1.push_back(getElements(vec_x, q1_inds));
	q1.push_back(getElements(vec_y, q1_inds));
	vector<double> q1_x{ thresh_x, bounds_x[1] };
	vector<double> q1_y{ thresh_y, bounds_y[1] };
	q1.push_back(q1_x); q1.push_back(q1_y);

	// 2nd quadrant (top right)
	vector<int> q2_inds;
	set_difference(lInds_x.begin(), lInds_x.end(), lInds_y.begin(), lInds_y.end(),
		back_inserter(q2_inds));
	vector<vector<double>> q2;
	q2.push_back(getElements(vec_x, q2_inds));
	q2.push_back(getElements(vec_y, q2_inds));
	vector<double> q2_x{ bounds_x[0], thresh_x };
	vector<double> q2_y{ thresh_y, bounds_y[1] };
	q2.push_back(q2_x); q2.push_back(q2_y);

	// 3rd quadrant (bot left)
	vector<int> q3_inds;
	set_intersection(lInds_x.begin(), lInds_x.end(), lInds_y.begin(), lInds_y.end(),
		back_inserter(q3_inds));
	vector<vector<double>> q3;
	q3.push_back(getElements(vec_x, q3_inds));
	q3.push_back(getElements(vec_y, q3_inds));
	vector<double> q3_x{ bounds_x[0], thresh_x };
	vector<double> q3_y{ bounds_y[0], thresh_y };
	q3.push_back(q3_x); q3.push_back(q3_y);

	// 4th quadrant (bot right)
	vector<int> q4_inds;
	set_difference(lInds_y.begin(), lInds_y.end(), lInds_x.begin(), lInds_x.end(),
		back_inserter(q4_inds));
	vector<vector<double>> q4;
	q4.push_back(getElements(vec_x, q4_inds));
	q4.push_back(getElements(vec_y, q4_inds));
	vector<double> q4_x{ thresh_x, bounds_x[1] };
	vector<double> q4_y{ bounds_y[0], thresh_y };
	q4.push_back(q4_x); q4.push_back(q4_y);

	// build final vector
	vector<vector<vector<double>>> split_vector;
	split_vector.push_back(q1);
	split_vector.push_back(q2);
	split_vector.push_back(q3);
	split_vector.push_back(q4);
	return(split_vector);
}

// takes in a 3D vector from vectorSplit and a quantile threshol 
// returns a boolean true/false for quadrant independence
// [[Rcpp::export]]
bool independenceTest(vector<vector<vector<double>>> vec_split, double q_thresh) {
	// count the number of samples in each quadrant
	vector<double> count_vec;
	for (size_t index = 0; index < vec_split.size(); ++index) {
		count_vec.push_back(vec_split[index][0].size());
	}
	double count_mean = accumulate(count_vec.begin(), count_vec.end(), 0.0) / count_vec.size();

	// calculate q statistic and test
	double chi_q;
	vector<double> count_dif(count_vec.size());
	for (size_t index = 0; index < count_vec.size(); ++index) {
		count_dif[index] = pow(count_vec[index] - count_mean, 2);
	}
	chi_q = accumulate(count_dif.begin(), count_dif.end(), 0.0) / count_mean;
	return chi_q < q_thresh;
}

// takes in the range of the subsample (vector), the total range (vector), the number of sub samples (int), total number of samples (int)
// returns the point-wise MI
// [[Rcpp::export]]
double calcMI(int num_samps, int tot_samps, double sub_range) {
	// calculate
	double sample_fraction = double(num_samps) / tot_samps;
	double mi = sample_fraction * log(sample_fraction / pow(sub_range, 2));
	// return if finite
	if (isfinite(mi)) {
		return mi;
	}
	else {
		return 0.0;
	}
}

// iterator function for recursive AP
// [[Rcpp::export]]
vector<double> apmiIterator(vector<vector<double>> quadrant_obj, double q_thresh, int size_thresh, int tot_samps) {	
	vector<double> mi_vals;
	int num_samps = quadrant_obj[0].size();

	// return if not enough samples in the partition
	if (num_samps < size_thresh) {
		mi_vals.push_back(calcMI(num_samps, tot_samps, quadrant_obj[2][1] - quadrant_obj[2][0]));
		return mi_vals;
	}

	// make split
	vector<vector<vector<double>>> vec_split = vectorSplit(quadrant_obj);
	// test for independence
	if (independenceTest(vec_split, q_thresh)) {
		mi_vals.push_back(calcMI(num_samps, tot_samps, quadrant_obj[2][1] - quadrant_obj[2][0]));
		return mi_vals;
	}

	// iterate for each quadrant
	vector<double> mi_vec;
	for (size_t index = 0; index < vec_split.size(); ++index) {
		mi_vec = apmiIterator(vec_split[index], q_thresh, size_thresh, tot_samps);
		mi_vals.insert(mi_vals.end(), mi_vec.begin(), mi_vec.end());
	}
	return mi_vals;
}

// WRAPPER for external use
// takes in two vectors, p-value and size thresholds
// returns vector of MI values for each region
// [[Rcpp::export]]
vector<double> APMI(vector<double> vec_x, vector<double> vec_y, double q_thresh = 7.814728, int size_thresh = 4) {
	int tot_samps = vec_x.size();
	// make initial quad_obj
	vector<vector<double>> quad_obj;
	quad_obj.push_back(vec_x); quad_obj.push_back(vec_y);
	vector<double> bounds{ 0.0, 1.0 };
	quad_obj.push_back(bounds); quad_obj.push_back(bounds);
	// make initial split
	vector<vector<vector<double>>> vec_split = vectorSplit(quad_obj);

	// iterate through initial four quadrants
	vector<double> mi_vals;
	vector<double> mi_vec;
	for (size_t index = 0; index < vec_split.size(); ++index) {
		mi_vec = apmiIterator(vec_split[index], q_thresh, size_thresh, tot_samps);
		mi_vals.insert(mi_vals.end(), mi_vec.begin(), mi_vec.end());
	}
	return mi_vals;
}