#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
// Count number of matching sequences between pairs
int count_sequential_matches(std::string str1, std::string str2) {
  int count = 0;
  int minLength = std::min(str1.length(), str2.length());
  for (int i = 0; i < minLength; i++) {
	if (str1[i] == '0' || str2[i] == '0'){
		break;
	} 
    else if (str1[i] == str2[i]) {
      count++;
    } else {
      break;
    }
  }
  
  return count;
}

 
// [[Rcpp::export]]
IntegerMatrix get_distance_sequential(ListOf<StringVector> strList) {
	// return distance matrix 
	// input: list of n mx1 barcodes
	// return: mxm distance matrix (sum pairwise dists over all barcodes)
  	int m = strList[0].size(); // num tips
  	int n = strList.size(); // num barcodes
	IntegerMatrix result(m, m);
	// iterate over pairs
  	for (int i = 0; i < m; i++) {
    	for (int j = i; j < m; j++) {
      		int count = 0;
      		for (int k = 0; k < n; k++) {
        		count += count_sequential_matches(as<std::string>(strList[k][i]), as<std::string>(strList[k][j]));
      		}
			result(i,j) = count;
			result(j,i) = count;
		}
	}
	return result;
}

 
int find_outgrp(int i1, int i2, int i3, NumericMatrix D) {
  int id_out = -1;

  // Check D(a, b) < min(D(a, c), D(b, c))
  if (D(i1, i2) < std::min(D(i1, i3), D(i2, i3))) {
    id_out = i3;
  } else if (D(i1, i3) < std::min(D(i1,i2), D(i3,i2))) {
	id_out = i2;
  } else if (D(i2,i3) < std::min(D(i2,i1), D(i3,i1))){
	id_out = i1;
  } else {
	id_out = -1;
  }
  
  return id_out;

}

// Get prop triplets which agree
// [[Rcpp::export]]
double proportion_triplets(NumericMatrix D_true, NumericMatrix D) {
  int n = D_true.nrow();
  int true_count = 0;  
  for (int i1 = 0; i1 < n - 2; ++i1) {
    for (int i2 = i1 + 1; i2 < n - 1; ++i2) {
      for (int i3 = i2 + 1; i3 < n; ++i3) {
        // Find outgroup out of i1,i2,i3 based on D_true
        int outgrp = find_outgrp(i1,i2,i3,D_true);
		int outgrp_test = find_outgrp(i1,i2,i3,D);
		// check
		if (outgrp_test == -1 | outgrp != outgrp_test) {
		  true_count++; 
        }
      }
    }
  }  
  return double(true_count); 
}



//[[Rcpp::export]]
// Remove nonsequential identical edits
// Return updated barcodes
CharacterVector remove_homo_edits(
	CharacterVector bars){
	CharacterVector newbars = clone(bars); //copy	
	int m = bars.size(); // num samples
	int k = bars[0].size(); // length bar
	for (int i = 0; i < m; ++i) {
		for (int j = i+1; j < m; ++j) {
			std::string b1 = Rcpp::as<std::string>(bars[i]);
		    std::string b2 = Rcpp::as<std::string>(bars[j]);
			int s = count_sequential_matches(b1,b2);
	     // change additional matching characters
      	for (int p = s; p < k; ++p) {
        	if (b1[p] == b2[p]) {
          		b1[p] = b2[p] = '*';
        	}
      	}
      newbars[i] = b1;
      newbars[j] = b2;
    }
  }
  return newbars;
}


// [[Rcpp::export]]
// censored Poisson 
// max_n > 0
double trunc_poi(int n, int max_n, double rate) {
  double p_x;
  if (n > max_n) {
    p_x = 0.0;
  } else if (n == max_n) {
    double sum = 0.0;
    for (int i = 0; i < max_n; ++i) {
      sum += R::dpois(i, rate,false);
    }
    p_x = 1.0 - sum;
  } else {
    p_x = R::dpois(n, rate,false);
  }
  return p_x;
}


//[[Rcpp::export]]
// Probability Homoplasy
// Conditioned on number of remaining sites k
double corr_prob_homo(int x, double lambda, int k, 
	double d, double q) {
	double p = 0.0;
	// i = x
	p += trunc_poi(x,k, lambda*(1.-d)) * 
		trunc_poi(x,k,lambda*(1.-d)) * 
		pow(q,x);
	for (int j = x+1; j <= k; ++j){
		p += 2.0 * trunc_poi(x,k, lambda*(1.-d)) *
			trunc_poi(j, k, lambda*(1.-d)) *
			pow(q,x);
	}
	// i > x
	if (q < 1.0) {
		for (int i = x+1; i <= k; ++i){
			p += pow(trunc_poi(i,k,lambda*(1.-d)),2) *
				pow(q,x) * (1.-q); 
			for (int j = i+1; j <= k; ++j){
				p += 2.0 * trunc_poi(i,k,lambda*(1.-d)) *
				trunc_poi(j, k, lambda*(1.-d)) *
				pow(q,x) * (1.-q);
			}
		}
	}
	return p;
}

// Probability x edits on branch of depth d and length ell
//[[Rcpp::export]]
double prob_on_branch(int x, int k, double lambda, double ell, double d){
   double p = 0.0;
   for (int k1=0; k1 <= k; ++k1){
        p += trunc_poi(x, k-k1, lambda*ell)*trunc_poi(k1, k, lambda*d);
   }
   return p;
}

// Probability x homo edits on branch of depth d (length 1-d)
//[[Rcpp::export]]
double prob_on_branch_homo(int x, double d, double lambda, int k, double q){
   double p = 0.0;
   for (int k1=0; k1 <= k; ++k1){
        p += corr_prob_homo(x,lambda,k-k1,d,q)*trunc_poi(k1, k, lambda*d);
   }
   return p;
}

//--------------------------//
// --- DIFFERENCE-TRIPLETS--//
//--------------------------//
//[[Rcpp::export]] 
// Prob x shared edits ingrp 
// given k sites left
double prob_x_in(int x, double du, double dv, double lambda, int k,
	double q){
	double p_in = 0.0;
	double p_x1 = 0.0;
	double p_x2 = 0.0;
	//x=x1+x2
	for (int x1 = 0; x1 <= x; ++x1){
		int x2 = x-x1;
		// edits from du to dv
		p_x1 = trunc_poi(x1, k, lambda*(dv-du));	
		// prob homoplasy in group
   		p_x2 = corr_prob_homo(x2,lambda, k-x1,dv,q);
		p_in += p_x1*p_x2;	
    }
	return p_in;
}


// du: depth of LCA(a,b,c)
// dv: depth of LCA of (a,b)
//[[Rcpp::export]]
double prob_seq_edits_in_grp(int x, double du, double dv, double lambda,int k, double q){
	double p_x = 0.0;
	double p_k1 = 0.0;
	double p_x2 = 0.0;
	double p_x1 = 0.0;
	for (int k1 = 0; k1 <= k; ++k1){	
		for (int x1 = 0; x1 <= x; ++x1){
			int x2 = x-x1;
			// prob x1 edits
			p_x1 = trunc_poi(x1, k-k1, lambda*(dv-du));	
		    // p k1 edits
			p_k1 = trunc_poi(k1, k, lambda*du);
			// prob homoplasy
        	p_x2 = corr_prob_homo(x2,lambda, k-x1-k1,dv,q);
        	p_x += p_x1*p_x2*p_k1;
		}
	}	
	return p_x;
}


//[[Rcpp::export]]
// d: depth, lambda: rate, k: max # sites
// q: collision prob
double expected_trips(double du, double dv, double lambda, int k,double q){ 
    double expect_in = 0.0;
    double expect_out = 0.0;
    for (int x = 0; x <= k; ++x){
        expect_in += x * prob_seq_edits_in_grp(x,du,dv,lambda,k,q);
        expect_out += x * prob_on_branch_homo(x,du,lambda,k,q);
    }
	return expect_in - expect_out;	
}



///*** DOES NOT WORK ***///


//[[Rcpp::export]]
// prob difference given k sites left
// d: depth, lambda: rate, k: max # sites
// q: collision prob
// x in [-k,k]
double p_trip_0(int x, double du, double dv, double lambda, int k, double q){ 
	double p_diff = 0.0;
   	for (int y = 0; y <= k; ++y){
		int y2 = y+x;
		if (y2 >= 0 && y2 <= k){
			// prob y2 shared edits in-grp
			double p_in = prob_x_in(y2,du,dv,lambda,k,q);
			// prob y shared edits out-grp
       		double p_out = corr_prob_homo(y,lambda,k,du,q);
			p_diff += p_in*p_out;	
		}
	} 
	return p_diff;	
}

//[[Rcpp::export]]
// Prob distribution of difference ingrp-outgrp
double p_trip00(int x, double du, double dv, double lambda,
	int k, double q){
	double p = 0.0;
	double p_k1 = 0.0;
	// process from root to du
	for (int k1 = 0; k1 <= k; ++k1){
		p_k1 = trunc_poi(k1, k, lambda*du);
		// they will share k1 edits
		p += p_trip_0(x, du, dv, lambda, k-k1,q)*p_k1;
	}	
	return p;
}
//[[Rcpp::export]]
// d: depth, lambda: rate, k: max # sites
// q: collision prob
// NOTE: same as other one, just slower
double expected_trips00(double du, double dv, double lambda, int k,
	double q){ 
    double expect_diff = 0.0;
    for (int x = -k; x <= k; ++x){
        expect_diff += x * p_trip00(x,du,dv,lambda,k,q); 
    }
	return expect_diff;	
}


