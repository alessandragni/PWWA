#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// --- Helper Functions with Bounds Checking ---

// [[Rcpp::export]]
int which_max_safe(IntegerVector v) {
  if (v.size() == 0) return 0;
  int max_val = v[0];
  int idx = 0;
  for (int i = 1; i < v.size(); i++) {
    if (v[i] > max_val) {
      max_val = v[i];
      idx = i;
    }
  }
  return idx;
}

// [[Rcpp::export]]
NumericVector cumsum_cpp(NumericVector x) {
  int n = x.size();
  NumericVector res(n);
  double acc = 0;
  for (int i = 0; i < n; i++) {
    acc += x[i];
    res[i] = acc;
  }
  return res;
}

// Safer version of index lookup to avoid which_max issues
// [[Rcpp::export]]
double look_up_f_cpp(double tim, NumericVector stim, NumericVector fval) {
  int n = stim.size();
  if (n == 0 || fval.size() == 0) return 0.0;
  
  int best_idx = -1;
  for (int i = 0; i < n; i++) {
    if (stim[i] <= tim) {
      best_idx = i; // Assumes stim is sorted; finds the last (largest) index <= tim
    }
  }
  
  if (best_idx == -1) return 0.0;
  // Ensure we don't exceed fval bounds
  return fval[std::min(best_idx, (int)fval.size() - 1)];
}

// [[Rcpp::export]]
NumericVector apply_look_up_f_cpp(NumericVector tim, NumericVector stim, NumericVector fval) {
  int n = tim.size();
  NumericVector res(n);
  for (int i = 0; i < n; i++) {
    res[i] = look_up_f_cpp(tim[i], stim, fval);
  }
  return res;
}

// [[Rcpp::export]]
NumericMatrix vvmult(NumericVector V1, NumericVector V2) {
  int n1 = V2.size();
  int n2 = V1.size();
  NumericMatrix res(n1, n2);
  for (int i = 0; i < n1; i++) {
    res(i, _) = V1 * V2[i];
  }
  return res;
}

// [[Rcpp::export]]
int max_lesseq_cpp(NumericVector tau, double r) {
  int k = tau.size();
  int l1 = 0;
  bool found = false;
  for (int i = 0; i < k; ++i) {
    if (tau[i] <= r) {
      l1 = i;
      found = true;
    } else {
      break; 
    }
  }
  return found ? l1 : 0;
}

// [[Rcpp::export]]
int min_greq_cpp(NumericVector tau, double r) {
  int k = tau.size();
  for (int i = 0; i < k; ++i) {
    if (tau[i] >= r) return i;
  }
  return 0; 
}

// [[Rcpp::export]]
NumericVector pmin_cpp(NumericVector vec, double value) {
  int n = vec.size();
  NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = std::min(vec[i], value);
  }
  return result;
}

// [[Rcpp::export]]
NumericVector cumsum_right_cpp(NumericVector x) {
  int n = x.size();
  if (n == 0) return x;
  NumericVector result(n);
  result[n - 1] = x[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    result[i] = result[i + 1] + x[i];
  }
  return result;
}

// [[Rcpp::export]]
NumericVector subset_by_indices(NumericVector x, int low, int high) {
  int n = x.size();
  if (n == 0 || low > high || low >= n) return NumericVector(0);
  int safe_high = std::min(high, n - 1);
  int safe_low = std::max(0, low);
  Rcpp::Range range(safe_low, safe_high);
  return x[range];
}

// --- Main Computational Functions ---

// [[Rcpp::export]] 
NumericVector Heff_cpp(double t, double k,
                       NumericMatrix Lam1_tau1, NumericMatrix Lam2_tau1, 
                       NumericVector tau1, NumericMatrix dLam1,
                       NumericVector tau12, NumericMatrix dLam12, 
                       NumericMatrix Lam12_tau1, NumericMatrix Lam12_tau12) {
  int n = dLam1.ncol();
  int k12 = tau12.size();
  int indt = max_lesseq_cpp(tau1, t);
  
  NumericMatrix H(indt + 1, n);
  NumericMatrix Kt2mint_tau12(k12, n);
  NumericMatrix Kt2mint_tau1(indt + 1, n);
  
  for(int i = 0; i < n; i++) {
    NumericVector inner = pow(pmin_cpp(tau12, t), -1./k) * exp(-Lam12_tau12(_, i)) * dLam12(_, i);
    NumericVector csr = cumsum_right_cpp(inner);
    
    // Safety check for index 1
    if (csr.size() > 1) csr[0] = csr[1];
    Kt2mint_tau12(_, i) = csr;
    
    Kt2mint_tau1(_, i) = apply_look_up_f_cpp(tau1[Rcpp::Range(0, indt)], tau12, Kt2mint_tau12(_, i));
    H(_, i) = cumsum_cpp( Kt2mint_tau1(_, i) * exp(Lam12_tau1(_, i)) * exp(-Lam1_tau1(_, i) - Lam2_tau1(_, i)) * dLam1(_, i) );
  }
  
  return H(indt, _);
} 

// [[Rcpp::export]] 
NumericVector CM1_cpp(double t, double k, NumericVector Ttilde2, NumericVector Ttilde1, 
                      NumericVector delta1, 
                      NumericVector tau12, NumericMatrix dLam12,
                      NumericVector tau1, NumericMatrix dLam1,
                      NumericMatrix Lam1_tau1, 
                      NumericMatrix Lam2_tau1, 
                      NumericMatrix Lam1_tauc, 
                      NumericMatrix Lam2_tauc, NumericMatrix Lam12_tau1, NumericMatrix Lam12_tau12,
                      NumericVector tauc, NumericMatrix KC_tauc, NumericMatrix dLamC) {
  
  int n = dLamC.ncol();
  int kc = tauc.size();
  int max_integral1_tau2 = max_lesseq_cpp(tau12, t);
  int max_integral1_tau1 = max_lesseq_cpp(tau1, t);
  int min_integral2_tau2 = min_greq_cpp(tau12, t);
  
  NumericVector CM1(n);
  
  for(int i = 0; i < n; i++) {
    // 1. Martingale calculation
    NumericVector dMc(kc);
    for(int j = 0; j < kc; j++) {
      double dNc = (Ttilde2[i] == tauc[j]) ? 1.0 : 0.0;
      double Y = (tauc[j] <= Ttilde2[i]) ? 1.0 : 0.0;
      dMc[j] = dNc - Y * dLamC(j, i);
    }
    
    // 2. Inner integrals
    NumericVector Kt2_v = subset_by_indices(pow(tau12, -1./k) * exp(-Lam12_tau12(_, i)) * dLam12(_, i), 0, max_integral1_tau2);
    NumericVector Kt2_tau12 = cumsum_right_cpp(Kt2_v);
    if (Kt2_tau12.size() > 1) Kt2_tau12[0] = Kt2_tau12[1];
    
    NumericVector Kt_in = exp(-Lam12_tau12(_, i)) * dLam12(_, i);
    // Set zero before min_integral2_tau2
    for(int j = 0; j < std::min(min_integral2_tau2, (int)Kt_in.size()); j++) Kt_in[j] = 0;
    NumericVector Kt_tau12 = cumsum_right_cpp(Kt_in);
    
    NumericVector Kt2_tau1 = apply_look_up_f_cpp(tau1, tau12, Kt2_tau12);
    NumericVector Kt_tau1 = apply_look_up_f_cpp(tau1, tau12, Kt_tau12);
    
    NumericVector H1_v = subset_by_indices(Kt2_tau1 * exp(Lam12_tau1(_, i)) * exp(-Lam1_tau1(_, i) - Lam2_tau1(_, i)) * dLam1(_, i), 0, max_integral1_tau1);
    NumericVector H1_tau1 = cumsum_right_cpp(H1_v);
    if (H1_tau1.size() > 1) H1_tau1[0] = H1_tau1[1];
    
    NumericVector H2_tau1 = pow(t, -1./k) * cumsum_right_cpp(Kt_tau1 * exp(Lam12_tau1(_, i)) * exp(-Lam1_tau1(_, i) - Lam2_tau1(_, i)) * dLam1(_, i));
    if (H2_tau1.size() > 1) H2_tau1[0] = H2_tau1[1];
    
    NumericVector H1_tauc = apply_look_up_f_cpp(tauc, tau1, H1_tau1);
    NumericVector H2_tauc = apply_look_up_f_cpp(tauc, tau1, H2_tau1);
    
    int ind = max_lesseq_cpp(tauc, std::min(Ttilde1[i], t));
    
    // Final accumulation (replaces inefficient cumsum_cpp(...)(ind))
    NumericVector final_integrand = exp(Lam1_tauc(_, i) + Lam2_tauc(_, i)) * (H1_tauc + H2_tauc) * (1.0 / KC_tauc(_, i)) * dMc;
    double res = 0;
    for(int j = 0; j <= ind && j < final_integrand.size(); j++) res += final_integrand[j];
    CM1[i] = res;
  }
  return CM1; 
} 

// [[Rcpp::export]] 
NumericVector CM2_cpp(double t, double k, NumericVector Ttilde2, NumericVector Ttilde1, 
                      NumericVector delta1, 
                      NumericVector tau12, NumericMatrix dLam12, NumericMatrix Lam12_tauc, NumericMatrix Lam12_tau12,
                      NumericVector tauc, NumericMatrix KC_tauc, NumericMatrix dLamC) {
  
  int n = dLamC.ncol();
  int kc = tauc.size();
  int min_integral2 = min_greq_cpp(tauc, t);
  NumericVector CM2(n);
  
  for(int i = 0; i < n; i++) {
    if ((delta1[i] == 1) && (Ttilde1[i] <= t)) {
      NumericVector dMc(kc);
      for(int j = 0; j < kc; j++) {
        dMc[j] = (Ttilde2[i] == tauc[j] ? 1.0 : 0.0) - (tauc[j] <= Ttilde2[i] ? 1.0 : 0.0) * dLamC(j, i);
      }
      
      NumericVector Kt2_csr = cumsum_right_cpp(pow(pmin_cpp(tau12, t), -1./k) * exp(-Lam12_tau12(_, i)) * dLam12(_, i));
      if (Kt2_csr.size() > 1) Kt2_csr[0] = Kt2_csr[1];
      
      NumericVector Kt_csr = pow(t, -1./k) * cumsum_right_cpp(exp(-Lam12_tau12(_, i)) * dLam12(_, i));
      if (Kt_csr.size() > 1) Kt_csr[0] = Kt_csr[1];
      
      NumericVector Kt2mint_tauc = apply_look_up_f_cpp(tauc, tau12, Kt2_csr) * exp(Lam12_tauc(_, i));
      NumericVector Kt_tauc = apply_look_up_f_cpp(tauc, tau12, Kt_csr) * exp(Lam12_tauc(_, i));
      
      int min1 = min_greq_cpp(tauc, Ttilde1[i]);
      int max1 = max_lesseq_cpp(tauc, std::min(Ttilde2[i], t));
      
      double sum1 = 0;
      if (min1 <= max1) {
        NumericVector sub = subset_by_indices(Kt2mint_tauc * (1.0 / KC_tauc(_, i)) * dMc, min1, max1);
        for(int j=0; j<sub.size(); j++) sum1 += sub[j];
      }
      
      double sum2 = 0;
      if (t < Ttilde2[i]) {
        int max2 = min_greq_cpp(tauc, Ttilde2[i]);
        if (min_integral2 <= max2) {
          NumericVector sub = subset_by_indices(Kt_tauc * (1.0 / KC_tauc(_, i)) * dMc, min_integral2, max2);
          for(int j=0; j<sub.size(); j++) sum2 += sub[j];
        }
      }
      CM2[i] = sum1 + sum2;
    } else {
      CM2[i] = 0;
    }
  }
  return CM2; 
}