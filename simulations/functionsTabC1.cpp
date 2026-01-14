#include <Rcpp.h>

using namespace Rcpp;
using namespace sugar;

// [[Rcpp::export]]
int which_maxCpp(IntegerVector v) {
  int z = which_max(v);
  return z;
}



// [[Rcpp::export]]
NumericVector cumsum_cpp(NumericVector x){
  double acc = 0;
  NumericVector res(x.size());
  for(int i = 0; i < x.size(); i++){
    acc += x[i];
    res[i] = acc;
  }
  return res;
}  



// [[Rcpp::export]]
double look_up_f_cpp(double tim, NumericVector stim, NumericVector fval){
  int n = stim.size();
  double res = 0;
  if (min(stim)<=tim){
    IntegerVector x(n);
    for(int i = 0; i < n; i++){
      if(stim[i]<=tim){
        x[i] = i;
      }
      else{
        x[i] = 0;
      }
    }
    int index = which_maxCpp(x);
    res = fval[index];
  }
  return res;
}


// [[Rcpp::export]]
NumericVector apply_look_up_f_cpp(NumericVector tim, NumericVector stim, NumericVector fval){
  int n = tim.size();
  NumericVector res(n);
  for(int i = 0; i < n; i++){
    res[i] = look_up_f_cpp(tim[i], stim, fval);
  }
  return res;
}


// [[Rcpp::export]]
NumericMatrix vvmult(NumericVector V1, NumericVector V2){
  int n1 = V2.size();
  int n2 = V1.size();
  NumericMatrix res(n1,n2);
  for(int i = 0; i < n1; i++){
    res(i,_) = V1*V2(i);
  }
  return res;
}



// [[Rcpp::export]]
int max_lesseq_cpp(NumericVector tau, double r) {
  int k = tau.size();
  
  int l1 = k - 1; // Initialize to the last index to handle case where all values less or equal than r
  bool all_greater_r = true;
  
  // Check if all values are greater than r
  for (int i = 0; i < k; ++i) {
    if (tau[i] <= r) {
      l1 = i; // Update the index whenever tau[i] <= r
      all_greater_r = false;
    } else {
      break;
    }
  }
  
  // Return 0 if all values are greater than r
  if (all_greater_r)
    return 0;
  
  return l1;
}


// [[Rcpp::export]]
int min_greq_cpp(NumericVector tau, double r) {
  int k = tau.size();
  
  int l1 = 0; // Initialize to k-1 to handle case where all values are greater or equal to r
  bool all_less_r = true;
  
  // Check if all values are less than r
  for (int i = 0; i < k; ++i) {
    if (tau[i] >= r) {
      l1 = i; // Update the index whenever tau[i] >= r
      all_less_r = false;
      break; // Break the loop as we've found the first index where tau[i] >= r
    }
  }
  
  if(all_less_r)
    return 0;
  
  return l1;
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
  NumericVector result(n);
  
  result[n - 1] = x[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    result[i] = result[i + 1] + x[i];
  }
  return result;
}


// [[Rcpp::export]] 
NumericVector Heff_cpp(double t, double k,
                       NumericMatrix Lam1_tau1, NumericMatrix Lam2_tau1, 
                       NumericVector tau1, NumericMatrix dLam1,
                       NumericVector tau12, NumericMatrix dLam12, 
                       NumericMatrix Lam12_tau1, NumericMatrix Lam12_tau12){
  int n = dLam1.ncol();
  // int k1 = tau1.size();
  int k12 = tau12.size();
  
  int indt = max_lesseq_cpp(tau1, t);
  // Rcpp::Rcout << "The value of indt is: " << indt << std::endl;
  // Rcpp::Rcout << "The size of tau1 is: " << tau1.size() << std::endl;
  // Rcpp::Rcout << "The tau1[indt] is: " << tau1[indt] << std::endl;
  // Rcpp::Rcout << "The tau1[indt+1] is: " << tau1[indt+1] << std::endl;
  
  NumericMatrix H(indt + 1, n);
  
  // Torben's code
  // NumericMatrix Kt2mint_tau12(k12, n);
  // NumericMatrix d_inner_int(k12, n);
  // 
  // d_inner_int(k12-1, _) = pow(std::min(tau12(k12-1), t), -1./k) * exp(-Lam12_tau12(k12-1, _)) * dLam12(k12-1, _);
  // Kt2mint_tau12(k12-1, _) = d_inner_int(k12-1, _);
  // for (int j = k12 - 2; j >= 1; --j) {
  //     d_inner_int(j, _) = pow(std::min(tau12(j), t), -1./k) * exp(-Lam12_tau12(j, _)) * dLam12(j, _);
  //     Kt2mint_tau12(j, _) = Kt2mint_tau12(j + 1, _) + d_inner_int(j, _);
  //   }
  // Kt2mint_tau12(0, _) = Kt2mint_tau12(1, _);
  
  // More compact
  NumericMatrix Kt2mint_tau12(k12, n);
  NumericMatrix Kt2mint_tau1(indt + 1, n);
  
  for(int i = 0; i < n; i++){
    
    Kt2mint_tau12(_, i) = cumsum_right_cpp(pow(pmin_cpp(tau12, t), -1./k) * exp(-Lam12_tau12(_, i)) * dLam12(_, i));
    Kt2mint_tau12(0, i) = Kt2mint_tau12(1, i);
    // Rcpp::Rcout << "The value of Kt2mint_tau12(0, i) is:" << Kt2mint_tau12(0, i) << std::endl;
    
    // Change support
    Kt2mint_tau1(_, i) = apply_look_up_f_cpp(tau1[Rcpp::Range(0, indt)], tau12, Kt2mint_tau12(_, i));
    
    H(_, i) = cumsum_cpp( Kt2mint_tau1(_, i) * exp(Lam12_tau1(_, i)) * exp( - Lam1_tau1(_, i) - Lam2_tau1(_, i) ) * dLam1(_, i) );
  }
  
  return H(indt, _);
} 



// [[Rcpp::export]]
NumericVector subset_by_indices(NumericVector x, int low, int high) {
  Rcpp::Range range(low, high); // Create a Range object for the desired indices
  
  return x[range]; // Subset the vector using the Range object
}



// [[Rcpp::export]]
double sum_vector(NumericVector vec) {
  double sum = 0.0; // Initialize sum
  
  for (int i = 0; i < vec.size(); ++i) {
    sum += vec[i]; // Accumulate sum
  }
  return sum;
}


// [[Rcpp::export]]
NumericVector setZeroBefore(NumericVector vec, int j) {
  // Set elements before position j to zero
  std::fill(vec.begin(), vec.begin() + j, 0);
  return vec;
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
                      NumericVector tauc, NumericMatrix KC_tauc, NumericMatrix dLamC){
  
  int n = dLamC.ncol();
  int kc = tauc.size();
  int k12 = tau12.size();
  int k1 = tau1.size();
  int min_integral2_tau2 = min_greq_cpp(tau12, t);
  int max_integral1_tau2 = max_lesseq_cpp(tau12, t);
  //int min_integral2_tau1 = min_greq_cpp(tau1, t);
  int max_integral1_tau1 = max_lesseq_cpp(tau1, t);
  
  NumericVector CM1(n);
  
  // NumericMatrix Kt2_tau12(k12, n);
  // NumericMatrix Kt_tau12(k12, n);
  
  for(int i = 0; i < n; i++){
    // Rcpp::Rcout << "The value of i is: " << i << std::endl;
    
    // We compute the censoring martingale
    NumericVector dNc(kc);
    dNc = (Ttilde2(i) == tauc);
    
    NumericVector Y(kc);
    Y = (tauc <= Ttilde2(i));
    
    NumericVector dMc(kc);
    dMc = dNc - Y * dLamC(_, i);
    
    // We compute the two inner integrals in t2
    
    // For the first integral
    
    NumericVector Kt2_tau12(max_integral1_tau2 + 1);
    Kt2_tau12 = cumsum_right_cpp(subset_by_indices(pow(tau12, -1./k) * exp(-Lam12_tau12(_, i)) * dLam12(_, i), 
                                                   0, max_integral1_tau2));
    Kt2_tau12(0) = Kt2_tau12(1);                           
    
    
    // For the second integral we split it in two int_r^t int_t^tau dt2 dt1 and int_t^tau int_t1^tau dt2 dt1
    
    // t2 from t to tau
    // double Kt_tau12_part1;
    // Kt_tau12_part1 = sum_vector(subset_by_indices(exp(-Lam12_tau12(_, i)) * dLam12(_, i), 
    //                                               min_integral2_tau2, k12-1));
    // // t2 from t1 to tau
    // NumericVector Kt_tau12_part2(k12);
    // Kt_tau12_part2 = cumsum_right_cpp(exp(-Lam12_tau12(_, i)) * dLam12(_, i));
    // Kt_tau12_part2(0) = Kt_tau12_part2(1);
    NumericVector Kt_tau12(k12);
    Kt_tau12 = cumsum_right_cpp(setZeroBefore(exp(-Lam12_tau12(_, i)) * dLam12(_, i), 
                                              min_integral2_tau2));
    
    // We change the support from tau12 to tau1
    
    NumericVector Kt2_tau1(k1);
    Kt2_tau1 = apply_look_up_f_cpp(tau1, tau12, Kt2_tau12);
    
    //NumericVector Kt_tau1_part2(k1);
    //Kt_tau1_part2 = apply_look_up_f_cpp(tau1, tau12, Kt_tau12_part2);
    NumericVector Kt_tau1(k1);
    Kt_tau1 = apply_look_up_f_cpp(tau1, tau12, Kt_tau12);
    
    
    // We compute the two inner integrals in t1
    NumericVector H1_tau1(max_integral1_tau1 + 1);
    //NumericVector H2_tau1(max_integral1_tau1 + 1);
    //NumericVector H2_tau1_part1(max_integral1_tau1 + 1);
    NumericVector H2_tau1(k1);
    
    // For the first integral
    H1_tau1 = cumsum_right_cpp( subset_by_indices(Kt2_tau1 * exp(Lam12_tau1(_, i)) * exp( - Lam1_tau1(_, i) - Lam2_tau1(_, i) ) * dLam1(_, i), 
                                                  0, max_integral1_tau1));
    H1_tau1(0) = H1_tau1(1);
    
    // For the second integral
    //H2_tau1_part1 = cumsum_right_cpp( subset_by_indices(Kt_tau12_part1 * exp(Lam12_tau1(_, i)) * exp( - Lam1_tau1(_, i) - Lam2_tau1(_, i) ) * dLam1(_, i), 
    //                                                    0, max_integral1_tau1));
    //H2_tau1_part1(0) = H2_tau1_part1(1);
    //H2_tau1 = pow(t, -1./k) * (H2_tau1_part1 + 
    //          sum_vector(subset_by_indices(Kt_tau1_part2 * exp(Lam12_tau1(_, i)) * exp( - Lam1_tau1(_, i) - Lam2_tau1(_, i) ) * dLam1(_, i), 
    //                                       min_integral2_tau1, k1-1)));
    H2_tau1 = pow(t, -1./k) * cumsum_right_cpp(Kt_tau1 * exp(Lam12_tau1(_, i)) * exp( - Lam1_tau1(_, i) - Lam2_tau1(_, i) ) * dLam1(_, i));
    H2_tau1(0) = H2_tau1(1);
    
    // We change the support from tau1 to tauc for both terms
    NumericVector H1_tauc(kc);
    H1_tauc = apply_look_up_f_cpp(tauc, tau1, H1_tau1);
    NumericVector H2_tauc(kc);
    H2_tauc = apply_look_up_f_cpp(tauc, tau1, H2_tau1);
    
    // We compute the outer integral summing the two contributions
    // and extracting the correct value (ind)
    int ind = max_lesseq_cpp(tauc, std::min(Ttilde1(i), t));
    
    CM1(i) = cumsum_cpp(exp(Lam1_tauc(_, i) + Lam2_tauc(_, i)) * (H1_tauc + H2_tauc) * pow(KC_tauc(_, i), -1) * dMc)(ind);
    
  }
  return CM1; 
} 




// [[Rcpp::export]] 
NumericVector CM2_cpp(double t, double k, NumericVector Ttilde2, NumericVector Ttilde1, 
                      NumericVector delta1, 
                      NumericVector tau12, NumericMatrix dLam12, NumericMatrix Lam12_tauc, NumericMatrix Lam12_tau12,
                      NumericVector tauc, NumericMatrix KC_tauc, NumericMatrix dLamC){
  
  int n = dLamC.ncol();
  int kc = tauc.size();
  int k12 = tau12.size();
  
  int min_integral2 = min_greq_cpp(tauc, t);
  
  NumericVector CM2_1(n);
  NumericVector CM2_2(n);
  NumericVector CM2(n);
  
  for(int i = 0; i < n; i++){
    
    if ((delta1(i) == 1) && (Ttilde1(i) <= t)) {
      
      NumericVector dNc(kc);
      dNc = (Ttilde2(i) == tauc);
      
      NumericVector Y(kc);
      Y = (tauc <= Ttilde2(i));
      
      NumericVector dMc(kc);
      dMc = dNc - Y * dLamC(_, i);
      
      NumericVector Kt2mint_tau12(k12);
      NumericVector Kt_tau12(k12);
      
      // For the first integral
      
      Kt2mint_tau12 = cumsum_right_cpp(pow(pmin_cpp(tau12, t), -1./k) * exp(-Lam12_tau12(_, i)) * dLam12(_, i));
      Kt2mint_tau12(0) = Kt2mint_tau12(1);
      
      
      // For the second integral
      Kt_tau12 = pow(t, -1./k) * cumsum_right_cpp( exp(-Lam12_tau12(_, i)) * dLam12(_, i)); 
      Kt_tau12(0) = Kt_tau12(1);
      
      
      // We change the support
      NumericVector Kt2mint_tauc(kc);
      Kt2mint_tauc = apply_look_up_f_cpp(tauc, tau12, Kt2mint_tau12) * exp(Lam12_tauc(_, i));
      
      NumericVector Kt_tauc(kc);
      Kt_tauc = apply_look_up_f_cpp(tauc, tau12, Kt_tau12) * exp(Lam12_tauc(_, i));
      
      
      // We compute the first integral
      
      int min_integral1 = min_greq_cpp(tauc, Ttilde1(i));
      int max_integral1 = max_lesseq_cpp(tauc, std::min(Ttilde2(i), t));
      
      if(min_integral1 <= max_integral1){
        CM2_1(i) = sum_vector(subset_by_indices(Kt2mint_tauc * pow(KC_tauc(_, i), -1) * dMc,
                              min_integral1, max_integral1));
      } else {
        CM2_1(i) = 0;
      }
      
      
      // We compute the second integral
      
      if (t < Ttilde2(i)) {
        int max_integral2 = min_greq_cpp(tauc, Ttilde2(i));
        
        if(min_integral2 <= max_integral2){
          CM2_2(i) = sum_vector(subset_by_indices(Kt_tauc * pow(KC_tauc(_,i), -1) * dMc,
                                min_integral2, max_integral2));
        } else {
          CM2_2(i) = 0; 
        }
      } else {
        CM2_2(i) = 0;
      }
      
      // We sum the two integrals
      CM2(i) = CM2_1(i) + CM2_2(i);
      
    } else {
      CM2(i) = 0;
    }
  }
  return CM2; 
} 





