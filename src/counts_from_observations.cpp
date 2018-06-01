#include <Rcpp.h>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

//' @title
//' Count number of occurences
//'
//' @param features A Matrix with doubles.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix  counts_from_observations(NumericMatrix features)
{
    std::unordered_map<double,int> counter;
    // Generate a map with one key for every unique value in the features vector
    // and a zero count for each value. If the elements of the feature vector
    // are floating point numbers, a new key will be created in the map unless
    // the elements are exactly equal.
    for (int i=0; i<features.nrow(); i++) {
      counter[features[i]]=0;
    }
    // Accumulate the counts for each unique element in the feature vector
    for (int i=0;i<features.nrow();i++) {
      counter[features[i]]+=1;
    }
    // Convert the unordered map into a vector so we can sort it by the 
    // feature elements
    std::vector<std::pair<double,int> > keys(counter.begin(),counter.end());
    std::sort(keys.begin(),keys.end());
    // Create an output matrix to return to R
    NumericMatrix result(keys.size(),2);
    for (int i=0;i<keys.size();i++)
    {
        result(i,0)=keys[i].first;
        result(i,1)=keys[i].second;
    }
    return result;
}
