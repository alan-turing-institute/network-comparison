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
    int i;
    // Silly solution to the float comparison problem
    for (i=0;i<features.nrow();i++) {
      counter[features[i]]=0;
    }
    for (i=0;i<features.nrow();i++) {
      counter[features[i]]+=1;
    }
    std::vector<std::pair<double,int> > keys(counter.begin(),counter.end());
    std::sort(keys.begin(),keys.end());
    NumericMatrix result(keys.size(),2);
    for (i=0;i<keys.size();i++)
    {
        result(i,0)=keys[i].first;
        result(i,1)=keys[i].second;
    }
    return result;
}
