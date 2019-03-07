// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// NOTE: We don't use Doxygen //' comments here, otherwise the function  
// docstring is added to RcppExports.R with NULL content and this generates a
// "missing name" error when Doxygen generated documentation. The two fixes 
// for this seem to be (1) Add an @name annotation, which allows Doxygen to
// generate and .Rd file from the docstring in RcppExports.R or (2) Do not
// provide a Doxygen docstring for functions like this that we are not exporting
//

//' @title
//' Compute Earth Mover's Distance (EMD) between two Empirical Cumulative 
//' Density Functions (ECDFs)
//'
//' @param locations1 Locations for ECDF 1
//' @param values1 Cumulative masses for ECDF 1
//' @param locations2 Locations for ECDF 2
//' @param values2 Cumulative masses for ECDF 2
//'
//' @export
// [[Rcpp::export]]
double emd_fast_no_smoothingYdir(NumericVector locations1, NumericVector values1, 
                      NumericVector locations2, NumericVector values2)
{
  double segment_y_end   = 0;
  double segment_y_begin = 0; 
  double segment_area = 0;
  // Set start location of sweep below the minimum location across both ECDFs so
  // that we accumulate the full area between the two ECDFs
  double segment_loc1_start = locations1[0];
  double segment_loc2_start = locations2[0];
  if (values1[0] < values2[0])
  {
  	segment_y_end=values1[0];
  }
  else
  {
  	segment_y_end=values2[0];
  }
  double segment_value1 = 0;
  double segment_value2 = 0;
  int index1 = 0;
  int index2 = 0;
  double emd = 0;
  double max_value_ecdf = 1.0;
  double compensation = 0.0;
  double temp1 = 0;
  // We scan across all locations in both ECDFs, calculting the area of
  // each rectangular segment between the two ECDFs.
  double loc1 = locations1[0]-1;
  double loc2 = locations2[0]-1;
  int i;
  double curval = 0; 
  for (i=0;i< 1 + locations1.size()+locations2.size();i++)
  {
    // Calculate the area of the next rectangular segment...
    if ((values1[index1] < values2[index2]) || (index2==values2.size()))
    {
      // need to update the 
      segment_area += (values1[index1] - curval)*std::abs(locations1[index1]-locations2[index2]);
      temp1=(values1[index1] - curval)*std::abs(locations1[index1]-locations2[index2]);
      curval = values1[index1];
//      std::cout << i << " i1=" << index1 << " i2=" << index2 << " l1="<< loc1 << " l2=" << loc2 << " v1=" << temp1 << "\n";
      index1+=1;
    }
    else
    {
      segment_area += (values2[index2] - curval)*std::abs(locations1[index1]-locations2[index2]);
      temp1=(values2[index2] - curval)*std::abs(locations1[index1]-locations2[index2]);
      curval = values2[index2];
 //     std::cout << i << " i1=" << index1 << " i2=" << index2 << " l1="<< loc1 << " l2=" << loc2 << " v1=" << temp1 << "\n";
      index2+=1;
    }
	}
   return segment_area;
}




//' @title
//' Compute Earth Mover's Distance (EMD) between two Empirical Cumulative 
//' Density Functions (ECDFs)
//'
//' @param locations1 Locations for ECDF 1
//' @param values1 Cumulative masses for ECDF 1
//' @param locations2 Locations for ECDF 2
//' @param values2 Cumulative masses for ECDF 2
//'
//' @export
// [[Rcpp::export]]
double median_fast_no_smoothing(NumericVector locations1, NumericVector values1, 
                      NumericVector locations2, NumericVector values2,double limit)
{
  double segment_area = 0;
  // Set start location of sweep below the minimum location across both ECDFs so
  // that we accumulate the full area between the two ECDFs
  int index1 = 0;
  int index2 = 0;
  // We scan across all locations in both ECDFs, calculting the area of
  // each rectangular segment between the two ECDFs.
  int i;
  double curval = 0; 
  std::vector<std::pair< double, double> > offsetsGroup;
  std::vector<std::pair< double, double> > alloffsetsGroup;
  double temp1;
  double temp2;
  for (i=0;i< locations1.size()+locations2.size();i++)
  {
    // Calculate the area of the next rectangular segment...
    if ((values1[index1] < values2[index2]) || (index2==values2.size()))
    {
      // need to update the 
      temp1 = (values1[index1] - curval);
      if (std::abs(temp1)<0.0000000001)
      {
        curval = values1[index1];
        index1+=1;
        continue;
      }
      temp2 = (locations1[index1]-locations2[index2]);
      std::pair<double, double> tempDouble;
      tempDouble.first=temp2;
      tempDouble.second=temp1;
      offsetsGroup.push_back(tempDouble);
 //     alloffsetsGroup.push_back(tempDouble);
      curval = values1[index1];
      index1+=1;
    }
    else
    {
      temp1 = (values2[index2] - curval);
      if (std::abs(temp1)<0.0000000001)
      {
        curval = values2[index2];
        index2+=1;
        continue;
      }
      temp2 = (locations1[index1]-locations2[index2]);
      std::pair<double, double> tempDouble;
      tempDouble.first=temp2;
      tempDouble.second=temp1;
      offsetsGroup.push_back(tempDouble);
    //  alloffsetsGroup.push_back(tempDouble);
      curval = values2[index2];
      index2+=1;
    }
	}
  std::sort(offsetsGroup.begin(),offsetsGroup.end());
  double curTotal = 0;
  double offset;
  for (i=0;i<offsetsGroup.size();i++)
  {
    curTotal+=offsetsGroup[i].second; 
    if (curTotal >0.5)
    {
      offset = offsetsGroup[i].first;
      break;
    }
  }
  double emd=0;
  for (i=0;i<offsetsGroup.size();i++)
  {
    emd += offsetsGroup[i].second*std::abs(offsetsGroup[i].first-offset);
  }
  return emd;
}