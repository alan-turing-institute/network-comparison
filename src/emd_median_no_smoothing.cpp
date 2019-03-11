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
                      NumericVector locations2, NumericVector values2)
{
  double segment_area = 0;
  int index1 = 0;
  int index2 = 0;
  int i;
  double curval = 0; 
  std::vector<std::pair< double, double> > offsetsGroup;
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
  int split1;
  for (i=0;i<offsetsGroup.size();i++)
  {
    curTotal+=offsetsGroup[i].second; 
    if (curTotal >0.5)
    {
      offset = offsetsGroup[i].first;
      split1=i;
      break;
    }
  }
  double emd=0;
  for (i=0;i<split1;i++)
  {
    emd += offsetsGroup[i].second*(offset - offsetsGroup[i].first);
  }
  for (i=split1+1;i<offsetsGroup.size();i++)
  {
    emd += offsetsGroup[i].second*(offsetsGroup[i].first - offset);
  }
  return emd;
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
double median_fast_no_smoothingLimit(NumericVector locations1, NumericVector values1, 
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
  double curTotal = 0;
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
      alloffsetsGroup.push_back(tempDouble);
      if (temp2<-limit)
      {
        curTotal += temp1;
      }
      else if (temp2<limit)
      {
        offsetsGroup.push_back(tempDouble);
      }
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
      alloffsetsGroup.push_back(tempDouble);
      if (temp2<-limit)
      {
        curTotal += temp1;
      }
      else if (temp2<limit)
      {
        offsetsGroup.push_back(tempDouble);
      }
      curval = values2[index2];
      index2+=1;
    }
	}
  std::sort(offsetsGroup.begin(),offsetsGroup.end());
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
  for (i=0;i<alloffsetsGroup.size();i++)
  {
    emd += alloffsetsGroup[i].second*std::abs(alloffsetsGroup[i].first-offset);
  }
  std::cout << offsetsGroup.size() << " " << alloffsetsGroup.size() << "\n"; 
  return emd;
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
NumericVector median_fast_no_smoothingSlow3(NumericVector locations1, NumericVector values1, 
                      NumericVector locations2, NumericVector values2,
                      NumericVector locations3, NumericVector values3)
{
 NumericVector result(3);  
  result[0] =   median_fast_no_smoothing(locations1,values1,locations2,values2);
  result[1] =   median_fast_no_smoothing(locations1,values1,locations3,values3);
  result[2] =   median_fast_no_smoothing(locations2,values2,locations3,values3);
  return result;
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
NumericVector median_fast_no_smoothing3(NumericVector locations1, NumericVector values1, 
                      NumericVector locations2, NumericVector values2,
                      NumericVector locations3, NumericVector values3)
{
  double segment_area = 0;
  int index1 = 0;
  int index2 = 0;
  int index3 = 0;
  int i;
  double curval = 0; 
  std::vector<std::pair< double, double> > offsetsGroup12;
  std::vector<std::pair< double, double> > offsetsGroup13;
  std::vector<std::pair< double, double> > offsetsGroup23;
  double temp1;
  double temp2;
  double temp12; 
  double temp23; 
  double temp13; 
  temp12 = (locations1[index1]-locations2[index2]);
  temp13 = (locations1[index1]-locations3[index3]);
  temp23 = (locations2[index2]-locations3[index3]);
  double temp12_new = temp12; 
  double temp23_new = temp13; 
  double temp13_new = temp23; 
  for (i=0;i< locations1.size()+locations2.size()+locations3.size();i++)
  {
    
    // Calculate the area of the next rectangular segment...
    if ((values1[index1] < values2[index2]) || (index2==values2.size()))
    {
      if ((values1[index1] < values3[index3]) || (index3==values3.size()))
      {
        temp1 = (values1[index1] - curval);
        curval = values1[index1];
        index1+=1;
        temp12_new = (locations1[index1]-locations2[index2]);
        temp13_new = (locations1[index1]-locations3[index3]);
      }
      else
      {
        temp1 = (values3[index3] - curval);
        curval = values3[index3];
        index3+=1;
        temp13_new = (locations1[index1]-locations3[index3]);
        temp23_new = (locations2[index2]-locations3[index3]);
      }
      
    }
    else 
    {
      if ((values2[index2] < values3[index3]) || (index3==values3.size()))
      {
          temp1 = (values2[index2] - curval);
          curval = values2[index2];
          index2+=1;
          temp12_new = (locations1[index1]-locations2[index2]);
          temp23_new = (locations2[index2]-locations3[index3]);
      }
      else
      {
        temp1 = (values3[index3] - curval);
        curval = values3[index3];
        index3+=1;
        temp13_new = (locations1[index1]-locations3[index3]);
        temp23_new = (locations2[index2]-locations3[index3]);
      }
    }
      if (temp1>0.0000001)
      {
        
      // temp12 
      std::pair<double, double> tempDouble12;
      tempDouble12.first=temp12;
      tempDouble12.second=temp1;
      offsetsGroup12.push_back(tempDouble12);
      
      // temp13 
      std::pair<double, double> tempDouble13;
      tempDouble13.first=temp13;
      tempDouble13.second=temp1;
      offsetsGroup13.push_back(tempDouble13);
      
      // temp23 
      std::pair<double, double> tempDouble23;
      tempDouble23.first=temp23;
      tempDouble23.second=temp1;
      offsetsGroup23.push_back(tempDouble23);
      
      }
      temp12 = temp12_new; 
      temp23 = temp13_new; 
      temp13 = temp23_new; 
	}
  
  //test 1 v 2
  std::sort(offsetsGroup12.begin(),offsetsGroup12.end());
  double curTotal = 0;
  double offset;
  int split1=0;
  for (i=0;i<offsetsGroup12.size();i++)
  {
    curTotal+=offsetsGroup12[i].second; 
    if (curTotal >0.5)
    {
      offset = offsetsGroup12[i].first;
      split1=i;
      break;
    }
  }
  double emd12=0;
  for (i=0;i<split1;i++)
  {
    emd12 += offsetsGroup12[i].second*(offset - offsetsGroup12[i].first);
  }
  for (i=split1;i<offsetsGroup12.size();i++)
  {
    emd12 += offsetsGroup12[i].second*(offsetsGroup12[i].first-offset);
  }
  
 
  //test 1 v 3
  std::sort(offsetsGroup13.begin(),offsetsGroup13.end());
  curTotal = 0;
  for (i=0;i<offsetsGroup13.size();i++)
  {
    curTotal+=offsetsGroup13[i].second; 
    if (curTotal >0.5)
    {
      offset = offsetsGroup13[i].first;
      break;
    }
  }
  double emd13=0;
  for (i=0;i<offsetsGroup13.size();i++)
  {
    emd13 += offsetsGroup13[i].second*std::abs(offsetsGroup13[i].first-offset);
  }
  
  
 
  //test 2 v 3
  std::sort(offsetsGroup23.begin(),offsetsGroup23.end());
  curTotal = 0;
  for (i=0;i<offsetsGroup23.size();i++)
  {
    curTotal+=offsetsGroup23[i].second; 
    if (curTotal >0.5)
    {
      offset = offsetsGroup23[i].first;
      break;
    }
  }
  double emd23=0;
  for (i=0;i<offsetsGroup23.size();i++)
  {
    emd23 += offsetsGroup23[i].second*std::abs(offsetsGroup23[i].first-offset);
  }
  NumericVector result(3);  
  result[0] =  emd12;
  result[1] =  emd13; 
  result[2] =  emd23; 
  return result;
}