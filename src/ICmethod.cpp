#include <Rcpp.h>

#define LOG_LEVEL 0

#if LOG_LEVEL >= 2

#define LOG_DEBUG( msg ) Rcpp::Rcout << msg << "\n";

#else

#define LOG_DEBUG( msg )

#endif

typedef double go_dist_func_t( double mica, double ic1, double ic2, double mic );

double go_dist_Resnik( double mica, double ic1, double ic2, double mic )
{
  // Resnik does not consider how distant the terms are from their common ancestor.
  return mica;
}

double go_dist_Lin( double mica, double ic1, double ic2, double mic )
{
  // Lin takes common ancestor distances into account.
  return 2 * mica / (ic1+ic2);
}

double go_dist_Jiang( double mica, double ic1, double ic2, double mic )
{
  // Jiang takes common ancestor distances into account.
  return std::max( 0.0, 1.0 - ( -2 * mica + ic1 + ic2 ) );
}

double go_dist_Rel( double mica, double ic1, double ic2, double mic )
{
  // mica*mic equals to the original IC value.
  // and exp(-mica*mic) equals to the probability of the term's occurence.
  return 2 * mica/(ic1+ic2) * (1-exp(-mica*mic));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix infoContentMethod_cpp(
  Rcpp::StringVector&  id1_,
  Rcpp::StringVector&  id2_,
  Rcpp::List&          anc_,
  Rcpp::NumericVector& ic_,
  const std::string&   method_,
  const std::string&   ont_
){
  go_dist_func_t* go_dist;
  // Resnik does not consider how distant the terms are from their common ancestor.
  //  Lin and Jiang take that distance into account.
  if (method_ == "Resnik") {
    go_dist = &go_dist_Resnik;
  }
  else if (method_ == "Lin") {
    go_dist = &go_dist_Lin;
  }
  else if (method_ == "Jiang") {
    go_dist = &go_dist_Jiang;
  }
  else if (method_ == "Rel") {
    go_dist = &go_dist_Rel;
  }
  else {
    throw std::runtime_error( "Unknown GO distance method" );
  }
  
  typedef std::string term_id_t;
  typedef std::set<term_id_t> term_set_t;

  // calculate the maximum IC and build the map of normalized IC
  typedef std::map<term_id_t, double> ic_map_t;
  ic_map_t normIcMap;
  // more specific term, larger IC value.
  // Normalized, all divide the most informative IC.
  // all IC values range from 0(root node) to 1(most specific node)
  double mic = NA_REAL;
  {
    Rcpp::StringVector icNames( ic_.names() );
    for (std::size_t i=0; i < ic_.size(); i++ ) {
      const double cic = ic_[i];
      if ( Rcpp::NumericVector::is_na( cic ) || cic == R_PosInf ) continue;
      if ( Rcpp::NumericVector::is_na( mic ) || mic < cic ) mic = cic;
    }
    LOG_DEBUG( "mic=" << mic );
    for (std::size_t i=0; i < ic_.size(); i++ ) {
      const double cic = ic_[i];
      if ( Rcpp::NumericVector::is_na( cic ) || cic == R_PosInf ) continue;
      normIcMap.insert( std::make_pair( (std::string) icNames[i], cic / mic ) );
    }
  }

  // set root node IC to 0
  if(ont_ == "DO") {
    normIcMap["DOID:4"] = 0;
  } else {
    normIcMap["all"] = 0;
  }

  // convert anc_ into map of sets
  typedef std::map<term_id_t, term_set_t> anc_map_t;
  anc_map_t ancMap;
  {
    Rcpp::StringVector goTerms( anc_.names() );
    for (std::size_t i=0; i < anc_.size(); i++ ) {
      const std::vector<std::string> ancVec = Rcpp::as<std::vector<std::string> >( anc_[i] );
      term_set_t ancestors( ancVec.begin(), ancVec.end() );
      // term itself is also considered an ancestor
      ancestors.insert( (std::string)goTerms[i] );
      ancMap.insert( std::make_pair( (std::string) goTerms[i], ancestors ) );
    }
  }

  Rcpp::NumericMatrix res( id1_.size(), id2_.size() );
  res.attr("dimnames") = Rcpp::Rcpp_list2( id1_, id2_ );
  for ( std::size_t i = 0; i < id1_.size(); i++ ) {
    const std::string id1_term = (std::string)id1_[i];
    const ic_map_t::const_iterator iIcIt = normIcMap.find( id1_term );
    if ( iIcIt != normIcMap.end() && iIcIt->second != 0 ) {
      const double iIc = iIcIt->second;
      LOG_DEBUG( "ic[" << id1_term << "]=" << iIc );
      const anc_map_t::const_iterator iAncsIt = ancMap.find( id1_term );
      for ( std::size_t j = 0; j < id2_.size(); j++ ) {
        const std::string id2_term = (std::string)id2_[j];
        const ic_map_t::const_iterator jIcIt = normIcMap.find( id2_term );
        if ( jIcIt != normIcMap.end() && jIcIt->second != 0 ) {
          const anc_map_t::const_iterator jAncsIt = ancMap.find( id2_term );
          // find common ancestors
          term_set_t commonAncs;
          if ( iAncsIt != ancMap.end() && jAncsIt != ancMap.end() ) {
            std::set_intersection( iAncsIt->second.begin(), iAncsIt->second.end(),
                            jAncsIt->second.begin(), jAncsIt->second.end(),
                            std::inserter( commonAncs, commonAncs.end() ) );
          }
          LOG_DEBUG( "n(commonAncs(" << id1_term << "," << id2_term << "))=" << commonAncs.size() );

          // Information Content of the most informative common ancestor (MICA)
          double mica = 0;
          for ( term_set_t::const_iterator termIt = commonAncs.begin(); termIt != commonAncs.end(); ++termIt ) {
            ic_map_t::const_iterator ancIcIt = normIcMap.find( *termIt );
            if ( ancIcIt != normIcMap.end() && mica < ancIcIt->second ) mica = ancIcIt->second; 
          }
          LOG_DEBUG( "mica(" << id1_term << "," << id2_term << ")=" << mica );
          res(i,j) = go_dist( mica, iIc, jIcIt->second, mic );
        } else {
          res(i,j) = NA_REAL;
        }
      }
    } else {
      for ( std::size_t j = 0; j < id2_.size(); j++ ) {
        res(i,j) = NA_REAL;
      }
    }
  }
  return ( res );
}
