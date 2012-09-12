#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP infoContentMethod_cpp(SEXP id1_, SEXP id2_, SEXP anc1_, SEXP anc2_, SEXP icname_, SEXP ic_, SEXP method_, SEXP ont_) {
  std::string ID1 = as<std::string>(id1_);
  std::string ID2 = as<std::string>(id2_);
  std::vector<std::string> anc1 = as<std::vector<std::string> >(anc1_);
  std::vector<std::string> anc2 = as<std::vector<std::string> >(anc2_);
  std::vector<std::string> ICName = as<std::vector<std::string> >(icname_);
  std::vector<double> IC = as<std::vector<double> >(ic_);
  std::string method = as<std::string>(method_);
  std::string ont = as<std::string>(ont_);

  // more specific term, larger IC value.
  // Normalized, all divide the most informative IC.
  // all IC values range from 0(root node) to 1(most specific node)
  double mic = 0;
  for(int i=0; i < IC.size(); i++) {
    if (IC[i] != R_PosInf && IC[i] > mic)
      mic = IC[i];
  }

  std::string topNode;
  if(ont == "DO") {
    topNode = "DOID:4";
  } else {
    topNode = "all";
  }

  if (topNode == "all") {
    // for GO
    ICName.insert(ICName.end(), topNode);
    IC.insert(IC.end(), 0);
  } else {
    // for DO
    for(int i=0; i < ICName.size(); i++) {
      if (ICName[i] == topNode) {
	IC[i] = 0;
	break;
      }
    }
  }

  double ic1, ic2;
  bool f1=false, f2=false;
  for (int i=0; i < ICName.size(); i++) {
    if (ICName[i] == ID1) {
      ic1 = IC[i]/mic;
      f1 = true;
    }
    if (ICName[i] == ID2) {
      ic2 = IC[i]/mic;
      f2 = true;
    }
    if (f1 && f2)
      break;
  }
  if ( ! (f1 && f2) )
    return wrap(NA_REAL);

  if ( ic1 == 0 || ic2 == 0)
    return wrap(NA_REAL);

  std::vector<std::string> commonAnc;
  anc1.insert(anc1.begin(), ID1);
  anc2.insert(anc2.begin(), ID2);

  for (int i=0; i < anc1.size(); i++) {
    for (int j=0; j < anc2.size(); j++) {
      if (anc1[i] == anc2[j])
	commonAnc.insert(commonAnc.end(), anc1[i]);
    }
  }

  double mica = 0;
  for (int i=0; i<commonAnc.size(); i++) {
    for (int j=0; j < ICName.size(); j++) {
      if (commonAnc[i] == ICName[j]) {
	double ica = IC[j];
	if (mica < ica)
	  mica = ica;
      }
    }
  }

  // Information Content of the most informative common ancestor (MICA)
  mica /= mic;

  double sim;

  // Resnik does not consider how distant the terms are from their common ancestor.
  //  Lin and Jiang take that distance into account.
  if (method == "Resnik") {
    sim = mica;
  }
  if (method == "Lin") {
    sim = 2 * mica/(ic1+ic2);
  }
  if (method == "Jiang") {
    double d = -2 * mica + ic1 + ic2;
    if ( d >= 1)
      d = 1;
    sim = 1- d;
  }
  if (method == "Rel") {
    sim = 2 * mica/(ic1+ic2) * (1-exp(-mica*mic));
    // mica*mic equals to the original IC value.
    // and exp(-mica*mic) equals to the probability of the term's occurence.
  }
  return wrap(sim);
}
