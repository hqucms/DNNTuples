/*
 * TreeWriter.h
 *
 * Class to help with tree writing.
 *
 */

#ifndef NTUPLECOMMONS_INTERFACE_TREEWRITER_H_
#define NTUPLECOMMONS_INTERFACE_TREEWRITER_H_

#include <string>
#include <vector>
#include <TString.h>
#include <TTree.h>

namespace deepntuples {

class TreeWriter {

public :
  TreeWriter(TTree *tree, const char *treename="tree"): fTreeName(treename), fTree(tree) {
    fTree->SetName(fTreeName.Data());
  }

  ~TreeWriter() { delete fTree; }

  TString getTreeName() const { return fTreeName; }
  void    setTreeName(TString n) { fTreeName = n;  fTree->SetName(fTreeName.Data()); }
  void    fill() { fTree->Fill();  }
  TTree   *getTree() { return fTree; }

  template<class T>
  void    book(const char *name, T& var, const char *type)  { fTree->Branch(name, &var, TString(name).Append("/").Append(type).Data()); }

  template<class T>
  void    book(const char *name, std::vector<T>& varv)    { fTree->Branch(name, &varv); }

protected :
  TString fTreeName;
  TTree   *fTree;

}; // TreeWriter



} /* namespace deepntuples */

#endif /* NTUPLECOMMONS_INTERFACE_TREEWRITER_H_ */
