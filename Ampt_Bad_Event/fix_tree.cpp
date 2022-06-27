
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

using namespace std;


void fix_tree(string bad_path, string fix_path, string tree_name, vector<int> bad_events);

struct ampt_tree_branches {
	int event;
	int refmult;
	int refmult2;
	int refmult3;

	float imp;
	float qx, qy;

	int npp, npt, nesp, ninesp, nest, ninest;

	vector<int>* pid = 0;
	vector<float>* px = 0;
	vector<float>* py = 0;
	vector<float>* pz = 0;

	TBranch* branch_pid = 0;
	TBranch* branch_px = 0;
	TBranch* branch_py = 0;
	TBranch* branch_pz = 0;
};


void set_ampt_tree_branches(TTree* tree, ampt_tree_branches& branches) {
	tree->SetBranchAddress("event", &branches.event);

	tree->SetBranchAddress("refmult", &branches.refmult);
	tree->SetBranchAddress("refmult2", &branches.refmult2);
	tree->SetBranchAddress("refmult3", &branches.refmult3);

	tree->SetBranchAddress("imp", &branches.imp);
	tree->SetBranchAddress("qx", &branches.qx);
	tree->SetBranchAddress("qy", &branches.qy);

	tree->SetBranchAddress("npp", &branches.npp);
	tree->SetBranchAddress("npt", &branches.npt);
	tree->SetBranchAddress("nesp", &branches.nesp);
	tree->SetBranchAddress("ninesp", &branches.ninesp);
	tree->SetBranchAddress("nest", &branches.nest);
	tree->SetBranchAddress("ninest", &branches.ninest);

	tree->SetBranchAddress("pid", &branches.pid, &branches.branch_pid);
	tree->SetBranchAddress("px", &branches.px, &branches.branch_px);
	tree->SetBranchAddress("py", &branches.py, &branches.branch_py);
	tree->SetBranchAddress("pz", &branches.pz, &branches.branch_pz);
}


void fix_tree(string bad_path, string fix_path, string tree_name, vector<int> bad_events) {
	TFile *f_in = new TFile(bad_path.data(), "READ");
	TTree *tree_in = (TTree*)f_in->Get(tree_name.data());
	
	ampt_tree_branches branches;
	set_ampt_tree_branches(tree_in, branches);

	TFile* f_out = new TFile(fix_path.data(), "RECREATE");
	TTree* tree_out = new TTree(tree_name.data(), "AMPT Data");

	cout << "Fixing " << bad_path << " to " << fix_path << endl;

	int buffer_size = 5000000;
	int split_level = 1;

	//Define event branches:-------------------------------------------
	tree_out->Branch("event", &branches.event, "event/I");
	tree_out->Branch("refmult", &branches.refmult, "refmult/I");
	tree_out->Branch("refmult2", &branches.refmult2, "refmult2/I");
	tree_out->Branch("refmult3", &branches.refmult3, "refmult3/I");
	tree_out->Branch("qx", &branches.qx, "qx/F");
	tree_out->Branch("qy", &branches.qy, "qy/F");
	tree_out->Branch("imp", &branches.imp, "imp/F");

	tree_out->Branch("npp", &branches.npp, "npp/I");
	tree_out->Branch("npt", &branches.npt, "npt/I");
	tree_out->Branch("nesp", &branches.nesp, "nesp/I");
	tree_out->Branch("ninesp", &branches.ninesp, "ninesp/I");
	tree_out->Branch("nest", &branches.nest, "nest/I");
	tree_out->Branch("ninest", &branches.ninest, "ninest/I");

	//particle branches:
	tree_out->Branch("pid", &branches.pid, buffer_size, split_level);
	tree_out->Branch("px", &branches.px, buffer_size, split_level);
	tree_out->Branch("py", &branches.py, buffer_size, split_level);
	tree_out->Branch("pz", &branches.pz, buffer_size, split_level);

	int event_index = 0;
	while (tree_in->GetEvent(event_index++)) {
		if (find(bad_events.begin(), bad_events.end(), event_index) == bad_events.end()) {
			tree_out->Fill();  // Only fill event if not in bad_event list
		}
	}

	f_out->cd();
	tree_out->Write();
	f_out->Close();

	cout << bad_path << " fixed successfully to " << fix_path << endl;

	tree_in->ResetBranchAddresses();
	f_in->Close();
}