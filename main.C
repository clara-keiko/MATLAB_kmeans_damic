#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>
#include <typeinfo>
#include "Riostream.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <fstream>
#include <string>

#include <cstdio>

#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"

#include "ConfigValue.C"
#include "Methods.C"


int rootTocsv(char *str, int index)
{
 
	ofstream myfile;
	myfile.open ("out.csv", ofstream::app);
	TFile *f2 = new TFile(str,"READ");
	TTree *tree2 = (TTree*)f2->Get("clusters_tree");

	const Double_t k2 = 2.906e-4; // calibration factors
	const Double_t k4 = 2.635e-4;
	const Double_t k5 = 2.725e-4;
	const Double_t k6 = 2.687e-4;
	const Double_t thr = 0.05; // threshold in keV
	const Int_t n = 100000; // arbitrary initial size of the arrays

	Int_t j = 0; // index for number of clusters with E > threshold 
	Int_t cluster_id2 = NULL;
	Double_t ocuparea;
	TParameter<double>* center_x2 = NULL;
	TParameter<double>* center_y2 = NULL;
	TParameter<double>* charge_total2 = NULL; // total charge of the cluster
	TParameter<double>* npix = NULL; // number of pixels in the cluster
	TParameter<double>* charge_maximum2 = NULL; // maximum pixel value in the track
	TParameter<double>* charge_max_x2 = NULL; // x position of maximum value pixel 
	TParameter<double>* charge_max_y2 = NULL; // y position of maximum value pixel
	TParameter<double>* track_startx2 = NULL; // x initial position of the track
	TParameter<double>* track_starty2 = NULL; // y initial position of the track
	TParameter<double>* track_endx2 = NULL; // x final position of the track
	TParameter<double>* track_endy2 = NULL; // y final position of the track
	TParameter<double>* track_length2 = NULL; // track length in pixels
	TParameter<double>* charge_rms_x2 = NULL; // x RMS in pixels
	TParameter<double>* charge_rms_y2 = NULL; // y RMS in pixels
	TParameter<double>* raw_nsat2 = NULL;
	TParameter<double>* curve_track_rms2 = NULL;
	TParameter<double>* curve_correlation_factor2 = NULL;
	TParameter<double>* size_aspect_ratio2 = NULL;
	TParameter<double>* sigma_mean2 = NULL;

	Int_t veccluster[n];
	//Double_t NCLUSTERS;
	//Double_t veccarea[n];
	Double_t vecccenterx[n];
	Double_t vecccentery[n];
	Double_t vecct2[n];
	Double_t veccnpix[n];
	Double_t veccm2[n];
	Double_t veccmx2[n];
	Double_t veccmy2[n];
	Double_t vectrackstartx2[n];
	Double_t vectrackstarty2[n];
	Double_t vectrackendx2[n];
	Double_t vectrackendy2[n];
	Double_t vectracklength2[n];
	Double_t veccharge_rms_x2[n];
	Double_t veccharge_rms_y2[n];

	clusters_tree->SetBranchAddress("cluster_id",&cluster_id2);
	clusters_tree->SetBranchAddress("charge_total",&charge_total2);
	clusters_tree->SetBranchAddress("size_npixels",&npix);
	clusters_tree->SetBranchAddress("charge_maximum",&charge_maximum2);
	clusters_tree->SetBranchAddress("charge_max_x",&charge_max_x2);
	clusters_tree->SetBranchAddress("charge_max_y",&charge_max_y2);
	clusters_tree->SetBranchAddress("center_x",&center_x2);
	clusters_tree->SetBranchAddress("center_y",&center_y2);
	clusters_tree->SetBranchAddress("track_startx",&track_startx2);
	clusters_tree->SetBranchAddress("track_starty",&track_starty2);
	clusters_tree->SetBranchAddress("track_endx",&track_endx2);
	clusters_tree->SetBranchAddress("track_endy",&track_endy2);
	clusters_tree->SetBranchAddress("size_slength",&track_length2);
	clusters_tree->SetBranchAddress("charge_rms_x",&charge_rms_x2);
	clusters_tree->SetBranchAddress("charge_rms_y",&charge_rms_y2);
	clusters_tree->SetBranchAddress("raw_nsat",&raw_nsat2);
	clusters_tree->SetBranchAddress("curve_track_rms",&curve_track_rms2);
	clusters_tree->SetBranchAddress("curve_correlation_factor",&curve_correlation_factor2);
	clusters_tree->SetBranchAddress("size_aspect_ratio",&size_aspect_ratio2);
	clusters_tree->SetBranchAddress("sigma_mean",&sigma_mean2);
	for (int i=0; i<tree2->GetEntries(); i++)
	{
		tree2->GetEntry(i);

		if (k2*charge_total2->GetVal() > thr) // cuts on events above the threshold
		{

			if(npix->GetVal() > 3)
			{
			char *png_filec; 
			string cluster_id_s = Form("%d",cluster_id2);
			string ids = Form("%d",index);
			string png_file = ids + ".png";
			png_filec = (char *)malloc(png_file.size() + 1);
						
			memcpy(png_filec, png_file.c_str(), png_file.size() + 1);
			
			
			veccluster[j] = cluster_id2;
			vecct2[j] = k2*charge_total2->GetVal();
			veccnpix[j] = npix->GetVal();
			//veccarea[n] = veccnpix[j]/area;
			veccm2[j] = k2*charge_maximum2->GetVal();
			veccmx2[j] = charge_max_x2->GetVal();
			veccmy2[j] = charge_max_y2->GetVal();
			vectrackstartx2[j] = track_startx2->GetVal();
			vectrackstarty2[j] = track_starty2->GetVal();
			vectrackendx2[j] = track_endx2->GetVal();
			vectrackendy2[j] = track_endy2->GetVal();
			vectracklength2[j] = track_length2->GetVal();
			veccharge_rms_x2[j] = charge_rms_x2->GetVal();
			veccharge_rms_y2[j] = charge_rms_y2->GetVal();
			
			double area1 = Viewer(str,cluster_id2, png_filec);
			ocuparea = veccnpix[j]/area1;
			double center_parameter = sqrt(((center_x2->GetVal() - veccmx2[j])*(center_x2->GetVal() - veccmx2[j]))+((center_y2->GetVal() - veccmy2[j])*(center_y2->GetVal() - veccmy2[j])));
			double area2;
			double raio,ratio2;
			raio=sqrt((((track_endx2->GetVal())-(center_x2->GetVal()))*((track_endx2->GetVal())-(center_x2->GetVal())))+(((track_endy2->GetVal())-(center_y2->GetVal()))*((track_endy2->GetVal())-(center_y2->GetVal()))));
			area2=3.14*(raio*raio);
			ratio2=(npix->GetVal())/area2;
			double ratio3 = raw_nsat2->GetVal()/npix->GetVal();
			//cout << "String " << png_filec << endl;
			//cout << j << endl;
			//cout << "" << endl;
			myfile << vecct2[j] << ';' << veccnpix[j] << ';' << vectrackstartx2[j] << ';' << vectrackstarty2[j] << ';' << vectrackendx2[j] << ';' << vectrackendy2[j] << ';' << vectracklength2[j] << ';' << ocuparea << ';' << veccluster[j] << ';' << center_parameter << ';' << ratio3 << ';' << curve_track_rms2->GetVal() << ';' << curve_correlation_factor2->GetVal() << ';' << size_aspect_ratio2->GetVal() << ';' << sigma_mean2->GetVal() << '\n';
			j += 1;
			index += 1;
			free(png_filec);
			}
		}     
	}

	myfile.close();
	f2->Close();
	const Long64_t jj = j;
	Double_t Vecct2[jj];
	Double_t Veccnpix[jj];
	Double_t Veccm2[jj];
	Double_t Veccmx2[jj];
	Double_t Veccmy2[jj];
	Double_t Vectrackstartx2[jj];
	Double_t Vectrackstarty2[jj];
	Double_t Vectrackendx2[jj];
	Double_t Vectrackendy2[jj];
	Double_t Vectracklength2[jj];
	Double_t Veccharge_rms_x2[jj];
	Double_t Veccharge_rms_y2[jj];

	for (int i=0; i<jj; i++) // loop to reduce the size of the array
	{
		Vecct2[i] = vecct2[i];
		Veccnpix[i] = veccnpix[i];
		Veccm2[i] = veccm2[i];
		Veccmx2[i] = veccmx2[i];
		Veccmy2[i] = veccmy2[i];
		Vectrackstartx2[i] = vectrackstartx2[i];
		Vectrackstarty2[i] = vectrackstarty2[i];
		Vectrackendx2[i] = vectrackendx2[i];
		Vectrackendy2[i] = vectrackendy2[i];
		Vectracklength2[i] = vectracklength2[i];
		Veccharge_rms_x2[i] = veccharge_rms_x2[i];
		Veccharge_rms_y2[i] = veccharge_rms_y2[i];
	}

	//cout << "" << endl;
	//cout << "Number of candidates: " << sizeof(Veccm2)/sizeof(Long64_t) << endl;
	return index;


}

double Viewer(TString fname, Int_t cnum, char *str_file)
{
    Double_t area;
    LoadConfigFile("ConfigValue.C");
    Int_t zoomx = ConfigValue("canvas_xsize");
    Int_t zoomy = ConfigValue("canvas_ysize");
    
    TFile* f = new TFile(fname, "READ");
    TFile* canvas_file = new TFile("template.root", "READ");
    
    if(f->IsZombie()){
        std::cout << "Could not open " << fname << std::endl;
        return;
    }
    
    TTree* t = NULL;
    t = (TTree*) f->Get("clusters_tree");
    if(t==NULL){
        std::cout << "Could not open clusters_tree" << std::endl;
        std::cout << "Run FindClusters first!" << std::endl;
        f->Close();
        return;
    }
    
    TCanvas* c = NULL;
    c = (TCanvas*) canvas_file->Get("clusters_canvas");
    /*if(c==NULL){
        std::cout << "Could not open clusters_canvas" << std::endl;
        std::cout << "Run FindClusters first!" << std::endl;
        f->Close();
        return;
    }
    */

    //to know if simulated file
    TParameter<double>* hwid = NULL;
    t->SetBranchAddress("HWID", &hwid);
    t->GetEntry(0);
    bool is_sim = false;
    if(hwid!=NULL && hwid->GetVal()<0) is_sim = true;
    
    //arrays to store object
    TArrayD* pixel_x = NULL;
    TArrayD* pixel_y = NULL;
    TArrayD* pixel_val = NULL;
    TArrayD* pixel_sigma = NULL;
    TArrayD* pixel_ped = NULL;
    
    TArrayD* pixel_sime = NULL;
    TArrayD* pixel_simz = NULL;
    TArrayD* pixel_simx = NULL;
    TArrayD* pixel_simy = NULL;
    
    t->SetBranchAddress("pixel_x", &pixel_x);
    t->SetBranchAddress("pixel_y", &pixel_y);
    t->SetBranchAddress("pixel_val", &pixel_val);
    t->SetBranchAddress("pixel_sigma", &pixel_sigma);
    t->SetBranchAddress("pixel_ped", &pixel_ped);
    
    
    //for simulated
    /*if(is_sim){
     
        t->SetBranchAddress("pixel_sime", &pixel_sime);
        t->SetBranchAddress("pixel_simz", &pixel_simz);
        t->SetBranchAddress("pixel_simx", &pixel_simx);
        t->SetBranchAddress("pixel_simy", &pixel_simy);
    }*/
    

    //list to store processed cluster values
    TList* vals = new TList();
    
    Int_t size = t->GetEntries();
    
    if(cnum>=size || cnum<0)
        std::cout << "Cluster number out of range" << std::endl;
    
    else{
        
        t->GetEntry(cnum);
        

        TH2F* cluster = ArraysToTH2F(pixel_x, pixel_y, pixel_val);
        cluster->SetNameTitle("cluster","cluster");
        TH2F* sigma = ArraysToTH2F(pixel_x, pixel_y, pixel_sigma);
        sigma->SetNameTitle("sigma","sigma");
        TH2F* pedestal = ArraysToTH2F(pixel_x, pixel_y, pixel_ped);
        pedestal->SetNameTitle("pedestal","pedestal");
        
        TString ctitle = "Cluster #";
        ctitle+=cnum;
        cluster->SetTitle(ctitle);
        sigma->SetTitle(ctitle);
        pedestal->SetTitle(ctitle);
        
        Int_t xmin = cluster->GetXaxis()->GetXmin();
        Int_t xmax = cluster->GetXaxis()->GetXmax();
        Int_t ymin = cluster->GetYaxis()->GetXmin();
        Int_t ymax = cluster->GetYaxis()->GetXmax();
        
        Double_t yrange = ymax-ymin;
        Double_t xrange = xmax-xmin;
        Double_t area = yrange*xrange;

	//Double_t raio;
	//raio=sqrt(((xrange/2)*(xrange/2))+((yrange/2)*(yrange/2)));
	//Double_t carea=3.14*(raio*raio)/2;
	/*int countpx = 0;
	len = (sizeof(pixel_x)/sizeof(*pixel_x));
	for(int i=0; i<len; i++)
	{
		if(srqt(((pixel_x - (xrange/2))*(pixel_x - (xrange/2)))+((pixel_y - (yrange/2))*(pixel_y - (yrange/2))))<raio)
		{
			cout << 'teste: ' << pixel_x << '\n';			
			countpx ++;
		}		
	}
	*/
	//Double_t ratio=(countpx/carea*(1/area));
	
	
        xmin -= (zoomx - xrange)/2;
        xmax += (zoomx - xrange)/2;
        ymin -= (zoomy - yrange)/2;
        ymax += (zoomy - yrange)/2;
        
        ((TH2F*) c->FindObject("image"))->SetAxisRange(xmin,xmax,"X");
        ((TH2F*) c->FindObject("image"))->SetAxisRange(ymin,ymax,"Y");
        c->Draw();
        	

        //Draw cluster
        TCanvas* myc = new TCanvas();
        cluster->Draw("COLZ");
        myc->SaveAs(str_file);
	myc->Close();
       
        //Print vals
        for(Int_t i=0; i<vals->GetSize(); i++)
        {    
		std::cout << vals->At(i)->GetName() << " " << ( (TParameter<double>*) vals->At(i))->GetVal() << std::endl;
        }
        vals->Clear();
    }

    
    return area;
}

void main()
{
	std::ifstream infile("rootfiles.csv");
	std::string line;
	int index = 1;
	while (getline(infile, line))
	{
		const char * line_c = line.c_str();
		index = rootTocsv(line_c,index);
		
		cout << line << '\n';
	}
	

}




