#include "Tables.h"

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>



bool stringCompare(TString left_,TString right_ ){
  const string left(left_, left_.Sizeof());
  const string right(right_, right_.Sizeof());
  for( string::const_iterator lit = left.begin(), rit = right.begin(); lit != left.end() && rit != right.end(); ++lit, ++rit )
      if( tolower( *lit ) < tolower( *rit ) )
         return true;
      else if( tolower( *lit ) > tolower( *rit ) )
         return false;
   if( left.size() < right.size() )
      return true;
   return false;
}


Tables::Tables(TString name){
  Name=name;
}

Tables::~Tables(){
}
  
void  Tables::MakeEffTable(std::vector<TH1D> histo, std::vector<TString>  names,float Lumi,std::vector<float> CrossSectionandAcceptance,
			       std::vector<float> nevents){


  ofstream output;
  output.open(Name+".tex", ios::out);
  (output) << "\\documentclass[11pt]{article}" << std::endl;
  (output) << "\\usepackage[width=6.5in, height=8.5in]{geometry}" << std::endl;
  (output) << "\\usepackage{epsfig}" << std::endl;
  (output) << "\\usepackage{cite}" << std::endl;
  (output) << "\\usepackage{fancyhdr}" << std::endl;
  (output) << "\\usepackage{lscape,graphicx}" << std::endl;
  TString title=Name;
  title.ReplaceAll("_"," ");
  (output) << "\\title{Analysis Results for: " << title <<"}" << std::endl;
  (output) << "\\begin{document}" << std::endl; 
  (output) << "\\maketitle" << std::endl;
  (output) << "\\tableofcontents" << std::endl;
  (output) << "\\newpage" << std::endl;
  (output) << "\\listoftables" << std::endl;
  (output) << "\\newpage" << std::endl;
  (output) << "\\listoffigures" << std::endl; 
  (output) << "\\clearpage" << std::endl; 


  (output) << "\\section{Event Cut Flow Tables}" << std::endl; 
  //Makes the Selection Table
  if(histo.size()>0){
    int nhist=histo.size();
    int ncuts=histo[0].GetNbinsX()-1;

    bool flag=true;
    int ncol=6;
    int k=0;
    int l=0;
    int t=0;
    (output) << "\\begin{landscape}" << std::endl;
    while(flag){
      (output) << "\\begin{table}[t]" << std::endl;
      (output) << "\\tiny " << std::endl;
      (output) << "\\begin{center}" << std::endl;
      (output) << "\\begin{tabular}{|p{5cm}|";
    
      for(Int_t i=0; i<ncol && i<nhist-k;i++){
	(output) << "p{1.75cm}|";
      }
      (output) << "} \\hline" << endl;
      t=0;
      for(Int_t j=0; j<nhist;j++){
	if(k<=j && j<nhist && j<k+ncol && t<histo.size()){
	  TString id=histo[t].GetTitle();
	  id.ReplaceAll("N Passed ","");
	  (output) << " &  $" << id  << "$";
	}
	t++;
      }
      (output) << " \\\\  \\hline" << std::endl;

      t=0;
      (output) << " Before Skim";
      for(Int_t j=0; j<nhist;j++){
        if(k<=j && j<nhist && j<k+ncol && t<histo.size()){
	  (output) << " &  $" << nevents[j]  << "$";
        }
        t++;
      }
      (output) << " \\\\  " << std::endl;

      
      for(int i=-1; i<ncuts;i++){
	t=0;
	if(i>=0){
	  TString title=names[i];
	  (output) << title;
	}
	else{
	  (output) << "Before Cuts";
	}

	for(Int_t j=0; j<nhist;j++){
	  if(k<=j  && j<nhist && j<k+ncol && t<histo.size() ){
	    (output) << " &   "<< histo[t].GetBinContent(i+2);
	  }
	  t++;
	}
	(output) << " \\\\  " << std::endl;
      }
      (output) << " \\hline" << endl;
      
      (output) << "\\end{tabular}"                                << std::endl;
      (output) << "\\caption[Raw Event Selection ]{Raw Event Selection }" << std::endl;
      (output) << "\\end{center}"                                 << std::endl;
      (output) << "\\end{table}"                               << std::endl;
      (output) << "\\normalsize" << std::endl;
      k+=ncol;
      if(k>=nhist){
	flag=false;
      }
      l++;
    }
    (output) << "\\end{landscape}" << std::endl;
  }


 //Makes the Lumi Norm Selection Table
  if(histo.size()>0){
    int nhist=histo.size();
    int ncuts=histo[0].GetNbinsX()-1;

    bool flag=true;
    int ncol=6;
    int k=0;
    int l=0;
    int t=0;
    (output) << "\\begin{landscape}" << std::endl;
    while(flag){
      (output) << "\\begin{table}[t]" << std::endl;
      (output) << "\\tiny " << std::endl;
      (output) << "\\begin{center}" << std::endl;
      (output) << "\\begin{tabular}{|p{5cm}|";
    
      for(Int_t i=0; i<ncol && i<nhist-k;i++){
	(output) << "p{1.75cm}|";
      }
      (output) << "} \\hline" << endl;
      t=0;
      for(Int_t j=0; j<nhist;j++){
	if(k<=j && j<nhist && j<k+ncol && t<histo.size()){
	  TString id=histo[t].GetTitle();
	  id.ReplaceAll("N Passed ","");
	  (output) << " &  $" << id  << "$";
	}
	t++;
      }
      
      (output) << " \\\\  \\hline" << std::endl;

      t=0;
      (output) << "Cross Section" << std::endl;
      for(Int_t j=0; j<nhist;j++){
	if(k<=j && j<nhist && j<k+ncol && t<histo.size()){
	  (output) << " &  " << CrossSectionandAcceptance[j] << "$pb$";
	}
	t++;
      }
      (output) << " \\\\  " << std::endl;

      t=0;
      (output) << "Lumi" << std::endl;
      for(Int_t j=0; j<nhist;j++){
	if(k<=j && j<nhist && j<k+ncol && t<histo.size()){
	  (output) << " &  " << Lumi << "$pb^{-1}$ ";
	}
	t++;
      }
      
      (output) << " \\\\  \\hline" << std::endl;
      t=0;
      (output) << " Before Skim";
      for(Int_t j=0; j<nhist;j++){
        if(k<=j && j<nhist && j<k+ncol && t<histo.size()){
          if(CrossSectionandAcceptance[j]>0)(output) << "  &  $" << nevents[j]*Lumi*CrossSectionandAcceptance[j]/nevents[j]  << "$";
	  else (output) << "  &  $" << nevents[j]  << "$";
	}
        t++;
      }
      (output) << " \\\\  " << std::endl;


      for(int i=-1; i<ncuts;i++){
	t=0;
	if(i>=0){
	  TString title=names[i];
	  (output) << title;
	}
	else{
	  (output) << "Before Cuts";
	}

	for(Int_t j=0; j<nhist;j++){
	  if(k<=j  && j<nhist && j<k+ncol && t<histo.size() ){
	  if(nevents[t]>0){
	    if(CrossSectionandAcceptance[t]>0){
	      (output) << " &   "<< histo[t].GetBinContent(i+2)*(Lumi*CrossSectionandAcceptance[t])/nevents[t];
	    }
	    else{
	      (output) << " &   "<< histo[t].GetBinContent(i+2);
	    }
	  }
	  else{
	    (output) << " &   NA";
	  }
	}
	t++;
	}
	(output) << " \\\\  " << std::endl;
      }
      (output) << " \\hline" << endl;
      
      (output) << "\\end{tabular}"                                << std::endl;
      (output) << "\\caption[Luminosity Normilized Event Selection ]{Luminosity Normilized Event Selection }" << std::endl;
      (output) << "\\end{center}"                                 << std::endl;
      (output) << "\\end{table}"                               << std::endl;
      (output) << "\\normalsize" << std::endl;
      k+=ncol;
      if(k>=nhist){
	flag=false;
      }
      l++;
    }
    (output) << "\\end{landscape}" << std::endl;
  }


  // Makes Accumulative Efficiency table
  if(histo.size()>0){
    int nhist=histo.size();
    int ncuts=histo[0].GetNbinsX()-1;

    bool flag=true;
    int ncol=6;
    int k=0;
    int l=0;
    int t=0;
    (output) << "\\begin{landscape}" << std::endl;
    while(flag){
      (output) << "\\begin{table}[t]" << std::endl;
      (output) << "\\tiny " << std::endl;
      (output) << "\\begin{center}" << std::endl;
      (output) << "\\begin{tabular}{|p{5cm}|";
      
      for(Int_t i=0; i<ncol && i<nhist-k;i++){
	(output) << "p{1.75cm}|";
      }
      (output) << "} \\hline" << endl;
      t=0;
      for(Int_t j=0; j<nhist;j++){
	if(k<=j && j<nhist && j<k+ncol && t<histo.size()){
	  TString id=histo[t].GetTitle();
	  id.ReplaceAll("N Passed ","");
	  (output) << " &  $" << id  << "$";
	}
	t++;
      }
      
      (output) << " \\\\  \\hline" << std::endl;
      
      for(int i=0; i<=ncuts;i++){
	t=0;
	TString title="Before Skim";
	if(i>0) title=names[i-1];
	(output) << title;
	
	for(Int_t j=0; j<nhist;j++){
	  if(k<=j  && j<nhist && j<k+ncol && t<histo.size() ){
	  if(nevents[t]>0){
	    (output) << " &   " << histo[t].GetBinContent(i+1)/nevents[t];
	  }
	  else{
	    (output) << " &   0";
	  }
	}
	t++;
	}
	(output) << " \\\\  " << std::endl;
      }
      (output) << " \\hline" << endl;
      
      (output) << "\\end{tabular}"                                << std::endl;
      (output) << "\\caption[Accumulative Efficiency ]{Accumulative Efficiency }" << std::endl;
      (output) << "\\end{center}"                                 << std::endl;
      (output) << "\\end{table}"                               << std::endl;
      (output) << "\\normalsize" << std::endl;
      k+=ncol;
      if(k>=nhist){
	flag=false;
      }
      l++;
    }
    (output) << "\\end{landscape}" << std::endl;
  }




  // Makes the Relative Efficiency Table
  if(histo.size()>0){
    int nhist=histo.size();
    int ncuts=histo[0].GetNbinsX()-1;

    bool flag=true;
    int ncol=6;
    int k=0;
    int l=0;
    int t=0;
    (output) << "\\begin{landscape}" << std::endl;
    while(flag){
      (output) << "\\begin{table}[t]" << std::endl;
      (output) << "\\tiny " << std::endl;
      (output) << "\\begin{center}" << std::endl;
      (output) << "\\begin{tabular}{|p{5cm}|";
      
      for(Int_t i=0; i<ncol && i<nhist-k;i++){
	(output) << "p{1.75cm}|";
      }
      (output) << "} \\hline" << endl;
      t=0;
      for(Int_t j=0; j<nhist;j++){
	if(k<=j && j<nhist && j<k+ncol && t<histo.size()){
	  TString id=histo[t].GetTitle();
	  id.ReplaceAll("N Passed ","");
	  (output) << " &  $" << id  << "$";
	}
	t++;
      }
      
      (output) << " \\\\  \\hline" << std::endl;
      
      for(int i=0; i<ncuts;i++){
	t=0;
	TString title=names[i];
	(output) << title;
	
	for(Int_t j=0; j<nhist;j++){
	  if(k<=j  && j<nhist && j<k+ncol && t<histo.size() ){
	  if(nevents[t]>0){
	    (output) << " &   "<< histo[t].GetBinContent(i+2)/histo[t].GetBinContent(i+1);
	  }
	  else{
	    (output) << " &   0";
	  }
	}
	t++;
	}
	(output) << " \\\\  " << std::endl;
      }
      (output) << " \\hline" << endl;
      
      (output) << "\\end{tabular}"                                << std::endl;
      (output) << "\\caption[Relative Efficiency ]{Relative Efficiency }" << std::endl;
      (output) << "\\end{center}"                                 << std::endl;
      (output) << "\\end{table}"                               << std::endl;
      (output) << "\\normalsize" << std::endl;
      k+=ncol;
      if(k>=nhist){
	flag=false;
      }
      l++;
    }
    (output) << "\\end{landscape}" << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Add Plots

  //(output) << "\\clearpage" << std::endl;

  vector<TString> files;
  string dir="./EPS/";
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    cout << "Error(" << errno << ") opening " << dir << endl;
  }
  else{
    while ((dirp = readdir(dp)) != NULL) {
      files.push_back(string(dirp->d_name));
    }
    closedir(dp);
   
  }
  std::sort(files.begin(),files.end(),stringCompare);
  for(int j=0;j<2;j++){
    if(j==0)  (output) << "\\clearpage \n\\section{Nminus-1 Plots}" << std::endl; 
    if(j==1)  (output) << "\\clearpage \n\\section{Assorted Plots}" << std::endl;
    for(int l=0;l<names.size();l++){
    for(int i=0;i<files.size();i++){
      if(files[i].Contains(".eps") && files[i].Contains(Name)){
	//cout << files[i] << endl;
	TString index="_index_";
	index+=l;
	if(j==0 && files[i].Contains("Nminus1_") && files[i].Contains(index)){
	  //cout << files[i] << endl;
	  TString EPSName0=files[i];
	  TString EPSName1=EPSName0;
	  EPSName1.ReplaceAll("Nminus1","Nminus0");
	  TString EPSName2=EPSName0;
	  EPSName2.ReplaceAll("Nminus1","Nminus1dist");
	  TString EPSName3=EPSName0;
	  EPSName3.ReplaceAll("Nminus1","Accumdist");

	  bool is1(false),is2(false),is3(false);
	  for(int k=0;k<files.size();k++){
	    if(files[k].Contains(EPSName1)) is1=true;
	    if(files[k].Contains(EPSName2)) is2=true;
	    if(files[k].Contains(EPSName3)) is3=true;
	  }

	  TString name=files[i];
	  name.ReplaceAll(Name,"");
	  name.ReplaceAll(index,"");
	  name.ReplaceAll("_"," ");
	  name.ReplaceAll("Data.eps","");
	  name.ReplaceAll("Nminus1","");

	  (output) << "\\begin{figure}[p]" << std::endl;
	  (output) << "  \\begin{center} " << std::endl;
	  (output) << "    \\begin{tabular}{c}" << std::endl;
	  (output) << "      \\begin{minipage}[h!]{404pt}" << std::endl;
	  (output) << "        \\begin{center}" << std::endl;
	  (output) << "          \\begin{minipage}[h!]{200pt}" << std::endl;
	  (output) << "            \\begin{center}" << std::endl;
	  (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName0 << "}" << std::endl;
	  if(is1 && (is2 || is3) ) (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName1 << "}" << std::endl;
	  (output) << "            \\end{center}" << std::endl;
	  (output) << "          \\end{minipage}" << std::endl;
	  (output) << "          \\begin{minipage}[h!]{200pt}" << std::endl;
	  (output) << "            \\begin{center}" << std::endl;
	  if(is1 && !(is2 || is3) ) (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName1 << "}" << std::endl;
	  if(is2)(output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName2 << "}" << std::endl;
	  if(is3)(output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName3 << "}" << std::endl;
	  (output) << "            \\end{center}" << std::endl; 
	  (output) << "          \\end{minipage}" << std::endl;
	  (output) << "          \\caption[ "<< name <<"  ]{"<< std::endl;
	  if(is2 || is3)            (output) << "(Upper Left) The N minus 1 plot for the "<< name << " cut. " << std::endl;
	  if(!(is2 || is3))         (output) << "(Left) The N minus 1 plot for the "<< name << " cut. " << std::endl;
	  if(is1 &&  (is2 || is3))  (output) << "(Lower Left) The N minus 0 plot for the "<< name << " cut. " << std::endl;
	  if(is1 && !(is2 || is3))  (output) << "(Right) The N minus 0 plot for the "<< name << " cut. " << std::endl;
	  if(is2)                   (output) << "(Upper Right) The associated plot for the "<< name << " cut, with all cuts except the "<< name << " cut applied. " << std::endl;
	  if(is3)                   (output) << "(Lower Right) The accumulative associated plot for the "<< name << " cut. This plot is  accumulative, meaning that all cuts from the cut flow before the "<< name << " cut have been applied. " << std::endl;
	  (output) << " }" << std::endl; 
	  (output) << "        \\end{center}" << std::endl; 
	  (output) << "      \\end{minipage}" << std::endl;
	  (output) << "    \\end{tabular}" << std::endl; 
	  (output) << "  \\end{center}" << std::endl;
	  (output) << "\\end{figure}" << std::endl; 
	  (output) << "\\clearpage" << std::endl;
	}
	else if(j==1 && files[i].Contains(index) && !files[i].Contains("Nminus1") && !files[i].Contains("Nminus0") && !files[i].Contains("Accumdist") && !files[i].Contains("Nminus1dist") && !files[i].Contains("_log_") && !files[i].Contains("_sig") && !files[i].Contains("_sigtobkg")){
	  cout << "in File loop" << endl;
	  TString EPSName1=files[i];
	  TString EPSName2=EPSName1;
	  EPSName2.ReplaceAll("_index_","_log_index_");
	  TString EPSName3=EPSName2;
	  EPSName3.ReplaceAll("_log_index_","_sig_index_");
	  TString EPSName4=EPSName3;
	  EPSName4.ReplaceAll("_sig_index_","_sigtobkg_index_");
	  TString EPSName5=EPSName3;
          EPSName5.ReplaceAll("_sig_index_","_siglt_index_");
	  TString EPSName6=EPSName3;
          EPSName6.ReplaceAll("_sig_index_","_sigtobkglt_index_");
	  TString EPSName7=EPSName3;
          EPSName7.ReplaceAll("_sig_index_","_siggt_index_");
	  TString EPSName8=EPSName3;
          EPSName8.ReplaceAll("_sig_index_","_sigtobkggt_index_");

	  TString name=files[i];
	  name.ReplaceAll(Name,"");
	  name.ReplaceAll(index,"");
	  name.ReplaceAll("_"," ");
	  name.ReplaceAll("Data.eps","");

	  bool f1(false),f2(false),f3(false),f4(false),f5(false),f6(false),f7(false),f8(false);
	  for(int i=0;i<files.size();i++){
	    if(files[i].Contains(EPSName1)) f1=true;
	    if(files[i].Contains(EPSName2)) f2=true;
	    if(files[i].Contains(EPSName3)) f3=true;
	    if(files[i].Contains(EPSName4)) f4=true;
	    if(files[i].Contains(EPSName5)) f5=true;
	    if(files[i].Contains(EPSName6)) f6=true;
	    if(files[i].Contains(EPSName7)) f7=true;
	    if(files[i].Contains(EPSName8)) f8=true;
	  }

	  if(f1 && !f2 && !f3 && !f4 && !f5 && !f6 && !f7 && !f8){
	    cout << "A" << endl;
	    (output) << "\\begin{figure}[p]" << std::endl;
	    (output) << "\\begin{center} " << std::endl;
	    (output) << "\\begin{tabular}{c}" << std::endl;
	    (output) << "\\begin{minipage}[h!]{400pt}" << std::endl;
	    (output) << "\\begin{center}" << std::endl;
	    (output) << "\\includegraphics*[width=400pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName1 << "}" << std::endl;
	    (output) << "\\caption[ "<< name <<"  ]{ A plot of "<< name << ".}" << std::endl;
	    (output) << "\\end{center}" << std::endl;
	    (output) << "\\end{minipage}" << std::endl;
	    (output) << "\\end{tabular}" << std::endl;
	    (output) << "\\end{center}" << std::endl;
	    (output) << "\\end{figure}" << std::endl;
	    (output) << "\\clearpage" << std::endl;
	  }
	  else{
	    cout << "B" << endl;
	    (output) << "\\begin{figure}[p]" << std::endl;
	    (output) << "  \\begin{center} " << std::endl;
	    (output) << "    \\begin{tabular}{c}" << std::endl;
	    (output) << "      \\begin{minipage}[h!]{404pt}" << std::endl;
	    (output) << "        \\begin{center}" << std::endl;
	    (output) << "          \\begin{minipage}[h!]{200pt}" << std::endl;
	    (output) << "            \\begin{center}" << std::endl;
	    (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName1 << "}" << std::endl;
	    if(f3)   (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName3 << "}" 
			      << std::endl;
	    (output) << "            \\end{center}" << std::endl;
	    (output) << "          \\end{minipage}" << std::endl;
	    (output) << "          \\begin{minipage}[h!]{200pt}" << std::endl;
	    (output) << "            \\begin{center}" << std::endl;
	    (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName2 << "}" << std::endl;
	    if(f4)(output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName4 << "}" << std::endl;
	    (output) << "            \\end{center}" << std::endl; 
	    (output) << "          \\end{minipage}" << std::endl;
	    (output) << "          \\caption[ "<< name <<"  ]{"<< std::endl;
	    if((!f3 && !f4))          (output) << "(Left) The plot of "<< name << ". (Right) A log plot of "<< name << "." << std::endl;
	    if((f3 || f4))            (output) << "(Upper Left) The plot of "<< name << ". (Upper Right) The log plot of "<< name 
					       << "." << "(Lower Left) The significance plot of "<< name 
					       << ". (Lower Right) The signal to background plot of "<< name << "." <<  std::endl;
	    (output) << " }" << std::endl; 
	    (output) << "        \\end{center}" << std::endl; 
	    (output) << "      \\end{minipage}" << std::endl;
	    (output) << "    \\end{tabular}" << std::endl; 
	    (output) << "  \\end{center}" << std::endl;
	    (output) << "\\end{figure}" << std::endl; 
	    (output) << "\\clearpage" << std::endl;
	  }
	  if(f5 && f6 && f7 && f8){
            cout << "B" << endl;
            (output) << "\\begin{figure}[p]" << std::endl;
            (output) << "  \\begin{center} " << std::endl;
            (output) << "    \\begin{tabular}{c}" << std::endl;
            (output) << "      \\begin{minipage}[h!]{404pt}" << std::endl;
            (output) << "        \\begin{center}" << std::endl;
            (output) << "          \\begin{minipage}[h!]{200pt}" << std::endl;
            (output) << "            \\begin{center}" << std::endl;
            (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName5 << "}" << std::endl;
	    (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName7 << "}"
		     << std::endl;
            (output) << "            \\end{center}" << std::endl;
            (output) << "          \\end{minipage}" << std::endl;
            (output) << "          \\begin{minipage}[h!]{200pt}" << std::endl;
            (output) << "            \\begin{center}" << std::endl;
            (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName6 << "}" << std::endl;
            (output) << "              \\includegraphics*[width=200pt,bb=0pt 0pt 567pt 550pt]{./EPS/" << EPSName8 << "}" << std::endl;
            (output) << "            \\end{center}" << std::endl;
            (output) << "          \\end{minipage}" << std::endl;
            (output) << "          \\caption[ "<< name <<"  ]{"<< std::endl;
	    (output) << "(Upper Left) The significance plot (lt) of "<< name << ". (Upper Right) The Purity plot (lt) of "<< name
		     << "." << "(Lower Left) The significance plot (gt) of "<< name
		     << ". (Lower Right) The purity (gt) plot of "<< name << "." <<  std::endl;
            (output) << " }" << std::endl;
            (output) << "        \\end{center}" << std::endl;
            (output) << "      \\end{minipage}" << std::endl;
            (output) << "    \\end{tabular}" << std::endl;
            (output) << "  \\end{center}" << std::endl;
            (output) << "\\end{figure}" << std::endl;
            (output) << "\\clearpage" << std::endl;
	  }
	}
      }
    }
    }
  }

  (output) << "\\end{document}" << std::endl;
  output.close();

  TString cmd1="latex "+Name+".tex -o "+Name+".dvi";
  TString cmd2="dvips "+Name+".dvi -o "+Name+".ps";
  TString cmd3="ps2pdf "+Name+".ps "+Name+".pdf";
  system(cmd1.Data());
  system(cmd1.Data());
  system(cmd1.Data());
  system(cmd2.Data());
  system(cmd3.Data());
}



