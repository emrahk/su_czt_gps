//#include <iostream>
//#include <string>

//int main()
{	
	int nbins = 200;
	float fmin = 0;
	float fmax = 1000;
	int npix;
	float eng, px, py, pz;
	float tote;

	ifstream ifs;	
	ifs.open("output/imaging_output.dat");

	TH1F * old = (TH1F*)gDirectory->Get("eCh");
        if (old) old->Delete();
        string tstr = "Energy Spectrum for Channel ";
        TH1F * eCh = new TH1F("eCh", tstr.c_str(), nbins, fmin, fmax);

	while (!ifs.eof())
	{
		tote = 0;
		ifs >> npix;
		ifs >> eng;
		ifs >> px;
		ifs >> py;
		ifs >> pz;
		tote = eng;
		for (int i=0; i< npix-1; i++)
		{
			ifs >> eng;
			ifs >> eng;
        	        ifs >> px;
        	        ifs >> py;
        	        ifs >> pz;
			tote += eng;
		}
		//cout << tote << endl;
		eCh->Fill(tote);
	}
	eCh->Draw();
	ifs.close();
	
	return 0;
}
