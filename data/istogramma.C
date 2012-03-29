void plot(TString input_file, TString histo_title, TString output_file, double min, double max)
{

    TCanvas c1;
	TH1F histo ("histo", histo_title, 100, min, max);
		std::ifstream leggo (input_file,std::ios::in);
		while (!leggo.eof())
		{
			double var;
			leggo >> var;
			histo.Fill (var);
		}
		leggo.close();

	TF1 *func = new TF1 ("func", "[0]*exp(-x)*sqrt(x)", min, max);
	func->SetParameter(0, 5000);
	histo.SetFillColor(0);
	histo.Fit("func", "R");
	histo.Draw();
	/*TCanvas c1;
	histo.SetFillColor(2);
	histo.Draw ();
	c1.Print("isto.gif","gif");
	histo.Fit ("gauss","L");
	histo.Draw ();
	c1.Print("istogauss.gif","gif");
	*/
	
	c1.Print(output_file,"gif");
}
