
{

	TChain* collisionsChain = new TChain("l1NtupleProducer/L1Tree");

	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_1.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_10.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_11.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_12.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_13.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_2.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_3.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_4.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_5.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_6.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_7.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_8.root");
	collisionsChain->AddFile("/data/0a/digiovan/CMSSW_3_3_3/root/HaloRunRegistry/123893//CollisionsData-run123893_9.root");

	cout << "collisionsChain.C finished building the TChain\n";
}
