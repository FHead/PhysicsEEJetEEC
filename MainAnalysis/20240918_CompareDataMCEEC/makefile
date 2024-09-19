default: TestRun

TestRun: Execute
	./Execute --InputData ../20240313_SimpleHistogramFiller/Output/LEP1Data1994_Plots.root \
		  --InputMCs  ../20240313_SimpleHistogramFiller/Output/LEP1MC1994_Reco_Plots.root,../20240313_SimpleHistogramFiller/Output/Herwig_Plots.root,../20240313_SimpleHistogramFiller/Output/Sherpa_Plots.root \
		  --InputMCLegs  "Archived MC",Herwig,Sherpa \
		  --InputSysts ../20240313_SimpleHistogramFiller/Output/Sherpa_Plots.root

Execute: CompareDataMCEEC.cpp
	g++ CompareDataMCEEC.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
