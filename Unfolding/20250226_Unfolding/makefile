default: TestRunScheme2


TestRunScheme2: Execute ExeUnfoldingHist ExeMatchingEffCorr
	./Execute --Input /data/ALEPH/MC/LEP1MC/LEP1MCMerged.root \
		--Output LEP1MCMerged_Matched.root \
		--Gen tgen --Reco t --Fraction 1.0

	./ExeMatchingEffCorr --Input matchingScheme2/LEP1MCMerged_Matched.root \
			     --MatchingEffName matchingScheme2/MatchingEff.root  --FakeCorrFileName matchingScheme2/FakeCorr.root --MakeMatchingEffCorrFactor true --MakeFakeCorrFactor true

	./ExeMatchingEffCorr --Input matchingScheme2/LEP1MCMerged_Matched.root \
			     --MatchingEffName matchingScheme2/MatchingEff.root --FakeCorrFileName matchingScheme2/FakeCorr.root


TestRun: Execute ExeUnfoldingHist ExeMatchingEffCorr
	./Execute --Input /data/ALEPH/MC/LEP1MC/LEP1MC1994_recons_aftercut-001.root \
		--Output LEP1MC1994_recons_aftercut-001_Matched.root \
		--Gen tgen --Reco t --Fraction 1.00
	./Execute --Input /data/ALEPH/MC/LEP1MC/LEP1MC1994_recons_aftercut-001.root \
		--Output LEP1MC1994_recons_aftercut-001_Matched.root \
		--Gen tgen --Reco t --Fraction 1.00 --MatchingSchemeChoice 1
	./Execute --Input /data/ALEPH/MC/LEP1MC/LEP1MC1994_recons_aftercut-001.root \
		--Output LEP1MC1994_recons_aftercut-001_Matched.root \
		--Gen tgen --Reco t --Fraction 1.00 --MatchingSchemeChoice 3
	./Execute --Input /data/ALEPH/MC/LEP1MC/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
		--Output LEP1MC1994_recons_aftercut-001_Matched.root \
		--Gen tgen --Reco t --Fraction 1.00 --MatchingSchemeChoice 4

	./ExeUnfoldingHist

	./ExeMatchingEffCorr --Input matchingScheme1/LEP1MC1994_recons_aftercut-001_Matched.root \
			     --MatchingEffName matchingScheme1/MatchingEff.root --FakeCorrFileName matchingScheme1/FakeCorr.root --MakeMatchingEffCorrFactor true --MakeFakeCorrFactor true

	./ExeMatchingEffCorr --Input matchingScheme1/LEP1MC1994_recons_aftercut-001_Matched.root \
			     --MatchingEffName matchingScheme1/MatchingEff.root --FakeCorrFileName matchingScheme1/FakeCorr.root

	./ExeMatchingEffCorr --Input matchingScheme2/LEP1MC1994_recons_aftercut-001_Matched.root \
			     --MatchingEffName matchingScheme2/MatchingEff.root  --FakeCorrFileName matchingScheme2/FakeCorr.root --MakeMatchingEffCorrFactor true --MakeFakeCorrFactor true

	./ExeMatchingEffCorr --Input matchingScheme2/LEP1MC1994_recons_aftercut-001_Matched.root \
			     --MatchingEffName matchingScheme2/MatchingEff.root --FakeCorrFileName matchingScheme2/FakeCorr.root

	./ExeMatchingEffCorr --Input matchingScheme3/LEP1MC1994_recons_aftercut-001_Matched.root \
			     --MatchingEffName matchingScheme3/MatchingEff.root --FakeCorrFileName matchingScheme3/FakeCorr.root --MakeMatchingEffCorrFactor true --MakeFakeCorrFactor true

	./ExeMatchingEffCorr --Input matchingScheme3/LEP1MC1994_recons_aftercut-001_Matched.root \
			     --MatchingEffName matchingScheme3/MatchingEff.root --FakeCorrFileName matchingScheme3/FakeCorr.root

	./ExeMatchingEffCorr --Input matchingScheme4/LEP1MC1994_recons_aftercut-001_Matched.root \
			     --MatchingEffName matchingScheme4/MatchingEff.root --FakeCorrFileName matchingScheme4/FakeCorr.root --MakeMatchingEffCorrFactor true --MakeFakeCorrFactor true

	./ExeMatchingEffCorr --Input matchingScheme4/LEP1MC1994_recons_aftercut-001_Matched.root \
			     --MatchingEffName matchingScheme4/MatchingEff.root --FakeCorrFileName matchingScheme4/FakeCorr.root

FullRun: Execute
	bash MatchSome.sh

Execute: MatchEEC.cpp
	g++ MatchEEC.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

ExeUnfoldingHist: createUnfoldingHistograms.C
	g++ createUnfoldingHistograms.C -o ExeUnfoldingHist \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

ExeMatchingEffCorr: matchingEffCorr.cpp
	g++ matchingEffCorr.cpp -o ExeMatchingEffCorr \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o