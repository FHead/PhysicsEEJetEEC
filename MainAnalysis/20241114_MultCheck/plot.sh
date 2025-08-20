# ------ Plot just the mc ----------
# ./Execute --Input Files/MatchingTest_v5.root,Files/PYTHIA8.root,Files/PYTHIA8_DIRE.root,Files/PYTHIA8_VINCIA.root,Files/Sherpa.root,Files/HERWIG.root \
#         --Output ModelComparisons \
#         --Label "Archived MC Gen Level","PYTHIA8","PYTHIA8 DIRE","PYTHIA8 VINCIA","SHERPA","HERWIG" \
#         --Prefix "0807024" \
#         --DoRatio true \
#         --DoWeight false \
#         --DoReflection false \


./Execute --Input ArchivedMC_E2C_MultCut_0.root,ArchivedMC_E2C_MultCut_10.root,ArchivedMC_E2C_MultCut_20.root,ArchivedMC_E2C_MultCut_30.root,ArchivedMC_E2C_MultCut_40.root \
        --Output ModelComparisons \
        --Label "Archived MC Gen Level","nChargedHadronsHP > 10","nChargedHadronsHP > 20","nChargedHadronsHP > 30","nChargedHadronsHP > 40" \
        --Prefix "11122024" \
        --DoRatio true \
        --DoWeight false \
        --DoReflection false \


# ------ Make the reflected plot ----------
# ./Execute --Input Files/MatchingTest_v5.root \
#         --Output ArchivedMCGen \
#         --Label "Archived MC Gen Level" \
#         --Prefix "0807024" \
#         --DoRatio false \
#         --DoWeight false \
#         --DoReflection true \

# ------ Make the mc data comparison plot ---------
# ./Execute --Input Files/MatchingTest_v5.root,Files/DataTree.root \
#         --Output DataComp2 \
#         --Label "Archived MC Gen Level","Uncorrected Data" \
#         --Prefix "0717024" \
#         --DoRatio true \
#         --DoWeight false \




# ------ Make the multiplicity plot ---------
# ./Execute --Input Files/MatchingTest_v5.root,Files/MatchingTest_v6_nCh10.root,Files/MatchingTest_v6_nCh20.root,Files/MatchingTest_v6_nCh30.root,Files/MatchingTest_v6_nCh35.root,Files/MatchingTest_v6_nCh40.root\
#         --Output MultComp \
#         --Label "Archived MC Gen Level","nChargedHadronsHP > 10","nChargedHadronsHP > 20","nChargedHadronsHP > 30","nChargedHadronsHP > 35","nChargedHadronsHP > 40"\
#         --Prefix "0717024" \
#         --DoRatio true \
#         --DoWeight false \
