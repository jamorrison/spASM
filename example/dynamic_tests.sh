#-------------------------------------------------------------------------------
#
# Dynamic testing of the GVSU CIS 611 end of semester project
#
# Note, this testing script only works on UNIX systems. It's possible it could
# be made to work on Windows, but effort would need to be made to ensure this.
#
# Created by:
#     Jacob Morrison
#
# Date created:
#     Nov 2022
#
# Release notes:
#
#    Nov 2022
#        - Initial creation
#
#-------------------------------------------------------------------------------

# Interact with via command line interface
echo -e "Testing command line interface works"
echo -e "\t../target/release/cis_611_project test.epibed.gz"
../target/release/cis_611_project test.epibed.gz

# Choose a region to focus on
echo -e "Testing a specific region can be chosen (whole chromosome)"
echo -e "\t../target/release/cis_611_project -r \"chr1\" test.epibed.gz"
../target/release/cis_611_project --region "chr1" test.epibed.gz

echo -e "Testing a specific region can be chosen (chromosome and region)"
echo -e "\t../target/release/cis_611_project -r \"chr1:0-100\" test.epibed.gz"
../target/release/cis_611_project --region "chr1:0-100" test.epibed.gz

echo -e "Testing user can combine reads into individual DNA fragments"
echo -e "\t../target/release/cis_611_project --fragment test.epibed.gz"
../target/release/cis_611_project --fragment test.epibed.gz

echo -e "Testing user can output biscuit ASM file format"
echo -e "\t../target/release/cis_611_project --biscuit test.epibed.gz"
../target/release/cis_611_project --biscuit test.epibed.gz

echo -e "Testing various false discovery correction implementations"
echo -e "\t../target/release/cis_611_project --fdr BH test.epibed.gz"
../target/release/cis_611_project --fdr BH test.epibed.gz
echo -e "\t../target/release/cis_611_project --fdr BY test.epibed.gz"
../target/release/cis_611_project --fdr BY test.epibed.gz
echo -e "\t../target/release/cis_611_project --fdr Bonferroni test.epibed.gz"
../target/release/cis_611_project --fdr Bonferroni test.epibed.gz
echo -e "\t../target/release/cis_611_project --fdr Holm test.epibed.gz"
../target/release/cis_611_project --fdr Holm test.epibed.gz
echo -e "\t../target/release/cis_611_project --fdr Hochberg test.epibed.gz"
../target/release/cis_611_project --fdr Hochberg test.epibed.gz
echo -e "\t../target/release/cis_611_project --fdr None test.epibed.gz"
../target/release/cis_611_project --fdr None test.epibed.gz

echo -e "Testing p-value cutoff change"
echo -e "\t../target/release/cis_611_project --pcutoff 0.4 --fdr None test.epibed.gz"
../target/release/cis_611_project --pcutoff 0.4 --fdr None test.epibed.gz
