#-------------------------------------------------------------------------------
#
# Examples of how to run spASM
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
#    Feb 2023
#        - Update to release examples
#
#-------------------------------------------------------------------------------

# Interact with via command line interface
echo -e "Testing command line interface works"
echo -e "\t../target/release/spasm example.fa.gz test.epibed.gz"
../target/release/spasm example.fa.gz test.epibed.gz

# Choose a region to focus on
echo -e "Testing a specific region can be chosen (whole chromosome)"
echo -e "\t../target/release/spasm -r chr1 example.fa.gz test.epibed.gz"
../target/release/spasm --region chr1 example.fa.gz test.epibed.gz

echo -e "Testing a specific region can be chosen (chromosome and region)"
echo -e "\t../target/release/spasm -r chr1:0-100 example.fa.gz test.epibed.gz"
../target/release/spasm --region chr1:0-100 example.fa.gz test.epibed.gz

echo -e "Testing user can combine reads into individual DNA fragments"
echo -e "\t../target/release/spasm --fragment example.fa.gz test.epibed.gz"
../target/release/spasm --merge-mates example.fa.gz test.epibed.gz

echo -e "Testing user can output biscuit ASM file format"
echo -e "\t../target/release/spasm --biscuit example.fa.gz test.epibed.gz"
../target/release/spasm --biscuit example.fa.gz test.epibed.gz

echo -e "Testing various false discovery correction implementations"
echo -e "\t../target/release/spasm --fdr BH example.fa.gz test.epibed.gz"
../target/release/spasm --fdr BH example.fa.gz test.epibed.gz
echo -e "\t../target/release/spasm --fdr BY example.fa.gz test.epibed.gz"
../target/release/spasm --fdr BY example.fa.gz test.epibed.gz
echo -e "\t../target/release/spasm --fdr Bonferroni example.fa.gz test.epibed.gz"
../target/release/spasm --fdr Bonferroni example.fa.gz test.epibed.gz
echo -e "\t../target/release/spasm --fdr Holm example.fa.gz test.epibed.gz"
../target/release/spasm --fdr Holm example.fa.gz test.epibed.gz
echo -e "\t../target/release/spasm --fdr Hochberg example.fa.gz test.epibed.gz"
../target/release/spasm --fdr Hochberg example.fa.gz test.epibed.gz
echo -e "\t../target/release/spasm --fdr None example.fa.gz test.epibed.gz"
../target/release/spasm --fdr None example.fa.gz test.epibed.gz

echo -e "Testing p-value cutoff change"
echo -e "\t../target/release/spasm --pcutoff 0.4 --fdr None example.fa.gz test.epibed.gz"
../target/release/spasm --pcutoff 0.4 --fdr None example.fa.gz test.epibed.gz
