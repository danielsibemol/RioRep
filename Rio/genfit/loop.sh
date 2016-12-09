#!/bin/bash 

#
#script to genarate and fit with Rio++ several times to test how the fit is going
#

data='../GenFit/FastMC_0.root'
input_toy='../GenFit/Model1_pars3K_toy.txt'
input_ll='../GenFit/Model1_pars3K_fit.txt'
sample_root='Toy.root'


cd ../toy/.
make clean && make
cd ../fit/.
make clean && make 
cd ../GenFit/.

# loop of generation and fit
for i in `seq 501 1000`;

do
  echo 'loop ' $i

#generate
cd ../toy/.
./ToyMCGenerator $input_toy
mv FastMC_0.root ../GenFit/.

echo " Finish Generation "

cd ../GenFit/.

#root -l makeRandomInput.C << EOF
#.q
#EOF

##fit
cd ../fit/
ln -s $data $sample_root
echo $input_ll
./Log_Likelihood_Fitter $input_ll
rm -rf $sample_root
mv fitparameters.txt ../GenFit/FitPars/fitparameters$i.txt

echo " Finish Fit "

done

echo "done"
