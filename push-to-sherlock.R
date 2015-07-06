
# this pushes the files from local GitHub repo to Sherlock
# note: will give lots of permissions errors, but it does work

# push just R files
system( "scp -r ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/for-sherlock/r/* mmathur@sherlock:/share/PI/manishad/naEst/r" )


# push all files (slow)
system( "scp -r ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/for-sherlock/* mmathur@sherlock:/share/PI/manishad/naEst" )



sbatch /share/PI/manishad/naEst/test-sbatch.sbatch -p manishad