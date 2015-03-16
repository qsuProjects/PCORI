# this pushes the files from local GitHub repo to Sherlock

# push all the individual files
system( "scp ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/* mmathur@sherlock:/share/PI/manishad/genCov" )

# push the sbatch folder
system( "scp -r ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/sbatch_files mmathur@sherlock:/share/PI/manishad/genCov" )

# run them
# sbatch 1.sbatch -p manishad