########################### FUNCTION: STITCH FILES ###########################

# given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

stitch_files = function(.results.singles.path, .results.stitched.write.path=.results.singles.path,
                        .name.prefix, .stitch.file.name="stitched_model_fit_results.csv") {
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose names start with .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # grab variable names from first file
  #################
  #### BOOKMARK - THIS LINE GIVES "INCOMPLETE FINAL LINE" ERROR UNLESS YOU MANUALLY INSERT HARD RETURN AT END OF FILE!!
  #################
  names = names( read.csv(keepers[1]) )[-1]
  
  # initialize stitched dataframe
  s = as.data.frame( matrix(nrow=1, ncol=length(names)) )
  names(s) = names
  
  # stitch the files
  for ( i in 1:length(keepers) ) {
    row = read.csv(keepers[i])[1,-1]
    s[i,] = row
  }
  
  write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )
}

# TEST - WORKS!
#.results.singles.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/results"
#.name.prefix = "right_results" 
#stitch_files(.results.singles.path, .results.singles.path, .name.prefix)

