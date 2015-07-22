########################### FUNCTION: STITCH FILES ###########################

# given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

stitch_files = function(.results.singles.path, .results.stitched.write.path=.results.singles.path,
                        .name.prefix, .stitch.file.name="stitched_model_fit_results.csv") {
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose names includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # grab variable names from first file
  # SPECIFIC TO THIS PROJECT - RENAME THE FIRST COLUMN "STAT"
  names = c("stat", names( read.csv(keepers[1]) )[-1])
  
  # initialize stitched dataframe
  s = as.data.frame( matrix(nrow=1, ncol=length(names)) )
  names(s) = names

  # stitch the files
  for ( i in 1:length(keepers) ) {
    new.chunk = read.csv(keepers[i])
    names(new.chunk)[1] = "stat"  # SPECIFIC TO THIS PROJECT
    s = rbind(s, new.chunk)
    
    #row = read.csv(keepers[i])[1,-1]
    #s[i,] = row
  }
  
  s = s[-1,]  # delete annoying NA row
  write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )
}

# TEST - WORKS!
#.results.singles.path = "~/Dropbox/QSU/Mathur/MY_PAPERS/TVC/Code/git_repo/time-varying-coeffs/2015-03-25_local_test/results"
#.name.prefix = "right_results" 
#stitch_files(.results.singles.path, .results.singles.path, .name.prefix)

# # LOCAL TEST
# stitch_files( "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test",
#               "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test",
#               .name.prefix="results",
#               .stitch.file.name="stitched.csv"
# )


# stitch files from Sherlock
stitch_files( "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/unstitched",
              "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/stitched",
                    .name.prefix="NA_frailty",
              .stitch.file.name="NA_frailty_stitched.csv"
              )

stitch_files( "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/unstitched",
              "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/stitched",
              .name.prefix="NA_naive",
              .stitch.file.name="NA_naive_stitched.csv"
)

stitch_files( "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/unstitched",
              "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/stitched",
              .name.prefix="NA_log-t",
              .stitch.file.name="NA_log-t_stitched.csv"
)

stitch_files( "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/unstitched",
              "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/stitched",
              .name.prefix="complete.case",
              .stitch.file.name="CC_stitched.csv"
)

stitch_files( "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/unstitched",
              "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/stitched",
              .name.prefix="full",
              .stitch.file.name="full_stitched.csv"
)






