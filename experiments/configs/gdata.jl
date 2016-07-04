#  Parameter settings for the program gendata.jl
#  The values in N_list and mu_list are not used.  
N_list   = [100]   # Should be of length 1.  The value is not used in gendata.jl
#cp_list = [-0.10, -0.06, 0.00, 0.06, 0.10 ]  # The position of 0.0 is used in hyptest.jl
cp_list = [-0.06, 0.00, 0.06 ]   
mu_list    = [ 6.0 ]   # Should be of length 1.  The value is not used in gendata.jl
mean_list = [25.0, 15.0, 5.0]  # The means of the populations corresponding to values in cp_list.
#  The length of mean_list must be the same as the length of cp_list
stddv = 2.5   # The standard deviation of all generated populations
T = 1000    # Number of trials
