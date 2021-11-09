simcoal=function(N, i){
	n = i # number of sampled haploid individuals
	coal_times = NULL # record the time of coalescences from 'i' lineages to the MRCA
	n_lineages = NULL # record the number of lineages over times

	# loop over the coalescent tree
	while(i>=2){
		P_one_coal = 2*N / (i * (i-1)) # expected time to have one coalescence among i
		coal_time_tmp = rexp(n = 1, rate=1/P_one_coal) # random sampling of this expected time
		coal_times = c(coal_times, coal_time_tmp) # put the sampled coalescence time into the vector
		n_lineages = c(n_lineages, i) # put the current number of lineages into the vector
		i = i-1 # after the coalescence event : one lineage is lost
	}
	T_MRCA = sum(coal_times) # the age of the Most Recent Common Ancestor is equal to the sum of all times of coalescence
	T_Total = sum(coal_times * n_lineages) # the total length of the tree (in generations) is equal to the sum of i*T_i

	# format the final output	
	res = c(T_MRCA,2*N*(1-1/n), T_Total, 2*N*sum(1/(1:(n-1))), cumsum(coal_times))
	names(res) = c('T_MRCA (sim)', 'T_MRCA (exp)', 'T_Total (sim)', 'T_Total (exp)', paste('T', n:2, sep="_"))
	return(res)
}
