### Method one : exact likelihood following the Bayes formula
bayes = function(observations, dices){
	# dices = c(4, 6, 8, 12, 20)
	for(obs in 1:length(observations)){
		if(obs==1){
			P_dice = rep(1/length(dices), length(dices))
		}else{
			P_dice = P_dice_given_observations
		}
		
		P_observations_given_dice = NULL
		for(i in dices){ # loop over the list of dices
			if(observations[obs] > i){
				# if the observed value is greater than the largest number on the dice, then the proba to observe this value for the dice 'i' is equal to zero
				P_observations_given_dice = c(P_observations_given_dice, 0)
			}else{
				# if the observed value is smaller or equal the largest number on the dice, then the proba to observe this value for the dice 'i' is equal to 1/i
				P_observations_given_dice = c(P_observations_given_dice, 1/i)
			}
		}

		P_observation = 0
		for(i in 1:length(dices)){
			P_observation = P_observation +  (P_dice[i] * P_observations_given_dice[i])
		}

		P_dice_given_observations = ( P_observations_given_dice * P_dice ) / P_observation
	}
	
	# format the output for easier interpretations	
	names(P_dice_given_observations) = dices
	best_model = paste('D', names(which(P_dice_given_observations == max(P_dice_given_observations))), sep='')
	proba_best_model = as.numeric(round(max(P_dice_given_observations),5))
	
	res = NULL
	for(i in P_dice_given_observations){
		res = c(res, round(i,5))
	}
	
	res = c(res, best_model, proba_best_model)
	res = as.data.frame(t(res))

	names(res) = c(names(P_dice_given_observations), 'best_model', 'post_proba')
	
	return(res) # posterior probabilities of all of the dices given the two observations

}


### Method two : approximate Bayesian computation 
abc = function(observations, nSimulations){
	# observations is a data.frame with: each row=a set of observations; each row=the observation for each of the nRolls
	# nSimulations = number os simulations for the ABC inference
	library(abcrf)

	nRolls = ncol(observations)

	## random simulations of each of the 5 alternative models
	D4 = D6 = D8 = D12 = D20 = NULL

	for(i in 1:nSimulations){
		D4 = rbind(D4, sample(1:4, nRolls, replace=T))
		D6 = rbind(D6, sample(1:6, nRolls, replace=T))
		D8 = rbind(D8, sample(1:8, nRolls, replace=T))
		D12 = rbind(D12, sample(1:12, nRolls, replace=T))
		D20 = rbind(D20, sample(1:20, nRolls, replace=T))
	}

	colnames(D4) = paste(rep("obs", nRolls), 1:nRolls, sep="")
	colnames(D6) = paste(rep("obs", nRolls), 1:nRolls, sep="")
	colnames(D8) = paste(rep("obs", nRolls), 1:nRolls, sep="")
	colnames(D12) = paste(rep("obs", nRolls), 1:nRolls, sep="")
	colnames(D20) = paste(rep("obs", nRolls), 1:nRolls, sep="")

	## gathering all simulations in one table names sumsta (for summary statistics)
	sumsta = rbind(D4, D6, D8, D12, D20)

	## producing a vector of model indicators, to associate a row of sumsta to any of the 5 compared models
	modIndex = rep(c('D4', 'D6', 'D8', 'D12', 'D20'), each = nSimulations)

	## model comparison
	mod = abcrf(modIndex~., data = data.frame(modIndex, sumsta), ntree = 1000, paral = T, ncores = 3)

	## model prediction
	predicted_dice = predict(mod, observations, training=data.frame(modIndex, sumsta), ntree = 1000, paral = T, ncores=3)
	return(data.frame(best_model=predicted_dice$allocation, post_proba=predicted_dice$post.prob))
}


# ABC: plot the joint distribution for the dice example
dice_sumStats_plot = function(nSampling, nSimulations){
	# ABC: plot the joint distribution for the dice example
	library(MASS)
	obs4 = obs6 = obs8 = obs12 = obs20 = matrix(NA, ncol=2, nrow=nSimulations)
	
	for(i in 1:nSimulations){
		tmp4 = sample(1:4, nSampling, replace=T)
		tmp6 = sample(1:6, nSampling, replace=T)
		tmp8 = sample(1:8, nSampling, replace=T)
		tmp12 = sample(1:12, nSampling, replace=T)
		tmp20 = sample(1:20, nSampling, replace=T)

		obs4[i,] = c(mean(tmp4), sd(tmp4))
		obs6[i,] = c(mean(tmp6), sd(tmp6))
		obs8[i,] = c(mean(tmp8), sd(tmp8))
		obs12[i,] = c(mean(tmp12), sd(tmp12))
		obs20[i,] = c(mean(tmp20), sd(tmp20))
	}

	tmp4 = kde2d(obs4[,1], obs4[,2])
	tmp6 = kde2d(obs6[,1], obs6[,2])
	tmp8 = kde2d(obs8[,1], obs8[,2])
	tmp12 = kde2d(obs12[,1], obs12[,2])
	tmp20 = kde2d(obs20[,1], obs20[,2])

	plot(range(obs4[,1], obs6[,1], obs8[,1], obs12[,1], obs20[,1]), range(obs4[,2], obs6[,2], obs8[,2], obs12[,2], obs20[,2]), col="white", xlab = "mean observation", ylab = "standard deviation")

	contour(tmp20, add=T, col=rainbow(5)[5], lwd=2)
	contour(tmp12, add=T, col=rainbow(5)[4], lwd=2)
	contour(tmp8, add=T, col=rainbow(5)[3], lwd=2)
	contour(tmp6, add=T, col=rainbow(5)[2], lwd=2)
	contour(tmp4, add=T, col=rainbow(5)[1], lwd=2)
}




# plot the percentage of inferences among nReplicates for which we do not support the exact correct model
plot_bayes_error = function(good_dice=8, nReplicates=1000){
	# plot the percentage of inferences among nReplicates for which we do not support the exact correct model
	## example:
	# dev.new(width=13.4, height=4.2)
	# par(mfrow=c(1,5), mar=c(4.5,3.5,3.5,1.5))
	# plot_bayes_error(4, 1000)
	# plot_bayes_error(6, 1000)
	# plot_bayes_error(8, 1000)
	# plot_bayes_error(12, 1000)
	# plot_bayes_error(20, 1000)
	dices = c(4,6,8,12,20)
	
	best_N_2 = best_N_5 = best_N_10 = best_N_100 = best_N_1000 = NULL
	for(i in 1:nReplicates){
		N=2
		best_N_2 = c(best_N_2, as.vector(bayes(observations=sample(1:good_dice, N, replace=T), dices=dices)$best_model))
		N=5
		best_N_5 = c(best_N_5, as.vector(bayes(observations=sample(1:good_dice, N, replace=T), dices=dices)$best_model))
		N=10
		best_N_10 = c(best_N_10, as.vector(bayes(observations=sample(1:good_dice, N, replace=T), dices=dices)$best_model))
		N=100
		best_N_100 = c(best_N_100, as.vector(bayes(observations=sample(1:good_dice, N, replace=T), dices=dices)$best_model))
		N=1000
		best_N_1000 = c(best_N_1000, as.vector(bayes(observations=sample(1:good_dice, N, replace=T), dices=dices)$best_model))
	}

	error_N_2 = length(which(best_N_2 != paste('D', good_dice, sep=''))) / nReplicates 
	error_N_5 = length(which(best_N_5 != paste('D', good_dice, sep=''))) / nReplicates 
	error_N_10 = length(which(best_N_10 != paste('D', good_dice, sep=''))) / nReplicates 
	error_N_100 = length(which(best_N_100 != paste('D', good_dice, sep=''))) / nReplicates 
	error_N_1000 = length(which(best_N_1000 != paste('D', good_dice, sep=''))) / nReplicates 

	barplot(c(error_N_2, error_N_5, error_N_10, error_N_100, error_N_1000), names = c('#obs=2', '#obs=5', '#obs=10', '#obs=100', '#obs=1000'), main = paste('good dice = ', good_dice,sep=''), ylim=c(0, 0.6), ylab='error rate', cex.lab=1.1)
	abline(h=0.05, col="red", lty=2)
	
	summary = matrix(NA, ncol=5, nrow=5)
	for(i in 1:5){ summary[1,i] = length(which(best_N_2 == paste('D', dices[i], sep=''))) }
	for(i in 1:5){ summary[2,i] = length(which(best_N_5 == paste('D', dices[i], sep=''))) }
	for(i in 1:5){ summary[3,i] = length(which(best_N_10 == paste('D', dices[i], sep=''))) }
	for(i in 1:5){ summary[4,i] = length(which(best_N_100 == paste('D', dices[i], sep=''))) }
	for(i in 1:5){ summary[5,i] = length(which(best_N_1000 == paste('D', dices[i], sep=''))) }
	
	colnames(summary) = paste('D', dices, sep="")
	rownames(summary) = c('N2', 'N5', 'N10', 'N100', 'N1000')

	cat('\n')
	cat('Distribution of supported models\n')	
	cat(paste('(made over ', nReplicates, ' replicates)\n', sep=''))
	cat(paste('The correct dice is : ', good_dice, '\n', sep=''))
	cat('\n')
	print(summary)
	cat('\n')
	return(summary)
}

## comparisons between the two methods
#real_dice = 20
#nObservations = 1000
#nRolls = 10
#
#observations = NULL
#for(i in 1:nRolls){
#	observations = cbind(observations, sample(1:real_dice, nObservations, replace=T))
#}
#
#observations = as.data.frame(observations)
#colnames(observations) = paste(rep("obs", nRolls), 1:nRolls, sep="")
#
#
#res_bayes = NULL # store the relative posterior probabilities of the 5 dice
## loop over the observations: compute the relative posterior probabilities by using the Bayes formula
#for(i in 1:nrow(observations)){
#	res_bayes = rbind(res_bayes, bayes(observations[i,], dices=c(4,6,8,12,20)))
#}
#
#res_abc = abc(observations, nSimulations=5000)
#
# table(res_bayes$best_model)
# table(res_abc$best_model)
#
#



example_ABC=function(){
	# build the reference table
	N = 10 # number of observed values from the dice to guess, rolled N times

	nSimulations = 10000 # number of random simulations to perform under each model

	simulations_die_4 = as.data.frame(matrix(NA, ncol=N, nrow=nSimulations)) # N columns (corresponding to the number of observations) and nSimulations rows (the number of independent simulations)
	simulations_die_6 = as.data.frame(matrix(NA, ncol=N, nrow=nSimulations))
	simulations_die_8 = as.data.frame(matrix(NA, ncol=N, nrow=nSimulations))
	simulations_die_12 = as.data.frame(matrix(NA, ncol=N, nrow=nSimulations))
	simulations_die_20 = as.data.frame(matrix(NA, ncol=N, nrow=nSimulations))

	# loop in order to perform nSimulations simulations
	for(i in 1:nSimulations){
		simulations_die_4[i,] = sample(1:4, N, replace=T)
		simulations_die_6[i,] = sample(1:6, N, replace=T)
		simulations_die_8[i,] = sample(1:8, N, replace=T)
		simulations_die_12[i,] = sample(1:12, N, replace=T)
		simulations_die_20[i,] = sample(1:20, N, replace=T)
	}

	# just rename the columns, because I like it
	colnames(simulations_die_4) = paste( rep('obs', N), 1:N, sep='')
	colnames(simulations_die_6) = paste( rep('obs', N), 1:N, sep='')
	colnames(simulations_die_8) = paste( rep('obs', N), 1:N, sep='')
	colnames(simulations_die_12) = paste( rep('obs', N), 1:N, sep='')
	colnames(simulations_die_20) = paste( rep('obs', N), 1:N, sep='')

	# combine all simulations in a single large reference table
	all_simulations = rbind(simulations_die_4, simulations_die_6, simulations_die_8, simulations_die_12, simulations_die_20) 

	# a vector used to labellize each of the simulations, usefull to help the random forest to recognize the different models
	modIndexes = as.factor(rep(c('D4', 'D6', 'D8', 'D12', 'D20'), each = nSimulations))

	# training the forest
	mod = abcrf(modIndexes~., data = data.frame(modIndexes, all_simulations), ntree = 1000, paral = T, ncores = 3)

	# lets try the forest now, to see if it works
	nSimulations = 10000 # here, we'll try it nSimulations times, in order to do some stats about the error rate of the forest

	simulations_die_8 = as.data.frame(matrix(NA, ncol=N, nrow=nSimulations))

	for(i in 1:nSimulations){
		observations = matrix(sample(1:8, N, replace=T), nrow=1, ncol=N)
		observations = as.data.frame(observations)
		simulations_die_8[i,] = observations
	}
		
	colnames(simulations_die_8) = paste( rep('obs', N), 1:N, sep='')

	# now, use the previously trained forest in order to predict which model is the best for the nSimulations simulations
	predicted_dice = predict(mod, simulations_die_8, training=data.frame(modIndexes, all_simulations), ntree = 1000, paral = T, ncores=2)
}

