###

connectedComponents <- function( X, tol.epsilon=1e-3 ){

	N <- ncol( X )
	
	## initialize components matrix
	C <- X
	
	## initialize components vector
	labels <- vector( mode="integer", length=N )
	
	K <- 1 
	labels[1] <- 1
	C[,1] <- X[,1]
	
	# pb <- txtProgressBar( min=0, max=N, style=3 )
	
	## efficient connected component algorithm
	for( n in 2:N ){
		
		assigned <- FALSE
				
		for( k in 1:K ){
			
			distance <- distanceFunction( X[,n], C[,k] )
			
			if( distance < tol.epsilon ){
				
				labels[n] <- k
				assigned <- TRUE
				break
				
			}
			
		}
		
		if( !assigned ){
			
			K <- K + 1
			labels[n] <- K
			C[,K] <- X[,n]
			
		}
		
		# setTxtProgressBar( pb, n )
		
	}
	
	C <- as.matrix( C[,1:K] )
	colnames( C ) <- paste( "mode", 1:K, sep="" )
	
	labels <- as.integer( labels )
	
	output <- list( components=C, labels=labels )
	
	# close( pb )
	
	message( "\nThe algorithm found ", as.character( K ),
	" clusters.\n")
	
	return( output )
		
}
