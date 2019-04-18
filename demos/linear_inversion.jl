
#define a function that takes a contaminated matrix, A, uses SVD to decompose it, and then gives the solution to Ax=b
#Requires matrix A, vector b,  k is the no. of iterations, and lambda is your dampening effect ,

function linear_invert(A, b, k, l)
    U, Σ, V = svd(A, full=false)
	#number of `columns` in the solution s, or length of diagnol
	c1=length(Σ)
	Si=Diagonal(Σ)
	for i=1:c1
		Si[i,i]= 1/Si[i,i] - 1/Si[i,i]*(l/(l+Si[i,i]^2))^k
    end
	#solution to x using noisy b and iterative T approach
	x1=V*Si*U'*b

end
