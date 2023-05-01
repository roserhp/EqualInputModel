# Uncomment the next two lines to install the packages LinearAlgebra and Combinatorics.
#import Pkg; Pkg.add("Combinatorics")
#import Pkg; Pkg.add("LinearAlgebra")
using LinearAlgebra
using Combinatorics
B = [[0 0 0 0 0 0 -1 -1 0]
     [1 1 1 0 0 0 0 0 -1]
     [0 0 0 0 0 0 0 -2 -2]
     [0 0 0 0 0 0 -2 0 0]
     #-------------
     [0 0 -1 0 0 0 0 0 0]
     [0 -1 0 0 0 0 0 0 0]
     [-1 0 0 0 0 0 0 0 0]
     #-------------
     [0 0 0 0 0 -1 0 1 0]
     [0 0 0 0 -1 0 0 1 0]
     [0 0 0 -1 0 0 0 1 0]
     #-------------
     [0 0 1 0 0 1 1 0 0]
     [0 1 0 0 1 0 1 0 0] 
     [1 0 0 1 0 0 1 0 0]
     #-------------
     [0 0 0 0 0 1 0 0 1]
     [0 0 0 0 1 0 0 0 1]
     [0 0 0 1 0 0 0 0 1]
     #-------------
     [0 0 -1 0 0 -1 0 0 0]
     [0 -1 0 0 -1 0 0 0 0]
     [-1 0 0 -1 0 0 0 0 0]]

R=collect(powerset([1:1:size(B)[1];],2,8));
C=collect(powerset([1:1:size(B)[2];],3,9));

MSM=0;
MSMInd=[];
for j in C
    for i in R
        if length(i)<length(j)
            L=0;
            for k in j 
                aux = sign.(B[i,k])
                if 1 in aux && -1 in aux 
                    L=L+1;
                end
            if L == length(j) 
                MSM+=1;
                MSMInd=vcat(MSMInd,[[i,j]]);
            end
            end
        end
    end
end
print("There are ", MSM , " submatrices with more columns than rows");
        
                
                
            
            
            
