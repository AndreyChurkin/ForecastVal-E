### FOR JULIA v1.1
function Shapley_matrix(n) #this function returns "u" Matrix of scenarios, not Shapley values

# n=3 #add number of players
B=zeros(2^n,n)

        #composing Basic Shapley Matrix B
        v=1
        ρ1=1
        ρ2=0
        P=zeros(n)
        for k=1:n #constructing the first column
            ρ=2^(n-k)
            ρ2=ρ2+ρ #global
            B[ρ1:ρ2,1].=v
            ρ1=ρ1+ρ #global
            v=v+1 #global
            P[k]=ρ2
        end

        for f=2:n #constructing the remaining columns
            ρρ1=2
            ρρ2=1
            for k=1:n
            ρρ2=floor(Int,P[k]+1)
            B[ρρ1:ρρ2,f]=view(B[ρρ2:2^n,f-1],1:ρρ2-ρρ1+1)
            ρρ1=floor(Int,P[k]+2)
            end
        end

        #repositioning of the elements
        B2=zeros(2^n,n)
        for i=1:2^n
            for j=1:n
                if floor(Int,B[i,j]) != 0
                B2[i,floor(Int,B[i,j])]=B[i,j]
                end
            end
        end
        # convert(Array{Int64}, B)

#constructing vectors with number of factors
Σ1=zeros(2^n,1)
Σ2=zeros(2^n,1)
gl=zeros(2^n,2)
gl0=zeros(2^n,2)
gl1=zeros(2^n,2)
gl2=zeros(2^n,2)
for i=1:2^n
    gl[i,1]=count(!iszero,B2[i,1:n])
    gl0[i,1]=gl[i,1]
    gl1[i,1]=gl0[i,1]-1
    gl2[i,1]=n-gl[i,1]-1
    Σ1[i,1]=maximum(gl1[i,:])
    Σ2[i,1]=maximum(gl2[i,:])
end

#composing a matrix of a non-missing elements
u=zeros(2^n,n)
for j=1:n
    for i=1:2^n
        if Int(B2[i,j]) != 0
                u[i,j]=1
        end
    end
end
u1=u.-1

#composing a matrix with coeficients M
# M=zeros(2^n,n)
# for i=1:2^n
#     for j=1:n
#         M[i,j]=(u[i,j]*factorial(Σ1[i])*factorial(n-Σ1[i]-1) + u1[i,j]*factorial(Σ2[i])*factorial(n-Σ2[i]-1))/factorial(n)
#     end
# end

# M=view(M,1:2^n-1,1:n) #throwing away the last raw

                    #  ------utility functions examples-------
# V=[5036.47, 15648.38, 33842.62, 28860.13, 10416.08, 23516.25, 18798.95] #X.Tan paper example
# V=[1, 10, 2] #two players gloves example
# V=[290.0, 211.026, 151.078, 211.026, 290.0, 309.026, 290.0] #my-3bus example
# V=[0, 1, 1, 1, 0, 0, 0] #3 playes veto example from coursera
# V=[0,60,80,80,80,40,80,20,0,60,80,40,0,20,0] #game#1 example from 1971 paper (4 players)
# V=[-40,10,50,0,0,0,-10,-50,10,50,40,0,0,-10,-50] #game#2 example from 1971 paper (4 players)


#here is the Shapley vector
# E=V'*M


return u
end
