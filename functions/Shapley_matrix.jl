"""
This function creates a coalitional structure for a cooperative game with `n` players.
That is, it generates a matrix of all possible coalitions, which can later be used to compute the Shapley value and other solution concepts.

Andrey Churkin https://andreychurkin.ru/

"""

function Shapley_matrix(n)

    B = zeros(2^n,n)
            v = 1
            ρ1 = 1
            ρ2 = 0
            P = zeros(n)
            for k = 1:n
                ρ = 2^(n-k)
                ρ2 = ρ2+ρ 
                B[ρ1:ρ2,1] .= v
                ρ1 = ρ1+ρ 
                v = v+1 
                P[k] = ρ2
            end

            for f = 2:n 
                ρρ1 = 2
                ρρ2 = 1
                for k = 1:n
                ρρ2 = floor(Int,P[k]+1)
                B[ρρ1:ρρ2,f] = view(B[ρρ2:2^n,f-1],1:ρρ2-ρρ1+1)
                ρρ1 = floor(Int,P[k]+2)
                end
            end

            B2 = zeros(2^n,n)
            for i = 1:2^n
                for j = 1:n
                    if floor(Int,B[i,j]) != 0
                    B2[i,floor(Int,B[i,j])] = B[i,j]
                    end
                end
            end

    Σ1 = zeros(2^n,1)
    Σ2 = zeros(2^n,1)
    gl = zeros(2^n,2)
    gl0 = zeros(2^n,2)
    gl1 = zeros(2^n,2)
    gl2 = zeros(2^n,2)
    for i = 1:2^n
        gl[i,1] = count(!iszero,B2[i,1:n])
        gl0[i,1] = gl[i,1]
        gl1[i,1] = gl0[i,1]-1
        gl2[i,1] = n-gl[i,1]-1
        Σ1[i,1] = maximum(gl1[i,:])
        Σ2[i,1] = maximum(gl2[i,:])
    end

    u = zeros(2^n,n)
    for j = 1:n
        for i = 1:2^n
            if Int(B2[i,j]) !=  0
                    u[i,j] = 1
            end
        end
    end

    u1 = u.-1

    return u
end
