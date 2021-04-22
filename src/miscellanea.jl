function real_eigvectors(A)
    n = size(A, 1)

    vals, vecs = eigen(A)
    ix = imag.(vals) .>= 0

    real.(vals[ix]), real.(vecs[:, ix]), isreal.(vals[ix])
end


function rosenbrok(A, B, C)
    if length(B) == 0
        [
            A
            C
        ]
    elseif length(C) == 0
        [
            A B
        ]
    else
        [
            A B
            C zeros(size(C, 1), size(B, 2))
        ]
    end
end


function fixed_modes(A, B, C)
    ra = eigvals(A)
    rn = length(ra)
    n = size(A, 1)
    N = length(B)
    
    minsvd = Inf .* ones(rn)
    
    for i in 1:rn
        for j in 1:2^N
            mask = digits(Bool, j - 1, base=2, pad=N)
            BB = reduce(hcat, B[mask], init=zeros(n, 0))
            CC = reduce(vcat, C[.! mask], init=zeros(0, n))
            s = svdvals(DistributedControl.rosenbrok(A - ra[i] * I, BB, CC))
            if minsvd[i] > s[n]
                minsvd[i] = s[n]
            end
        end
    end
    minsvd
end