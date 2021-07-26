# note: indexing in Julia starts at 1 => array indices of 1 correspond to occupation numbers of 0, etc.
"Prefactor from applying a creation operator"
function m_up(k, n_cut)
    if k-1 >= n_cut
        0
    else
        sqrt(k)
    end
end

"Prefactor from applying an annihilation operator"
function m_down(k)
    if k-1 <= 0
        0
    else
        sqrt(k-1)
    end
end

"Prefactor from applying the S^z operator"
function m_z(k)
    k-1
end

"Convenience function to avoid out-of-bounds indices"
function lower(k, n_cut)
    if k-1 <= 0
        n_cut+2
    else 
        k-1
    end
end

"Adjusted sign-function that returns 1 for m=0"
function sign_new(m)
    if sign(m) == -1
        -1
    else
        1
    end
end


"""
Definition of the state preparation protocol using the approximated Hamiltonian (eq. )

Steps:
1. squeezing on one BEC
2. pi-pulse
3. splitting
4. local squeezing on both BECs
5. pi/8 rotation on b => measure in setting (pi/8, 0)
6. pi/4 rotation on b => measure in setting (3pi/8, 0)
7. pi/4 rotation on a => measure in setting (3pi/8, pi/4)
8. -pi/4 rotation on b => measure in setting (pi/8, pi/4)
"""
function H_approx!(du, u, p, t)
    r, n_cut = p

    @inbounds begin
        for b2 in 1:n_cut+1
            for b1 in 1:n_cut+1
                for a2 in 1:n_cut+1
                    for a1 in 1:n_cut+1
                        if t <= 1  # squeezing on one BEC
                            du[a1, a2, b1, b2] = -1im * 2*r * 
                                (m_up(a1, n_cut) * m_up(a2, n_cut) * u[a1+1, a2+1, b1, b2] + 
                                 m_down(a1) * m_down(a2) * u[lower(a1, n_cut), lower(a2, n_cut), b1, b2])

                        elseif t <= 2  # pi-pulse
                            du[a1, a2, b1, b2] = -1im * pi * 
                                (m_z(a1) * u[a1, a2, b1, b2] + m_z(b1) * u[a1, a2, b1, b2])
                        
                        elseif t <= 3  # splitting
                            du[a1, a2, b1, b2] = -1im * pi/4 *
                            (-1im * m_up(a1, n_cut) * m_down(b1) * u[a1+1, a2, lower(b1, n_cut), b2] + 
                             1im * m_down(a1) * m_up(b1, n_cut) * u[lower(a1, n_cut), a2, b1+1, b2]-
                             1im * m_up(a2, n_cut) * m_down(b2) * u[a1, a2+1, b1, lower(b2, n_cut)] + 
                             1im * m_down(a2) * m_up(b2, n_cut) * u[a1, lower(a2, n_cut), b1, b2+1])
                        
                        elseif t <= 4  # local squeezing on both BECs
                            du[a1, a2, b1, b2] = -1im * r * 
                                (m_up(a1, n_cut) * m_up(a2, n_cut) * u[a1+1, a2+1, b1, b2] + 
                                 m_down(a1) * m_down(a2) * u[lower(a1, n_cut), lower(a2, n_cut), b1, b2] +
                                 m_up(b1, n_cut) * m_up(b2, n_cut) * u[a1, a2, b1+1, b2+1] + 
                                 m_down(b1) * m_down(b2) * u[a1, a2, lower(b1, n_cut), lower(b2, n_cut)])
                        
                        elseif t <= 5  # pi/8 rotation on site B
                            du[a1, a2, b1, b2] = -1im * pi/8 *
                                (-1im * m_up(b1, n_cut) * m_down(b2) * u[a1, a2, b1+1, lower(b2, n_cut)] + 
                                 1im * m_down(b1) * m_up(b2, n_cut) * u[a1, a2, lower(b1, n_cut), b2+1])
                            
                        elseif t <= 6  # pi/4 rotation on B
                            du[a1, a2, b1, b2] = -1im * pi/4 *
                                (-1im * m_up(b1, n_cut) * m_down(b2) * u[a1, a2, b1+1, lower(b2, n_cut)] + 
                                 1im * m_down(b1) * m_up(b2, n_cut) * u[a1, a2, lower(b1, n_cut), b2+1])
                        
                        elseif t <= 7  # pi/4 rotation on A
                            du[a1, a2, b1, b2] = -1im * pi/4 *
                                (-1im * m_up(a1, n_cut) * m_down(a2) * u[a1+1, lower(a2, n_cut), b1, b2] + 
                                 1im * m_down(a1) * m_up(a2, n_cut) * u[lower(a1, n_cut), a2+1, b1, b2])
                        
                        elseif t <= 8  # -pi/4 rotation on B
                            du[a1, a2, b1, b2] = 1im * pi/4 *
                                (-1im * m_up(b1, n_cut) * m_down(b2) * u[a1, a2, b1+1, lower(b2, n_cut)] + 
                                 1im * m_down(b1) * m_up(b2, n_cut) * u[a1, a2, lower(b1, n_cut), b2+1])
                        
                        end

                    end
                end
            end
        end

    end

end

"Modified version of m_up() for the full Hamiltonian"
function m_up0(k, N)
    if k-1 >= N  # the zero-modes are constrained by N, not by n_cut
        0
    else
        sqrt(k)
    end
end

"Modified version of lower() for the full Hamiltonian"
function lower0(k, N)
    if k-1 <= 0
        N+2
    else
        k-1
    end
end

"""
Definition of the state preparation protocol using the full Hamiltonian (eq. )
- the zero-modes a0 and b0 are included, affecting steps 1, 3, and 4 of the state preparation
"""
function H_exact!(du, u, p, t)
    r, n_cut, N = p
    
    @inbounds begin
        for b0 in 1:N+1
            for b2 in 1:n_cut+1
                for b1 in 1:n_cut+1
                    for a0 in 1:N+1
                        for a2 in 1:n_cut+1
                            for a1 in 1:n_cut+1
                                if t <= 1  # squeezing on one BEC
                                    du[a1, a2, a0, b1, b2, b0] = -1im * 2*r/N * 
                                        (m_up(a1, n_cut) * m_up(a2, n_cut) * m_down(a0) * m_down(a0-1) * u[a1+1, a2+1, lower0(lower0(a0, N), N), b1, b2, b0] + 
                                         m_down(a1) * m_down(a2) * m_up0(a0, N) * m_up0(a0+1, N) * u[lower(a1, n_cut), lower(a2, n_cut), a0+2, b1, b2, b0])
                                    
                                elseif t <= 2  # pi-pulse
                                    du[a1, a2, a0, b1, b2, b0] = -1im * pi * 
                                        (m_z(a1) * u[a1, a2, a0, b1, b2, b0] + 
                                         m_z(b1) * u[a1, a2, a0, b1, b2, b0])
                                    
                                elseif t <= 3  # splitting
                                    du[a1, a2, a0, b1, b2, b0] = -1im * pi/4 *
                                        (-1im * m_up(a1, n_cut) * m_down(b1) * u[a1+1, a2, a0, lower(b1, n_cut), b2, b0] + 
                                         1im * m_down(a1) * m_up(b1, n_cut) * u[lower(a1, n_cut), a2, a0, b1+1, b2, b0] -
                                         1im * m_up(a2, n_cut) * m_down(b2) * u[a1, a2+1, a0, b1, lower(b2, n_cut), b0] + 
                                         1im * m_down(a2) * m_up(b2, n_cut) * u[a1, lower(a2, n_cut), a0, b1, b2+1, b0] -
                                         1im * m_up0(a0, N) * m_down(b0) * u[a1, a2, a0+1, b1, b2, lower0(b0, N)] + 
                                         1im * m_down(a0) * m_up0(b0, N) * u[a1, a2, lower0(a0, N), b1, b2, b0+1])
                                    
                                elseif t <= 4  # local squeezing on both BECs
                                    du[a1, a2, a0, b1, b2, b0] = -1im * r/N * 
                                        (m_up(a1, n_cut) * m_up(a2, n_cut) * m_down(a0) * m_down(a0-1) * u[a1+1, a2+1, lower0(lower0(a0, N), N), b1, b2, b0] + 
                                         m_down(a1) * m_down(a2) * m_up0(a0, N) * m_up0(a0+1, N) * u[lower(a1, n_cut), lower(a2, n_cut), a0+2, b1, b2, b0] +
                                         m_up(b1, n_cut) * m_up(b2, n_cut) * m_down(b0) * m_down(b0-1) * u[a1, a2, a0, b1+1, b2+1, lower0(lower0(b0, N), N)] + 
                                         m_down(b1) * m_down(b2) * m_up0(b0, N) * m_up0(b0+1, N) * u[a1, a2, a0, lower(b1, n_cut), lower(b2, n_cut), b0+2])
                                    
                                elseif t <= 5  # pi/8 rotation on B
                                    du[a1, a2, a0, b1, b2, b0] = -1im * pi/8 *
                                        (-1im * m_up(b1, n_cut) * m_down(b2) * u[a1, a2, a0, b1+1, lower(b2, n_cut), b0] + 
                                         1im * m_down(b1) * m_up(b2, n_cut) * u[a1, a2, a0, lower(b1, n_cut), b2+1, b0])
                                    
                                elseif t <= 6  # pi/4 rotation on B
                                    du[a1, a2, a0, b1, b2, b0] = -1im * pi/4 *
                                        (-1im * m_up(b1, n_cut) * m_down(b2) * u[a1, a2, a0, b1+1, lower(b2, n_cut), b0] + 
                                         1im * m_down(b1) * m_up(b2, n_cut) * u[a1, a2, a0, lower(b1, n_cut), b2+1, b0])
                                    
                                elseif t <= 7  # pi/4 rotation on A
                                    du[a1, a2, a0, b1, b2, b0] = -1im * pi/4 *
                                        (-1im * m_up(a1, n_cut) * m_down(a2) * u[a1+1, lower(a2, n_cut), a0, b1, b2, b0] + 
                                         1im * m_down(a1) * m_up(a2, n_cut) * u[lower(a1, n_cut), a2+1, a0, b1, b2, b0])
                                    
                                else  # -pi/4 rotation on B
                                    du[a1, a2, a0, b1, b2, b0] = 1im * pi/4 *
                                        (-1im * m_up(b1, n_cut) * m_down(b2) * u[a1, a2, a0, b1+1, lower(b2, n_cut), b0] + 
                                         1im * m_down(b1) * m_up(b2, n_cut) * u[a1, a2, a0, lower(b1, n_cut), b2+1, b0])
                                end
                            end
                        end
                    end
                end
            end
        end
    end                                                                                                        
    
end