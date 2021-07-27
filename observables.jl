function approach1(state, n_cut)  
    chsh1 = 0.0;
    chsh1_norm = 0.0;
    chsh2 = 0.0;
    chsh2_norm = 0.0;
    chsh3 = 0.0;
    chsh3_norm = 0.0;
    chsh4 = 0.0;
    chsh4_norm = 0.0;

    prob1 = 0.0;
    prob2 = 0.0;
    prob3 = 0.0;
    prob4 = 0.0;

    @inbounds begin
        for b2 in 1:n_cut+1
            for b1 in 1:n_cut+1
                for a2 in 1:n_cut+1
                    for a1 in 1:n_cut+1
                        prob1 = abs2(state[1][a1, a2, b1, b2]);
                        chsh1 += prob1 * (m_z(a1) - m_z(a2)) * (m_z(b1) - m_z(b2));
                        chsh1_norm += prob1 * (m_z(a1) + m_z(a2)) * (m_z(b1) + m_z(b2));
                        
                        prob2 = abs2(state[2][a1, a2, b1, b2]);
                        chsh2 += prob2 * (m_z(a1) - m_z(a2)) * (m_z(b1) - m_z(b2));
                        chsh2_norm += prob2 * (m_z(a1) + m_z(a2)) * (m_z(b1) + m_z(b2));
                        
                        prob3 = abs2(state[3][a1, a2, b1, b2]);
                        chsh3 += prob3 * (m_z(a1) - m_z(a2)) * (m_z(b1) - m_z(b2));
                        chsh3_norm += prob3 * (m_z(a1) + m_z(a2)) * (m_z(b1) + m_z(b2));
                        
                        prob4 = abs2(state[4][a1, a2, b1, b2]);
                        chsh4 += prob4 * (m_z(a1) - m_z(a2)) * (m_z(b1) - m_z(b2));
                        chsh4_norm += prob4 * (m_z(a1) + m_z(a2)) * (m_z(b1) + m_z(b2));
                    end
                end
            end
        end
    end
    return abs(-chsh1/chsh1_norm + chsh2/chsh2_norm + chsh3/chsh3_norm + chsh4/chsh4_norm)
end


function approach23(state, n_cut) 
    chsh1_approach3 = 0.0;
    vacuumAB1 = abs2(state[1][1,1,1,1]);
    chsh2_approach3 = 0.0;
    vacuumAB2 = abs2(state[2][1,1,1,1]);
    chsh3_approach3 = 0.0;
    vacuumAB3 = abs2(state[3][1,1,1,1]);
    chsh4_approach3 = 0.0;
    vacuumAB4 = abs2(state[4][1,1,1,1]);
    
    chsh1_approach2 = 0.0;
    chsh2_approach2 = 0.0;
    chsh3_approach2 = 0.0;
    chsh4_approach2 = 0.0;
    
    vacuumA1 = 0.0;
    vacuumA2 = 0.0;
    vacuumA3 = 0.0;
    vacuumA4 = 0.0;
    vacuumB1 = 0.0;
    vacuumB2 = 0.0;
    vacuumB3 = 0.0;
    vacuumB4 = 0.0;

    prob1 = 0.0;
    prob2 = 0.0;
    prob3 = 0.0;
    prob4 = 0.0;
    
    @inbounds begin
        for b2 in 1:n_cut+1
            for b1 in 1:n_cut+1
                for a2 in 1:n_cut+1
                    for a1 in 1:n_cut+1
                        
                        prob1 = abs2(state[1][a1, a2, b1, b2]);
                        prob2 = abs2(state[2][a1, a2, b1, b2]);
                        prob3 = abs2(state[3][a1, a2, b1, b2]);
                        prob4 = abs2(state[4][a1, a2, b1, b2]);
                        
                        chsh1_approach2 += prob1 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                        chsh2_approach2 += prob2 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                        chsh3_approach2 += prob3 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                        chsh4_approach2 += prob4 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                        
                        if !(a1 == a2 == 1) && !(b1 == b2 == 1)
                            chsh1_approach3 += prob1 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                            chsh2_approach3 += prob2 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                            chsh3_approach3 += prob3 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                            chsh4_approach3 += prob4 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                        end
                        
                        if (b1 == 1) && (b2 == 1)
                            vacuumB1 += prob1;
                            vacuumB2 += prob2;
                            vacuumB3 += prob3;
                            vacuumB4 += prob4;
                        end
                        
                        if (a1 == 1) && (a2 == 1)
                            vacuumA1 += prob1;
                            vacuumA2 += prob2;
                            vacuumA3 += prob3;
                            vacuumA4 += prob4;
                        end
                    end
                end
            end
        end
        
    end
    
    chsh_approach3 = abs(- chsh1_approach3 / (1 - vacuumA1 - vacuumB1 + vacuumAB1) 
                         + chsh2_approach3 / (1 - vacuumA2 - vacuumB2 + vacuumAB2) 
                         + chsh3_approach3 / (1 - vacuumA3 - vacuumB3 + vacuumAB3) 
                         + chsh4_approach3 / (1 - vacuumA4 - vacuumB4 + vacuumAB4))
    chsh_approach2 = abs(-chsh1_approach2 + chsh2_approach2 + chsh3_approach2 + chsh4_approach2)
    
    return chsh_approach2, chsh_approach3
end


function G_loss(k, l, n_cut, γ)
    sum = 0
    for i in 0:n_cut
        for j in 0:i-1
            sum += binomial(k, k-j) * binomial(l, l-i) * γ^(i+j) * (1 - γ)^(k-j+l-i)
        end
    end
    return 1 - 2*sum
end


function approach23_loss!(state, n_cut, γ, loss_array)
    chsh1 = 0.0;
    vacuumAsgnB1 = 0.0;
    vacuumBsgnA1 = 0.0;
    vacuumA1 = 0.0;
    vacuumB1 = 0.0;
    vacuumAB1 = 0.0;
    chsh2 = 0.0;
    vacuumAsgnB2 = 0.0;
    vacuumBsgnA2 = 0.0;
    vacuumA2 = 0.0;
    vacuumB2 = 0.0;
    vacuumAB2 = 0.0;
    chsh3 = 0.0;
    vacuumAsgnB3 = 0.0;
    vacuumBsgnA3 = 0.0;
    vacuumA3 = 0.0;
    vacuumB3 = 0.0;
    vacuumAB3 = 0.0;
    chsh4 = 0.0;
    vacuumAsgnB4 = 0.0;
    vacuumBsgnA4 = 0.0;
    vacuumA4 = 0.0;
    vacuumB4 = 0.0;
    vacuumAB4 = 0.0;

    prob1 = 0.0;
    prob2 = 0.0;
    prob3 = 0.0;
    prob4 = 0.0;
    
    lossA = 0.0;
    lossB = 0.0;
    
    @inbounds begin
        for i in 1:n_cut+1
            for j in 1:n_cut+1
                loss_array[i, j] = G_loss(i-1, j-1, n_cut, γ);
            end
        end

        for b2 in 1:n_cut+1
            for b1 in 1:n_cut+1
                for a2 in 1:n_cut+1
                    for a1 in 1:n_cut+1
                        
                        lossA = loss_array[a1-1, a2-1];
                        lossB = loss_array[b1-1, b2-1];
                        
                        prob1 = abs2(state[1][a1, a2, b1, b2]);
                        chsh1 += prob1 * lossA * lossB;
                        vacuumAsgnB1 += prob1 * (1-γ)^(a1-1+a2-1) * lossB;
                        vacuumBsgnA1 += prob1 * lossA * (1-γ)^(b1-1+b2-1);
                        vacuumA1 += prob1 * (1-γ)^(a1-1+a2-1);
                        vacuumB1 += prob1 * (1-γ)^(b1-1+b2-1);
                        vacuumAB1 += prob1 * (1-γ)^(a1-1+a2-1) * (1-γ)^(b1-1+b2-1)
                        
                        prob2 = abs2(state[2][a1, a2, b1, b2]);
                        chsh2 += prob2 * lossA * lossB;
                        vacuumAsgnB2 += prob2 * (1-γ)^(a1-1+a2-1) * lossB;
                        vacuumBsgnA2 += prob2 * lossA * (1-γ)^(b1-1+b2-1);
                        vacuumA2 += prob2 * (1-γ)^(a1-1+a2-1);
                        vacuumB2 += prob2 * (1-γ)^(b1-1+b2-1);
                        vacuumAB2 += prob2 * (1-γ)^(a1-1+a2-1) * (1-γ)^(b1-1+b2-1);
                        
                        prob3 = abs2(state[3][a1, a2, b1, b2]);
                        chsh3 += prob3 * lossA * lossB;
                        vacuumAsgnB3 += prob3 * (1-γ)^(a1-1+a2-1) * lossB;
                        vacuumBsgnA3 += prob3 * lossA * (1-γ)^(b1-1+b2-1);
                        vacuumA3 += prob3 * (1-γ)^(a1-1+a2-1);
                        vacuumB3 += prob3 * (1-γ)^(b1-1+b2-1);
                        vacuumAB3 += prob3 * (1-γ)^(a1-1+a2-1) * (1-γ)^(b1-1+b2-1);
                        
                        prob4 = abs2(state[4][a1, a2, b1, b2]);
                        chsh4 += prob4 * lossA * lossB;
                        vacuumAsgnB4 += prob4 * (1-γ)^(a1-1+a2-1) * lossB;
                        vacuumBsgnA4 += prob4 * lossA * (1-γ)^(b1-1+b2-1);
                        vacuumA4 += prob4 * (1-γ)^(a1-1+a2-1);
                        vacuumB4 += prob4 * (1-γ)^(b1-1+b2-1);
                        vacuumAB4 += prob4 * (1-γ)^(a1-1+a2-1) * (1-γ)^(b1-1+b2-1);
                    end
                end
            end
        end
    end
    chsh_approach3 = abs(-(chsh1 - vacuumAsgnB1 - vacuumBsgnA1 + vacuumAB1)/(1 - vacuumA1 - vacuumB1 + vacuumAB1) + 
                          (chsh2 - vacuumAsgnB2 - vacuumBsgnA2 + vacuumAB2)/(1 - vacuumA2 - vacuumB2 + vacuumAB2) + 
                          (chsh3 - vacuumAsgnB3 - vacuumBsgnA3 + vacuumAB3)/(1 - vacuumA3 - vacuumB3 + vacuumAB3) + 
                          (chsh4 - vacuumAsgnB4 - vacuumBsgnA4 + vacuumAB4)/(1 - vacuumA4 - vacuumB4 + vacuumAB4))
    chsh_approach2 = abs(-chsh1 + chsh2 + chsh3 + chsh4)
    return chsh_approach2, chsh_approach3
end

function approach1_exact(state, n_cut, N)  
    chsh1 = 0.0;
    chsh1_norm = 0.0;
    chsh2 = 0.0;
    chsh2_norm = 0.0;
    chsh3 = 0.0;
    chsh3_norm = 0.0;
    chsh4 = 0.0;
    chsh4_norm = 0.0;

    prob1 = 0.0;
    prob2 = 0.0;
    prob3 = 0.0;
    prob4 = 0.0;

    @inbounds begin
        for b0 in 1:N+1
            for b2 in 1:n_cut+1
                for b1 in 1:n_cut+1
                    for a0 in 1:N+1
                        for a2 in 1:n_cut+1
                            for a1 in 1:n_cut+1
                                prob1 = abs2(state[1][a1, a2, a0, b1, b2, b0]);
                                chsh1 += prob1 * (m_z(a1) - m_z(a2)) * (m_z(b1) - m_z(b2));
                                chsh1_norm += prob1 * (m_z(a1) + m_z(a2)) * (m_z(b1) + m_z(b2));

                                prob2 = abs2(state[2][a1, a2, a0, b1, b2, b0]);
                                chsh2 += prob2 * (m_z(a1) - m_z(a2)) * (m_z(b1) - m_z(b2));
                                chsh2_norm += prob2 * (m_z(a1) + m_z(a2)) * (m_z(b1) + m_z(b2));

                                prob3 = abs2(state[3][a1, a2, a0, b1, b2, b0]);
                                chsh3 += prob3 * (m_z(a1) - m_z(a2)) * (m_z(b1) - m_z(b2));
                                chsh3_norm += prob3 * (m_z(a1) + m_z(a2)) * (m_z(b1) + m_z(b2));

                                prob4 = abs2(state[4][a1, a2, a0, b1, b2, b0]);
                                chsh4 += prob4 * (m_z(a1) - m_z(a2)) * (m_z(b1) - m_z(b2));
                                chsh4_norm += prob4 * (m_z(a1) + m_z(a2)) * (m_z(b1) + m_z(b2));
                            end
                        end
                    end
                end
            end
        end
    end
    return abs(-chsh1/chsh1_norm + chsh2/chsh2_norm + chsh3/chsh3_norm + chsh4/chsh4_norm)
end


function approach23_exact(state, n_cut, N) 
   
    chsh1_approach3 = 0.0;
    chsh2_approach3 = 0.0;
    chsh3_approach3 = 0.0;
    chsh4_approach3 = 0.0;
    
    
    chsh1_approach2 = 0.0;
    chsh2_approach2 = 0.0;
    chsh3_approach2 = 0.0;
    chsh4_approach2 = 0.0;
    
    vacuumAB1 = 0;
    vacuumAB2 = 0;
    vacuumAB3 = 0;
    vacuumAB4 = 0;
    
    vacuumA1 = 0.0;
    vacuumA2 = 0.0;
    vacuumA3 = 0.0;
    vacuumA4 = 0.0;
    
    vacuumB1 = 0.0;
    vacuumB2 = 0.0;
    vacuumB3 = 0.0;
    vacuumB4 = 0.0;

    prob1 = 0.0;
    prob2 = 0.0;
    prob3 = 0.0;
    prob4 = 0.0;
    
    @inbounds begin
        for b0 in 1:N+1
            for b2 in 1:n_cut+1
                for b1 in 1:n_cut+1
                    for a0 in 1:N+1
                        for a2 in 1:n_cut+1
                            for a1 in 1:n_cut+1

                                prob1 = abs2(state[1][a1, a2, a0, b1, b2, b0]);
                                prob2 = abs2(state[2][a1, a2, a0, b1, b2, b0]);
                                prob3 = abs2(state[3][a1, a2, a0, b1, b2, b0]);
                                prob4 = abs2(state[4][a1, a2, a0, b1, b2, b0]);
                                

                                chsh1_approach2 += prob1 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                                chsh2_approach2 += prob2 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                                chsh3_approach2 += prob3 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                                chsh4_approach2 += prob4 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                                
                                if (a1 == a2 == b1 == b2 == 1)
                                    vacuumAB1 += prob1;
                                    vacuumAB2 += prob2;
                                    vacuumAB3 += prob3;
                                    vacuumAB4 += prob4;
                                end

                                if !(a1 == a2 == 1) && !(b1 == b2 == 1)
                                    chsh1_approach3 += prob1 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                                    chsh2_approach3 += prob2 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                                    chsh3_approach3 += prob3 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                                    chsh4_approach3 += prob4 * sign_new(m_z(a1) - m_z(a2)) * sign_new(m_z(b1) - m_z(b2));
                                end

                                if (b1 == 1) && (b2 == 1)
                                    vacuumB1 += prob1;
                                    vacuumB2 += prob2;
                                    vacuumB3 += prob3;
                                    vacuumB4 += prob4;
                                end

                                if (a1 == 1) && (a2 == 1)
                                    vacuumA1 += prob1;
                                    vacuumA2 += prob2;
                                    vacuumA3 += prob3;
                                    vacuumA4 += prob4;
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    chsh_approach3 = abs(- chsh1_approach3 / (1 - vacuumA1 - vacuumB1 + vacuumAB1) 
                         + chsh2_approach3 / (1 - vacuumA2 - vacuumB2 + vacuumAB2) 
                         + chsh3_approach3 / (1 - vacuumA3 - vacuumB3 + vacuumAB3) 
                         + chsh4_approach3 / (1 - vacuumA4 - vacuumB4 + vacuumAB4))
    chsh_approach2 = abs(-chsh1_approach2 + chsh2_approach2 + chsh3_approach2 + chsh4_approach2)
    
    return chsh_approach2, chsh_approach3
end


function G_loss(k, l, n_cut, γ)
    sum = 0
    for i in 0:n_cut
        for j in 0:i-1
            sum += binomial(k, k-j) * binomial(l, l-i) * γ^(i+j) * (1 - γ)^(k-j+l-i)
        end
    end
    return 1 - 2*sum
end


function approach23_loss_exact!(state, n_cut, N, γ, loss_array)
    chsh1 = 0.0;
    vacuumAsgnB1 = 0.0;
    vacuumBsgnA1 = 0.0;
    vacuumA1 = 0.0;
    vacuumB1 = 0.0;
    vacuumAB1 = 0.0;
    chsh2 = 0.0;
    vacuumAsgnB2 = 0.0;
    vacuumBsgnA2 = 0.0;
    vacuumA2 = 0.0;
    vacuumB2 = 0.0;
    vacuumAB2 = 0.0;
    chsh3 = 0.0;
    vacuumAsgnB3 = 0.0;
    vacuumBsgnA3 = 0.0;
    vacuumA3 = 0.0;
    vacuumB3 = 0.0;
    vacuumAB3 = 0.0;
    chsh4 = 0.0;
    vacuumAsgnB4 = 0.0;
    vacuumBsgnA4 = 0.0;
    vacuumA4 = 0.0;
    vacuumB4 = 0.0;
    vacuumAB4 = 0.0;

    prob1 = 0.0;
    prob2 = 0.0;
    prob3 = 0.0;
    prob4 = 0.0;
    
    lossA = 0.0;
    lossB = 0.0;
    
    @inbounds begin
        for i in 1:n_cut+1
            for j in 1:n_cut+1
                loss_array[i, j] = G_loss(i-1, j-1, n_cut, γ);
            end
        end

        for b0 in 1:N+1
            for b2 in 1:n_cut+1
                for b1 in 1:n_cut+1
                    for a0 in 1:N+1
                        for a2 in 1:n_cut+1
                            for a1 in 1:n_cut+1

                                lossA = loss_array[a1-1, a2-1];
                                lossB = loss_array[b1-1, b2-1];

                                prob1 = abs2(state[1][a1, a2, a0, b1, b2, b0]);
                                chsh1 += prob1 * lossA * lossB;
                                vacuumAsgnB1 += prob1 * (1-γ)^(a1-1+a2-1) * lossB;
                                vacuumBsgnA1 += prob1 * lossA * (1-γ)^(b1-1+b2-1);
                                vacuumA1 += prob1 * (1-γ)^(a1-1+a2-1);
                                vacuumB1 += prob1 * (1-γ)^(b1-1+b2-1);
                                vacuumAB1 += prob1 * (1-γ)^(a1-1+a2-1) * (1-γ)^(b1-1+b2-1)

                                prob2 = abs2(state[2][a1, a2, a0, b1, b2, b0]);
                                chsh2 += prob2 * lossA * lossB;
                                vacuumAsgnB2 += prob2 * (1-γ)^(a1-1+a2-1) * lossB;
                                vacuumBsgnA2 += prob2 * lossA * (1-γ)^(b1-1+b2-1);
                                vacuumA2 += prob2 * (1-γ)^(a1-1+a2-1);
                                vacuumB2 += prob2 * (1-γ)^(b1-1+b2-1);
                                vacuumAB2 += prob2 * (1-γ)^(a1-1+a2-1) * (1-γ)^(b1-1+b2-1);

                                prob3 = abs2(state[3][a1, a2, a0, b1, b2, b0]);
                                chsh3 += prob3 * lossA * lossB;
                                vacuumAsgnB3 += prob3 * (1-γ)^(a1-1+a2-1) * lossB;
                                vacuumBsgnA3 += prob3 * lossA * (1-γ)^(b1-1+b2-1);
                                vacuumA3 += prob3 * (1-γ)^(a1-1+a2-1);
                                vacuumB3 += prob3 * (1-γ)^(b1-1+b2-1);
                                vacuumAB3 += prob3 * (1-γ)^(a1-1+a2-1) * (1-γ)^(b1-1+b2-1);

                                prob4 = abs2(state[4][a1, a2, a0, b1, b2, b0]);
                                chsh4 += prob4 * lossA * lossB;
                                vacuumAsgnB4 += prob4 * (1-γ)^(a1-1+a2-1) * lossB;
                                vacuumBsgnA4 += prob4 * lossA * (1-γ)^(b1-1+b2-1);
                                vacuumA4 += prob4 * (1-γ)^(a1-1+a2-1);
                                vacuumB4 += prob4 * (1-γ)^(b1-1+b2-1);
                                vacuumAB4 += prob4 * (1-γ)^(a1-1+a2-1) * (1-γ)^(b1-1+b2-1);
                            end
                        end
                    end
                end
            end
        end
    end
    chsh_approach3 = abs(-(chsh1 - vacuumAsgnB1 - vacuumBsgnA1 + vacuumAB1)/(1 - vacuumA1 - vacuumB1 + vacuumAB1) + 
                          (chsh2 - vacuumAsgnB2 - vacuumBsgnA2 + vacuumAB2)/(1 - vacuumA2 - vacuumB2 + vacuumAB2) + 
                          (chsh3 - vacuumAsgnB3 - vacuumBsgnA3 + vacuumAB3)/(1 - vacuumA3 - vacuumB3 + vacuumAB3) + 
                          (chsh4 - vacuumAsgnB4 - vacuumBsgnA4 + vacuumAB4)/(1 - vacuumA4 - vacuumB4 + vacuumAB4))
    chsh_approach2 = abs(-chsh1 + chsh2 + chsh3 + chsh4)
    return chsh_approach2, chsh_approach3
end