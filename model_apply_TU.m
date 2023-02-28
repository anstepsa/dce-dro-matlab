


function Cv = model_apply_TU(k, Cp, t)
% Apply the tissue-uptake model (TU)...
% Implementation (adapted slightly) from ROCKETSHIP code

    % N.B. Time unit is 'min' here throughout...
    
    PS = k(1); vp = k(2); Fp = k(3);
    
    E  = PS / (Fp + PS);
    Tp = vp / (Fp + PS);

    Cv = zeros(1, numel(t));
    
    N = numel(t);
    
    for u = 1:N

        T  = t(1:u);
        CP = Cp(1:u);

        F = CP .* Fp .* ( ( (1 - E) .* exp(-(T(end)-T)./Tp) ) + E );
        
        if(numel(T) == 1)
            M = 0;
        else
            M = trapz(T,F); 
        end

        Cv(u) = M;
        
    end
    
    Cv = Cv';

end

