% To grab non-zero parts of PSD
function [x,C,P] = get_valids(x_raw,C_raw,P_raw)
    x = [];
    C = [];
    P = [];
    for i = 1:length(x_raw)
        if (P_raw(i) > 0.0) && (C_raw(i) < 1.0)
            x = [x x_raw(i)];
            C = [C C_raw(i)];
            P = [P P_raw(i)];
        end
    end
end