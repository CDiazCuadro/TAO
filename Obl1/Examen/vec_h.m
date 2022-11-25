%%%%% Funcion para calcular h(ê)
function h = vec_h(e)
    n = length(e);
    h = zeros(n,1);
    for idx = 1 : n
        ei = e(idx);
        h(idx) = ei^2/(abs(ei)+1);
    end
end