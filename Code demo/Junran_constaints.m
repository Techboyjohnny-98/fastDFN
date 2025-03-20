function [c,ceq] = Junran_constaints(x)
    run param/params_Samsung30T.m
    p.L_p           = x(1);
    p.L_n           = x(2);
    p.c_s_p_max     = x(3);
    p.c_s_n_max     = x(4);
    p.epsilon_s_p   = x(5);
    p.epsilon_s_n   = x(6);
    p.Area          = x(7);
    p.R_s_p         = x(8);
    p.R_s_n         = x(9);
    p.D_s_p0        = x(10);
    p.D_s_n0        = x(11);
    p.R_f_n         = x(12);
    p.epsilon_e_p   = x(13);
    p.epsilon_e_n   = x(14);
    p.c_e           = x(15);
    p.L_s           = x(16);
    p.t_plus        = x(17);
    p.n_Li_s = 0.008*(p.c_s_n_max*p.epsilon_s_n*p.L_n*p.Area) ...
        +p.epsilon_s_p*p.L_p*p.Area*0.915*p.c_s_p_max;
%     [cn_low,cp_low] = init_cs(p,p.volt_min);
%     [cn_high,cp_high] = init_cs(p,p.volt_max);
%     Delta_cn = cn_high-cn_low;
%     Delta_cp = cp_low-cp_high;
%     p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);
%     c = [];
%     ceq = p.OneC - 3;
    c = [];
    ceq(1) = 0.661*p.c_s_p_max*p.epsilon_s_p*p.L_p*p.Faraday/3600 - p.Capacity;
    ceq(2) = 0.847*p.c_s_n_max*p.epsilon_s_n*p.L_n*p.Faraday/3600 - p.Capacity;
    
end