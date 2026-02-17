
using MAGEMin_C, Plots

data        =   Initialize_MAGEMin("mp", verbose=false);

X           = [0.5922, 0.1813, 0.006, 0.0223, 0.0633, 0.0365, 0.0127, 0.0084, 0.0016, 0.0007, 0.075]
Xoxides     = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "MnO", "H2O"]
sys_unit    = "wt"
n = 50
P = 10.0
T = collect(range(400.0,800.0,n))


out = Vector{out_struct}(undef,n)
for i=1:n
    T_calc  = T[i]       # retrieves the temperature from the temperature array we just defined  
    out[i]  = single_point_minimization(P, T_calc, data, X=X, Xoxides=Xoxides, sys_in=sys_unit, name_solvus = true)
    if "H2O" in out[i].ph
        id_h2o      = findfirst(out[i].ph .== "H2O")
        h2o_wt      = out[i].ph_frac_wt[id_h2o]
        h2o_comp_wt = out[i].PP_vec[id_h2o - out[i].n_SS].Comp_wt

        X           = X .- (h2o_wt .* h2o_comp_wt)
    end
end

q_vol = zeros(Float64,n)

for i=1:n
    if "q" in out[i].ph     #here we check if out[i].ph contains "q"
        id_q = findfirst(out[i].ph .== "q")
        q_vol[i] = out[i].ph_frac_vol[id_q]
    end
end

plot(   T,  q_vol .* 100, 
        label = "qtz",
        xlabel = "T°C",
        ylabel = "vol%")


