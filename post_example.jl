#example script for interpreting results of trio-GCTA models. CBCL attention used in example
#author: Espen Eilertsen, Laura Hegemann 

using LinearAlgebra
using VCModels
using JLD
using StatsModels, StatsBase, Statistics
using DataFrames, CSV

cd("Z:/projects/Neurodev_trioGCTA/results")

vari = "att_ADHD_scale_std"
m_full = JLD.load(vari * "_full_64599.jld")["model"]

m_nocov = JLD.load(vari * "_nocov_64599.jld")["model"]
m_direct = JLD.load(vari * "_direct_64599.jld")["model"]
m_null = JLD.load(vari * "_null_64599.jld")["model"]

# Compare fits
# -------------------------------------------
lrtest(m_full, m_nocov, m_direct, m_null)
ms = [m_full, m_nocov, m_direct, m_null]
round.(aic.(ms), digits = 2)
round.(deviance.(ms), digits = 2)

aic_m = round.(aic.(ms), digits = 2)
bic_m = round.(bic.(ms), digits = 2)
deviance_m = round.(deviance.(ms), digits = 2)

test_1 = lrtest(m_full, m_nocov, m_direct, m_null)
test_2 = hcat(collect(test_1.dof),collect(test_1.deviance), collect(test_1.pval))
test = DataFrame(test_2, ["DOF", "Deviance", "p(>Chisq)"])
test.ΔDOF = [NaN; diff(test.DOF)]
test.ΔDeviance = [NaN; diff(test.Deviance)]
CSV.write(vari * "_lrtest.csv", test)


fit_1 = hcat(aic_m, bic_m, deviance_m)
fit_2 = DataFrame(fit_1, ["aic", "bic", "deviance"])
CSV.write(vari * "_fit.csv", fit_2)

# Full model
# --------------------------------------------
# rs = [A_m, A_f, A_c, D_fm, D_cm, D_cf, R]
function VCModels.transform(m::VCModel) 
    θ = m.θ
    δ = Vector{eltype(θ)}(undef, length(θ))
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[4]^2
    δ[3] = θ[3]^2 + θ[5]^2 + θ[6]^2
    δ[4] = θ[2] * θ[1]
    δ[5] = θ[3] * θ[1]
    δ[6] = θ[5] * θ[4] + θ[3] * θ[2]
    δ[7] = θ[7]
    println(δ)
    println("worked")
    δ
end
δ = m_full.δ
VCModels.transform(m_full)
println(δ)
Δ = [δ[1] δ[4] δ[5] 0;
     δ[4] δ[2] δ[6] 0;
     δ[5] δ[6] δ[3] 0;
     0 0 0 δ[7]]
# Rescale to % variance
v_pos = [1, 2, 3, 5, 6, 7]
var_tot = sum(δ[v_pos])
sc = diagm(fill(1 / sqrt(var_tot), 4))
Δ_std = sc * Δ * sc


par_full = [Δ_std[3, 3], Δ_std[1, 1], Δ_std[2, 2],
Δ_std[3, 1], Δ_std[3, 2], Δ_std[2, 1], Δ_std[4, 4]]
round.(par_full, digits = 3)'
sum(par_full[[1, 2, 3, 4, 5, 7]])


# Correlations
s = sqrt.(diag(Δ))
cor_mat = diagm(1 ./ s) * Δ * diagm(1 ./ s)

#[mf, cm, cf]
cor_1 = hcat( ["mf", "cm", "cf"], [cor_mat[1,2], cor_mat[1,3],cor_mat[2,3]]) 

cor = DataFrame(cor_1, ["var", "cor"])
CSV.write(vari * "_cor.csv", cor)
# No covariances model
# --------------------------------------------

function VCModels.transform(m::VCModel)
    θ = m.θ
    δ = Vector{eltype(θ)}(undef, length(θ))
    δ[1] = θ[1]^2    
    δ[2] = θ[2]^2    
    δ[3] = θ[3]^2  
    δ[4] = θ[4]
    δ
end

δ_nocov = m_nocov.δ
#VCModels.transform(m_nocov)
Δ_nocov = [δ_nocov[1] 0 0 0; 
           0 δ_nocov[2] 0 0;
           0 0 δ_nocov[3] 0;
           0 0 0 δ_nocov[4]]

var_tot_nocov = sum(δ_nocov)
sc_nocov = diagm(fill(1 / sqrt(var_tot_nocov), size(Δ_nocov, 1)))
Δ_std_nocov = sc_nocov * Δ_nocov * sc_nocov
par_nocov = [Δ_std_nocov[3, 3],Δ_std_nocov[1, 1], Δ_std_nocov[2, 2],0,0,0, Δ_std_nocov[4, 4]]
round.(par_nocov, digits = 3)'
sum(par_nocov)


# Direct model
# --------------------------------------------
function VCModels.transform(m::VCModel)
    θ = m.θ
    δ = Vector{eltype(θ)}(undef, length(θ))
    δ[1] = θ[1]^2
    δ[2] = θ[2]
    δ
end


δ_direct = m_direct.δ
Δ_direct = [δ_direct[1] 0;
            0 δ_direct[2]]

var_tot_direct = sum(δ_direct)
sc_direct = diagm(fill(1 / sqrt(var_tot_direct), size(Δ_direct, 1)))
Δ_std_direct = sc_direct * Δ_direct * sc_direct
par_direct = [Δ_std_direct[1, 1], 0, 0, 0, 0, 0, Δ_std_direct[2, 2]]
round.(Δ_std_direct, digits = 3)
var_tot_direct = sum(Δ_std_direct)

#std errors
function VCModels.transform(θ::Vector)
    δ = similar(θ)
    δ[1] = θ[1]^2
    δ[2] = θ[2]
     
    δ ./ sum(δ)
end
    
VCModels.transform(m_direct.θ)


sqrt.(diag(vcovvc(m_direct)))

J_direct= jacobian(m_direct)
C_direct = vcovvc(m_direct)
Ct_direct = J_direct * C_direct * J_direct'
set_direct = sqrt.(diag(Ct_direct))


# Null model
# --------------------------------------------
function VCModels.transform(m::VCModel)
    m.θ
end

δ_null = VCModels.transform(m_null)
sc_null = diagm(fill(1 / sqrt(δ_null[1]), 1))
std_null = sc_null * δ_null * sc_null
par_null = [0, 0, 0, 0, 0, 0, 1.0]


# Write results
# ----------------------------------------------
par = vcat(par_full')
datres = DataFrame(par, ["o", "m", "p", "om", "op", "mp", "e"])
datres[!, :model] = [ "full", "nocov", "direct", "null"]
datresl = stack(datres)
CSV.write(vari * ".csv", datresl)
round.(datres[:, 1:7], digits = 2)
