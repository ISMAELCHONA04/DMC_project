# Domain: Trial wavefunction for importance sampling and node handling

"""
    TrialWF(logpsi, gradlogpsi, lapllogpsi)
    TrialWF(logpsi, gradlogpsi, lapllogpsi, signpsi)

Callback bundle describing a real trial wavefunction through `log(abs(psi_T))`
and its derivatives. `signpsi` defaults to `+1` everywhere, which is suitable
for nodeless bosonic trials but not for fixed-node calculations.
"""
struct TrialWF{LP,GL,LL,SP}
    logpsi::LP # R -> log|psi_T(R)|
    gradlogpsi::GL # R -> ∇ log|psi_T(R)| (Vector)
    lapllogpsi::LL # R -> ∇² log|psi_T(R)| (Scalar)
    signpsi::SP # R -> sign(psi_T(R)) ∈ {-1, 0, 1}
end

TrialWF(logpsi, gradlogpsi, lapllogpsi) = TrialWF(logpsi, gradlogpsi, lapllogpsi, R -> 1.0)
