# Domain: Trial wavefunction for importance sampling and node handling

struct TrialWF{LP,GL,LL,SP}
    logpsi::LP # R -> log|psi_T(R)|
    gradlogpsi::GL # R -> ∇ log|psi_T(R)| (Vector)
    lapllogpsi::LL # R -> ∇² log|psi_T(R)| (Scalar)
    signpsi::SP # R -> sign(psi_T(R)) ∈ {-1, 0, 1}
end

TrialWF(logpsi, gradlogpsi, lapllogpsi) = TrialWF(logpsi, gradlogpsi, lapllogpsi, R -> 1.0)
