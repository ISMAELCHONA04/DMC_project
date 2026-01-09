# Use-cases: runtime parameters for DMC

struct DMCParams
    dt::Float64
    nsteps::Int
    nequil::Int
    targetN::Int   # target walker population   
    ET0::Float64            # initial reference-energy guess
    pop_gain::Float64       # how aggressively to adjust ET (optional)
    branch_cap::Int        # max branching factor (optional)
    nblocks::Int          # number of blocks for averaging (optional)
end