module Field_Induced_Hall_Effect
    using Optics_in_the_length_gauge
    using Brillouin
    using Bravais
    # using BZpath
    using Arpack
    using LinearAlgebra
    using Cubature
    using PhysicalConstants
    using PhysicalConstants.CODATA2018
    using Unitful
    using SparseArrays
    using StaticArrays
    using Parameters
    using CSV
    using JLD2
    using Serialization
    include("model.jl")
    include("wrapper.jl")
    # include("cluster/cluster_evaluate_subcube.jl")
    include("cluster/cluster_tools.jl")
    include("cluster/create_bashfile.jl")

    export sigma_abc_TaRh2B2_wrapper, hamiltonian_TaRh2B2, d_hamiltonian_TaRh2B2, d_d_hamiltonian_TaRh2B2
    export create_bashfile, slurm_conductivities, modelpresets

end
