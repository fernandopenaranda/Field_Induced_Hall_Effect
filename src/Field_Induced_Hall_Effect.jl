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
    include("model_TaRh2B2.jl")
    # include("model_ferroaxial_graphite.jl")
    include("three_band_ferroaxial_model.jl")
    include("wrappers.jl")
    # include("cluster/cluster_evaluate_subcube.jl")
    include("cluster/cluster_tools.jl")
    include("cluster/create_bashfile.jl")
    export sigma_abc_TaRh2B2_wrapper, hamiltonian_TaRh2B2, d_hamiltonian_TaRh2B2, d_d_hamiltonian_TaRh2B2, sigma_abc_ferroaxial_wrapper_3d, sigma_abc_ferroaxial_wrapper_2d
    export ferroaxial_ham2d, d_ferroaxial_ham2d, ferroaxial_ham3d, d_ferroaxial_ham3d, d2_ferroaxial_ham2d, d2_ferroaxial_ham3d, lattice_vectors, threeband_ferroaxial_ham3d, threeband_d_ferroaxial_ham3d, threeband_d2_ferroaxial_ham3d
    export create_bashfile, create_bashfile_FA, slurm_conductivities, slurm_conductivities_FA, modelpresets, modelpresets_FA, FerroAxialHam, TaRh2B2Params, processing, find_folder, modelpresets_3bandFA, ThreeBandsFerroAxialHam
end
