project: CoarseMC
summary: Flexible hybrid-parallel Monte Carlo simulation program for coarse-grained molecular models.
author: Jouni Karjalainen
year: 2016
email: jouni.m.karjalainen@gmail.com
project_dir: ./src
output_dir: ./ford-doc
display: public
source: false
exclude: m_fileunit.f90
exclude: json_module.F90
exclude: diamagnetic.f90
exclude: m_rodsphere_interaction.f90
exclude: printforces.f90
exclude: printpotentials.f90
exclude: ljwall_cephes.f90
warn: true
graph: true

Here is the API documentation for the Monte Carlo molecular simulation
program CoarseMC. The documentation is aimed at
developers interested in extending CoarseMC. The purpose is to provide
sufficient information about the interfaces and derived types that
provide points of extension. These are mostly concentrated in the
m_particle module, which contains the particle, pair_interaction and
single_interaction types.

For those who want to modify the behaviour of the simulation, the
mc_engine module is the place to begin and m_nvt_update, m_npt_update
and beta_exchange modules contain the specifics of the Monte Carlo
updates.

Notes on usage and input/output of the program are given in the
User Manual, included in the program package.