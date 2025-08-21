using DataFrames, Statistics, CSV,CairoMakie, Revise
using Hadron
df = DataFrame(CSV.File("data/simulation_log_224911.txt"))
bs = Bootstrap(x->exp(-x), df, :delta_H)
res = boot(bs,mean, skip=499)