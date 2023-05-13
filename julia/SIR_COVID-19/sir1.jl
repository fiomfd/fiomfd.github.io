# sir_tokyo.jl
######################
# Simulation of SIR model for COVID-19 in the next D days
#
# t: time [day]
# A: area
# N: total population
# d0: starting date
# d1: ending date
# D=d1-d0 [days]
# S=u(1): the number of susceptible individuals
# I=u(2): active cases
# R=u(3): the number of dischrged individuals
# S+I+R=N always holds.
# I+R: total cases
# T0: total cases at t=0
# T1: daily new cases at t=0
# R0: total dischrged at t=0
# R1: daily discharged at t=0
# beta: infection rate
# gamma: recovery rate
#
# Input A, d0, d1, T0, T1, R0, R1 and N.
# S(0), I(0), R(0), gamma and D are determined 
# by T0, T1, R0, R1, d0, D and N. 
# Run the program in Julia. 
# Julia will show you the result.
#######################

using DifferentialEquations
using Plots
gr()
using Dates

# input data
#A=string(Tokyo);
d0=Date(2021,1,2);
T0=61774;
T1=814;
R0=51295;
R1=51295-50861;
D=240;
N=13999568;
d1=d0*Day(D);

# computing
beta=T1/(T0-R0)/(N-T0);
gamma=R1/(T0-R0);
tspan=(0,D);
f(u,p,t) = [beta*u[3]*(N-u[1]); gamma*u[3]; u[3]*(beta*(N-u[1])-gamma)];
u0=[T0;R0;T0-R0];
prob = ODEProblem(f,u0,tspan);

# solving
sol = solve(prob);

# plotting
l0=string(d0);
l1=string(d0+Day(100));
l2=string(d0+Day(200));
plot(sol, 
    grid=false,
    linewidth=3, 
    title="SIR model for COVID-19 in Tokyo", 
    xlabel="",
    yaxis="population",
    label=["total cases" "discharged" "active cases"], 
    legend = :topleft)
plot!(xticks = ([0:100:200;], [l0, l1, l2]))
