# sir2.jl
######################
# Simulation of SIR model for COVID-19 in the next D days
#
# t: time [day]
# A: area
# N: total population
# d0: starting date t=0
# d1: ending date t=D
# d2: starting date of taking average t=-d
# D=d1-d0 [days]
# S=u(1): the number of susceptible individuals
# I=u(2): active cases
# R=u(3): the number of dischrged individuals
# S+I+R=N always holds.
# I+R: total cases
# T0: total cases at t=0
# T2: total cases at t=-d
# R0: total dischrged at t=0
# R2: total discharged at t=-d
# beta: infection rate
# gamma: recovery rate
#
# Input A, d0, D, d2, T0, R0, T2, R2 and N.
# S(0), I(0), R(0), beta, gamma, D and d are determined 
# by T0, R0, T2, R2, d0, d1, d2. 
# Run the program in Julia
# Julia will show you the result.
#######################

using DifferentialEquations
using Plots
gr()
using Dates

# input data
d0=Date(2021,1,3);
T0=62590;
R0=51657;
#
D=240
d1=d0+Day(D);
N=13999568;

# computing
d2=Date(2020,12,27);
T2=56559;
R2=48442;
T1=(T0-T2)/Dates.value(d0-d2);
R1=(R0-R2)/Dates.value(d0-d2);
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
    title="SIR model for COVID-19 in Tokyo (14M)", 
    xlabel="",
    yaxis="population",
    label=["total cases" "discharged" "active cases"], 
    legend = :topleft)
plot!(xticks = ([0:100:200;], [l0, l1, l2]))