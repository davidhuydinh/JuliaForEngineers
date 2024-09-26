using ModelingToolkit
using OrdinaryDiffEq #DifferentialEquations
using Plots

# ------------------------------------------------
# Part 1: Steady State Modeling ------------------
# ------------------------------------------------
pars = @parameters A=0.1 ẋ=1 c=1000 pₛ=300e5 pᵣ=0 ρ=1000 Cₒ=2.7 m=100 ẍ=0
vars = @variables p₁=300e5 p₂=0e5 Aₒ=0.001

# symbolic expressions
u = ẋ * (A/Aₒ)

# equations
eqs = [
    pₛ - p₁ ~ (1/2)*ρ*u^2*Cₒ
    p₂ - pᵣ ~ (1/2)*ρ*u^2*Cₒ
    m*ẍ ~ (p₂ - p₁)*A - c*ẋ
]

@named nlsys = NonlinearSystem(eqs, vars, pars)
sys = structural_simplify(nlsys)
prob = NonlinearProblem(sys, [], []) # [initial conditions], [parameters] 
sol = solve(prob)

sol.resid

sol[Aₒ] #<-- solution!



# how to quickly make a new solution 
orifices = []
velocity_limits = 1.0:0.1:2.0
for velocity_limit in velocity_limits
    prob′ = remake(prob; p=[ẋ => velocity_limit])
    sol′ = solve(prob′; abstol=1e-9)
    push!(orifices, sol′[Aₒ])
end
plot(velocity_limits, orifices; xlabel="velocity limit [m/s]", ylabel="orifice size [m^2]")





# ------------------------------------------------
# Part 2: Dynamic Modeling (DAEs) ----------------
# ------------------------------------------------

using ModelingToolkit: t_nounits as t, D_nounits as D
# @parameters t
# D = Differential(t)



pars = @parameters A=0.1 pₛ=300e5 pᵣ=0 ρ=1000 C₀=2.7 m=100 Aₒ=0.00094 c=1000
vars = @variables begin
    x(t)=0, [guess=0]
    ẋ(t)=0, [guess=0]
    p₁(t), [guess=300e5]
    p₂(t), [guess=0e5]
    ẍ(t), [guess=(p₂-p₁)*A/m]
end


# symbolic expressions
u = ẋ * (A/Aₒ)

# equations
eqs = [
    D(x) ~ ẋ
    D(ẋ) ~ ẍ

    pₛ - p₁ ~ (1/2)*ρ*u^2*C₀
    p₂ - pᵣ ~ (1/2)*ρ*u^2*C₀

    m*ẍ ~ (p₂-p₁)*A - c*ẋ
]

@named odesys = ODESystem(eqs, t, vars, pars)
sys = structural_simplify(odesys)

initsys = ModelingToolkit.generate_initializesystem(sys)
initsys = structural_simplify(initsys)
initprob = NonlinearProblem(initsys, [t=>0], [])
initsol = solve(initprob)

us = unknowns(sys)
u0 = us .=> initsol[us]

prob = ODEProblem(sys, u0, (0.0, 0.0001), [])
sol = solve(prob)

# explain sol object...
plot(sol.t, sol[x]; marker=:circle, ylabel="position [m]")
plot(sol, idxs=[x]; ylabel="position [m]")
plot(sol, idxs=[ẋ]; ylabel="velocity [m/s]")
plot(sol, idxs=[ẍ]; ylabel="acceleration [m/s^2]")
plot(sol, idxs=[p₁, p₂]; ylabel="pressure [Pa]")

# for comparison with compressible system
prob′ = remake(prob, tspan=(0, 0.1))
sol_ic = solve(prob′)




# ------------------------------------------------
# Part 3: Component Based Modeling ---------------
# ------------------------------------------------

# Connectors ----
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/connectors/connections/
@connector Port begin
    p(t), [guess=0]
    ṁ(t), [guess=0, connect = Flow]
end

@connector Flange begin
    ẋ(t), [guess=0]
    f(t), [guess=0, connect = Flow]
end

regPow(x, a, delta = 0.01) = x * (x * x + delta * delta)^((a - 1) / 2);

x = -1e-1:1e-5:1e-1
plot(x, regPow.(x,2))
plot!(x, x.^2)

x = -1e-2:1e-5:1e-2
plot(x, regPow.(x,2))
plot!(x, x.^2)



# Components ----
@mtkmodel Orifice begin
    @parameters begin
        Cₒ=2.7
        Aₒ=0.00094
        ρ₀=1000
    end
    @variables begin
        ṁ(t), [guess=0]
        p₁(t), [guess=0]
        p₂(t), [guess=0]
    end
    @components begin
        port₁ = Port()
        port₂ = Port()
    end
    begin
        u = ṁ/(ρ₀*Aₒ)
    end
    @equations begin
        ṁ ~ +port₁.ṁ
        ṁ ~ -port₂.ṁ
        p₁ ~ port₁.p
        p₂ ~ port₂.p
        
        p₁ - p₂ ~ (1/2)*ρ₀*regPow(u,2)*Cₒ
    end
end

@mtkmodel Volume begin
    @parameters begin
        A=0.1
        ρ₀=1000
        β=2e9
        direction=+1
    end
    @variables begin
        p(t), [guess=0]
        x(t), [guess=0]
        ṁ(t), [guess=0]
        f(t), [guess=0]
        ẋ(t), [guess=0]
        r(t), [guess=0]
        ṙ(t), [guess=0]
    end
    @components begin
        port = Port()
        flange = Flange()
    end
    @equations begin
        D(x) ~ ẋ
        D(r) ~ ṙ
        
        p ~ +port.p
        ṁ ~ +port.ṁ # mass is entering
        f ~ -flange.f * direction # force is leaving
        ẋ ~ flange.ẋ * direction

        r ~ ρ₀*(1 + p/β)
        ṁ ~ (r*ẋ*A) + (ṙ*x*A)
        f ~ p * A
    end
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
    end
    @variables begin
        f(t), [guess=0]
        x(t), [guess=0]
        ẋ(t), [guess=0]
        ẍ(t), [guess=0]
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        D(x) ~ ẋ
        D(ẋ) ~ ẍ

        f ~ flange.f
        ẋ ~ flange.ẋ

        m*ẍ ~ f
    end
end

@mtkmodel Actuator begin
    @components begin
        port₁ = Port()
        port₂ = Port()
        vol₁ = Volume(direction=-1, x=0.5)
        vol₂ = Volume(direction=+1, x=0.5)
        mass = Mass()
        flange = Flange()
    end
    @equations begin
        connect(port₁, vol₁.port)
        connect(port₂, vol₂.port)
        connect(vol₁.flange, vol₂.flange, mass.flange, flange)
    end
end

@mtkmodel Source begin
    @parameters begin
        p
    end
    @components begin
        port = Port()
    end    
    @equations begin
        port.p ~ p
    end
end

@mtkmodel Damper begin
    @parameters begin
        c = 1000
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        flange.f ~ c*flange.ẋ
    end
end


@mtkmodel System begin
    @components begin
        res₁ = Orifice()
        res₂ = Orifice()
        act = Actuator()
        src = Source(p=300e5)
        snk = Source(p=0)
        dmp = Damper()
    end
    @equations begin
        connect(src.port, res₁.port₁)
        connect(res₁.port₂, act.port₁)
        connect(act.port₂, res₂.port₁)
        connect(res₂.port₂, snk.port)
        connect(dmp.flange, act.flange)
    end
end

@mtkbuild sys = System()

initialization_eqs = [
    sys.act.mass.x ~ 0
    sys.act.mass.ẋ ~ 0
    sys.act.vol₁.p ~ 300e5
    sys.act.vol₂.p ~ 0
]

initsys = ModelingToolkit.generate_initializesystem(sys; initialization_eqs)
initsys = structural_simplify(initsys)
initprob = NonlinearProblem(initsys, [t=>0], [])
initsol = solve(initprob)

us = unknowns(sys)
u0 = us .=> initsol[us]


prob = ODEProblem(sys, u0, (0, 0.1), [])

# Solving with ImplicitEuler ---------------------------------------------------------
sol = solve(prob)
sol_r = solve(prob, Rodas5P())


# velocity comparison (incompressible vs. compressible)
plot(sol_r, idxs=[sys.act.mass.ẋ]; ylabel="velocity [m/s]", label="Compressible (Rodas5P)")
plot!(sol_ic, idxs=[ẋ], label="Incompressible")


# What's Next --> Using the ModelingToolkitStandardLibrary
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/
# RC Circuit
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/rc_circuit/
# DC Motor
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/dc_motor_pi/
