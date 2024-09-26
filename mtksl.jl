using ModelingToolkit
using OrdinaryDiffEq #DifferentialEquations
using Plots

using ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible: Valve, HydraulicPort, HydraulicFluid, FixedPressure, Pressure, liquid_density
using ModelingToolkitStandardLibrary.Blocks: Constant, Step, Ramp
using ModelingToolkitStandardLibrary.Mechanical.Translational: Mass, MechanicalPort, Damper, Fixed
using ModelingToolkit: t_nounits as t, D_nounits as D

@component function Volume(;
        #parameters
        area,
        direction = +1, name)
    pars = @parameters begin
        area = area
    end

    vars = @variables begin
        x(t), [guess=0]
        dx(t), [guess=0]
        p(t), [guess=0]
        f(t), [guess=0]
        rho(t), [guess=0]
        drho(t), [guess=0]
        dm(t), [guess=0]
    end

    systems = @named begin
        port = HydraulicPort()
        flange = MechanicalPort()
    end

    eqs = [
        # connectors
        port.p ~ p
        port.dm ~ dm
        flange.v * direction ~ dx
        flange.f * direction ~ -f

        # differentials
        D(x) ~ dx
        D(rho) ~ drho

        # physics
        rho ~ liquid_density(port, p)
        f ~ p * area
        dm ~ drho * x * area + rho * dx * area]

    ODESystem(eqs, t, vars, pars; name, systems)
end

@mtkmodel Orifice begin
    @parameters begin
        Aₒ=0.00094
        Cₒ=2.7
    end
    @components begin
        area = Constant(k=Aₒ)
        valve = Valve(Cd=Cₒ)
        port₁ = HydraulicPort()
        port₂ = HydraulicPort()
    end
    @equations begin
        connect(valve.area, area.output)
        connect(valve.port_a, port₁)
        connect(valve.port_b, port₂)
    end
end

@mtkmodel Actuator begin
    # @parameters begin
    #     L=1.0
    #     area=0.1
    # end
    @components begin
        port₁ = HydraulicPort()
        port₂ = HydraulicPort()
        vol₁ = Volume(direction=-1, area=0.1)
        vol₂ = Volume(direction=+1, area=0.1)
        mass = Mass(m=100)
        flange = MechanicalPort()
    end
    @equations begin
        connect(port₁, vol₁.port)
        connect(port₂, vol₂.port)
        connect(vol₁.flange, vol₂.flange, mass.flange, flange)
    end
end

@mtkmodel System begin
    @components begin
        res₁ = Orifice()
        res₂ = Orifice()
        act = Actuator()
        src = Pressure()
        ramp = Ramp(height = 300e5, start_time = 0.1, duration = 0.01, smooth=false)
        snk = FixedPressure(p=0)
        dmp = Damper(d=1000)
        wall = Fixed()

        fluid = HydraulicFluid(density = 1000, bulk_modulus = 2.0e9)
    end
    @equations begin
        connect(src.p, ramp.output)
        connect(src.port, res₁.port₁)
        connect(res₁.port₂, act.port₁)
        connect(act.port₂, res₂.port₁)
        connect(res₂.port₂, snk.port)
        connect(dmp.flange_a, act.flange)
        connect(dmp.flange_b, wall.flange)

        connect(fluid, src.port, snk.port)
    end
end

@mtkbuild sys = System()

initialization_eqs = [
    sys.act.mass.s ~ 0
    sys.act.mass.v ~ 0
    sys.act.vol₁.x ~ 0.5
    sys.act.vol₂.x ~ 0.5
    sys.act.vol₁.drho ~ 0.0
    sys.act.vol₂.drho ~ 0.0
]

initsys = ModelingToolkit.generate_initializesystem(sys; initialization_eqs)
initsys = structural_simplify(initsys)
initprob = NonlinearProblem(initsys, [t=>0], [])
initsol = solve(initprob)

us = unknowns(sys)
u0 = us .=> initsol[us]


prob = ODEProblem(sys, u0, (0, 0.3), [])

# Solving with ImplicitEuler ---------------------------------------------------------
sol = solve(prob, Rodas5P())

plot(sol, idxs=sys.src.port.p)
plot(sol, idxs=sys.act.mass.v; ylabel="velocity [m/s]", label="Compressible (Rodas5P)")


