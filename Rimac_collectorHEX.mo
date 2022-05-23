model Rimac_collectorHEX
  
  record Coldplate "geometry of the coldplate"
    Real length;
    Real height;
    Real width;
  end Coldplate;
  
  record Collector "geometry of circular collector pipe"
    Real length;
    Real diameter;
  end Collector;
  
  record ThermalLayer "geometry and thermal properties of battery/coldplate layers"
    Real thickness;
    Real width;
    Real length;
    Real lambda;
  end ThermalLayer;
  
  function hydraulicD "Hydraulic diameter for non circular pipes"
    input Coldplate cp;
    output Real dHyd;
  algorithm
    dHyd := (4*cp.width*cp.height/(cp.width*2+cp.height*2));
  end hydraulicD;
  
  function RectPipeRe "Reynolds number for non circular pipe"
    input Real width;
    input Real height;
    input Real mu;
    input Real ro;
    input Real velocity;
    output Real Re;
  algorithm
    Re := velocity * (4*width*height/(width*2+height*2))*ro/mu;
  end RectPipeRe;
  
  /*Boundary conditions and parts properties--------->>>>*/
  parameter Real batteryHeatGen = 120;
  parameter Real massFlowrate = 0.6 "kg/s";
  replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Medium.ThermodynamicState state = Medium.setState_pT(p=1e5, T=20+273.15) "state of the water";
  
  parameter Coldplate cp1(length = 0.36, height = 0.008, width = 0.063);
  parameter Coldplate cp2(length = 0.36, height = 0.010, width = 0.063);
  parameter Coldplate cp3(length = 0.36, height = 0.012, width = 0.063);
  
  parameter Collector collector(length = 0.2, diameter = 0.036);
  
  parameter ThermalLayer tim(thickness=0.001, width=0.063, length=0.3, lambda=1.6);
  parameter ThermalLayer alWall(thickness=0.002, width=0.063, length=0.3, lambda=160);
  
  parameter Real Rth_tim = tim.thickness/(tim.lambda * tim.width*tim.length);
  parameter Real Rth_alWall = alWall.thickness/(alWall.lambda * alWall.width*alWall.length);
  parameter Real Rth_total = Rth_tim + Rth_alWall;
  /*<<<<<<<<----------Boundary conditions and parts properties*/
  
  
  /*****Coldplate HTC calculation----------->>>>>>>>***********/
  
  Medium.Density ro = Medium.density(state);
  Medium.DynamicViscosity mu = Medium.dynamicViscosity(state);
  Medium.ThermalConductivity lambdaWater = Medium.thermalConductivity(state);
  
  Real cp1Velocity = massFlowrate/ro/(cp1.height*cp1.width);
  Real cp2Velocity = massFlowrate/ro/(cp2.height*cp2.width);
  Real cp3Velocity = massFlowrate/ro/(cp3.height*cp3.width);
  
  Real cp1Re = RectPipeRe(width = cp1.width, height = cp1.height, velocity = cp1Velocity, mu = mu, ro = ro);
  Real cp2Re = RectPipeRe(width = cp2.width, height = cp2.height, velocity = cp2Velocity, mu = mu, ro = ro);
  Real cp3Re = RectPipeRe(width = cp3.width, height = cp3.height, velocity = cp3Velocity, mu = mu, ro = ro);
  
  Real cpContactArea = cp1.length*cp1.width;
  
  Real cp1Nu = 0.021*cp1Re^(0.8)*0.71^(0.4);
  Real cp1HTC_info = cp1Nu / hydraulicD(cp1) * lambdaWater *cpContactArea;
  Real cp2Nu = 0.021*cp2Re^(0.8)*0.71^(0.4);
  Real cp2HTC_info = cp2Nu / hydraulicD(cp2) * lambdaWater *cpContactArea;
  Real cp3Nu = 0.021*cp3Re^(0.8)*0.71^(0.4);
  Real cp3HTC_info = cp3Nu / hydraulicD(cp3) * lambdaWater *cpContactArea;
  
  parameter Real cp1HTC = 42 "Attention!!! W/K; this value must correspond to cp1HTC_info";
  parameter Real cp2HTC = 34 "Attention!!! W/K; this value must correspond to cp2HTC_info";
  parameter Real cp3HTC = 28 "Attention!!! W/K; this value must correspond to cp3HTC_info";
  /*****<<<<<<<<--------Coldplate HTC calculation ***********/
  
  inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, momentumDynamics = Modelica.Fluid.Types.Dynamics.SteadyState)  annotation(
    Placement(visible = true, transformation(origin = {-226, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = Medium, T = 293.15, m_flow = massFlowrate, nPorts = 1, use_m_flow_in = false) annotation(
    Placement(visible = true, transformation(origin = {-272, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, nPorts = 1, p = 100000) annotation(
    Placement(visible = true, transformation(origin = {-268, -104}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Pipes.DynamicPipe pipe(redeclare package Medium = Medium, allowFlowReversal = false, crossArea = cp1.height * cp1.width, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, isCircular = false, length = cp1.length, perimeter = (cp1.width + cp1.height)*2, use_HeatTransfer = true)  annotation(
    Placement(visible = true, transformation(origin = {-136, -42}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow = batteryHeatGen)  annotation(
    Placement(visible = true, transformation(origin = {-30, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  Modelica.Thermal.HeatTransfer.Components.Convection convection annotation(
    Placement(visible = true, transformation(origin = {-104, -42}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant convectionConstant(k = cp1HTC) annotation(
    Placement(visible = true, transformation(origin = {-102, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistor(R = Rth_total)  annotation(
    Placement(visible = true, transformation(origin = {-68, -42}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Fluid.Pipes.StaticPipe staticPipe(redeclare package Medium = Medium, allowFlowReversal = true, diameter = collector.diameter, length = collector.length) annotation(
    Placement(visible = true, transformation(origin = {-198, -104}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Fluid.Pipes.StaticPipe pipe1 (redeclare package Medium = Medium, allowFlowReversal = false, diameter = collector.diameter, length = collector.length) annotation(
    Placement(visible = true, transformation(origin = {-202, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Pipes.StaticPipe pipe2(redeclare package Medium = Medium, diameter = collector.diameter, length = collector.length)  annotation(
    Placement(visible = true, transformation(origin = {-68, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Pipes.StaticPipe pipe4(redeclare package Medium = Medium, allowFlowReversal = true, diameter = collector.diameter, length = collector.length)  annotation(
    Placement(visible = true, transformation(origin = {-60, -104}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
  Modelica.Fluid.Fittings.TeeJunctionIdeal teeJunctionIdeal (redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-136, 26}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Fluid.Fittings.TeeJunctionIdeal teeJunctionIdeal1 (redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {166, -106}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Fluid.Fittings.TeeJunctionIdeal teeJunctionIdeal2 (redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {12, -104}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Fittings.TeeJunctionIdeal teeJunctionIdeal3 (redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {12, 26}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Fluid.Pipes.StaticPipe pipe3(redeclare package Medium = Medium, diameter = collector.diameter, length = collector.length)  annotation(
    Placement(visible = true, transformation(origin = {90, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Pipes.StaticPipe pipe5(redeclare package Medium = Medium, diameter = collector.diameter, length = collector.length)  annotation(
    Placement(visible = true, transformation(origin = {90, -104}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.Convection convection1 annotation(
    Placement(visible = true, transformation(origin = {46, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Fluid.Pipes.DynamicPipe pipe6(redeclare package Medium = Medium, allowFlowReversal = false, crossArea = cp2.height * cp2.width, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, isCircular = false, length = cp2.length, perimeter = (cp2.width + cp2.height)*2, use_HeatTransfer = true) annotation(
    Placement(visible = true, transformation(origin = {12, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Sources.Constant const(k = cp2HTC) annotation(
    Placement(visible = true, transformation(origin = {46, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistor1(R = Rth_total) annotation(
    Placement(visible = true, transformation(origin = {82, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow1(Q_flow = batteryHeatGen) annotation(
    Placement(visible = true, transformation(origin = {120, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  Modelica.Fluid.Pipes.DynamicPipe pipe7(redeclare package Medium = Medium, allowFlowReversal = false, crossArea = cp3.height*cp3.width, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, isCircular = false, length = cp3.length, perimeter = (cp3.width+cp3.height)*2, use_HeatTransfer = true) annotation(
    Placement(visible = true, transformation(origin = {168, -38}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Components.Convection convection11 annotation(
    Placement(visible = true, transformation(origin = {202, -38}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = cp3HTC) annotation(
    Placement(visible = true, transformation(origin = {202, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistor11(R = Rth_total) annotation(
    Placement(visible = true, transformation(origin = {238, -38}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow11(Q_flow = batteryHeatGen) annotation(
    Placement(visible = true, transformation(origin = {280, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  Modelica.Fluid.Fittings.TeeJunctionIdeal teeJunctionIdeal21 (redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {168, 26}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Fluid.Fittings.TeeJunctionIdeal teeJunctionIdeal4 (redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-134, -104}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
equation
  connect(convectionConstant.y, convection.Gc) annotation(
    Line(points = {{-102, -17}, {-103, -17}, {-103, -32}, {-104, -32}}, color = {0, 0, 127}));
  connect(convection.fluid, pipe.heatPorts[1]) annotation(
    Line(points = {{-114, -42}, {-114, -43}, {-132, -43}, {-132, -42}}, color = {191, 0, 0}));
  connect(fixedHeatFlow.port, thermalResistor.port_a) annotation(
    Line(points = {{-40, -42}, {-40, -43}, {-58, -43}, {-58, -42}}, color = {191, 0, 0}));
  connect(thermalResistor.port_b, convection.solid) annotation(
    Line(points = {{-78, -42}, {-94, -42}}, color = {191, 0, 0}));
  connect(pipe2.port_b, teeJunctionIdeal3.port_1) annotation(
    Line(points = {{-58, 26}, {2, 26}}, color = {0, 127, 255}));
  connect(pipe4.port_a, teeJunctionIdeal2.port_1) annotation(
    Line(points = {{-50, -104}, {2, -104}}, color = {0, 127, 255}));
  connect(teeJunctionIdeal2.port_2, pipe5.port_a) annotation(
    Line(points = {{22, -104}, {80, -104}}, color = {0, 127, 255}));
  connect(teeJunctionIdeal3.port_2, pipe3.port_a) annotation(
    Line(points = {{22, 26}, {80, 26}}, color = {0, 127, 255}));
  connect(teeJunctionIdeal3.port_3, pipe6.port_a) annotation(
    Line(points = {{12, 16}, {12, -30}}, color = {0, 127, 255}));
  connect(pipe6.port_b, teeJunctionIdeal2.port_3) annotation(
    Line(points = {{12, -50}, {12, -94}}, color = {0, 127, 255}));
  connect(convection1.fluid, pipe6.heatPorts[1]) annotation(
    Line(points = {{36, -40}, {16, -40}}, color = {191, 0, 0}));
  connect(thermalResistor1.port_b, convection1.solid) annotation(
    Line(points = {{72, -40}, {56, -40}}, color = {191, 0, 0}));
  connect(fixedHeatFlow1.port, thermalResistor1.port_a) annotation(
    Line(points = {{110, -38}, {92, -38}, {92, -40}}, color = {191, 0, 0}));
  connect(const.y, convection1.Gc) annotation(
    Line(points = {{46, -17}, {46, -31}}, color = {0, 0, 127}));
  connect(convection11.fluid, pipe7.heatPorts[1]) annotation(
    Line(points = {{192, -38}, {172, -38}}, color = {191, 0, 0}));
  connect(thermalResistor11.port_b, convection11.solid) annotation(
    Line(points = {{228, -38}, {212, -38}}, color = {191, 0, 0}));
  connect(fixedHeatFlow11.port, thermalResistor11.port_a) annotation(
    Line(points = {{270, -38}, {248, -38}}, color = {191, 0, 0}));
  connect(const1.y, convection11.Gc) annotation(
    Line(points = {{202, -15}, {202, -28}}, color = {0, 0, 127}));
  connect(boundary.ports[1], pipe1.port_a) annotation(
    Line(points = {{-262, 26}, {-212, 26}}, color = {0, 127, 255}));
  connect(pipe1.port_b, teeJunctionIdeal.port_1) annotation(
    Line(points = {{-192, 26}, {-146, 26}}, color = {0, 127, 255}));
  connect(teeJunctionIdeal.port_3, pipe.port_a) annotation(
    Line(points = {{-136, 16}, {-136, -32}}, color = {0, 127, 255}));
  connect(teeJunctionIdeal.port_2, pipe2.port_a) annotation(
    Line(points = {{-126, 26}, {-78, 26}}, color = {0, 127, 255}));
  connect(pipe5.port_b, teeJunctionIdeal1.port_2) annotation(
    Line(points = {{100, -104}, {156, -104}, {156, -106}}, color = {0, 127, 255}));
  connect(pipe7.port_b, teeJunctionIdeal1.port_3) annotation(
    Line(points = {{168, -48}, {166, -48}, {166, -96}}, color = {0, 127, 255}));
  connect(pipe3.port_b, teeJunctionIdeal21.port_1) annotation(
    Line(points = {{100, 26}, {158, 26}}, color = {0, 127, 255}));
  connect(teeJunctionIdeal21.port_3, pipe7.port_a) annotation(
    Line(points = {{168, 16}, {168, -28}}, color = {0, 127, 255}));
  connect(boundary1.ports[1], staticPipe.port_a) annotation(
    Line(points = {{-258, -104}, {-208, -104}}, color = {0, 127, 255}));
  connect(staticPipe.port_b, teeJunctionIdeal4.port_2) annotation(
    Line(points = {{-188, -104}, {-144, -104}}));
  connect(teeJunctionIdeal4.port_1, pipe4.port_b) annotation(
    Line(points = {{-124, -104}, {-70, -104}}, color = {0, 127, 255}));
  connect(pipe.port_b, teeJunctionIdeal4.port_3) annotation(
    Line(points = {{-136, -52}, {-134, -52}, {-134, -94}}, color = {0, 127, 255}));
  annotation(
    uses(Modelica(version = "4.0.0")));
end Rimac_collectorHEX;
