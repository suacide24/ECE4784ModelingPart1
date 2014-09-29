function ModelingPhase1 (probNo) 
%Ramon Sua
%ECE 4784 Modeling Project Part 1
%Due date 09/29/14
%This function models the Hodgkins Huxley's Model of the cell membrane
%Input - ProbNo
%Output - Corresponding Simulation

%-------------------------Simulation Constants------------------------------------%
%Biological Constants (K for potassium, Na for sodium)
%Conductances in [mS/cm^2]
 gK_given = 36;
 gNa_given = 120;
 gL_given = 0.3;
%Nerst Potentials in [mV]
 VK = -12;
 VNa = 115;
 VL = 10.6;
% Other Constants
 Vrest = -70; %[mV]
 Cm = 1;     %   [microF/cm^2]
 
%Simulation Parameters
simtime = 100 ; %[time in ms]
delta = 0.0001 ; %smaller delta produces smoother results

%---------------------------End of Constants--------

%Initialization
t = 0:delta:simtime ; 
l = length(t) ; 
I = zeros(1,l) ;

%Initialize the variables needed to implement Euler's method
V = zeros(1, l);
Vm = zeros(1, l);
m = zeros(1, l);
n = zeros(1, l);
h = zeros(1, l);
gna = zeros(1, l);
gk = zeros(1, l);
dm = 0;
dn = 0;
dh = 0;
dVm = 0; 

%If statement to choose which problem # and assign the value of I
%accordingly.
if probNo == 1
    settitle = 'response when no current supplied' ; 
    %No changes needed to the current vector as it is already initialized
    %to 0.
elseif probNo ==2
    settitle = 'response when 5microA/cm^2 within 0.5ms of current supplied' ;
    I(1:(1+ 0.5/delta)) = 5 ;   %Supply the 0.5 ms of current at the very beginning (index1)
elseif probNo ==3
    settitle = 'response when constant of 5microA/cm^2 is supplied' ; 
    I = ones(1,l) * 5 ;
end

%Setting Initial Values
%This will be helpful so we can index x-1 without getting an error for 
%indexing 0

%Initial votlage values
V(1) = Vrest ; 
Vm(1) = 0 ; 

%Calculate Eulers method cosntants using the equations given for ind=1
%A for alpha B for beta
Am = 0.1*((25-Vm(1))/(exp((25-Vm(1))/10) - 1));
Bm = 4*exp(-Vm(1)/18);
An = .01*((10-Vm(1))/(exp((10-Vm(1))/10) - 1));
Bn = .125*exp(-Vm(1)/80);
Ah = .07*exp(-Vm(1)/20);
Bh = 1/(exp((30-Vm(1))/10) + 1);

%calculate m, n and h. each are define as A/(A+B) for ind=1
m(1) = Am/(Am+Bm) ; 
n(1) = An/(An+Bn) ; 
h(1) = Ah/(Ah+Bh) ; 


%Currents (with exception of the input current) 
%have a value of zero at the beginning of the simulation

%Initial conductance values 
gna(1) = gNa_given*(m(1)^3)*h(1) ; 
gk(1) = gK_given*(n(1)^4) ; 

%Iteration from 2 till the end of simulation
for ind = 2:l   %we skip 1 since we already defined it in the previous lines
    Vm(ind) = Vm(ind-1) +dVm * delta ; %membrane voltage w/out Vrest
    V(ind) = Vm(ind) + Vrest; %took into account Vrest
    
    %Euler's Formula (the value of the next depends on the value before)
    m(ind) = m(ind-1) + dm * delta;
    n(ind) = n(ind-1) + dn * delta;
    h(ind) = h(ind-1) + dh * delta;


    %Find the value of the components through curve fitting
    gna(ind) = (m(ind)^3)*gNa_given*(h(ind));
    gk(ind) = (n(ind)^4)*gK_given;

    %Constants calculated using the formulas given (A for alpha, B for
    %beta)
    Am = 0.1*((25-Vm(ind))/(exp((25-Vm(ind))/10) - 1));
    Bm = 4*exp(-Vm(ind)/18);
    An = .01*((10-Vm(ind))/(exp((10-Vm(ind))/10) - 1));
    Bn = .125*exp(-Vm(ind)/80);
    Ah = .07*exp(-Vm(ind)/20);
    Bh = 1/(exp((30-Vm(ind))/10) + 1);
    
    %Calculate Resulting currents. Note that the current changes as the
    %conductances changes
    Ina = gna(ind)*(Vm(ind) - VNa);
    Ik = gk(ind)*(Vm(ind) - VK);
    Il = gL_given*(Vm(ind) - VL);
    Iion = I(ind) - Ina - Ik - Il;
    
    %Calculate new derivatives for next iteration
    dVm = Iion/Cm;
    dm = Am*(1-m(ind)) - Bm*m(ind);
    dn = An*(1-n(ind)) - Bn*n(ind);
    dh = Ah*(1-h(ind)) - Bh*h(ind);
end

    %------------------Plotting Results -------------------------------
    figure(1)
    plot(t,gNa_given*(m.^3).*h)
    hold on  %Set hold to on so you can plot 2 at the same graph
    plot(t,(n.^4).*gK_given, '--');  %-- makes it dashed
    legend('Na conductance', 'K conductance'); %sets the legend
    
    %axis labels and title
    title('sodium and potassium conductances');
    ylabel('conductance [mS/cm^2]');
    xlabel('time [ms]');
    hold off %sets hold back to off
    
    figure(2)
    plot(t,V);
    title(settitle) ;  %settitle was defined during the if statement
    xlabel('time [ms]');
    ylabel('membrane voltage [mV]');
    %change axis from 0 to the simtime. Y axis are arbitrary but make sure
    %to include -70mV which is the resting voltage
    axis([0 100 -90 90]);
end
