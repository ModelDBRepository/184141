%
%                ESTIMATION OF CONDUCTANCES 
%  2 DIMENSIONAL CONDUCTANCE-BASED MODELS OF QUADRATIC TYPE
%
%  This program estimates conductances in the subthreshold regime under the
%  presence of subthreshold-activated ionic currents. The method is
%  presented at: C.Vich, A.Guillamon (2015), preprint accepted to Journal 
%  of Computational Neuroscience.
%
%  The estimation can be done in models that exhibit a parabolic nullcline 
%  for the voltage and a linear nullcline for the recovery gating variable,
%  i.e. the resonant currents (e.g slow potassium) and amplifying currents 
%  (e.g., persistent sodium).
%  
%  CELL TYPE: medial entorhinal cortex layer II stellate cells and CA1 
%  pyramidal cells.
%
%
%  Model needs to be quadratized (see Horacio G. Rotstein (2015) for more
%  details about the quadratization procedure), that is, it has to be 
%  written as:
%
%       dv/dt = a*v^2 ? w + Isyn(t) + Iapp
%       dw/dt = eps*( alpha*v ? lambda - w)
%       
%       such that 
%           'v'     stands for the membrane potential.
%           'w'     stands for the set of gating variables.
%           'a'     controls the curvature of the v-nullcline.
%           'alpha' controls the slope of the w-nullcline.
%           'eps'   stands for the time scale separation between v and w,
%                   which tends to be small.
%           'lambda'controls the relative displacement between the two
%                   nullclines (the v one and the w one).



%                ESTIMATION PROCEDURE - Required information
%
%   We assume to have a quadratic approximation of a conductance-based neuron  
%   model, as in H.Rotstein (2015). Given the resulting membrane potential  
%   (v) and the course of the gating variable (w), this program estimates 
%   the synaptic current that the neuron is receiving at each time. 
%   Moreover, given the voltage traces for two different applied (steady) 
%   currents and the excitatory and inhibitory reversal potentials, the 
%   program estimates the excitatory and inhibitory conductances separately.
%   Finally, the program gives the option of estimating the synaptic
%   conductance. This conductance can be estimated in two different ways: 
%   (1)if only one voltage trace is given, the synaptic conductance
%   is estimated using the synaptic reversal potential;
%   (2) however, if two voltage traces are given (for two different applied   
%   currents), then the synaptic conductance can be either estimated using 
%   the synaptic reversal potential or the leak conductance.


%                INPUTS AND OUTPUTS
%
% * Input parameters:
%
%   --> Write on the command window: 
%
%           main_Estimation_Conductances(t,v,w,Iapplied)
%
%   Such that:
%
%         t: vector of length m containing the discretization of the time
%            sample.
%
%         v: m x n matrix containing n different samples of membrane
%            potential corresponding at n different values of applied
%            current.
%
%         w: m x n matrix containing n different samples of the set of
%            gating variables corresponding at n different values of applied
%            current.
%
%  Iapplied: vector of length n containing the different values of applied
%            current. If its length is greater than or equal to 2, method
%            discerns between the excitatory and the inhibitory conductances. 
%
%
%
% * Required parameters which are asked while running:
%
%         a: parameter that controls the curvature of the v-nullcline.
%     alpha: parameter that controls the slope of the w-nullcline.
%    lambda: parameter that controls the relative displacement between the 
%            two nullclines (the v one and the w one).
%        vE: excitatory reversal potential (if the Iapplied vector has
%            length greater than or equal to 2).
%        vI: inhibitory reversal potential (if the Iapplied vector has
%            length greater than or equal to 2).
%
%   and, if you want to estimate the synaptic conductance,
%   either
%      vsyn: synaptic reversal potential.
%   or
%        gL: synaptic conductance (only if the Iapplied vector has
%            length greater than or equal to 2).
%
%
% * Output parameters:
%
%   The output parameters are given inside a .mat file called 'estimation.mat'
%   such that
%
%      Isyn: an m x n matrix containing the total synaptic current. Each
%            column corresponds to a different value of applied current.
%
%      gsyn: an m x n matrix containing the total synaptic conductance.
%            Each column corresponds to a different value of applied
%            current.
%
%        gE: vector of length m containing the excitatory synaptic 
%            conductance. It is only computed if, at least, 2 applied
%            currents are considered.
%
%        gI: vector of length m containing the inhibitory synaptic 
%            conductance. It is only computed if, at least, 2 applied
%            currents are considered.
%
%


function [Isyn,gsyn,gE,gI] = main_Estimation_Conductances (t,v,w,Iapplied)


% Checking that the dimensions of the input matrices are consistent
[rv,cv]=size(v);
[rw,cw]=size(w);
[rI,cI]=size(Iapplied);

if (rI~=1) && (cI~=1)
    disp('The applied current must be a vector');
    stop;
end

if (cv~=length(Iapplied))
    disp('The number of columns of the voltage matrix has to be equal');
    disp('to the length of the applied current vector.');
    stop;
end

if (cw~=length(Iapplied))
    disp('The number of columns of the gating variable matrix has to be equal');
    disp('to the length of the applied current vector.');
    stop;
end

if (rv~=rw)
    disp('The size of the gating variable matrix has to be equal');
    disp('to the size of the voltage matrix.');
    stop;
end


disp(   '------ ESTIMATION OF CONDUCTANCES ------'  );

% values for the stellate cell considered in C.Vich and A.Guillamon (2015)
% a=0.1; alpha=0.4; lambda=-0.2;

disp('Write the parameter that controls the curvature of the v-nullcline');
a = input('a = ');

disp('Write the parameter that controls the slope of the w-nullcline.');
alpha = input('alpha = ');

disp('Write the parameter that controls the relative displacement between the two nullclines');
lambda = input('lambda = ');

disp('Estimating the synaptic current from v and w...');
for i=1:length(Iapplied)
    % ESTIMATION procedure
    
    c1=w(1,i);
    t=transpose(t);    
    Isyn(i,1)=-a*v(1,i)^2+exp(-t(1))*(c1+exp(t(1)).*(alpha*v(1,i)-lambda));
    for j=1:(length(t)-1)
        s=t(j) - t(j+1);
        c1e(i,j)=w(j,i)*exp(s);
        trapezi(i,j)=-s*(exp(s)*(alpha*v(j,i)-lambda) + (alpha*v(j+1,i)-lambda))/2;
        Isyn(i,j+1)=-a*v(j+1,i)^2+c1e(i,j)+trapezi(i,j);
    end
end

disp('The synaptic current has been estimated.');

% Estimate the conductances gE(t) and gI(t).
if length(Iapplied)>=2
    disp('Since the applied current vector has length greater than 1,');
    disp('we are going to estimate the excitatory and inhibitory conductances');
    disp('Write the excitatory reversal potential.');
    vE = input('vE = ');
    disp('Write the inhibitory reversal potential.');
    vI = input('vI = ');
    disp('Estimating the excitatory and inhibitory conductances...');
    for j=i:length(t)
        Ivector=[Isyn(1,j)-Iapplied(1); Isyn(2,j)-Iapplied(2)];
        Vmatrix=[v(j,1)-vE v(j,1)-vI ; v(j,2)-vE v(j,2)-vI];
        G=linsolve(-Vmatrix,Ivector);
        gE(j,1)=G(1);
        gI(j,1)=G(2);
    end
else
    gE=[];
    gI=[];
end

for i=1:length(Iapplied)
    Isyn(i,:)=Isyn(i,:)-Iapplied(i);
end

disp('Do you want to estimate the synaptic conductance? Yes (Y) or No (N)');
conditionGsyn = input('gSyn estimation = ');

if strcmp(conditionGsyn,'Y') || strcmp(conditionGsyn,'y')
    if length(Iapplied)>=2
        disp('Do you want to estimate the synaptic conductance using');
        disp('the synaptic reversal potential (vsyn) or the constant leak conductance (gL)');
        opGsyn = input('option to estimate the synaptic conductance = ');
        if strcmp(opGsyn,'gL') || strcmp(opGsyn,'gl')
            disp('Write the constant leak conductance');
            GL = input('gL = ');
            gsyn=gE+gI+GL;
        end
        
    end
    if (length(Iapplied)==1) || strcmp(opGsyn,'vsyn') 
        disp('Write the synaptic reversal potential');
        vsyn = input('vsyn = ');
        disp('Estimating the synaptic conductance...');
        for i=1:length(Iapplied)
            gsyn(:,i)=transpose(Isyn(i,:))./(v(:,i)-vsyn);
        end
    end
    if strcmp(opGsyn,'gL')==0 && strcmp(opGsyn,'gl')==0 && strcmp(opGsyn,'vsyn')==0
        disp('The option to estimate the synaptic conductance is not valid');
        stop;
    end
else
    gsyn=[];
end


save('estimation','Isyn', 'gsyn', 'gE', 'gI');

end

