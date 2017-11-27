clear;
clc;

parameters;
                                       
[f_initial_EAEP,Fs,Ts_EAEP] = Viscoelastic(time_end,delta_t,q,Gshear_eq,Kbulk_eq,...
                                           nVisco,Gshear_neq,Kbulk_neq,tau,...
                                           F0,deltaF);



times = linspace(0.0,time_end,n_steps+1);



figure(5)
set(gca,'FontSize',14)
plot(Fs(:,1,1),Ts_EAEP(:,1,1)-Ts_EAEP(:,3,3),plot_var)
xlabel('Stretch')
ylabel('Cauchy Stress Difference (MPa)')
set(gca,'FontSize',14);
hold on;



function [f_initial,Fs,Ts] = Viscoelastic(time_end,delta_t,q,Gshear_eq,Kbulk_eq,...
                                          nVisco,Gshear_neq,Kbulk_neq,tau,...
                                          F0,deltaF)
    n_steps = time_end / delta_t;
    
    %C = 1.0 / time_end;
    C = 1.e-3;
    % first initialize global arrays to return
    %
    fs = zeros(n_steps+1,1); % volume fractions
    f_initial = zeros(n_steps+1,1); % initial phase volume fractions
    Fs = zeros(n_steps+1,3,3); % deformation gradients from reference
    Ts = zeros(n_steps+1,3,3); % stresses
    
    % arrays for adding viscoelasticity
    %
    Fs_v = zeros(n_steps+1+1,nVisco,3,3); % viscous deformation gradients associated
                                 % with a particular phase
                                 % the first element is the viscous
                                 % deformation gradient associated with the
                                 % initial phase, the subsequent elements
                                 % are associated with the phases formed at
                                 % the i+1th time step
    
    % set initial conditions for viscous deformation gradients
    %
    for m=1:n_steps+1+1
        for n=1:nVisco
            Fs_v(m,n,:,:) = eye(3);
        end
    end
    
    % begin loop over time
    %
    delta_f = C*delta_t;
    fs(1,1) = 0.0;
    for m=2:n_steps+1 
        %
        % now initialize volume fraction array
        %
        if q==0
            fs(m,1) = delta_f;
        elseif q==1
            fs(m,1) = delta_f;
        elseif q==-1
            % this is a test case for no phase evolution
            %
            fs(m,1) = 0.0;
        else
            error('will never be programmed')
        end
    end
    %
    %
    %
    for m=1:n_steps+1
        f_m = sum(fs(1:m,1));
        %
        % new edition
        %
        if f_m >= 1
            f_m = 1;
        end
        %
        if f_m <= 0
            f_m = 0;
        end
        %%%%%%%%%
        if m==1
            F_tau = F0;
        else
            F_tau = deltaF^(m-1)*F0;
        end
        
        
        
        Teq_tau = Neohookean(Gshear_eq,Kbulk_eq,F_tau);
        
        T_total = (1-f_m)*Teq_tau;
        % get previous viscous deformation gradient for the initial
        % phase
        %
        for n=1:nVisco
            Fv_t = eye(3);
            Fv_t(:,:) = Fs_v(1,n,:,:);



            [Fv_tau,Tneq_tau] = ViscousUpdate(delta_t,Gshear_neq(n),Kbulk_neq(n),tau(n),...
                                              F_tau,Fv_t);

            % store the new viscous deformation gradient
            %
            Fs_v(1,n,:,:) = Fv_tau;
            T_total = T_total + (1-f_m)*Tneq_tau;
        end
        %
        % loop over evolving phases
        %
        for i=1:m
            if q==1 && i < m
                fs(i,1) = fs(i,1)*(1-delta_f);
            end
            %
            if fs(i,1) <= 0.0
                fs(i,1) = 0.0;
            end
            %
            %
            if fs(i,1) >= 1
                fs(i,1) = 1;
            end
            %
            F_temp = deltaF^(m-i-1);
            T_temp = Neohookean(Gshear_eq,Kbulk_eq,F_temp);
            
            T_total = T_total + fs(i,1)*T_temp;
            
            for n=1:nVisco
                Fv_temp_t = eye(3);
                Fv_temp_t(:,:) = Fs_v(i+1,n,:,:);

                [Fv_temp_tau,Tneq_temp] = ViscousUpdate(delta_t,Gshear_neq(n),Kbulk_neq(n),...
                                                        tau(n),F_temp,Fv_temp_t);

                Fs_v(i+1,n,:,:) = Fv_temp_tau;                                    

                T_total = T_total + fs(i,1)*Tneq_temp;
            end
        end
        %
        f_initial(m,1) = 1.0-f_m;
        Fs(m,:,:) = F_tau;
        Ts(m,:,:) = T_total;
        %
    end
end
%
%
%
function [Fv_tau,T] = ViscousUpdate(delta_t,G,K,tau,F_tau,Fv_t)
    %
    I = eye(3);
    % calculate elastic trial deformation gradient
    %
    Fe_trial = F_tau * Fv_t^-1;
    detFe_trial = det(Fe_trial);
    
    % perform polar decomposition of the trial elastic deformation
    % gradient
    %
    [Re_trial,Ue_trial,Ve_trial,Ee_trial] = PolarDecomposition(Fe_trial);
    
    trEe_trial = trace(Ee_trial);
    Ee_trial0 = Ee_trial - 1/3 * trEe_trial*I;
    
    % calculate elastic trial stress
    %
    Te_trial = (2*G*Ee_trial0 + K*trEe_trial*I)/detFe_trial;
    
    % Calculate Mandel stress
    %
    Me_trial = detFe_trial*Fe_trial.'*Te_trial*(Fe_trial.')^-1;
    
    Me_trial0 = Me_trial - 1/3 * trace(Me_trial)*I;
    
    %%%%Me_trial0 = Te_trial - 1/3 * trace(Te_trial)*I;
    
    Mbar_trial = 0.0;
    for i=1:3
        for j=1:3
            Mbar_trial = Mbar_trial + Me_trial0(i,j)*Me_trial0(i,j);
        end
    end
    
    Mbar_trial = sqrt(Mbar_trial/2);
    
    gamma_dot = Mbar_trial / (G*(tau+delta_t));
    
    if Mbar_trial <= 0.0
        Nv = zeros(3,3);
    else
        Nv = Me_trial0 / (sqrt(2)*Mbar_trial);
    end
    Dv = sqrt(1/2) * gamma_dot * Nv;
%     if tau < delta_t
%         Dv = zeros(3,3);
%     else
%         Dv = gamma_dot * Nv;
%     end
    
    Fe_tau = expm(-delta_t*Dv)*Fe_trial;
    
    detFe_tau = det(Fe_tau);
    
    [Re_tau,Ue_tau,Ve_tau,Ee_tau] = PolarDecomposition(Fe_tau);
    
    trEe_tau = trace(Ee_tau);
    Ee_tau0 = Ee_tau - 1/3 * trEe_tau * I;
    
    Fv_tau = Fe_tau^-1 * F_tau;
    
    T = (2*G*Ee_tau0 + K*trEe_tau*I)/detFe_tau;
end
%
%
%
function [R,U,V,E] = PolarDecomposition(F)
    

    B = F*F.';
    [B_vec,B_diag] = eig(B);
    
    V = zeros(3,3);
    for i=1:3
        V(i,i) = sqrt(B_diag(i,i));
    end
    V = B_vec*V*B_vec.';
    
    R = F*V^-1;
    
    U = F*R.';
    
    [V_vec,V_diag] = eig(V);
    
    E = zeros(3,3);
    for i=1:3
        E(i,i) = log(V_diag(i,i));
    end
    
    E = V_vec*E*V_vec.';
end
%
%
%
function T = Neohookean(G,K,F)
    J = det(F);
    B = F * F.';
    trB = trace(B);
    
    I = eye(3);
    
    T = J^(-2/3)*G*(B-1/3*trB*eye(3))+K*(J-1);
    %T = G*(B-I);
end