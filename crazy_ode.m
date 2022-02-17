function rhs = crazy_ode(x, u, pars)

    % Nonlinear ODE model for drone.
    %
    % The position and velocity refer to the inertial coordinate system.
    
    phi = x(7);   % Cardanian angles
    theta = x(8);
    psi = x(9);
    
    s_phi = sin(phi);
    c_phi = cos(phi);
    s_theta = sin(theta);
    c_theta = cos(theta);
    t_theta = s_theta / c_theta;
    s_psi = sin(psi);
    c_psi = cos(psi);
    
    % Rotation matrix
    R = [c_theta*c_psi s_phi*s_theta*c_psi-s_phi*s_psi c_phi*s_theta*c_psi+s_phi*s_psi; ...
         c_theta*s_psi s_phi*s_theta*s_psi+c_phi*c_psi c_phi*s_theta*s_psi-s_phi*c_psi; ...
         -s_theta      s_phi*c_theta                   c_phi*c_theta ];
     
    % Velocity in inertial coord. system
    v_x = x(4);
    v_y = x(5);
    v_z = x(6);
     
    % Velocitzy in body coord. system
    v_B = R'*x(4:6);
    
    % Drag force in body coord. system
    F_D_B = diag([pars.d_x pars.d_y pars.d_z])*v_B;
    
    % Drag force in inertial coord. system
    F_D = R*F_D_B;
    
    % Angular velocities
    om_x = x(10); 
    om_y = x(11);
    om_z = x(12);
    
    % Control variables
    F = u(1)+u(2)+u(3)+u(4);                % Total vertical force
    T_x = (u(1)-u(2)+u(3)-u(4))*pars.b;     % Torque about x-axis
    T_y = (-u(1)-u(2)+u(3)+u(4))*pars.b;    % Torque about y-axis
    T_z = (u(1)-u(2)-u(3)+u(4))*pars.kappa; % Torque about z-axis
    
    J_xx = pars.J_xx;
    J_yy = pars.J_yy;
    J_zz = pars.J_zz;
    
    v_x_dot = (F * (c_phi*s_theta*c_psi + s_phi*s_psi) - F_D(1)) / pars.m;
    v_y_dot = (F * (c_phi*s_theta*s_psi - s_phi*c_psi) - F_D(2)) / pars.m;
    v_z_dot = (F* c_phi*c_theta - pars.F_g - F_D(3)) / pars.m;
    
    om_x_dot = (T_x + (J_yy - J_zz) * om_y*om_z) / J_xx;
    om_y_dot = (T_y + (J_zz - J_xx) * om_x*om_z) / J_yy;
    om_z_dot = (T_z + (J_xx - J_yy) * om_x*om_y) / J_zz;
    
    phi_dot   = om_x + (s_phi*om_y + c_phi*om_z)*t_theta;
    theta_dot = c_phi*om_y - s_phi*om_z;
    psi_dot   = (s_phi*om_y + c_phi*om_z) / c_theta;
   
    rhs = [v_x; v_y; v_z; v_x_dot; v_y_dot; v_z_dot; ...
           phi_dot; theta_dot; psi_dot; om_x_dot; om_y_dot; om_z_dot];

end %function