%% Initialization

%% Array-matrix
row_elements = 8;
column_elements = 8;
r_prime = zeros(3,row_elements*column_elements);

%% Evenly spaced matrix initialization
%matrix array in the xy-plane
uni_distance = 2*pi/2;      % Give length in kd
%uni_distance = 2*pi*lambda_d;
c = 343;
frequency = 3805;
%uni_distance = 2*pi*frequency/c*0.02;
lambda = c/frequency;
lambda_rel = 0.02/lambda

element_index = 1;
for i = 1:row_elements
    for j = 1:column_elements
        r_prime(1,element_index) = i*uni_distance;
        r_prime(2,element_index) = j*uni_distance;
        element_index = element_index + 1;
    end
end

%set origin of the matrix to (0,0)
r_prime(1,:) = r_prime(1,:)-row_elements*uni_distance/2 - uni_distance/2;
r_prime(2,:) = r_prime(2,:)-column_elements*uni_distance/2 - uni_distance/2;

figure(1);
plot(r_prime(1,:),r_prime(2,:),'linestyle','none','marker','*');

%% Phase and amplitude correction
theta_listen_deg = 0;
phi_listen_deg = 0;

theta_l = theta_listen_deg/180*pi;
phi_l = phi_listen_deg/180*pi;
phase_corr = zeros(1,row_elements*column_elements);

element_index = 1;
for i = 1:row_elements*column_elements
    
    phase_corr(element_index) =  exp(-1i*(r_prime(1,i)*sin(theta_l)*cos(phi_l) + ...
            r_prime(2,i)*sin(theta_l)*sin(phi_l)));
        element_index = element_index + 1;
end

phase_amplitudes = zeros(1,row_elements*column_elements);


%Binomial distribution      Reduces the sidelobe levels drastically
element_index = 1;
for i = 0:row_elements-1
    row_amp = nchoosek(row_elements-1,i);
    for j = 0:column_elements-1
        column_amp = nchoosek(column_elements-1,j);
        phase_amplitudes(element_index) = row_amp*column_amp;
        element_index = element_index +1;
    end
end
%phase_corr = phase_amplitudes.*phase_corr;

%On-Off amplitude distribution
phase_amplitudes_OF = zeros(1,row_elements*column_elements);
mode = 1;
row_lim = ceil((row_elements)/mode);
column_lim = ceil((column_elements)/mode);
test = 0;
for i = 1:row_lim
    for j = 1:column_lim
        element_index = (mode*(i-1))*row_elements + mode*(j-1) +1;
        phase_amplitudes_OF(element_index) = 1;
    end
end

phase_corr = phase_amplitudes_OF.*phase_corr;
figure(39)
plot(r_prime(1,:).*phase_amplitudes_OF,r_prime(2,:).*phase_amplitudes_OF,'linestyle','none','marker','*');

%% Get far-field
resolution = 500;
[AF,theta,phi] = array_factor(r_prime,phase_corr,resolution);
%AF = AF_test;

%% Radiated power
antenna_element = 1;                
power = abs((AF.*antenna_element)).^2;    %P = |AF|^2
integrand = ((sin(theta))').*power;     %Integral to be calculated: S |AF|^2 sin(theta)dthetadphi

temp = cumtrapz(phi,cumtrapz(theta,integrand',2));
power_rad = (temp(end,end) - temp(1,1))/(2*pi);

AF_dbi = 10*log10(power/power_rad);     %Power in terms of isotropic radiation

% Check that normalized power is indeed NORMALIZED!
power_norm = power/power_rad;
integrand_check = ((sin(theta))').*power_norm;
temp_check = cumtrapz(phi,cumtrapz(theta,integrand_check',2));
check = (temp_check(end,end) - temp_check(1,1))/(2*pi);

%% Half power beamwidth
min_value = -40;     %Lowest dB value to plot
%Set listen-angle to [0,0]
max_dB = max(max(AF_dbi));
beamwidth = 0;
theta_half = zeros(2,1);
theta_half_db = zeros(2,1);
theta_memory_ind = 1;
db_limits_found = 1;

for theta_ind = 2:length(theta)
    if (AF_dbi(theta_memory_ind) < max_dB-3  && AF_dbi(theta_ind) > max_dB-3)
        theta_half(db_limits_found) = theta(theta_ind);
        theta_half_db(db_limits_found) = AF_dbi(theta_ind);
        db_limits_found = db_limits_found+1;
    end
    if(max_dB-3 < AF_dbi(theta_memory_ind,1) && AF_dbi(theta_ind) < max_dB-3)
        theta_half(db_limits_found) = theta(theta_memory_ind);
        theta_half_db(db_limits_found) = AF_dbi(theta_memory_ind,1);
        db_limits_found = db_limits_found+1;
    end
    theta_memory_ind = theta_memory_ind+1;
end



if db_limits_found == 2
    beamwidth = 2*theta_half(1)*180/pi;
    theta_x_lines = [theta_half(1) theta_half(1)]*180/pi;
    theta_y_lines = [min_value theta_half_db(1)];
else 
    beamwidth = abs(theta_half(1) - theta_half(2)) * 180/pi;
    theta_x_lines = [theta_half(1) theta_half(2);theta_half(1) theta_half(2) ]*180/pi;
    theta_y_lines = [min_value min_value ;theta_half_db(1) theta_half_db(2) ];
end

directivity = max_dB;

%% Plotting

figure(2)
plot(phi/pi*180,AF_dbi(250,:),'color','b');
xlim([0 360])
ylim([min_value directivity])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\varphi$ (deg)','Interpreter','latex','FontSize',20);
ylabel('$AF(45,\varphi)|^2$ (dBi)','Interpreter','latex','FontSize',20)
%legend('MATLAB−sim','CST−sim','Interpreter','latex','location','southeast')
grid on
set(gca,'GridALpha',0.1,'LineWidth',.3);

figure(3)
plot(theta/pi*180,AF_dbi(:,1),'color','b','linewidth',1);
hold on 
line(theta_x_lines, theta_y_lines,'color','k','linestyle','--'...
    ,'linewidth',1)
annotation('textbox', [0.5, 0.8, 0.35, 0.1], 'String', "Beamwidth: " +beamwidth...
    + "$^{\circ}$" + newline + "Directivity: " +directivity + " dB" ,'Interpreter',...
    'latex','FontSize',11,'BackgroundColor','w')
xlim([0 90])
ylim([min_value directivity])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\theta$ (deg)','Interpreter','latex','FontSize',20);
ylabel('$|AF(\theta,0)|^2$ (dBi)','Interpreter','latex','FontSize',20)
%legend('MATLAB−sim','CST−sim','Interpreter','latex','location','southeast')
set(gca,'GridALpha',.5,'LineWidth',.3);
grid on

% X = sin(theta)'*cos(phi);
% Y = sin(theta)'*sin(phi);
% 
% figure(4)
% surf(X,Y,abs(AF))

figure(5)
[x,y,z,x_dbi,y_dbi,z_dbi] = farfield(AF,theta,phi,-20);
linear3d_max = max(max([x,y,z]));
dB3d_max = max(max([x_dbi,y_dbi,z_dbi]));
s = surf(x,y,z);
s.EdgeColor = 'interp';
%s.FaceLighting = 'gouraud';
s.FaceColor = 'flat';
xlim([-linear3d_max linear3d_max])
ylim([-linear3d_max linear3d_max])
zlim([-linear3d_max linear3d_max])
hold on
s_matrix_array = surf(r_prime(1,1:row_elements:row_elements^2),r_prime(2,1:column_elements),zeros(row_elements,column_elements));

figure(6)
s = surfl(x_dbi,y_dbi,z_dbi);
s.EdgeColor = 'interp';
%s.FaceLighting = 'gouraud';
s.FaceColor = 'flat';
xlim([-dB3d_max dB3d_max])
ylim([-dB3d_max dB3d_max])
zlim([-dB3d_max dB3d_max])
hold on
s_matrix_array2 = surf(r_prime(1,1:row_elements:row_elements^2),r_prime(2,1:column_elements),zeros(row_elements,column_elements));