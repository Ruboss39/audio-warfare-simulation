%% Initialization

%% Array-matrix
row_elements = 8;
column_elements = 8;
r_prime = zeros(3,row_elements*column_elements);

%% Evenly spaced matrix initialization
%matrix array in the xy-plane
%uni_distance = 4*pi/2;      % Give length in kd
%uni_distance = 2*pi*lambda_d;
c = 343;
frequency = 2000;
uni_distance = 2*pi*frequency/c*0.02;
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

%% Array creation v2
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;0.3010 0.7450 0.9330];

uni_a_d = 120;
d_factor = c/(2*pi*frequency);

r_a1 = 2*pi*frequency/c*10^-3 * [uni_a_d;uni_a_d;0];
r_a2 = 2*pi*frequency/c*10^-3 * [-uni_a_d;uni_a_d;0];
r_a3 = 2*pi*frequency/c*10^-3 * [uni_a_d;-uni_a_d;0];
r_a4 = 2*pi*frequency/c*10^-3 * [-uni_a_d;-uni_a_d;0];

uni_distance_m = 2;

x_end = d_factor*(column_elements-4)*uni_distance;
y_end = d_factor*(row_elements-4)*uni_distance;

array_ends = [-x_end,x_end,x_end,-x_end,-x_end;...
    y_end,y_end,-y_end,-y_end,y_end];

array_shift = d_factor*[r_a1,r_a2,r_a3,r_a4];
x_min_ = -x_end - 10^-3 * uni_a_d;
x_max_ = x_end + 10^-3 * uni_a_d;
y_min_ = -y_end - 10^-3 * uni_a_d;
y_max_ = y_end + 10^-3 * uni_a_d;

array_matrix = matrix_array(r_a1,uni_distance,row_elements,column_elements);
array_matrix2 = matrix_array(r_a2,uni_distance,row_elements,column_elements);
array_matrix3 = matrix_array(r_a3,uni_distance,row_elements,column_elements);
array_matrix4 = matrix_array(r_a4,uni_distance,row_elements,column_elements);

array_matrices = [array_matrix,array_matrix2,array_matrix3,array_matrix4];
%r_prime = [array_matrix.r_prime , array_matrix2.r_prime , array_matrix3.r_prime, array_matrix4.r_prime];

sub_arrays = length(array_matrices);

figure(100);      %Plot the geometry in the xy-plane

%array_mode_shifts = [1,1;0,1;1,0;0,0];
array_mode_shifts = [0,0;1,0;0,1;1,1];
%array_mode_shifts = [0,0;0,0;0,0;0,0];
array_mode = 1;

r_prime_cm = [];

for array = 1:sub_arrays
    r_prime_mode = array_matrices(array).config_mode(array_mode,array_mode_shifts(array,1),...
        array_mode_shifts(array,2));
    r_prime_cm = [r_prime_cm,r_prime_mode];
    z_array_coord_mode = zeros(length(r_prime_mode(1,:)));
    plot3(d_factor *r_prime_mode(1,:),d_factor * r_prime_mode(2,:),...
        z_array_coord_mode,'linestyle','none','marker','o',...
        'MarkerSize',4,'color',colors(array,:),'MarkerFaceColor',colors(array,:));
    hold on

    z_array_coord = zeros(length(array_matrices(array).r_prime(1,:)));
    plot3(d_factor * array_matrices(array).r_prime(1,:),d_factor * array_matrices(array).r_prime(2,:),...
        z_array_coord,'linestyle','none','marker','o','MarkerSize',4,'color',colors(array,:));
    hold on
    plot3((array_ends(1,:) + array_shift(1,array)),(array_ends(2,:) + array_shift(2,array)),[0,0,0,0,0],'color',colors(array,:),'linewidth',0.6);
    axis square
    zlim([-0.5 0.5])
    hold on
end
line_1_x =  [x_end - 10^-3 *uni_a_d, x_end - 10^-3 *uni_a_d];
line_1_y =  [-y_end - 10^-3 *uni_a_d, -y_end - 10^-3 *uni_a_d];
line_1_z = [-0.5 0];

line_2_x =  [-x_end + 10^-3 *uni_a_d, -x_end + 10^-3 *uni_a_d];
line_2_y =  [-y_end - 10^-3 *uni_a_d, -y_end - 10^-3 *uni_a_d];
line_2_z = [-0.5 0];

plot3(line_1_x,line_1_y,line_1_z,'color','k','linestyle','--'...
    ,'linewidth',0.8)
hold on
plot3(line_2_x,line_2_y,line_2_z,'color','k','linestyle','--'...
    ,'linewidth',0.8)
hold on

% a_x = [0 0.4];
% a_y = [0 0];
% annotation('doublearrow',a_x,a_y,'String','y = x ')
e_factor = 0.00006;
c_factor = e_factor * (uni_a_d-80);
an = annotation('doublearrow');
an.X = [0.596-c_factor, 0.675 + c_factor];
an.Y =  [0.28-c_factor,  0.315 + c_factor];

dim = [.595 .1 .3 .3];
str = 'd_a';
annotation('textbox',dim,'String','$d_{array}$','Interpreter','latex','FontSize',12,'LineStyle','none')


set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlim([x_min_ x_max_])
ylim([y_min_ y_max_])
xlabel('x (m)','Interpreter','latex','FontSize',11);
ylabel('y (m)','Interpreter','latex','FontSize',11)
title(strcat('$d_{array} = $ ',{' '} ,int2str(2*(uni_a_d-80)),{''},'mm'),'Interpreter','latex','FontSize',15);
%title('Multiple arrays','Interpreter','latex','FontSize',15);
grid on
%exportgraphics(gcf,'multiple_sub_arrays.png','Resolution',300)

%r_prime = [array_matrix.r_prime , array_matrix2.r_prime , array_matrix3.r_prime, array_matrix4.r_prime];
r_prime = r_prime_cm;
%% Phase and amplitude correction
elements = length(r_prime(1,:));
theta_listen_deg = 0;
phi_listen_deg = 0;
anti_beam = 0;

theta_l = theta_listen_deg/180*pi;
phi_l = phi_listen_deg/180*pi;
phase_corr = zeros(1,elements);

element_index = 1;
for i = 1:elements
    
    phase_corr(element_index) = exp(-1i*(r_prime(1,i)*sin(theta_l)*cos(phi_l) + ...
            r_prime(2,i)*sin(theta_l)*sin(phi_l)));
        element_index = element_index + 1;
end

phase_amplitudes = zeros(1,elements);


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
%phase_corr = (phase_amplitudes).*phase_corr;

%On-Off amplitude distribution
phase_amplitudes_OF = zeros(1,elements);
mode = 1;
row_lim = ceil((row_elements)/mode);
column_lim = ceil((column_elements)/mode);
even_odd = -1;
element_size = 0;

anti_beam_mode = 1-mod(column_elements,2);

for i = 1:row_lim
    for j = 1:column_lim
        element_index = (mode*(i-1))*row_elements + mode*(j-1) +1;
        %phase_amplitudes_OF(element_index) = (-2*mod(i,1 + anti_beam_mode) + 1)*even_odd;
        phase_amplitudes_OF(element_index) = 1;
        element_size = element_size +1;

        if anti_beam == 1
            even_odd = -1*even_odd;
        end
    end
end


positive_phase = zeros(3,element_size);
negative_phase = zeros(3,element_size);

for i = 1:length(r_prime(1,:))
    if phase_amplitudes_OF(i) == -1
        negative_phase(1,i) = r_prime(1,i);
        negative_phase(2,i) = r_prime(2,i);
    elseif  phase_amplitudes_OF(i) == 1
        positive_phase(1,i) = r_prime(1,i);
        positive_phase(2,i) = r_prime(2,i);
    end
end

%phase_corr = phase_amplitudes_OF.*phase_corr;
figure(39)
plot(positive_phase(1,:),positive_phase(2,:),'linestyle','none','marker','*','color','b');
hold on
plot(negative_phase(1,:),negative_phase(2,:),'linestyle','none','marker','*','color','r');

%% Get far-field
resolution = 800;
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

directivity1 = max_dB;

%% Plotting
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;0.3010 0.7450 0.9330];

theta_index = 1 +ceil(theta_listen_deg/90 * resolution);

figure(2)
plot(phi/pi*180,AF_dbi(theta_index,:),'color','b');
xlim([0 360])
ylim([min_value directivity1])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\varphi$ (deg)','Interpreter','latex','FontSize',18);
ylabel('$AF(\theta_0,\varphi)|^2$ (dBi)','Interpreter','latex','FontSize',18)
%legend('MATLAB−sim','CST−sim','Interpreter','latex','location','southeast')
xticks([0 30 60 90 120 150 180 210 240 270 300 330 360])
grid on
set(gca,'GridALpha',0.1,'LineWidth',.3);

figure(3)
plot(theta/pi*180,AF_dbi(:,1),'color','b','linewidth',1);
hold on
plot(-theta/pi*180,AF_dbi(:,1),'color','b','linewidth',1);
hold on 
line(theta_x_lines, theta_y_lines,'color','k','linestyle','--'...
    ,'linewidth',0.8)
hold on
line(-theta_x_lines, theta_y_lines,'color','k','linestyle','--'...
    ,'linewidth',0.8)
annotation('textbox', [0.55, 0.2, 0.3, 0.1], 'String', "Beamwidth: " +beamwidth...
    + "$^{\circ}$" + newline + "Directivity: " +directivity1 + " dB" ,'Interpreter',...
    'latex','FontSize',11,'BackgroundColor','w')
% hold on
% x_ = [0.3 0.2];
% y_ = [0.85 0.89];
% annotation('textarrow',x_,y_,'String','grating lobe ')
% hold on
% x__ = [0.7 0.8];
% y__ = [0.85 0.89];
% annotation('textarrow',x__,y__,'String','grating lobe ')
xlim([-90 90])
ylim([min_value directivity1])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\theta$ (deg)','Interpreter','latex','FontSize',20);
ylabel('$|AF(\theta,0)|^2$ (dBi)','Interpreter','latex','FontSize',20)
%legend('MATLAB−sim','CST−sim','Interpreter','latex','location','southeast')
xticks([-90 -60 -30 0 30 60 90]);
grid on
set(gca,'GridALpha',0.1,'LineWidth',.3);
%exportgraphics(gcf,'4x4_array_example_2.png','Resolution',300)

%% la
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
zlim([0 linear3d_max])
axis equal
xticks([])
yticks([])
zticks([])
hold on
%s_matrix_array = surf(r_prime(1,1:row_elements:row_elements^2),r_prime(2,1:column_elements),zeros(row_elements,column_elements));
s_matrix_array = plot3(r_prime(1,:),r_prime(2,:),zeros(1,length(r_prime(1,:))),...
    'linestyle','none','marker','o','MarkerFaceColor',colors(2,:),'linewidth',...
    1.3,'MarkerSize',4,'color',colors(2,:));
%exportgraphics(gcf,'4x4_array_example_2_3d.png','Resolution',300)

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