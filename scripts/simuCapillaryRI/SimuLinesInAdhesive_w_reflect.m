%% Assumptions
% Neglect reflection (small angle & close RI)
% Fresnel equations with unpolarized  light sourse

%% Geommetry
R1=0.25e-3; % tapering
R2=0.5e-3;  % inner wall
R3=0.75e-3; % outer wall

% light source
R_src = 5e-2;
src_dis = 19e-2;

obj_dis = 0.6;   % object working distance
D = 0.05;    % pupil
F=0.3;      % focal lens     
camera_dis = obj_dis+1/(1/F-1/obj_dis)+0.000; % image plane
%% Optics
n1=1.56;
n2=1.57;
n3=1.536;
n4=1.56;
RI = [n1 n2 n3 n4];

incident_angle_range =[0 0];
incident_angle_h = pi/(10*180);
%%
incident_sample=linspace(-R_src,R_src,50000);
N=length(incident_sample);

% final result
final_intensity = zeros(1,1000);
M=length(final_intensity);
FoV = linspace(-1.1*R3,1.1*R3,M); %field of view
%%
for alpha =incident_angle_range(1):incident_angle_h:incident_angle_range(2)
    % Light meet the interface
    cross_point_x=nan(N,8);
    cross_point_x(:,1)=-src_dis;
    cross_point_y=nan(N,8);
    cross_point_y(:,1)=incident_sample;
    
    % reaching the lens
    lens_y = incident_sample;
    lens_angle = alpha+0*incident_sample;
    lens_Intensity = ones(size(incident_sample));
    image_y = nan(N,1);
    for k=1:1:N      
        c_angle = alpha; % current angle
        c_x = cross_point_x(k,1); % current x position
        c_y = cross_point_y(k,1); % current y position
        cross_num = 2;     
        %% hit R3?
        a=1; 
        b=-2*sin(c_angle)^2*c_x+2*sin(c_angle)*cos(c_angle)*c_y;
        c=sin(c_angle)^2*c_x^2+cos(c_angle)^2*c_y^2-2*sin(c_angle)*cos(c_angle)*c_x*c_y...
            -R3^2*cos(c_angle)^2;
        if b^2-4*a*c>0
            tmpx = (-b-(b^2-4*a*c)^0.5)/2;
            tmpy = tan(c_angle)*(tmpx-c_x)+c_y;
            c_x = tmpx;
            c_y = tmpy;
            cross_point_x(k,cross_num) = tmpx;
            cross_point_y(k,cross_num) = tmpy;
            Thetai= atan2(-tmpy,-tmpx)-c_angle;
            [Thetat, Transmission] = Fresnel(RI(4),RI(3),Thetai);
            if Transmission <0.2
                Thetat = pi-Thetai;
                if cos(atan2(-tmpy,-tmpx)-Thetat)<=0
                    lens_y(k) = inf;
                    lens_angle(k)=0;
                    lens_Intensity(k)=0;
                continue;
                end
                Transmission = 1-Transmission;
            end
            lens_Intensity(k)=lens_Intensity(k)*Transmission;
            c_angle = atan2(-tmpy,-tmpx)-Thetat;
            cross_num=cross_num+1;
            %% hit R2?
             a=1; 
             b=-2*sin(c_angle)^2*c_x+2*sin(c_angle)*cos(c_angle)*c_y;
             c=sin(c_angle)^2*c_x^2+cos(c_angle)^2*c_y^2-2*sin(c_angle)*cos(c_angle)*c_x*c_y...
                -R2^2*cos(c_angle)^2;
            if b^2-4*a*c>0
                tmpx = (-b-(b^2-4*a*c)^0.5)/2;
                tmpy = tan(c_angle)*(tmpx-c_x)+c_y;
                c_x = tmpx;
                c_y = tmpy;
                cross_point_x(k,cross_num) = tmpx;
                cross_point_y(k,cross_num) = tmpy;
                Thetai= atan2(-tmpy,-tmpx)-c_angle;
                [Thetat, Transmission] = Fresnel(RI(3),RI(2),Thetai);
                if Transmission <0.2
                    Thetat = pi-Thetai;
                    if cos(atan2(-tmpy,-tmpx)-Thetat)<=0
                        lens_y(k) = inf;
                        lens_angle(k)=0;
                        lens_Intensity(k)=0;
                    continue;
                    end
                    Transmission = 1-Transmission;
                end
                lens_Intensity(k)=lens_Intensity(k)*Transmission;
                c_angle = atan2(-tmpy,-tmpx)-Thetat;
                cross_num=cross_num+1;
                %% hit R1?
                a=1; 
                b=-2*sin(c_angle)^2*c_x+2*sin(c_angle)*cos(c_angle)*c_y;
                c=sin(c_angle)^2*c_x^2+cos(c_angle)^2*c_y^2-2*sin(c_angle)*cos(c_angle)*c_x*c_y...
                    -R1^2*cos(c_angle)^2;
                if b^2-4*a*c>0
                    tmpx = (-b-(b^2-4*a*c)^0.5)/2;
                    tmpy = tan(c_angle)*(tmpx-c_x)+c_y;
                    c_x = tmpx;
                    c_y = tmpy;
                    cross_point_x(k,cross_num) = tmpx;
                    cross_point_y(k,cross_num) = tmpy;
                    Thetai= atan2(-tmpy,-tmpx)-c_angle;
                    [Thetat, Transmission] = Fresnel(RI(2),RI(1),Thetai);
                    if Transmission <0.2
                        Thetat = pi-Thetai;
                        if cos(atan2(-tmpy,-tmpx)-Thetat)<=0
                            lens_y(k) = inf;
                            lens_angle(k)=0;
                            lens_Intensity(k)=0;
                        continue;
                        end
                        Transmission = 1-Transmission;
                    end
                    lens_Intensity(k)=lens_Intensity(k)*Transmission;
                    c_angle = atan2(-tmpy,-tmpx)-Thetat;
                    cross_num=cross_num+1;
                    %% Return to Region 2
                    a=1; 
                    b=-2*sin(c_angle)^2*c_x+2*sin(c_angle)*cos(c_angle)*c_y;
                    c=sin(c_angle)^2*c_x^2+cos(c_angle)^2*c_y^2-2*sin(c_angle)*cos(c_angle)*c_x*c_y...
                        -R1^2*cos(c_angle)^2;
                    tmpx = (-b+(b^2-4*a*c)^0.5)/2;
                    tmpy = tan(c_angle)*(tmpx-c_x)+c_y;
                    c_x = tmpx;
                    c_y = tmpy;
                    cross_point_x(k,cross_num) = tmpx;
                    cross_point_y(k,cross_num) = tmpy;
                    Thetai= atan2(tmpy,tmpx)-c_angle;
                    [Thetat, Transmission] = Fresnel(RI(1),RI(2),Thetai);
                    if Transmission <0.2
                        Thetat = pi-Thetai;
                        if cos(atan2(-tmpy,-tmpx)-Thetat)<=0
                            lens_y(k) = inf;
                            lens_angle(k)=0;
                            lens_Intensity(k)=0;
                        continue;
                        end
                        Transmission = 1-Transmission;
                    end
                    lens_Intensity(k)=lens_Intensity(k)*Transmission;
                    c_angle = atan2(tmpy,tmpx)-Thetat;
                    cross_num=cross_num+1;
                end
                %% Return to Region 3
                a=1; 
                b=-2*sin(c_angle)^2*c_x+2*sin(c_angle)*cos(c_angle)*c_y;
                c=sin(c_angle)^2*c_x^2+cos(c_angle)^2*c_y^2-2*sin(c_angle)*cos(c_angle)*c_x*c_y...
                    -R2^2*cos(c_angle)^2;
                tmpx = (-b+(b^2-4*a*c)^0.5)/2;
                tmpy = tan(c_angle)*(tmpx-c_x)+c_y;
                c_x = tmpx;
                c_y = tmpy;
                cross_point_x(k,cross_num) = tmpx;
                cross_point_y(k,cross_num) = tmpy;
                Thetai= atan2(tmpy,tmpx)-c_angle;
                [Thetat, Transmission] = Fresnel(RI(2),RI(3),Thetai);
                if Transmission <0.2
                    Thetat = pi-Thetai;
                    if cos(atan2(-tmpy,-tmpx)-Thetat)<=0
                        lens_y(k) = inf;
                        lens_angle(k)=0;
                        lens_Intensity(k)=0;
                    continue;
                    end
                    Transmission = 1-Transmission;
                end
                lens_Intensity(k)=lens_Intensity(k)*Transmission;
                c_angle = atan2(tmpy,tmpx)-Thetat;
                cross_num=cross_num+1;
            end
            %% Return to region 4
            a=1; 
            b=-2*sin(c_angle)^2*c_x+2*sin(c_angle)*cos(c_angle)*c_y;
            c=sin(c_angle)^2*c_x^2+cos(c_angle)^2*c_y^2-2*sin(c_angle)*cos(c_angle)*c_x*c_y...
                -R3^2*cos(c_angle)^2;
            tmpx = (-b+(b^2-4*a*c)^0.5)/2;
            tmpy = tan(c_angle)*(tmpx-c_x)+c_y;
            c_x = tmpx;
            c_y = tmpy;
            cross_point_x(k,cross_num) = tmpx;
            cross_point_y(k,cross_num) = tmpy;
            Thetai= atan2(tmpy,tmpx)-c_angle;
            [Thetat, Transmission] = Fresnel(RI(3),RI(4),Thetai);
            if Transmission <0.2
                Thetat = pi-Thetai;
                if cos(atan2(-tmpy,-tmpx)-Thetat)<=0
                    lens_y(k) = inf;
                    lens_angle(k)=0;
                    lens_Intensity(k)=0;
                continue;
                end
                Transmission = 1-Transmission;
            end
            lens_Intensity(k)=lens_Intensity(k)*Transmission;
            c_angle = atan2(tmpy,tmpx)-Thetat;
            cross_num=cross_num+1;
        end
        %% Hit the lens
        tmpx = obj_dis;
        tmpy = tan(c_angle)*(tmpx-c_x)+c_y;
        lens_Intensity(k)=1*lens_Intensity(k);
        cross_point_x(k,cross_num) = tmpx;
        cross_point_y(k,cross_num) = tmpy;
        lens_y(k)=tmpy;
        %% Imaging
        if abs(tmpy)<D/2
            image_y(k)=tmpy+(tan(c_angle)*F-tmpy)*(camera_dis-obj_dis)/F;
        end
    end
    %%
    temp_intensity = zeros(1,M);
    Smear = abs(diff(image_y))/mean(diff(incident_sample));
     lens_int = (lens_Intensity(1:end-1)+lens_Intensity(2:end))/2;
    for l=1:1:M
        Interval = image_y-FoV(l);
        IsBelong = Interval(1:(end-1)).*Interval(2:end);
        IsBelong = (IsBelong<=0);   
        temp_intensity(l)=sum(lens_int(IsBelong)'./Smear(IsBelong));
    end
    final_intensity = final_intensity+temp_intensity*cos(alpha);
end
final_intensity=final_intensity/length(incident_angle_range(1):incident_angle_h:incident_angle_range(2));
%%
subplot(2,1,1)
image(final_intensity*256)
colormap(gray(256))
subplot(2,1,2)
image(final_intensity*150)
colormap(gray(256))
  