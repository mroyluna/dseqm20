function schumann_all_zloop_sine(z, Lx, tmax, dx, dt, startind, period)
%SCHUMANN non-dimensional equations, only depend on z

dirpath  = ['Figures_z_', num2str(z), '_Lx_' num2str(Lx) '_tmax_', num2str(tmax), '_sine_', num2str(period)];
mkdir(dirpath);

nout     = 100;
nstreams = 10;

x    = [0:dx:Lx];
n    = length(x);

vv   = 1./z;

% initialize
T_solid = zeros(1, n); 
T_fluid = zeros(1, n);
T_f_try = T_fluid;
T_diff  = zeros(1, n);
T_solid_steps = [];
T_fluid_steps = [];
t_steps       = [];
legendtext    = [];
dT_fluid_dx   = gradient(T_fluid, dx); % starting gradient in fluid
delta_T       = T_fluid - T_solid;     % starting contrast

%Value of initial temperature prior to perturbation
T_in   = 0;
% Value of temp of fluid at left edge
T_0    = 0; %initial inlet temp

h1 = figure(1); grid on; box on; hold on
%h1 = subplot(211); grid on; box on; hold on
%h2 = subplot(212); grid on; box on; hold on
solid_plot = plot(x, T_solid);
fluid_plot = plot(x, T_fluid);
set(solid_plot, 'linewidth',[1]);
set(fluid_plot, 'linewidth',[1]);
ylim([-0.1,1]);
ylabel('Temperature/{\Delta}T'); xlabel('x');
set(gca,'fontname','helvetica','fontsize', [14])

diff_plot     = plot(x, T_diff);
[mdiff, indm] = max(T_diff);
mx(1,:)       = [0, x(indm)];
lm            = [mx 0; mx T_0];
max_plot      = plot(lm(:,1), lm(:,2));
set(max_plot,'linestyle','--','color','k')
legend('{T''}_s', '{T''}_f','{T''}_f - {T''}_s');

k    = 1;
step = 0;
maxiter = length([dt:dt:tmax]);
outiter = floor(maxiter/nout);

for iter = 1:maxiter 
    step  = step + 1;
    t     = dt*iter;
    
    front = vv*t; 
    in    = x <= front; % points that the fluid has passed so far
    
    % To eliminate numerical instabilities, impose a boundary condition on
    % all new points which have been reached since the last time step.
 
    new_points = in & x > vv * (t - dt);
    
    T_fluid(new_points) = T_0*exp(-front);
    
    % take an intermediate step to evolve the fluid temp 
    dtsub = dt*0.5;
    dT_fluid_dt(in) = - delta_T(in) - dT_fluid_dx(in);
    T_f_try(in) = T_fluid(in) + vv * dT_fluid_dt(in) * dtsub;
    % evolve solid temp at int step
    T_solid_try(in) = T_solid(in) + delta_T(in) * dtsub;
    % find fluid gradient at substep
    dT_fluid_dx_try(in) = gradient(T_f_try(in), dx);
    delta_T(in)         = T_f_try(in) - T_solid_try(in);
    
    % use values after intermediate step in order to take next step
    dT_fluid_dt(in) = - delta_T(in) - dT_fluid_dx_try(in);
    
    % take full step
    T_fluid(in) = T_fluid(in) + vv * dT_fluid_dt(in) * dt;
    T_solid(in) = T_solid(in) + delta_T(in) * dt;
    
    %T_fluid(in) = 0.15*T_f_try(in)+(1-0.15)*T_fluid(in);
    delta_T(in)     = T_fluid(in) - T_solid(in);
    dT_fluid_dx(in) = gradient(T_fluid(in), dx);
    
   
    T_fluid(1)  = 0.5*(1+sin(2*pi*t/period + pi/2));  
    % boundary condition on Tf only at left
    if iter == 1 
        T_fluid(2) = 0.5;
    end
    T_diff = T_fluid - T_solid;
    
    if (mod(iter, outiter) == 0)
        set(solid_plot, 'YDATA', T_solid);
        set(fluid_plot, 'YDATA', T_fluid);
        set(diff_plot,  'YDATA', T_diff);
        [mdiff, indm] = max(T_diff);
        mx        = [mx; t x(indm)];
        lm        = [x(indm) 0; x(indm) T_0];
        set(max_plot, 'XDATA', lm(:,1), 'YDATA', lm(:,2));
        drawnow
        tl = ['t = ' num2str(t)];
        title(tl);
        fn = [dirpath '/Temp_x_' num2str(t) '.png'];
        %print(fn, '-dpng')
	    fn = [dirpath '/Temp_x_' num2str(t) '.fig'];
	    savefig(fn)
        
        % save a growing array of all solutions through time
        Tf(k,:) = T_fluid;
        Ts(k,:) = T_solid;
        k       = k + 1;
    end
    
    % nstreams sets the level of decimation
    if mod(step, floor(1 / dt / nstreams)) == 0
        t;
        t_steps = [t_steps; t];
        % make a decimated sampling of Temperatures in the domain
        T_solid_steps = [T_solid_steps; T_solid(1:n/nstreams:end)];
        T_fluid_steps = [T_fluid_steps; T_fluid(1:n/nstreams:end)];
    end
    
    
end

x_points = x(1:n/nstreams:end); % make a decimated sampling of domain
time_of_contact = x_points / vv;
y = time_of_contact;
legendtext = strcat('x''=',string(num2cell(x_points)));
zz = (t_steps - time_of_contact); 

h2=figure(); hold on; grid on; box on; title('Solid');
xlabel('t'' = (t - x*z)'); ylabel('Ts');

for ind = 1:length(y)
    plot(zz(:,ind), T_solid_steps(:,ind))
end
xlim([0,zz(end,1)]);
axP = get(gca,'Position');
legend(legendtext,'Location','NorthEastOutside')
axP = get(gca,'Position');

h3=figure(); hold on; grid on; box on; title('Fluid');
xlabel('t'' = (t - x*z)'); ylabel('Tf');
for ind = 1:length(y)
    plot(zz(:,ind), T_fluid_steps(:,ind))
end
xlim([0,zz(end,1)]);
legend(legendtext,'Location','NorthEastOutside')

h4=figure();
hold on; grid on; box on; title('Location of Max(T_f - T_s) vs t');
plot(mx(:,1),mx(:,2),'k.');
%plot(mx(:,1),mx(:,2),'-');
p = polyfit(mx(startind:end,1), mx(startind:end,2), 1);
f = polyval(p, mx(:,1));
plot(mx(:,1),f,'k--');
xlabel('t'), ylabel('Location of max(T_f - T_s)')
title(['Slope = ' num2str(p(1)) ])

figure(2); set(gca,'fontname','helvetica','fontsize', [14])
f1 = [dirpath '/Tsolid.png'];
%print(f1, '-dpng')
f1 = [dirpath '/Tsolid.fig'];
savefig(f1);
figure(3); set(gca,'fontname','helvetica','fontsize', [14])
f2 = [dirpath '/Tfluid.png'];
%print(f2, '-dpng')
f2 = [dirpath '/Tfluid.fig'];
savefig(f2);
figure(4); set(gca,'fontname','helvetica','fontsize', [14])
f3 = [dirpath '/Maxt.png'];
%print(f3, '-dpng')
f3 = [dirpath '/Maxt.fig'];
savefig(f3);
f4 = [dirpath '/vars.mat'];
save(f4); % write all vars to a .mat file for later use in analysis


end




