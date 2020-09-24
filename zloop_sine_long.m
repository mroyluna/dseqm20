% 12/27/19 - run a whole matrix of runs using the schumann_all_zloop code, 
% storing the results in different folders - one for each z-value 
% and for short vs long runs.

clear all 
close all

 n        = 5000;                 % grid resolution
 nstreams = 10;                  % level of decimation when plotting
 nout     = 50;                  % number of files that are written for each case

 
 %zarr = [0.0096: 0.01: 0.2376];
 zarr = [ 0.1096];
 % do two sets of runs - one for short domain, looking at short time
 % evolution; one for longer domain, longer time evolution
 
for i = 1:length(zarr)
    
%         Lx   = 10;
%         tmax = 1.5;
%         dx   = Lx/n;
%         dt   = zarr(i)*1.0e-5*dx;          % small, for stability
%         startind = 50;
%         close all
%         schumann_all_zloop(zarr(i), Lx, tmax, dx, dt, startind);
% 
% 
%         Lx   = 100;
%         tmax = 80;
%         dx   = Lx/n;
%         dt   = zarr(i)*1.0e-4*dx;          % small, for stability
%         startind = 40;
%         close all
%         schumann_all_zloop(zarr(i), Lx, tmax, dx, dt, startind);

          Lx   = 5000;
          tmax = 15000;
          dx   = Lx/n;
          dt   = dx*1.0e-4;          % small, for stability
          startind = 40;
          period = [100]; % in dimensionless time
          close all
          for j = 1:length(period)
             schumann_all_zloop_sine(zarr(i), Lx, tmax, dx, dt, startind, period(j));
             close all
          end
end
