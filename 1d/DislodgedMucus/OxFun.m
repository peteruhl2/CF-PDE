%%% oxygen space function from cowley et al
%%% oxygen is just an exponential function over space
%%%
%%% started 11/7/2022

close all

%%% space scaling factor
xs = 500;

%%% oxygen scaling factor
ws = 200;

data7 = readmatrix("10^7 cells.csv");
data8 = readmatrix("10^8 cells.csv");
data9 = readmatrix("10^9 cells.csv");
data10 = readmatrix("10^10 cells.csv");

%%% scale everything
data7(:,1) = data7(:,1)/xs;
data7(:,2) = data7(:,2)/ws;

data8(:,1) = data8(:,1)/xs;
data8(:,2) = data8(:,2)/ws;

data9(:,1) = data9(:,1)/xs;
data9(:,2) = data9(:,2)/ws;

data10(:,1) = data10(:,1)/xs;
data10(:,2) = data10(:,2)/ws;

%%% space for curve
x = linspace(0,1);

a7 = 0.01;
y7 = exp(-a7*x);

a8 = 0.3;
y8 = exp(-a8*x);

a9 = 6;
y9 = exp(-a9*x);

a10 = 23;
y10 = exp(-a10*x);




%%% Plots =================================================================

% hold on; box on
% plot(data7(:,1),data7(:,2),'x','linewidth',2)
% plot(x,y7,'linewidth',2)
% 
% hold on; box on
% plot(data8(:,1),data8(:,2),'x','linewidth',2)
% plot(x,y8,'linewidth',2)
% 
% hold on; box on
% plot(data9(:,1),data9(:,2),'x','linewidth',2)
% plot(x,y9,'linewidth',2)

hold on; box on
plot(data10(:,1),data10(:,2),'x','linewidth',2)
plot(x,y10,'linewidth',2)
